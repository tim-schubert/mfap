#!/usr/bin/env python3
"""
analyze_patients_for_variant_specific_additional_TIS.py

For each patient in variant_list_with_mutated_sequences.csv, this script:
  1. Reconstructs the full transcript by combining the reference 5′ and 3′ UTRs
     (from Homo_sapiens_ENST00000672700_1_sequence_with_100_bp_flanks.fa in titer/data)
     with the patient’s mutated CDS (from the "Mutated_Sequence" column).
  2. Extracts candidate 203‑nt windows (with candidate codon in positions 100–102) from the full transcript.
  3. Uses an ensemble of 32 TITER model weights (from titer/model) to predict the probability 
     of translation initiation for each candidate.
  4. Computes each candidate’s location relative to the patient’s CDS (using one-based indexing) and 
     extracts the RNA sequence from the candidate TIS to the end of the CDS.
  5. Saves detailed candidate predictions for each patient in a dedicated subfolder under titer/data/individual_patients/<patient_ID>/.
  6. Summarizes for each patient:
       - Overall most likely non_canonical TIS (excluding the canonical start at CDS position 1),
       - Most likely non_canonical_in_frame TIS (i.e. candidate with CDS position ≠ 1 and in frame)
       along with the patient's DNA variant.
  7. Writes a summary CSV file (summary_patients_most_likely_additional_TIS.csv) in titer/data.
  
Assumptions:
 - The variant list CSV has columns "patient_ID", "Mutated_Sequence", and "DNA_variant".
 - The reference transcript FASTA (with 100 bp flanks) is at titer/data/Homo_sapiens_ENST00000672700_1_sequence_with_100_bp_flanks.fa.
 - Model weight files (bestmodel_0.hdf5 … bestmodel_31.hdf5) are in titer/model.
"""

import os
import sys
import csv
import numpy as np
from Bio import SeqIO

import tensorflow as tf
from tensorflow.keras.optimizers import Nadam
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, Conv1D, MaxPooling1D, LSTM
from tensorflow.keras.constraints import MaxNorm

# ----------------------------------------------------------------------
# Set up directories relative to this script file.
script_dir = os.path.dirname(os.path.realpath(__file__))
TITER_DIR = os.path.abspath(os.path.join(script_dir, os.pardir))
BASE_DIR = script_dir  # codes directory
DATA_DIR = os.path.join(TITER_DIR, "data")
MODEL_DIR = os.path.join(TITER_DIR, "model")
INDIVIDUAL_FOLDER = os.path.join(DATA_DIR, "individual_patients")
if not os.path.exists(INDIVIDUAL_FOLDER):
    os.makedirs(INDIVIDUAL_FOLDER)

# ----------------------------------------------------------------------
# Load the reference transcript to extract UTR sequences.
# Automatically detect the uploaded FASTA with flanks in DATA_DIR
flank_files = [f for f in os.listdir(DATA_DIR) if f.lower().endswith((".fa", ".fasta"))]
if not flank_files:
    sys.exit(f"No flanked FASTA found in {DATA_DIR}")
ref_fasta = os.path.join(DATA_DIR, flank_files[0])
ref_record = SeqIO.read(ref_fasta, "fasta")
ref_seq = str(ref_record.seq).upper()
# Assume first 100 bp and last 100 bp of this sequence are the UTRs
utr5 = ref_seq[:100]
utr3 = ref_seq[-100:]
print(f"Loaded reference transcript '{flank_files[0]}' (length {len(ref_seq)} bp) and extracted 5' and 3' UTRs.")

# ----------------------------------------------------------------------
# Define candidate codons and window parameters.
candidate_codons = ["ATG", "CTG", "GTG", "TTG", "ACG", "AAG", "AGG", "ATA", "ATT", "ATC"]
WINDOW_LENGTH = 203
CODON_START_IN_WINDOW = 100

def extract_candidate_windows(full_seq):
    """Extract candidate windows from a full transcript."""
    candidate_windows = []
    n = len(full_seq)
    for i in range(n - 2):
        codon = full_seq[i:i+3]
        if codon in candidate_codons:
            start = i - CODON_START_IN_WINDOW
            end = start + WINDOW_LENGTH
            if start >= 0 and end <= n:
                window_seq = full_seq[start:end]
                candidate_windows.append((i, codon, window_seq))
    return candidate_windows

def seq_matrix(seq_list):
    """Convert a list of sequences to a one-hot encoded numpy array."""
    num_seqs = len(seq_list)
    tensor = np.zeros((num_seqs, WINDOW_LENGTH, 8))
    for i, seq in enumerate(seq_list):
        for j, s in enumerate(seq):
            if j < CODON_START_IN_WINDOW or j > CODON_START_IN_WINDOW+2:
                if s == 'A':
                    tensor[i, j] = [1, 0, 0, 0, 0, 0, 0, 0]
                elif s == 'T':
                    tensor[i, j] = [0, 1, 0, 0, 0, 0, 0, 0]
                elif s == 'C':
                    tensor[i, j] = [0, 0, 1, 0, 0, 0, 0, 0]
                elif s == 'G':
                    tensor[i, j] = [0, 0, 0, 1, 0, 0, 0, 0]
                else:
                    tensor[i, j] = [0, 0, 0, 0, 0, 0, 0, 0]
            else:
                if s == 'A':
                    tensor[i, j] = [0, 0, 0, 0, 1, 0, 0, 0]
                elif s == 'T':
                    tensor[i, j] = [0, 0, 0, 0, 0, 1, 0, 0]
                elif s == 'C':
                    tensor[i, j] = [0, 0, 0, 0, 0, 0, 1, 0]
                elif s == 'G':
                    tensor[i, j] = [0, 0, 0, 0, 0, 0, 0, 1]
                else:
                    tensor[i, j] = [0, 0, 0, 0, 0, 0, 0, 0]
    return tensor

def build_model():
    """Build and return the TITER model architecture."""
    model = Sequential()
    model.add(Conv1D(
        filters=128,
        kernel_size=3,
        input_shape=(WINDOW_LENGTH, 8),
        padding='valid',
        kernel_constraint=MaxNorm(3),
        activation='relu',
        strides=1))
    model.add(MaxPooling1D(pool_size=3))
    model.add(Dropout(rate=0.21370950078747658))
    model.add(LSTM(units=256, return_sequences=True))
    model.add(Dropout(rate=0.7238091317104384))
    model.add(Flatten())
    model.add(Dense(1))
    model.add(Activation('sigmoid'))
    return model

# Load the codon prior dictionary.
DICT_PATH = os.path.join(BASE_DIR, "dict_piror_front_Gaotrain.npy")
codon_tis_prior = np.load(DICT_PATH, allow_pickle=True).item()

# ----------------------------------------------------------------------
# Open the CSV file with patient variant information.
variants_csv = os.path.join(DATA_DIR, "variant_list_with_mutated_sequences.csv")
with open(variants_csv, "r") as csvfile:
    reader = csv.DictReader(csvfile)
    variant_rows = list(reader)

if not variant_rows:
    sys.exit(f"No patient variant data found in {variants_csv}")

# Create a list to store summary information for each patient.
summary_rows = []

# Process each patient.
for idx, row in enumerate(variant_rows):
    patient_id = row.get("patient_ID", f"patient_{idx+1}").strip()
    dna_variant = row.get("DNA_variant", "").strip()
    mutated_sequence = row["Mutated_Sequence"].strip().upper()
    print(f"\nProcessing patient {patient_id} with mutated CDS of length {len(mutated_sequence)} bp.")
    
    full_patient_seq = utr5 + mutated_sequence + utr3
    cds_start = 100
    cds_end = 100 + len(mutated_sequence)
    
    candidate_windows_patient = extract_candidate_windows(full_patient_seq)
    print(f"  Found {len(candidate_windows_patient)} candidate windows for patient {patient_id}.")
    
    patient_folder = os.path.join(INDIVIDUAL_FOLDER, patient_id)
    if not os.path.exists(patient_folder):
        os.makedirs(patient_folder)
    
    candidate_windows_csv_patient = os.path.join(patient_folder, "candidate_windows.csv")
    with open(candidate_windows_csv_patient, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Candidate_Position", "Candidate_Codon", "Window_Sequence"])
        for pos, cand, win_seq in candidate_windows_patient:
            writer.writerow([pos, cand, win_seq])
    
    if not candidate_windows_patient:
        print(f"  No candidate windows found for patient {patient_id}; recording summary as 'none found'.")
        summary_rows.append({
            "patient_ID": patient_id,
            "DNA_variant": dna_variant,
            "Most_likely_non_canonical_TIS_CDS_Position": "none found",
            "RNA_sequence_most_likely_non_canonical_TIS": "none found",
            "Most_likely_non_canonical_in_frame_TIS_CDS_Position": "none found",
            "RNA_sequence_most_likely_non_canonical_in_frame_TIS": "none found"
        })
        continue
    
    window_seqs_patient = [entry[2] for entry in candidate_windows_patient]
    X_candidates_patient = seq_matrix(window_seqs_patient)
    
    model = build_model()
    model.compile(loss='binary_crossentropy', optimizer=Nadam(), metrics=['accuracy'])
    
    num_candidates = X_candidates_patient.shape[0]
    y_pred_n = np.zeros((num_candidates, 1))
    y_pred_p = np.zeros((num_candidates, 1))
    candidate_priors_patient = np.array([codon_tis_prior.get(cand, 1) for (_, cand, _) in candidate_windows_patient]).reshape((-1,1))
    
    print("  Predicting on candidate windows...")
    for i in range(32):
        weight_path = os.path.join(MODEL_DIR, f"bestmodel_{i}.hdf5")
        model.load_weights(weight_path)
        y_pred = model.predict(X_candidates_patient, verbose=0)
        y_pred_n += y_pred
        y_pred_p += y_pred * candidate_priors_patient
    y_pred_n /= 32
    y_pred_p /= 32
    
    candidate_cds_locations = []
    rna_sequence_from_pot_noncanonical_TIS_list = []
    for cp, cand, win_seq in candidate_windows_patient:
        if cp >= cds_start and cp <= cds_end - 3:
            candidate_cds_locations.append(cp - cds_start + 1)
        else:
            candidate_cds_locations.append("outside CDS")
        rna_seq = full_patient_seq[cp:cds_end]
        rna_sequence_from_pot_noncanonical_TIS_list.append(rna_seq)
    
    output_csv_patient = os.path.join(patient_folder, "candidate_predictions.csv")
    with open(output_csv_patient, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Candidate_Position_Full", "Candidate_Codon", "CDS_Position", "Window_Sequence",
                         "Probability_Without_Prior", "Probability_With_Prior", "RNA_sequence_from_pot_non_canonical_TIS"])
        for ((cp, cand, win_seq), cds_loc, pred_n, pred_p, rna_seq) in zip(
                candidate_windows_patient, candidate_cds_locations, y_pred_n, y_pred_p, rna_sequence_from_pot_noncanonical_TIS_list):
            writer.writerow([cp, cand, cds_loc, win_seq, pred_n[0], pred_p[0], rna_seq])
    
    print(f"  Predictions for patient {patient_id} saved to {output_csv_patient}.")
    
    best_candidate = None
    best_prob = -1.0
    best_candidate_in_frame = None
    best_prob_in_frame = -1.0
    for (cp, cand, win_seq), cds_loc, pred_p, rna_seq in zip(
            candidate_windows_patient, candidate_cds_locations, y_pred_p, rna_sequence_from_pot_noncanonical_TIS_list):
        try:
            prob = float(pred_p)
        except:
            prob = 0
        if prob > best_prob:
            best_prob = prob
            best_candidate = {"CDS_Position": cds_loc, "RNA_sequence_from_pot_non_canonical_TIS": rna_seq}
        if isinstance(cds_loc, int) and cds_loc != 1 and (cds_loc - 1) % 3 == 0:
            if prob > best_prob_in_frame:
                best_prob_in_frame = prob
                best_candidate_in_frame = {"CDS_Position": cds_loc, "RNA_sequence_from_pot_non_canonical_TIS": rna_seq}
    
    if best_candidate_in_frame is None:
        best_candidate_in_frame = {"CDS_Position": "none found", "RNA_sequence_from_pot_non_canonical_TIS": "none found"}
    
    summary_rows.append({
        "patient_ID": patient_id,
        "DNA_variant": dna_variant,
        "Most_likely_non_canonical_TIS_CDS_Position": best_candidate["CDS_Position"] if best_candidate else "none found",
        "RNA_sequence_most_likely_non_canonical_TIS": best_candidate["RNA_sequence_from_pot_non_canonical_TIS"] if best_candidate else "none found",
        "Most_likely_non_canonical_in_frame_TIS_CDS_Position": best_candidate_in_frame["CDS_Position"],
        "RNA_sequence_most_likely_non_canonical_in_frame_TIS": best_candidate_in_frame["RNA_sequence_from_pot_non_canonical_TIS"]
    })

# ----------------------------------------------------------------------
# Write summary CSV file for all patients.
summary_csv = os.path.join(DATA_DIR, "summary_patients_most_likely_additional_TIS.csv")
with open(summary_csv, "w", newline="") as csvfile:
    fieldnames = ["patient_ID", "DNA_variant", 
                  "Most_likely_non_canonical_TIS_CDS_Position", "RNA_sequence_most_likely_non_canonical_TIS",
                  "Most_likely_non_canonical_in_frame_TIS_CDS_Position", "RNA_sequence_most_likely_non_canonical_in_frame_TIS"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in summary_rows:
        writer.writerow(row)

print(f"\nSummary of patient TIS predictions saved to {summary_csv}.")
print("\nAll patient variants processed.")