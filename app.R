library(shiny)
library(shinyjs)
library(DT)
library(stringr)
library(dplyr)
library(Biostrings)
library(bslib)

titer_python_modules <- c("numpy", "tensorflow", "Bio")
titer_summary_columns <- c(
  "patient_ID",
  "Most_likely_non_canonical_TIS_CDS_Position",
  "RNA_sequence_most_likely_non_canonical_TIS",
  "Most_likely_non_canonical_in_frame_TIS_CDS_Position",
  "RNA_sequence_most_likely_non_canonical_in_frame_TIS"
)
titer_character_columns <- c(
  "Most_likely_non_canonical_TIS_CDS_Position",
  "RNA_sequence_most_likely_non_canonical_TIS",
  "Most_likely_non_canonical_in_frame_TIS_CDS_Position",
  "RNA_sequence_most_likely_non_canonical_in_frame_TIS",
  "Protein_from_most_likely_non_canonical_TIS",
  "Protein_from_most_likely_non_canonical_in_frame_TIS"
)
titer_integer_columns <- c(
  "Protein_Length_from_most_likely_non_canonical_TIS",
  "Protein_Length_from_most_likely_non_canonical_in_frame_TIS",
  "Unaltered_Length_from_most_likely_non_canonical_TIS",
  "Frameshift_Length_from_most_likely_non_canonical_TIS",
  "Unaltered_Length_from_most_likely_non_canonical_in_frame_TIS",
  "Frameshift_Length_from_most_likely_non_canonical_in_frame_TIS"
)

ensure_titer_python_ready <- function(required_modules = titer_python_modules) {
  py_ready <- tryCatch(
    reticulate::py_available(initialize = TRUE),
    error = function(e) FALSE
  )
  if (!isTRUE(py_ready)) {
    return(list(ok = FALSE, message = "Python not available for TITER in this session."))
  }

  missing <- required_modules[!vapply(required_modules, reticulate::py_module_available, logical(1))]
  if (length(missing) > 0) {
    return(list(
      ok = FALSE,
      message = paste0(
        "Missing Python modules for TITER: ",
        paste(missing, collapse = ", "),
        ". Ensure deployment installs dependencies from requirements.txt."
      )
    ))
  }

  list(ok = TRUE, message = NULL)
}

add_missing_titer_columns <- function(df) {
  for (nm in titer_character_columns) {
    if (!nm %in% names(df)) df[[nm]] <- NA_character_
  }
  for (nm in titer_integer_columns) {
    if (!nm %in% names(df)) df[[nm]] <- NA_integer_
  }
  df
}

# parse_mutation
parse_mutation <- function(mutation, wildtype_seq) {
  info <- list(
    Locus=NA_integer_, Mutation_Type=NA_character_, New_Base=NA_character_,
    Original_Base=NA_character_, Deleted_Bases=NA_character_,
    Duplicated_Bases=NA_character_, Inserted_Bases=NA_character_,
    Mutated_Sequence=NA_character_
  )
  if (is.na(mutation) || !nzchar(trimws(mutation))) return(info)
  if (is.na(wildtype_seq) || !nzchar(wildtype_seq)) return(info)
  mutation <- trimws(mutation)
  if (str_starts(mutation, "c.")) mutation <- substring(mutation, 3)
  if (str_detect(mutation, ">")) {
    parts <- str_split(mutation, ">")[[1]]
    if (length(parts)==2) {
      left <- parts[1]; nb <- parts[2]
      num <- str_sub(left, 1, -2)
      if (str_detect(num, "^[0-9]+$") && str_detect(nb, "^[ACGT]$")) {
        loc <- as.integer(num)
        if (loc>=1 && loc<=str_length(wildtype_seq)) {
          orig <- str_sub(wildtype_seq, loc, loc)
          mut  <- paste0(str_sub(wildtype_seq, 1, loc-1), nb, str_sub(wildtype_seq, loc+1))
          info <- modifyList(info, list(
            Locus=loc, Mutation_Type="Point", New_Base=nb, Original_Base=orig,
            Mutated_Sequence=mut
          ))
        }
      }
    }
  } else if (str_detect(mutation, "del") && !str_detect(mutation, "delins")) {
    parts <- str_split(mutation, "del")[[1]]
    rp <- parts[1]; extra <- if (length(parts)>1) parts[2] else ""
    if (str_detect(rp, "_")) {
      pos <- as.integer(str_split(rp, "_")[[1]])
      if (length(pos)==2 && !any(is.na(pos))) {
        db  <- str_sub(wildtype_seq, pos[1], pos[2])
        mut <- paste0(str_sub(wildtype_seq, 1, pos[1]-1), str_sub(wildtype_seq, pos[2]+1))
        info <- modifyList(info, list(Locus=pos[1], Mutation_Type="Deletion", Deleted_Bases=db, Mutated_Sequence=mut))
      }
    } else {
      p <- as.integer(rp)
      if (!is.na(p)) {
        db  <- if (extra!="" && str_detect(extra,"^[ACGT]$")) extra else str_sub(wildtype_seq, p, p)
        mut <- paste0(str_sub(wildtype_seq, 1, p-1), str_sub(wildtype_seq, p+1))
        info <- modifyList(info, list(Locus=p, Mutation_Type="Deletion", Deleted_Bases=db, Mutated_Sequence=mut))
      }
    }
  } else if (str_detect(mutation, "dup")) {
    parts <- str_split(mutation, "dup")[[1]]
    left <- parts[1]
    if (str_detect(left, "_")) {
      pos <- as.integer(str_split(left, "_")[[1]])
      if (length(pos)==2 && !any(is.na(pos))) {
        db  <- str_sub(wildtype_seq, pos[1], pos[2])
        mut <- paste0(str_sub(wildtype_seq, 1, pos[2]), db, str_sub(wildtype_seq, pos[2]+1))
        info <- modifyList(info, list(Locus=pos[1], Mutation_Type="Duplication", Duplicated_Bases=db, Mutated_Sequence=mut))
      }
    } else {
      p <- as.integer(left)
      if (!is.na(p)) {
        db  <- str_sub(wildtype_seq, p, p)
        mut <- paste0(str_sub(wildtype_seq, 1, p), db, str_sub(wildtype_seq, p+1))
        info <- modifyList(info, list(Locus=p, Mutation_Type="Duplication", Duplicated_Bases=db, Mutated_Sequence=mut))
      }
    }
  } else if (str_detect(mutation, "delins")) {
    parts <- str_split(mutation, "delins")[[1]]
    rp <- parts[1]; ins <- parts[2]
    pos <- as.integer(str_split(rp, "_")[[1]]); if (length(pos)==1) pos <- c(pos, pos)
    if (length(pos)==2 && !any(is.na(pos))) {
      db  <- str_sub(wildtype_seq, pos[1], pos[2])
      mut <- paste0(str_sub(wildtype_seq, 1, pos[1]-1), ins, str_sub(wildtype_seq, pos[2]+1))
      info <- modifyList(info, list(Locus=pos[1], Mutation_Type="Indel", Deleted_Bases=db, Inserted_Bases=ins, Mutated_Sequence=mut))
    }
  } else if (str_detect(mutation, "ins") && !str_detect(mutation, "delins")) {
    parts <- str_split(mutation, "ins")[[1]]
    if (length(parts)==2) {
      left <- parts[1]; ins <- parts[2]
      pos_vals <- as.integer(str_split(left, "_")[[1]])
      if (length(pos_vals)==1 && !is.na(pos_vals)) { start <- pos_vals; end <- pos_vals
      } else if (length(pos_vals)==2 && !any(is.na(pos_vals))) { start <- pos_vals[1]; end <- pos_vals[2]
      } else return(info)
      mut <- paste0(str_sub(wildtype_seq, 1, start), ins, str_sub(wildtype_seq, end+1))
      info <- modifyList(info, list(Locus=start, Mutation_Type="Insertion", Inserted_Bases=ins, Mutated_Sequence=mut))
    }
  }
  info
}

# translation
translate_prot <- function(dna_seq) {
  if (is.na(dna_seq) || nchar(dna_seq)<3) return(NA_character_)
  prot_full <- suppressWarnings(as.character(translate(DNAString(dna_seq), if.fuzzy.codon="X")))
  if (grepl("\\*", prot_full)) substr(prot_full, 1, regexpr("\\*", prot_full)[1]-1) else prot_full
}

# Refinement helpers
refine_point <- function(locus, seq, wt_seq, codon_tbl) {
  idx <- floor((locus-1)/3)*3 + 1
  if (is.na(locus) || (idx+2)>nchar(wt_seq) || (idx+2)>nchar(seq)) return("Point")
  wt_cod  <- toupper(substr(wt_seq, idx, idx+2))
  mut_cod <- toupper(substr(seq,    idx, idx+2))
  if (!wt_cod %in% names(codon_tbl) || !mut_cod %in% names(codon_tbl)) return("Point")
  wt_aa  <- codon_tbl[[wt_cod]]; mut_aa <- codon_tbl[[mut_cod]]
  if (wt_aa==mut_aa) "Silent" else if (mut_aa=="*") "Nonsense" else "Missense"
}
refine_indel <- function(mt, del, dup, ins) {
  net <- switch(mt, Deletion=-nchar(del), Duplication=nchar(dup), Insertion=nchar(ins), Indel=nchar(ins)-nchar(del), 0)
  if (length(net) != 1 || is.na(net)) return("Frameshifting indel")
  if (net %% 3 == 0) "In-Frame indel" else "Frameshifting indel"
}

# UI
mfap <- function(){
  tags$span(
    class="mfap-brand",
    tags$span("M"),
    tags$span(class="mfap-f","f"),
    tags$span("AP")
  )
}


ui <- tagList(
  bslib::page_navbar(
    fillable = FALSE,
    useShinyjs(),
    id = "topnav",
    title = tags$div(
      
      class = "brand",
      style = "display:flex;flex-direction:column;line-height:1.05;align-items:flex-start;padding-right:32px;",
      tags$span(
        class = "brand-title",
        tags$span("M"),
        tags$span(class="mfap-f","f"),
        tags$span("AP")
      ),
      
      tags$span(class="brand-subtitle","Molecular feature Association Pipeline")
    ),
    theme = bslib::bs_theme(version = 5, bootswatch = "darkly"),
    window_title = "MfAP — Molecular feature Association Pipeline",
    header = tags$head(
      tags$link(
        rel = "stylesheet",
        href = "https://fonts.googleapis.com/css2?family=Libre+Bodoni:ital,wght@0,700;0,800;1,700;1,800&display=swap"
      ),
      tags$script(HTML("
      Shiny.addCustomMessageHandler('mfap_intro_modal_class', function(msg){
        var m = document.getElementById('shiny-modal');
        if (!m) return;
      
        if (msg && msg.on) {
          m.classList.add('mfap-intro-modal');
      
          try {
            $(m).off('hidden.bs.modal.mfapCenterFix').on('hidden.bs.modal.mfapCenterFix', function(){
              m.classList.remove('mfap-intro-modal');
              $(m).off('hidden.bs.modal.mfapCenterFix');
            });
          } catch(e) {}
        }
      });

    ")),
      
      
      
      
      tags$style(HTML("
    
      .page-container{max-width:1180px;margin:0 auto;padding:0 14px;}
      .navbar .container-fluid{max-width:1180px;margin:0 auto;padding:0 14px;}
      .lit-item p{ margin:0 0 6px 0; }
      .navbar{background:var(--mfap-navbar-bg) !important;border-bottom:1px solid var(--mfap-navbar-border)!important;box-shadow:none !important;}
      .navbar-brand{margin-right:auto;}
      :root{
        --mfap-font: 'Libre Bodoni', Georgia, 'Times New Roman', serif;
        --mfap-navbar-bg: #1d1f21;
        --mfap-navbar-border: rgba(255,255,255,.08);
        --mfap-brand-title: #ffffff;
        --mfap-brand-subtitle: #e5e7eb;
        --mfap-link: #4dabf7;
        --mfap-link-hover: #74c0fc;
        --mfap-muted: #adb5bd;
        --mfap-nav-link: #e9ecef;
        --mfap-card-bg: #202326;
        --mfap-card-header-bg: #1b1e22;
        --mfap-card-border: rgba(255,255,255,0.08);
        --mfap-card-header-border: rgba(255,255,255,0.06);
        --mfap-action-dock-bg: rgba(32,35,38,0.35);
        --mfap-action-dock-border: rgba(255,255,255,0.12);
        --mfap-action-dock-shadow: 0 0 16px rgba(0,0,0,0.16);
        --mfap-details-bg: #262a2e;
        --mfap-details-bg-open: #2d3236;
        --mfap-details-border: rgba(255,255,255,0.08);
        --mfap-help-trigger: #6c757d;
        --mfap-help-trigger-hover: #adb5bd;
        --mfap-disclaimer-bg: rgba(255,255,255,0.03);
        --mfap-disclaimer-border: rgba(255,255,255,0.1);
        --mfap-disclaimer-text: #868e96;
        --mfap-disclaimer-sep: #495057;
        --mfap-about-lead: #cfd4da;
        --mfap-about-notes-bg: #23272b;
        --mfap-about-notes-border: rgba(255,255,255,.06);
        --mfap-bib-bg: #1f1f1f;
        --mfap-bib-border: #2a2a2a;
        --mfap-modal-bg: #24272b;
        --mfap-modal-border: rgba(255,255,255,0.14);
        --mfap-footer-color: #ffffff;
        --mfap-intro-bg: rgba(255,255,255,0.06);
        --mfap-intro-border: rgba(255,255,255,.14);
        --mfap-intro-text: #e9ecef;
        --mfap-intro-hover-text: #ffffff;
        --bs-body-bg: #2b3138;
        --mfap-form-bg: #25282c;
        --mfap-form-border: rgba(255,255,255,0.14);
        --mfap-form-color: #e9ecef;
        --mfap-form-focus-border: rgba(255,255,255,0.25);
        --mfap-btn-info-bg: var(--mfap-link);
        --mfap-btn-info-hover: var(--mfap-link-hover);
        --mfap-btn-info-active: var(--mfap-link-hover);
        --mfap-btn-info-color: #fff;
        --mfap-toast-message-bg: #2563eb;
        --mfap-toast-message-border: #1d4ed8;
        --mfap-toast-message-text: #ffffff;
      }

      html[data-mfap-theme='light']{
        --mfap-navbar-bg: #ffffff;
        --mfap-navbar-border: rgba(17,24,39,.10);
        --mfap-brand-title: #111827;
        --mfap-brand-subtitle: #4b5563;
        --mfap-link: #2563eb;
        --mfap-link-hover: #1d4ed8;
        --mfap-muted: #6b7280;
        --mfap-nav-link: #1f2937;
        --mfap-card-bg: #ffffff;
        --mfap-card-header-bg: #f7f9fc;
        --mfap-card-border: rgba(17,24,39,0.10);
        --mfap-card-header-border: rgba(17,24,39,0.08);
        --mfap-action-dock-bg: rgba(255,255,255,0.76);
        --mfap-action-dock-border: rgba(17,24,39,0.14);
        --mfap-action-dock-shadow: 0 0 16px rgba(15,23,42,0.07);
        --mfap-details-bg: #f8fafc;
        --mfap-details-bg-open: #f1f5f9;
        --mfap-details-border: rgba(17,24,39,0.12);
        --mfap-help-trigger: #4b5563;
        --mfap-help-trigger-hover: #1f2937;
        --mfap-disclaimer-bg: rgba(15,23,42,0.04);
        --mfap-disclaimer-border: rgba(15,23,42,0.18);
        --mfap-disclaimer-text: #4b5563;
        --mfap-disclaimer-sep: #94a3b8;
        --mfap-about-lead: #334155;
        --mfap-about-notes-bg: #f8fafc;
        --mfap-about-notes-border: rgba(17,24,39,0.10);
        --mfap-bib-bg: #f8fafc;
        --mfap-bib-border: #cbd5e1;
        --mfap-modal-bg: #eef2f7;
        --mfap-modal-border: rgba(17,24,39,0.18);
        --mfap-footer-color: #1f2937;
        --mfap-intro-bg: rgba(15,23,42,0.05);
        --mfap-intro-border: rgba(15,23,42,0.18);
        --mfap-intro-text: #1f2937;
        --mfap-intro-hover-text: #1f2937;
        --bs-body-bg: #e1e7ef;
        --bs-body-color: #1f2937;
        --bs-border-color: #dbe2ea;
        --mfap-btn-info-bg: #4dabf7;
        --mfap-btn-info-hover: #74c0fc;
        --mfap-btn-info-active: #74c0fc;
        --mfap-btn-info-color: #fff;
        --mfap-toast-message-bg: #2563eb;
        --mfap-toast-message-border: #1d4ed8;
        --mfap-toast-message-text: #ffffff;
      }
      
      #shiny-modal.mfap-intro-modal{
        display: block !important;
      }
      
      @media (min-height: 700px){
        #shiny-modal.mfap-intro-modal{
          display: flex !important;
          align-items: center !important;
          justify-content: center !important;
        }
      
        #shiny-modal.mfap-intro-modal .modal-dialog{
          margin: 0 !important;
        }
      }
      
      #shiny-modal.mfap-intro-modal .modal-dialog{
        max-width: 760px;
        width: calc(100% - 2rem);
        margin: 1rem auto !important; /* top/bottom breathing room when not centered */
      }
      
      #shiny-modal.mfap-intro-modal .modal-content{
        max-height: calc(100vh - 2rem);
        overflow: auto;
        -webkit-overflow-scrolling: touch;
      }
      .modal-content{
        background: var(--mfap-modal-bg);
        color: var(--bs-body-color);
        border: 1px solid var(--mfap-modal-border);
      }
      .modal-header{border-bottom:1px solid var(--mfap-modal-border);}
      .modal-footer{border-top:1px solid var(--mfap-modal-border);}
      .modal-title{color:var(--bs-body-color);}

      
      .mfap-brand{
        font-family: var(--mfap-font);
        font-weight: 700;
        letter-spacing: .12px;
        color: inherit;
      }
      
      .mfap-brand .mfap-f{
        font-style: italic;
      }

      .brand-title{
        font-family: var(--mfap-font);
        font-weight: 800;
        letter-spacing: .2px;
        color: var(--mfap-brand-title) !important;
        font-size: 2.2rem;
        line-height: 1.05;
      }
      
      @media (max-width: 768px){
        .brand-title{
          font-size: 1.75rem;
        }
      }
      
      .brand-subtitle{
        font-family: inherit;
        font-weight: 400;
        font-size: 1.05rem;
        color: var(--mfap-brand-subtitle) !important;
        margin-top: 2px;
        letter-spacing: 0;
      }

      .mfap-f{
        font-style: italic;
        font-weight: inherit;
        margin: 0 -0.03em;
      }



      a{
        color:var(--mfap-link);
        text-decoration: none;
      }
      a:hover, a:focus{
        color:var(--mfap-link-hover);
        text-decoration: none;
      }
      .site-footer a{
        color:var(--mfap-muted) !important;
        text-decoration:none;
      }
      
      .site-footer a:hover,
      .site-footer a:focus{
        color:var(--mfap-link) !important;
        text-decoration:none;
      }

      .navbar-nav{ margin-left:auto; gap: 30px; }
      .navbar .navbar-toggler,
      .navbar .navbar-toggle{
        border:0 !important;
        background:transparent !important;
        background-color:transparent !important;
        background-image:none !important;
        box-shadow:none !important;
        outline:none !important;
        padding:.3rem !important;
        appearance:none !important;
        -webkit-appearance:none !important;
      }
      .navbar .navbar-toggler:focus,
      .navbar .navbar-toggler:focus-visible,
      .navbar .navbar-toggler:hover,
      .navbar .navbar-toggler:active,
      .navbar .navbar-toggle:focus,
      .navbar .navbar-toggle:hover,
      .navbar .navbar-toggle:active{
        border:0 !important;
        background:transparent !important;
        background-color:transparent !important;
        background-image:none !important;
        box-shadow:none !important;
        outline:none !important;
      }
      .navbar .navbar-toggler-icon{
        width:1.15rem !important;
        height:1.15rem !important;
        background-color:transparent !important;
        background:none !important;
        background-repeat:no-repeat !important;
        background-position:center !important;
        background-size:100% 100% !important;
        background-image:url(\"data:image/svg+xml,%3csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 30 30'%3e%3cpath fill='none' stroke='rgba(107,114,128,0.95)' stroke-linecap='round' stroke-miterlimit='10' stroke-width='2.0' d='M4 7h22M4 15h22M4 23h22'/%3e%3c/svg%3e\") !important;
      }
      .navbar .navbar-toggle .icon-bar{
        display:block !important;
        position:static !important;
        width:20px !important;
        height:2px !important;
        margin:0 !important;
        border-radius:999px !important;
        background:rgba(107,114,128,0.95) !important;
        background-color:rgba(107,114,128,0.95) !important;
        box-shadow:none !important;
        transform:none !important;
      }
      .navbar .navbar-toggle .icon-bar + .icon-bar{
        margin-top:5px !important;
      }
      html[data-mfap-theme='light'] .navbar-toggler-icon,
      html[data-bs-theme='light'] .navbar-toggler-icon{
        background-image:url(\"data:image/svg+xml,%3csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 30 30'%3e%3cpath fill='none' stroke='rgba(107,114,128,0.95)' stroke-linecap='round' stroke-miterlimit='10' stroke-width='2.0' d='M4 7h22M4 15h22M4 23h22'/%3e%3c/svg%3e\") !important;
      }
      @media (max-width: 768px){
        .brand-title{   font-size: 1.55rem; }
        .brand-subtitle{font-size: 0.98rem; }
        .navbar-collapse{ margin-top:8px; }
        .navbar-nav{
          margin-left:0;
          gap:0 !important;
        }
        .navbar-nav .nav-link{
          padding:.4rem .35rem;
          font-size:.90rem;
        }
        .theme-toggle-link{padding:.25rem .45rem !important; line-height:1;}
      }
      .navbar-nav .nav-link{padding:.5rem .6rem;color:var(--mfap-nav-link) !important;border-radius:0 !important;border-bottom:2px solid transparent;}
      .navbar-nav .nav-link:hover{color:var(--mfap-link) !important;background:transparent !important;border-bottom-color:var(--mfap-link);}
      .navbar-nav .nav-link.active{color:var(--mfap-link) !important;background:transparent !important;border-bottom-color:var(--mfap-link);}

      .card{border:1px solid var(--mfap-card-border)!important;border-radius:10px!important;background:var(--mfap-card-bg);}
      .card-header{background:var(--mfap-card-header-bg)!important;border-bottom:1px solid var(--mfap-card-header-border)!important;color:var(--bs-body-color)!important;padding:10px 14px!important;font-weight:600;}
      .card-body{padding:14px!important;}
      .card-footer{background:var(--mfap-card-bg);border-top:1px solid var(--mfap-card-header-border);padding:12px 14px;}

      html:not([data-mfap-theme='light']) .form-control,
      html:not([data-mfap-theme='light']) .shiny-input-container input[type='text'].form-control,
      html:not([data-mfap-theme='light']) .shiny-input-container input[type='number'].form-control{
        background-color:var(--mfap-form-bg) !important;
        border-color:var(--mfap-form-border) !important;
        color:var(--mfap-form-color) !important;
      }
      html:not([data-mfap-theme='light']) .form-control:focus{
        background-color:var(--mfap-form-bg) !important;
        border-color:var(--mfap-form-focus-border) !important;
        color:var(--mfap-form-color) !important;
        box-shadow:0 0 0 0.2rem rgba(255,255,255,0.08) !important;
      }
      html:not([data-mfap-theme='light']) .form-control::placeholder{
        color:var(--mfap-muted) !important;
      }
      html:not([data-mfap-theme='light']) .input-group .form-control{
        background-color:var(--mfap-form-bg) !important;
        border-color:var(--mfap-form-border) !important;
        color:var(--mfap-form-color) !important;
      }

      .sidebar{position:sticky;top:16px;max-height:calc(100vh - 32px);overflow:auto;}
      .action-dock{
          position: sticky;
          bottom: 14px;
          background: var(--mfap-action-dock-bg);
          backdrop-filter: blur(8px) saturate(120%);
          -webkit-backdrop-filter: blur(8px) saturate(120%); /* Safari */
          border: 1px solid var(--mfap-action-dock-border);
          box-shadow: var(--mfap-action-dock-shadow);
          padding: 10px 12px;
          border-radius: 12px;
          margin-top: 10px;
          margin-bottom: 18px;
          z-index: 2;
        }
      .action-row{display:flex;gap:8px;}
      .action-dock .btn{width:100%;}
      .btn-mfap-primary{background-color:#2563eb !important;border-color:#2563eb !important;color:#fff !important;}
      .btn-mfap-primary:hover,.btn-mfap-primary:focus{background-color:#1d4ed8 !important;border-color:#1d4ed8 !important;color:#fff !important;}
      #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button,
      #res_tbl .dataTables_wrapper .dataTables_paginate .page-link,
      #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button.btn,
      #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button.btn-success{
        background:#495057 !important;
        border-color:#495057 !important;
        color:#fff !important;
        box-shadow:none !important;
      }
      #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button.current,
      #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button.current.btn-success,
      #res_tbl .dataTables_wrapper .dataTables_paginate .page-item.active .page-link{
        background:#4dabf7 !important;
        border-color:#4dabf7 !important;
        color:#fff !important;
        box-shadow:none !important;
      }
      #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button:hover:not(.disabled),
      #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button:focus:not(.disabled),
      #res_tbl .dataTables_wrapper .dataTables_paginate .page-item:not(.disabled):not(.active) .page-link:hover,
      #res_tbl .dataTables_wrapper .dataTables_paginate .page-item:not(.disabled):not(.active) .page-link:focus{
        background:#4dabf7 !important;
        border-color:#4dabf7 !important;
        color:#fff !important;
      }
      #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button.disabled,
      #res_tbl .dataTables_wrapper .dataTables_paginate .page-item.disabled .page-link{
        background:#343a40 !important;
        border-color:#343a40 !important;
        color:rgba(255,255,255,0.5) !important;
        cursor:default !important;
        opacity:1 !important;
      }
      html[data-mfap-theme='light'] #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button,
      html[data-mfap-theme='light'] #res_tbl .dataTables_wrapper .dataTables_paginate .page-link,
      html[data-mfap-theme='light'] #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button.btn,
      html[data-mfap-theme='light'] #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button.btn-success{
        background:#94a3b8 !important;
        border-color:#94a3b8 !important;
        color:#fff !important;
      }
      html[data-mfap-theme='light'] #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button.current,
      html[data-mfap-theme='light'] #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button.current.btn-success,
      html[data-mfap-theme='light'] #res_tbl .dataTables_wrapper .dataTables_paginate .page-item.active .page-link{
        background:#2563eb !important;
        border-color:#2563eb !important;
        color:#fff !important;
      }
      html[data-mfap-theme='light'] #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button:hover:not(.disabled),
      html[data-mfap-theme='light'] #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button:focus:not(.disabled),
      html[data-mfap-theme='light'] #res_tbl .dataTables_wrapper .dataTables_paginate .page-item:not(.disabled):not(.active) .page-link:hover,
      html[data-mfap-theme='light'] #res_tbl .dataTables_wrapper .dataTables_paginate .page-item:not(.disabled):not(.active) .page-link:focus{
        background:#3b82f6 !important;
        border-color:#3b82f6 !important;
        color:#fff !important;
      }
      html[data-mfap-theme='light'] #res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button.disabled,
      html[data-mfap-theme='light'] #res_tbl .dataTables_wrapper .dataTables_paginate .page-item.disabled .page-link{
        background:#cbd5e1 !important;
        border-color:#cbd5e1 !important;
        color:#64748b !important;
      }
      
      /* Results DT: header-driven sizing + ellipsis truncation */
      #res_tbl table.dataTable{
        table-layout: auto !important; 
        width: 100% !important;
      }
      
      #res_tbl table.dataTable th,
      #res_tbl table.dataTable td{
        white-space: nowrap;          
      }
      
      #res_tbl table.dataTable td{
        overflow: hidden;          
        text-overflow: ellipsis;
        max-width: 1px;              
      }


      details{background:var(--mfap-details-bg);border:1px solid var(--mfap-details-border);border-radius:8px;padding:8px 10px;}
      details>summary{cursor:pointer;font-weight:600;color:var(--bs-body-color);list-style:none;}
      details>summary::-webkit-details-marker{display:none;}
      details[open]{background:var(--mfap-details-bg-open);}
      .upload-card-body{padding-top:0.35rem;padding-bottom:0.35rem;}
      .upload-card-body > *{margin:0;margin-top:0.08rem;}
      .upload-card-body > *:first-child{margin-top:0;}
      .upload-card-body .form-group,.upload-card-body .shiny-input-container{margin-bottom:0 !important;}
      .upload-fields{display:flex;flex-direction:column;gap:0.05rem;}
      .upload-fields .form-group,.upload-fields .shiny-input-container{margin-bottom:0 !important;}
      .help-trigger{font-size:0.82rem;color:var(--mfap-help-trigger);text-decoration:none;display:inline-flex;align-items:center;gap:6px;}
      .help-trigger:hover{color:var(--mfap-help-trigger-hover);text-decoration:none;}
      .help-trigger-icon{opacity:0.85;font-size:0.95em;}
      .disclaimer-bar{background:var(--mfap-disclaimer-bg);border-left:3px solid var(--mfap-disclaimer-border);
        border-radius:6px;padding:10px 14px;margin-bottom:16px;font-size:0.88rem;
        display:flex;flex-wrap:wrap;align-items:center;gap:6px;}
      .disclaimer-one-line{margin:0;color:var(--mfap-disclaimer-text);}
      .disclaimer-link{color:var(--mfap-help-trigger);text-decoration:none;}
      .disclaimer-link:hover{color:var(--mfap-help-trigger-hover);text-decoration:underline;}
      .disclaimer-sep{color:var(--mfap-disclaimer-sep);}
      .disclaimer-cite{color:var(--mfap-disclaimer-text);}
      .disclaimer-bib{margin-left:2px;}
      .muted{color:var(--mfap-muted);}
      .download-row{display:flex;justify-content:flex-end;margin:0;}

      .about-mfap .card-body{padding:12px 14px !important;}
      .about-mfap .lead-line{font-size:1.15rem;color:var(--mfap-about-lead);margin-bottom:10px;}
      .about-mfap .notes{background:var(--mfap-about-notes-bg);border:1px solid var(--mfap-about-notes-border);
                         border-radius:10px;padding:12px 14px;margin-top:8px;}
      .about-mfap .notes p{margin:0 0 6px 0;font-size:.98rem;line-height:1.35;}
      .about-mfap .cite-row{display:flex;align-items:center;gap:.5rem;margin-top:10px;}
      .about-mfap .bib-link{cursor:pointer;text-decoration:underline;}

      
      #bibtex_block{white-space:pre-wrap;font-size:.86rem;background:var(--mfap-bib-bg);
      padding:10px;border-radius:8px;border:1px solid var(--mfap-bib-border);}
      .bib-block{
        white-space:pre-wrap;
        font-size:.86rem;
        background:var(--mfap-bib-bg);
        color:var(--bs-body-color);
        padding:10px;
        border-radius:8px;
        border:1px solid var(--mfap-bib-border);
        margin:0;
      }
      .intro-example-shell{
        border:1px solid var(--mfap-intro-border);
        border-radius:10px;
        padding:10px 12px;
        margin:12px 0;
        background:var(--mfap-intro-bg);
      }
      .intro-example-table{
        margin:0;
        width:auto;
        display:inline-table;
        font-size:0.85rem;
        line-height:1.1;
      }
      .intro-example-table th,
      .intro-example-table td{
        padding:2px 6px;
      }
      .modal-subheading{
        margin-top:0;
        color:var(--bs-body-color);
      }
      
      .sidebar.well,
      .well.sidebar {
        background: transparent !important;
        border: 0 !important;
        box-shadow: none !important;
        padding: 0 !important;
      }
      
      .site-footer__inner{
        max-width:1180px; margin:0 auto; padding:20px 14px;
        display:flex; align-items:center; justify-content:space-between; gap:14px;
      }
      .site-footer__container{
        max-width:1180px;
        margin:0 auto;
        padding:16px 14px;
        display:flex;
        align-items:center;
        justify-content:space-between;
        gap:16px;
      }
      .site-footer__logo{ height:60px; flex:0 0 auto; }
      html[data-mfap-theme='light'] .site-footer__logo{
        filter: invert(1);
      }
      .site-footer__center{ flex:1 1 auto; text-align:center; color:var(--mfap-muted); }
      
      .site-footer__link:hover,
      .site-footer__link:focus{
        color:var(--mfap-link);      
        text-decoration:none;
      }
      .site-footer{
  position: relative;
  z-index: 10;
  color: var(--mfap-footer-color) !important;
}
      
button.intro-cta{
  width:100%;
  font-weight:700;
  padding:12px 14px;
  border-radius:10px;

  background: var(--mfap-intro-bg);
  border:1px solid var(--mfap-intro-border);

  color:var(--mfap-intro-text) !important;
  text-align:center;

  position: relative;
  overflow: hidden;

  transition: transform .08s ease, box-shadow .12s ease, border-color .12s ease;
}

button.intro-cta::before{
  content:'';
  position:absolute;
  inset:0;
  background: rgba(77,171,247,0.85);
  transform: translateX(-101%);
  transition: transform .28s ease;
  z-index:0;
}

button.intro-cta > span{
  position:relative;
  z-index:1;
}

button.intro-cta:hover,
button.intro-cta:focus,
button.intro-cta:active{
  color:var(--mfap-intro-hover-text) !important;
  border-color: rgba(77,171,247,0.65);
  box-shadow: 0 10px 26px rgba(0,0,0,.35);
  background: var(--mfap-intro-bg) !important;
}

button.intro-cta:hover::before{
  transform: translateX(0%);
}

button.intro-cta:active{
  transform: translateY(1px);
}

      html[data-mfap-theme='light'] #shiny-modal [style*='color:#dee2e6']{
        color:#374151 !important;
      }
      html[data-mfap-theme='light'] [style*='color:#ccc']{
        color:#64748b !important;
      }
      .theme-toggle-link{
        cursor:pointer;
        user-select:none;
        font-size:1.15rem;
        line-height:1;
        min-width:auto;
        text-align:center;
        border-radius:0 !important;
        border-bottom-color:transparent !important;
        padding:.45rem .25rem !important;
        margin-left:0 !important;
        background:transparent !important;
        box-shadow:none !important;
      }
      .theme-toggle-link:hover,
      .theme-toggle-link:focus{
        background:transparent !important;
        border-bottom-color:transparent !important;
        box-shadow:none !important;
      }
      html[data-mfap-theme='light'] .form-check-input{
        background-color:#f8fafc;
        border-color:#64748b;
      }
      html[data-mfap-theme='light'] .form-check-input:focus{
        border-color:#475569;
        box-shadow:0 0 0 .2rem rgba(37,99,235,.18);
      }
      html[data-mfap-theme='light'] .form-check-input:checked{
        background-color:#2563eb;
        border-color:#1d4ed8;
      }
      .shiny-input-container .checkbox label{
        display:flex;
        align-items:flex-start;
        gap:.55rem;
      }
      .shiny-input-container .checkbox input[type='checkbox']{
        appearance:auto !important;
        -webkit-appearance:checkbox !important;
        opacity:1 !important;
        position:relative !important;
        margin:0;
        margin-top:.22rem;
        width:1rem;
        height:1rem;
        flex:0 0 auto;
      }
      html[data-mfap-theme='light'] .shiny-input-container .checkbox input[type='checkbox']{
        accent-color:#2563eb;
      }
      .btn-info,
      .shiny-input-container .input-group-btn .btn,
      .shiny-input-container .btn-file,
      .shiny-input-container label.btn-file{
        background-color:var(--mfap-btn-info-bg) !important;
        border-color:var(--mfap-btn-info-bg) !important;
        color:var(--mfap-btn-info-color) !important;
      }
      .btn-info:hover,.btn-info:focus,
      .shiny-input-container .input-group-btn .btn:hover,
      .shiny-input-container .input-group-btn .btn:focus,
      .shiny-input-container .btn-file:hover,
      .shiny-input-container .btn-file:focus{
        background-color:var(--mfap-btn-info-hover) !important;
        border-color:var(--mfap-btn-info-hover) !important;
        color:var(--mfap-btn-info-color) !important;
      }
      .btn-info:active,
      .shiny-input-container .input-group-btn .btn:active,
      .shiny-input-container .btn-file:active{
        background-color:var(--mfap-btn-info-active) !important;
        border-color:var(--mfap-btn-info-active) !important;
        color:var(--mfap-btn-info-color) !important;
      }
      .shiny-notification.shiny-notification-message{
        background:var(--mfap-toast-message-bg) !important;
        border:1px solid var(--mfap-toast-message-border) !important;
        color:var(--mfap-toast-message-text) !important;
      }
      .shiny-notification.shiny-notification-message .shiny-notification-content,
      .shiny-notification.shiny-notification-message .shiny-notification-close{
        color:var(--mfap-toast-message-text) !important;
      }
      .shiny-notification.shiny-progress-notification,
      .shiny-progress-container .shiny-notification,
      .shiny-progress-container .shiny-notification.shiny-notification-closeable{
        background-color:var(--mfap-toast-message-bg) !important;
        background-image:none !important;
        border:1px solid var(--mfap-toast-message-border) !important;
        color:var(--mfap-toast-message-text) !important;
        opacity:1 !important;
        backdrop-filter:none !important;
      }
      .shiny-notification.shiny-progress-notification .shiny-notification-content,
      .shiny-notification.shiny-progress-notification .shiny-notification-close{
        color:var(--mfap-toast-message-text) !important;
        opacity:1 !important;
      }
      .shiny-progress-container .progress-text .progress-message,
      .shiny-progress-container .progress-text .progress-detail,
      .shiny-progress-container .shiny-notification-close{
        color:var(--mfap-toast-message-text) !important;
      }
      .shiny-notification.shiny-progress-notification .progress,
      .shiny-progress-container .progress{
        background:rgba(255,255,255,0.25) !important;
      }
      .shiny-notification.shiny-progress-notification .progress .progress-bar,
      .shiny-progress-container .progress .progress-bar{
        background:var(--mfap-link-hover) !important;
      }
      html[data-mfap-theme='light'] #clear_uploads,
      html[data-mfap-theme='light'] #clr_dom,
      html[data-mfap-theme='light'] #clr_mot,
      html[data-mfap-theme='light'] #clear_example_data{
        background-color:#94a3b8 !important;
        border-color:#94a3b8 !important;
        color:#ffffff !important;
      }
      html[data-mfap-theme='light'] #clear_uploads:hover,
      html[data-mfap-theme='light'] #clear_uploads:focus,
      html[data-mfap-theme='light'] #clr_dom:hover,
      html[data-mfap-theme='light'] #clr_dom:focus,
      html[data-mfap-theme='light'] #clr_mot:hover,
      html[data-mfap-theme='light'] #clr_mot:focus,
      html[data-mfap-theme='light'] #clear_example_data:hover,
      html[data-mfap-theme='light'] #clear_example_data:focus{
        background-color:#64748b !important;
        border-color:#64748b !important;
        color:#ffffff !important;
      }
      html[data-mfap-theme='light'] #clear_uploads:active,
      html[data-mfap-theme='light'] #clr_dom:active,
      html[data-mfap-theme='light'] #clr_mot:active,
      html[data-mfap-theme='light'] #clear_example_data:active{
        background-color:#475569 !important;
        border-color:#475569 !important;
        color:#ffffff !important;
      }



      @media (max-width:768px){
        .site-footer__logo{ height:42px; }
        .site-footer__center{ font-size:0.88rem; }
      }
      .compact-btn-row{
        display:flex;
        flex-wrap:nowrap;
        justify-content:flex-start;
        gap:0;
        padding-left:0;
        padding-right:0;
      }
      .compact-btn-row > [class*='col-']{
        width:auto !important;
        max-width:none !important;
        flex:0 0 auto !important;
        padding-left:0;
        padding-right:0;
      }
      .compact-btn-row > [class*='col-']:first-child{
        padding-left:15px;
        padding-right:5px;
      }
      .compact-btn-row > [class*='col-']:last-child{
        padding-left:5px;
        padding-right:15px;
      }
    ")),
      tags$script(HTML("
      (function(){
        var storageKey = 'mfap_theme_preference_v1';
        var systemThemeQuery = null;

        function setToggleLabel(mode){
          var el = document.getElementById('theme_toggle');
          if (!el) return;
          if (el.parentElement) {
            el.parentElement.style.marginLeft =
              (window.matchMedia && window.matchMedia('(max-width: 768px)').matches) ? '0' : 'auto';
          }
          if (mode === 'light') {
            el.textContent = '\u263C';
            el.setAttribute('title', 'Switch to dark mode');
            el.setAttribute('aria-label', 'Switch to dark mode');
          } else {
            el.textContent = '\u263E';
            el.setAttribute('title', 'Switch to light mode');
            el.setAttribute('aria-label', 'Switch to light mode');
          }
        }

        function getSystemPreferredTheme(){
          return (window.matchMedia && window.matchMedia('(prefers-color-scheme: light)').matches) ? 'light' : 'dark';
        }

        function applyTheme(mode, persistPreference){
          var resolved = mode === 'light' ? 'light' : 'dark';
          document.documentElement.setAttribute('data-mfap-theme', resolved);
          document.documentElement.setAttribute('data-bs-theme', resolved);
          setToggleLabel(resolved);
          if (persistPreference) {
            try { localStorage.setItem(storageKey, resolved); } catch(e) {}
          }
        }

        window.mfapApplyTheme = function(mode){
          applyTheme(mode, true);
        };
        window.mfapToggleTheme = function(){
          var current = document.documentElement.getAttribute('data-mfap-theme') || 'dark';
          applyTheme(current === 'light' ? 'dark' : 'light', true);
        };
        document.addEventListener('DOMContentLoaded', function(){
          var saved = null;
          try { saved = localStorage.getItem(storageKey); } catch(e) {}
          if (saved === 'light' || saved === 'dark') {
            applyTheme(saved, false);
          } else {
            applyTheme(getSystemPreferredTheme(), false);
          }

          if (window.matchMedia) {
            systemThemeQuery = window.matchMedia('(prefers-color-scheme: light)');
            var onSystemThemeChange = function(e){
              var explicit = null;
              try { explicit = localStorage.getItem(storageKey); } catch(err) {}
              if (explicit === 'light' || explicit === 'dark') return;
              applyTheme(e.matches ? 'light' : 'dark', false);
            };
            if (systemThemeQuery.addEventListener) {
              systemThemeQuery.addEventListener('change', onSystemThemeChange);
            } else if (systemThemeQuery.addListener) {
              systemThemeQuery.addListener(onSystemThemeChange);
            }
          }
        });

        document.addEventListener('click', function(e){
          var btn = e.target.closest('#theme_toggle');
          if (!btn) return;
          e.preventDefault();
          window.mfapToggleTheme();
        });
      })();
    ")),
      tags$script(HTML("
  Shiny.addCustomMessageHandler('mfap_intro_once', function(msg){
    const key = (msg && msg.key) ? msg.key : 'mfap_intro_dismissed_v1';
    let dismissed = false;
    try { dismissed = sessionStorage.getItem(key) === '1'; } catch(e) {}
    if (dismissed) return;

    // trigger server to show the modal
    Shiny.setInputValue('intro_show', Date.now(), {priority: 'event'});
  });

  Shiny.addCustomMessageHandler('mfap_intro_store', function(msg){
    const key = (msg && msg.key) ? msg.key : 'mfap_intro_dismissed_v1';
    try { sessionStorage.setItem(key, '1'); } catch(e) {}
  });
")),
      
      tags$script(HTML("
      document.addEventListener('click', function(e){
        const btn = e.target.closest('[data-copy-target]');
        if (!btn) return;
        const sel = btn.getAttribute('data-copy-target');
        const pre = document.querySelector(sel);
        if (!pre) return;
        navigator.clipboard.writeText(pre.innerText).then(function(){
          const old = btn.textContent;
          btn.textContent = 'Copied';
          setTimeout(function(){ btn.textContent = old; }, 1200);
        });
      });
    "))
    ,
    tags$script(HTML("
      (function(){
        function cssVar(name, fallback){
          var v = getComputedStyle(document.documentElement).getPropertyValue(name);
          return (v && v.trim()) ? v.trim() : fallback;
        }
        function styleProgressToast(node){
          if (!node || !node.querySelector) return;
          var isProgress = node.classList.contains('shiny-progress-notification') ||
            !!node.querySelector('.progress') ||
            !!node.querySelector('.progress-message') ||
            /Running analysis/i.test(node.textContent || '');
          if (!isProgress) return;

          var bg = cssVar('--mfap-toast-message-bg', '#2563eb');
          var border = cssVar('--mfap-toast-message-border', '#1d4ed8');
          var text = cssVar('--mfap-toast-message-text', '#ffffff');
          var bar = cssVar('--mfap-link-hover', '#74c0fc');

          node.style.setProperty('background-color', bg, 'important');
          node.style.setProperty('background-image', 'none', 'important');
          node.style.setProperty('border-color', border, 'important');
          node.style.setProperty('color', text, 'important');
          node.style.setProperty('opacity', '1', 'important');
          node.style.setProperty('backdrop-filter', 'none', 'important');

          node.querySelectorAll('.shiny-notification-content, .shiny-notification-close, .progress-message, .progress-detail').forEach(function(el){
            el.style.setProperty('color', text, 'important');
            el.style.setProperty('opacity', '1', 'important');
          });
          node.querySelectorAll('.progress').forEach(function(el){
            el.style.setProperty('background', 'rgba(255,255,255,0.25)', 'important');
          });
          node.querySelectorAll('.progress-bar').forEach(function(el){
            el.style.setProperty('background', bar, 'important');
          });
        }

        function scan(){
          document.querySelectorAll('.shiny-notification, .shiny-progress-container .shiny-notification').forEach(styleProgressToast);
        }

        document.addEventListener('DOMContentLoaded', function(){
          scan();
          var observer = new MutationObserver(function(mutations){
            mutations.forEach(function(m){
              m.addedNodes.forEach(function(n){
                if (!n || n.nodeType !== 1) return;
                if (n.matches && (n.matches('.shiny-notification') || n.matches('.shiny-progress-container'))) {
                  styleProgressToast(n);
                }
                if (n.querySelectorAll) {
                  n.querySelectorAll('.shiny-notification, .shiny-progress-container .shiny-notification').forEach(styleProgressToast);
                }
              });
            });
            scan();
          });
          observer.observe(document.body, { childList: true, subtree: true });
        });
      })();
    "))
    ),
    
    # Analysis tab
    bslib::nav_panel(
      "Analysis",
      div(class="page-container",
          sidebarLayout(
            sidebarPanel(
              width = 3, class="sidebar",
              style = "max-height: calc(100vh - 220px); overflow-y: auto; overflow-x: hidden;",
              bslib::card(
                class = "upload-card",
                bslib::card_header("Upload"),
                bslib::card_body(
                  class = "upload-card-body",
                  tags$div(class = "help-trigger-wrap",
                    actionLink("upload_help_trigger",
                      tagList(icon("info-circle", class = "help-trigger-icon"), tags$span("Table format and FASTA download")),
                      class = "help-trigger")
                  ),
                  tags$div(class = "upload-fields",
                    uiOutput("csv_ui"),
                    uiOutput("csv_status"),
                    uiOutput("fasta_ui"),
                    uiOutput("fasta_status")
                  ),
                  checkboxInput("use_titer", "Include alternative translation initiation site prediction (TITER, developed by Zhang et al. 2017)", FALSE),
                  conditionalPanel(
                    condition = "input.use_titer == true",
                    uiOutput("fasta_flank_ui"),
                    uiOutput("fasta_flank_status")
                  ),
                  div(style="height:2px;"),
                  actionButton("clear_uploads", "Clear uploads", class="btn btn-secondary")
                )
              ),
              bslib::card(
                bslib::card_header("Domains"),
                bslib::card_body(
                  fluidRow(
                    column(4, textInput("dom_name","Name")),
                    column(4, numericInput("dom_start","Start AA", 1, min = 1)),
                    column(4, numericInput("dom_end","End AA", 1, min = 1))
                  ),
                  fluidRow(
                    class = "compact-btn-row",
                    column(6, class = "col-6", actionButton("add_dom","Add", class="btn btn-info")),
                    column(6, class = "col-6", actionButton("clr_dom","Clear", class="btn btn-secondary"))
                  ),
                  uiOutput("dom_tbl_wrapper")
                )
              ),
              
              bslib::card(
                bslib::card_header("Motifs"),
                bslib::card_body(
                  fluidRow(column(6, textInput("mot_name","Name")),
                           column(6, textInput("mot_pattern","AA sequence"))),
                  fluidRow(class = "compact-btn-row",
                           column(6, class = "col-6", actionButton("add_mot","Add", class="btn btn-info")),
                           column(6, class = "col-6", actionButton("clr_mot","Clear", class="btn btn-secondary"))),
                  uiOutput("mot_tbl_wrapper")
                )
              ),
              bslib::card(
                bslib::card_header("Demo"),
                bslib::card_body(
                  tags$p(class="muted", style="font-size:0.9em;",
                         "Load a fictional cohort with example variants, reference sequences, and example domains/motifs."),
                  fluidRow(
                    class = "compact-btn-row",
                    column(6, class = "col-6", actionButton("load_example_data", "Load demo", class="btn btn-info")),
                    column(6, class = "col-6", actionButton("clear_example_data", "Clear demo", class="btn btn-secondary"))
                  )
                )
              ),
              div(class="action-dock",
                  div(class="action-row",
                      actionButton("run", "Analyze", icon=icon("play"), class="btn btn-mfap-primary")
                  )
              )
            ),
            mainPanel(
              width = 9,
              style = "max-height: calc(100vh - 220px); overflow-y: auto; overflow-x: hidden;",
              bslib::card(
                bslib::card_header("Important notes"),
                bslib::card_body(
                  tags$p(
                    class = "muted",
                    style = "font-size:0.98rem; line-height:1.4; margin-bottom:10px;",
                    "This application is provided for research and educational purposes only. ",
                    "It is not intended for clinical use, medical decision-making, diagnosis, or patient management. ",
                    "No warranty is given as to the correctness, completeness, or suitability of the results. ",
                    "The authors and operators assume no liability for any use of the information generated by this application. ",
                    "Users are solely responsible for compliance with applicable data-protection, privacy, and ethical regulations."
                  ),
                  tags$div(
                    class = "muted",
                    tags$span("Please cite: Schubert et al. (2025) ·"),
                    actionLink("show_bib", label = "BibTeX", class = "bib-link")
                  ),
                  tags$div(
                    class = "muted",
                    tags$span("If using TITER, please additionally cite: Zhang et al. (2017)"),
                    " · ",
                    actionLink("show_titer_bibtex", label = "BibTeX", class = "bib-link")
                  )
                )
              ),
              
              bslib::card(
                bslib::card_header(uiOutput("results_header")),
                bslib::card_body(uiOutput("results_body"))
              )
              
              
            )
          )
      )
    ),
    
    # Background tab
    bslib::nav_panel(
      "Background",
      div(class="page-container",
          fluidRow(
            column(12,
                   bslib::card(
                     bslib::card_header("Overview"),
                     bslib::card_body(
                      p(mfap()," is a framework, which allows users to calculate and predict DNA and protein level features resulting from genetic variation in single-exon (i.e., intronless) genes. ",mfap()," requires minimal user input: a list of variants in HGVS cDNA notation and the reference DNA sequence of the gene of interest.")
                     )
                   ),
                   bslib::card(
                     bslib::card_header("Rationale"),
                     bslib::card_body(
                       p("Protein-truncating variants (PTVs) canonically lead to nonsense-mediated mRNA decay. This process necessitates an exon-exon junction, which is lacking in single-exon genes. Therefore, PTVs in single-exon genes can result in a truncated protein that may have dominant-negative, gain-of-function, or neomorphic effects. In our accompanying manuscript, we demonstrate that ",mfap(),"-predicted protein-level attributes of such truncated proteins provide a stronger causal link between genotype and phenotypic severity."),
                       tags$figure(
                         tags$img(
                           src = "https://lh3.googleusercontent.com/d/1gDNYo7HuDE4md9mvPLc4yiQXgOnl29Ly",
                           alt = "Rationale",
                           style = "display:block;margin:10px auto;max-width:60%;height:auto;border-radius:10px;border:1px solid rgba(255,255,255,.08);",
                           loading = "lazy"
                         ),
                         tags$figcaption(
                           "Figure: Evasion of nonsense-mediated decay hypothesis. Figure created with Biorender.com",
                           style="text-align:center;color:#cfd4da;font-size:.9em;margin-top:6px;"
                         )
                       )
                     )
                   ),
                   bslib::card(
                     bslib::card_header("Workflow"),
                     bslib::card_body(
                       p("",mfap()," outputs calculations and predictions of the effects of variants on the DNA and protein level, including predictions of non-canonical translation initiation sites (TISs) and resulting proteins. For more information and inspiration for downstream analyses, please read the accompanying paper."),
                       tags$figure(
                         tags$img(
                           src = "https://lh3.googleusercontent.com/d/1h-h2yfvsKnxscA3ko8VyeFXkPaa9iAkE",
                           alt = "MfAP workflow schematic",
                           style = "display:block;margin:10px auto;max-width:60%;height:auto;border-radius:10px;border:1px solid rgba(255,255,255,.08);",
                           loading = "lazy"
                         ),
                         tags$figcaption(
                           "Figure: ",mfap()," workflow. Figure created with Biorender.com",
                           style="text-align:center;color:#cfd4da;font-size:.9em;margin-top:6px;"
                         )
                       )
                     )
                   ),
                   bslib::card(
                     bslib::card_header("References"),
                     bslib::card_body(
                       tags$p(
                         style = "margin:0 0 8px 0;",
                         HTML("Schubert T, Tietzel A, Pottayil H, Caro P, Gilmore RB, Franke F, Althammer F, Schaaf CP. 2025. "),
                         em("A blueprint for protein-centric genotype-phenotype investigations in single-exon disease genes applied to MAGEL2 and Schaaf-Yang syndrome."),
                         HTML(" <i>Journal</i> <b>VOLUME</b>(ISSUE):PAGES · "),
                         tags$a(href="https://doi.org/10.xxxx/xxxxx", target="_blank", "doi:10.xxxx/xxxxx"),
                         HTML(" · "),
                         actionLink("bib_schubert_bg", "BibTeX", class = "bib-link")
                       ),
                       tags$p(
                         style = "margin:0;",
                         HTML("Zhang S, Hu H, Jiang T, Zhang L, Zeng J. 2017. "),
                         em("TITER: predicting translation initiation sites by deep learning."),
                         HTML(" <i>Bioinformatics</i> <b>33</b>(14):i234–i242 · "),
                         tags$a(href="https://doi.org/10.1093/bioinformatics/btx247", target="_blank", "doi:10.1093/bioinformatics/btx247"),
                         HTML(" · "),
                         actionLink("bib_titer_bg", "BibTeX", class = "bib-link")
                       )
                     )
                   ),
                  bslib::card(
                    bslib::card_header("Contributors"),
                     bslib::card_body(
                       tags$div(
                         style = "font-size:0.9em; line-height:1.4; margin-bottom:10px;",
                        tags$div("Tim Schubert", tags$sup("1")),
                         tags$div("Antonia Tietzel", tags$sup("1,2")),
                         tags$div("Hari Pottayil", tags$sup("1")),
                         tags$div("Pilar Caro", tags$sup("1")),
                         tags$div("Rachel B. Gilmore", tags$sup("1")),
                         tags$div("Felix Franke", tags$sup("1")),
                         tags$div("Ferdinand Althammer", tags$sup("1")),
                         tags$div("Christian P. Schaaf", tags$sup("1"))
                       ),
                       tags$div(
                         style = "font-size:0.88em; line-height:1.35;",
                         tags$div(tags$b("Affiliations")),
                         tags$div(tags$sup("1"), " Institute of Human Genetics, Heidelberg University, Heidelberg, Germany"),
                         tags$div(tags$sup("2"), " Clinical Cooperation Unit Neuropathology, German Cancer Research Center (DKFZ), Heidelberg, Germany")
                       )
                     )
                   ),
                   bslib::card(
                     bslib::card_header("Software"),
                     bslib::card_body(
                       tags$div(
                         style = "font-size:0.9em; line-height:1.4;",
                         tags$p(
                           style = "margin-bottom:10px;",
                           "R version 4.4.2"
                         ),
                         tags$p(
                           style = "margin-bottom:0;",
                           "Shiny package (version 1.11.1)"
                         )
                       )
                     )
                   ),
                   bslib::card(
                     bslib::card_header("IT Infrastructure"),
                     bslib::card_body(
                       tags$div(
                         style = "font-size:0.9em; line-height:1.4;",
                         tags$p(
                           style = "margin-bottom:10px;",
                           "de.NBI Cloud"
                         ),
                         tags$p(
                           style = "margin-bottom:0;",
                           "University Computing Centre Heidelberg"
                         )
                       )
                     )
                   )
            )
          )
      )
    ),  # end nav_panel Background
    bslib::nav_item(
      actionLink("open_intro", "Introduction", class = "nav-link")
    ),
    bslib::nav_item(
      tags$a(
        id = "theme_toggle",
        href = "#",
        class = "nav-link theme-toggle-link",
        "☾"
      )
    )
  )
  
  ,
  
  tags$footer(
    class = "site-footer",
    tags$div(
      class = "site-footer__inner",
      
      tags$img(
        src = "https://lh3.googleusercontent.com/d/12vCwvlmN8xrKnlU1GsDGJeDBgq_yl0sM",
        class = "site-footer__logo"
      ),
      
      tags$div(
        class = "site-footer__center",
        tags$p(
          class = "footer-links",
          style = "margin:0;",
          
          tags$span(class = "footer-item", HTML("&copy; 2026 "), tags$a(href = "https://tim-schubert.github.io", target = "_blank", rel = "noopener noreferrer", class = "footer-link", "Tim Schubert")),
          
          tags$span(class = "footer-sep", HTML("&middot;")),
          
          tags$a(
            href   = "https://www.apache.org/licenses/LICENSE-2.0.html",
            target = "_blank",
            "Apache License 2.0",
            class  = "footer-item footer-link"
          ),
          
          tags$span(class = "footer-sep", HTML("&middot;")),
          
          tags$a(
            href   = "#",
            "Legal notice (Impressum)",
            class  = "footer-item footer-link",
            onclick = "Shiny.setInputValue('open_legal', Date.now(), {priority:'event'}); return false;"
          )
        )
      ),
      
      tags$img(
        src = "https://lh3.googleusercontent.com/d/1BPl641wAs67xCTAKEEETVO0zjzOYUrda",
        class = "site-footer__logo"
      )
    )
  )
  
)  # end tagList




# Server
server <- function(input, output, session) {
  example_csv <- reactiveVal(NULL)
  example_fasta_seq <- reactiveVal(NULL)
  example_fasta_flank_seq <- reactiveVal(NULL)
  example_domains <- data.frame(Name=c("DomainA","DomainB"), Start=c(5,50), End=c(20,100), stringsAsFactors=FALSE)
  example_motifs  <- data.frame(Name=c("MotifX","MotifY"),  Pattern=c("ACDE","WXYZ"), stringsAsFactors=FALSE)
  csv_df_norm    <- reactiveVal(NULL)
  csv_ok         <- reactiveVal(FALSE)
  csv_msg        <- reactiveVal("")
  csv_renames    <- reactiveVal(character())
  
  fasta_ok       <- reactiveVal(FALSE)
  fasta_msg      <- reactiveVal("")
  flank_ok       <- reactiveVal(FALSE)
  flank_msg      <- reactiveVal("")
  
  canon_key <- function(x) tolower(gsub("[^a-z0-9]+", "", trimws(x)))
  
  normalize_columns <- function(df) {
    stopifnot(is.data.frame(df))
    cn <- colnames(df)
    ck <- canon_key(cn)
    
    aliases <- list(
      patient_ID  = c("patientid","patient_id","patient","sampleid","sample_id","individualid","individual_id","subjectid","subject_id","id"),
      DNA_variant = c("dnavariant","dna_variant","variant","variants","cdnvariant","cdna_variant","hgvs","hgvsc","hgvs_cdna","cvariant","c_variant","dnachange","dna_change")
    )
    
    find_best <- function(target) {
      hits <- which(ck %in% aliases[[target]] | ck == canon_key(target))
      if (length(hits) == 0) return(NA_integer_)
      if (length(hits) > 1) {
        stop(sprintf("Ambiguous column mapping for %s. Candidates: %s", target, paste(cn[hits], collapse = ", ")))
      }
      hits[1]
    }
    
    renames <- c()
    i_pat <- find_best("patient_ID")
    if (!is.na(i_pat) && cn[i_pat] != "patient_ID") {
      renames[cn[i_pat]] <- "patient_ID"
      cn[i_pat] <- "patient_ID"
    }
    i_var <- find_best("DNA_variant")
    if (!is.na(i_var) && cn[i_var] != "DNA_variant") {
      renames[cn[i_var]] <- "DNA_variant"
      cn[i_var] <- "DNA_variant"
    }
    
    colnames(df) <- cn
    list(df = df, renames = renames)
  }
  
  variant_issue <- function(x) {
    if (is.na(x) || !nzchar(trimws(x))) return("Missing/empty variant")
    x0 <- trimws(x)
    if (grepl("\\s", x0)) return("Contains whitespace (remove spaces)")
    if (grepl("[\\+\\-\\*\\(\\)\\?]", x0)) return("Contains intronic/UTR/uncertain HGVS characters (+, -, *, ?, parentheses)")
    x <- gsub("\\s+", "", x0)
    x2 <- sub("^(c\\.)", "", x, ignore.case = TRUE)
    
    if (grepl("^\\d+[ACGTNacgtn]>[ACGTNacgtn]$", x2)) return("")
    if (grepl("^\\d+(?:_\\d+)?del[ACGTNacgtn]*$", x2)) return("")
    if (grepl("^\\d+(?:_\\d+)?dup[ACGTNacgtn]*$", x2)) return("")
    if (grepl("^\\d+(?:_\\d+)?ins[ACGTNacgtn]+$", x2)) return("")
    if (grepl("^\\d+(?:_\\d+)?delins[ACGTNacgtn]+$", x2)) return("")
    
    
    "Unrecognized HGVS cDNA pattern."
  }
  
  parse_hgvs_minimal <- function(x) {
    x <- gsub("\\s+", "", trimws(x))
    x2 <- sub("^(c\\.)", "", x, ignore.case = TRUE)
    
    if (grepl(">", x2)) {
      m <- regexec("^(\\d+)([ACGTNacgtn])>([ACGTNacgtn])$", x2)
      r <- regmatches(x2, m)[[1]]
      if (length(r) == 4) return(list(type="snv", start=as.integer(r[2]), end=as.integer(r[2]), ref=toupper(r[3]), alt=toupper(r[4])))
      return(NULL)
    }
    
    m <- regexec("^(\\d+)(?:_(\\d+))?(delins)([ACGTNacgtn]+)$", x2)
    r <- regmatches(x2, m)[[1]]
    if (length(r) == 5) {
      s <- as.integer(r[2]); e <- ifelse(nzchar(r[3]), as.integer(r[3]), s)
      return(list(type="delins", start=s, end=e, ins=toupper(r[5])))
    }
    
    m <- regexec("^(\\d+)(?:_(\\d+))?(del)([ACGTNacgtn]*)$", x2)
    r <- regmatches(x2, m)[[1]]
    if (length(r) == 5) {
      s <- as.integer(r[2]); e <- ifelse(nzchar(r[3]), as.integer(r[3]), s)
      return(list(type="del", start=s, end=e))
    }
    
    m <- regexec("^(\\d+)(?:_(\\d+))?(dup)([ACGTNacgtn]*)$", x2)
    r <- regmatches(x2, m)[[1]]
    if (length(r) == 5) {
      s <- as.integer(r[2]); e <- ifelse(nzchar(r[3]), as.integer(r[3]), s)
      return(list(type="dup", start=s, end=e))
    }
    
    m <- regexec("^(\\d+)(?:_(\\d+))?(ins)([ACGTNacgtn]+)$", x2)
    r <- regmatches(x2, m)[[1]]
    if (length(r) == 5) {
      s <- as.integer(r[2]); e <- ifelse(nzchar(r[3]), as.integer(r[3]), s)
      return(list(type="ins", start=s, end=e, ins=toupper(r[5])))
    }
    
    NULL
  }
  
  validate_variants_against_ref <- function(v, ref_seq) {
    n <- nchar(ref_seq)
    probs <- character(length(v))
    
    for (i in seq_along(v)) {
      hg <- v[i]
      p <- parse_hgvs_minimal(hg)
      if (is.null(p)) {
        probs[i] <- paste0("Could not interpret variant for semantic checks: ", hg)
        next
      }
      
      if (is.na(p$start) || is.na(p$end)) {
        probs[i] <- paste0("Invalid coordinates: ", hg)
        next
      }
      
      if (p$start < 1 || p$end < 1) {
        probs[i] <- paste0("Position < 1: ", hg)
        next
      }
      
      if (p$start > p$end) {
        probs[i] <- paste0("Invalid range (start > end): ", hg)
        next
      }
      
      if (p$end > n) {
        probs[i] <- paste0("Position out of CDS bounds (", p$end, " > ", n, "): ", hg)
        next
      }
      
      if (identical(p$type, "snv")) {
        ref_base <- substr(ref_seq, p$start, p$start)
        if (toupper(ref_base) != toupper(p$ref)) {
          probs[i] <- paste0("Reference-base mismatch at c.", p$start, " (FASTA=", toupper(ref_base), ", variant=", toupper(p$ref), "): ", hg)
          next
        }
      }
      
      probs[i] <- ""
    }
    
    probs
  }
  
  parse_mutation_strict <- function(mutation, wildtype_seq) {
    info <- parse_mutation(mutation, wildtype_seq)
    ok <- !is.na(info$Locus) && !is.na(info$Mutation_Type) && !is.na(info$Mutated_Sequence) && nzchar(info$Mutated_Sequence)
    if (!ok) stop(paste0("Could not parse variant: ", mutation))
    info
  }
  
  translate_prot_strict <- function(dna_seq, variant_label = NA_character_) {
    if (is.na(dna_seq) || !nzchar(dna_seq) || nchar(dna_seq) < 3) {
      stop(paste0("Protein translation failed (empty/too short sequence) for variant: ", variant_label))
    }
    aa <- suppressWarnings(as.character(translate(DNAString(dna_seq), if.fuzzy.codon="X")))
    if (is.na(aa) || !nzchar(aa)) stop(paste0("Protein translation failed for variant: ", variant_label))
    if (grepl("\\*", aa)) substr(aa, 1, regexpr("\\*", aa)[1]-1) else aa
  }
  
  
  
  domains <- reactiveVal(data.frame(Name=character(),Start=integer(),End=numeric(),stringsAsFactors=FALSE))
  motifs  <- reactiveVal(data.frame(Name=character(),Pattern=character(),stringsAsFactors=FALSE))
  
  wt_aa_len <- reactiveVal(NULL)
  
  result_data <- reactiveVal(NULL)
  
  observe({
    has_csv   <- !is.null(input$csv)   || !is.null(example_csv())
    has_fasta <- !is.null(input$fasta) || !is.null(example_fasta_seq())
    has_flank_file <- if (isTRUE(input$use_titer)) {!is.null(input$fasta_flank) || !is.null(example_fasta_flank_seq())} else TRUE
    
    csv_ready   <- if (!is.null(example_csv())) TRUE else (has_csv && isTRUE(csv_ok()))
    fasta_ready <- if (!is.null(example_fasta_seq())) TRUE else (has_fasta && isTRUE(fasta_ok()))
    flank_ready <- if (!isTRUE(input$use_titer)) TRUE else (if (!is.null(example_fasta_flank_seq())) TRUE else (has_flank_file && isTRUE(flank_ok())))
    
    shinyjs::toggleState("run", csv_ready && fasta_ready && flank_ready)
  })
  
  
  observeEvent(input$csv, {
    csv_ok(FALSE); csv_msg(""); csv_renames(character()); csv_df_norm(NULL)
    
    if (is.null(input$csv)) return()
    
    vd_raw <- tryCatch(
      read.csv(input$csv$datapath, stringsAsFactors = FALSE),
      error = function(e) e
    )
    if (inherits(vd_raw, "error")) { csv_msg(paste0("Could not read CSV: ", vd_raw$message)); return() }
    if (nrow(vd_raw) == 0) { csv_msg("CSV has zero rows."); return() }
    
    norm <- tryCatch(normalize_columns(vd_raw), error = function(e) e)
    if (inherits(norm, "error")) { csv_msg(norm$message); return() }
    
    vd <- norm$df
    csv_df_norm(vd)
    
    missing <- setdiff(c("patient_ID", "DNA_variant"), colnames(vd))
    if (length(missing)) { csv_msg(paste0("Missing required column(s): ", paste(missing, collapse = ", "))); return() }
    
    miss_pid <- sum(is.na(vd$patient_ID) | !nzchar(trimws(vd$patient_ID)))
    if (miss_pid > 0) { csv_msg(paste0("Missing/empty values in patient_ID: ", miss_pid)); return() }
    
    miss_var <- sum(is.na(vd$DNA_variant) | !nzchar(trimws(vd$DNA_variant)))
    if (miss_var > 0) { csv_msg(paste0("Missing/empty values in DNA_variant: ", miss_var)); return() }
    
    dup_pid <- sum(duplicated(vd$patient_ID))
    if (dup_pid > 0) { csv_msg(paste0("Duplicated patient_ID values: ", dup_pid)); return() }
    
    issues <- vapply(vd$DNA_variant, variant_issue, character(1))
    bad <- which(!is.na(issues) & issues != "")
    if (length(bad)) {
      i <- bad[1]
      csv_msg(paste0("Variant issue in ", length(bad), " row(s). Example row ", i, " (", vd$DNA_variant[i], "): ", issues[i]))
      return()
    }
    
    if (length(norm$renames)) {
      csv_renames(norm$renames)
      showNotification(
        paste0("Mapped columns: ", paste0(names(norm$renames), " → ", unname(norm$renames), collapse = ", ")),
        type = "message"
      )
    }
    
    csv_ok(TRUE)
    csv_msg("File OK ✓")
  }, ignoreInit = TRUE)
  
  
  observeEvent(input$fasta, {
    fasta_ok(FALSE); fasta_msg("")
    if (is.null(input$fasta)) return()
    
    lines <- tryCatch(readLines(input$fasta$datapath, warn = FALSE), error = function(e) e)
    if (inherits(lines, "error")) { fasta_msg(paste0("Could not read FASTA: ", lines$message)); return() }
    if (!any(startsWith(lines, ">"))) { fasta_msg("Invalid FASTA: missing header line starting with '>'"); return() }
    
    seq_nt <- toupper(paste(lines[!startsWith(lines, ">")], collapse = ""))
    if (!nzchar(seq_nt)) { fasta_msg("Invalid FASTA: sequence is empty."); return() }
    if (nchar(seq_nt) %% 3 != 0) { fasta_msg("FASTA OK, but CDS length is not divisible by 3."); return() }
    
    fasta_ok(TRUE)
    fasta_msg("File OK ✓")
  }, ignoreInit = TRUE)
  
  observeEvent(input$fasta_flank, {
    flank_ok(FALSE); flank_msg("")
    if (is.null(input$fasta_flank)) return()
    
    lines <- tryCatch(readLines(input$fasta_flank$datapath, warn = FALSE), error = function(e) e)
    if (inherits(lines, "error")) { flank_msg(paste0("Could not read FASTA: ", lines$message)); return() }
    if (!any(startsWith(lines, ">"))) { flank_msg("Invalid FASTA: missing header line starting with '>'"); return() }
    
    seq_nt <- toupper(paste(lines[!startsWith(lines, ">")], collapse = ""))
    if (!nzchar(seq_nt)) { flank_msg("Invalid FASTA: sequence is empty."); return() }
    
    flank_ok(TRUE)
    flank_msg("File OK ✓")
  }, ignoreInit = TRUE)
  
  
  show_intro_modal <- function() {
    showModal(modalDialog(
      title = NULL,
      easyClose = TRUE,
      size = "l",
      footer = tagList(
        actionButton("intro_dismiss", "Let's start →", class = "intro-cta")
      ),
      tags$div(
        tags$p(tags$b("Welcome to ",mfap(),"!"),
               style = "font-size:1.35rem; margin-bottom:10px;"),
        
        tags$p(
          mfap()," helps you quantify protein-level consequences of cDNA variants in single-exon (i.e., intronless) genes.",
          mfap()," is built to support genotype–phenotype correlation studies. To use ",mfap(),", you will upload a table, where each row represents one individual (or sample), and columns contain variables such as a patient ID and the DNA variants in HGVS annotation. In addition, please upload the reference sequence for your gene of interest (e.g., from Ensembl). You may also include columns for phenotype or severity measures but ", mfap()," does not require those."
        ),
        
        tags$div(
          class = "intro-example-shell",
          
          tags$table(
            class = "table table-sm intro-example-table",
            tags$thead(
              tags$tr(
                tags$th("patient_ID"),
                tags$th("DNA_variant"),
                tags$th("severity"),
                tags$th("…")
              )
            ),
            tags$tbody(
              tags$tr(
                tags$td("P001"),
                tags$td("c.123delA"),
                tags$td("0.2"),
                tags$td("…")
              ),
              tags$tr(
                tags$td("P002"),
                tags$td("c.456C>T"),
                tags$td("0.8"),
                tags$td("…")
              ),
              tags$tr(
                tags$td("P003"),
                tags$td("c.89_90insG"),
                tags$td("0.1"),
                tags$td("…")
              ),
              tags$tr(
                tags$td("…"),
                tags$td("…"),
                tags$td("…"),
                tags$td("…")
              )
            )
          ),
          
          tags$div(
            class = "muted",
            style = "margin-top:6px; font-size:0.92rem;",
            "This is an example to explain the table structure. Add as many additional rows and columns as needed."
          )
        ),
        
        tags$p(
          tags$b("Imprinted genes: "),
          "If your gene of interest is imprinted, only include variants on the expressed allele, i.e. the allele that is transcriptionally active in the relevant biological context."
        )      )
    ))
    session$sendCustomMessage("mfap_intro_modal_class", list(on = TRUE))
    
  }
  
  # Show once per browser tab (sessionStorage lives on the client)
  observeEvent(TRUE, {
    session$sendCustomMessage("mfap_intro_once", list(key = "mfap_intro_dismissed_v1"))
  }, once = TRUE)
  
  # When the JS 'mfap_intro_once' handler triggers intro_show, show the modal
  observeEvent(input$intro_show, {
    show_intro_modal()
  }, ignoreInit = TRUE)
  
  # Manual reopen from navbar/link/button
  observeEvent(input$open_intro, {
    show_intro_modal()
  })
  
  # Dismiss: close + persist dismissal flag in sessionStorage
  observeEvent(input$intro_dismiss, {
    removeModal()
    session$sendCustomMessage("mfap_intro_store", list(key = "mfap_intro_dismissed_v1"))
  }, ignoreInit = TRUE)
  
  
  
  observeEvent(input$clear_uploads, {
    shinyjs::reset("csv"); shinyjs::reset("fasta"); shinyjs::reset("fasta_flank")
    example_csv(NULL); example_fasta_seq(NULL); example_fasta_flank_seq(NULL)
    csv_df_norm(NULL); csv_ok(FALSE); csv_msg(""); csv_renames(character())
    fasta_ok(FALSE); fasta_msg("")
    flank_ok(FALSE); flank_msg("")
    result_data(NULL)
    shinyjs::disable("run")
  })
  
  
  observeEvent(input$load_example_data, {
    err <- tryCatch({
      if (!file.exists("www/example_variants.csv"))
        stop("www/example_variants.csv not found. Run the app from the directory that contains app.R and www/.")
      example_csv(read.csv("www/example_variants.csv", stringsAsFactors=FALSE))
      if (!file.exists("www/example_reference.fasta"))
        stop("www/example_reference.fasta not found.")
      lines <- readLines("www/example_reference.fasta", warn=FALSE)
      example_fasta_seq(toupper(paste(lines[!startsWith(lines, ">")], collapse="")))
      if (!file.exists("www/example_reference_flank.fasta"))
        stop("www/example_reference_flank.fasta not found.")
      lines_flank <- readLines("www/example_reference_flank.fasta", warn=FALSE)
      example_fasta_flank_seq(toupper(paste(lines_flank[!startsWith(lines_flank, ">")], collapse="")))
      domains(example_domains); motifs(example_motifs)
      showNotification("Example data loaded.", type="message")
      NULL
    }, error = function(e) e)
    if (!is.null(err)) {
      showNotification(paste0("Could not load demo: ", conditionMessage(err)), type = "error", duration = 10)
    }
  })
  
  observeEvent(input$bib_schubert_bg, {
    showModal(modalDialog(
      title = "BibTeX — Schubert et al. (2025)",
      easyClose = TRUE, size = "m",
      footer = tagList(
        modalButton("Close"),
        tags$button(type="button", class="btn btn-primary", "Copy",
                    `data-copy-target` = "#bib_schubert_bg_text")
      ),
      tags$pre(
        id   = "bib_schubert_bg_text",
        class = "bib-block",
        "@article{Schubert2025,
  author  = {Schubert, T. and Tietzel, A. and Pottayil, H. and Caro, P. and Gilmore, R. B. and Franke, F. and Althammer, F. and Schaaf, C. P.},
  title   = {A blueprint for protein-centric genotype-phenotype investigations in single-exon disease genes applied to MAGEL2 and Schaaf-Yang syndrome},
  journal = {Journal},
  year    = {2025},
  volume  = {VOLUME},
  number  = {ISSUE},
  pages   = {PAGES},
  doi     = {10.xxxx/xxxxx}
}"
      )
    ))
  })
  
  observeEvent(input$bib_titer_bg, {
    showModal(modalDialog(
      title = "BibTeX — Zhang et al. (2017)",
      easyClose = TRUE, size = "m",
      footer = tagList(
        modalButton("Close"),
        tags$button(type="button", class="btn btn-primary", "Copy",
                    `data-copy-target` = "#bib_titer_bg_text")
      ),
      tags$pre(
        id   = "bib_titer_bg_text",
        class = "bib-block",
        "@article{TITER2017,
  author  = {Zhang, S. and Hu, H. and Jiang, T. and Zhang, L. and Zeng, J.},
  title   = {TITER: predicting translation initiation sites by deep learning},
  journal = {Bioinformatics},
  year    = {2017},
  volume  = {33},
  number  = {14},
  pages   = {i234--i242},
  doi     = {10.1093/bioinformatics/btx247}
}"
      )
    ))
  })
  
  
  observeEvent(input$add_dom, {
    req(input$dom_name, input$dom_start, input$dom_end)
    
    start_val <- as.integer(input$dom_start)
    end_val   <- as.integer(input$dom_end)
    
    if (is.na(start_val) || is.na(end_val)) {
      showNotification("Please provide valid numeric Start and End positions.", type = "error")
      return()
    }
    
    if (start_val < 1 || end_val < 1 || start_val > end_val) {
      showNotification("Domain coordinates must be within [1, end] and Start ≤ End.", type = "error")
      return()
    }
    
    max_aa <- wt_aa_len()
    if (!is.null(max_aa) && !is.na(max_aa) && is.finite(max_aa) && end_val > max_aa) {
      showNotification(paste0("Domain end (", end_val, ") exceeds protein length (", max_aa, ")."), type = "error")
      return()
    }
    
    df <- domains()
    if (nrow(df) > 0 && any(df$Name == input$dom_name & df$Start == start_val & df$End == end_val)) {
      showNotification("Exact duplicate domain exists.", type = "error")
      return()
    }
    
    domains(rbind(df, data.frame(
      Name  = as.character(input$dom_name),
      Start = start_val,
      End   = end_val,
      stringsAsFactors = FALSE
    )))
  })
  
  
  observeEvent(input$clr_dom, domains(data.frame(Name=character(),Start=integer(),End=numeric(),stringsAsFactors=FALSE)))
  output$dom_tbl <- renderTable(domains(), striped=TRUE, hover=TRUE)
  output$dom_tbl_wrapper <- renderUI({
    if (nrow(domains()) > 0) {
      tableOutput("dom_tbl") %>% tagAppendAttributes(style = "font-size:0.9em; margin-top:8px;")
    } else NULL
  })

  observeEvent(input$add_mot, {
    req(input$mot_name, input$mot_pattern)
    if (!str_detect(input$mot_pattern, "^[ACDEFGHIKLMNPQRSTVWYBZX]+$")) {showNotification("Pattern must use only IUPAC amino acid letters A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,B,Z,X", type="error"); return()}
    df <- motifs()
    if (any(df$Name==input$mot_name & df$Pattern==input$mot_pattern)) {showNotification("Exact duplicate motif exists.", type="error"); return()}
    motifs(rbind(df, data.frame(Name=input$mot_name, Pattern=input$mot_pattern, stringsAsFactors=FALSE)))
  })
  observeEvent(input$clr_mot, motifs(data.frame(Name=character(),Pattern=character(),stringsAsFactors=FALSE)))
  output$mot_tbl <- renderTable(motifs(), striped=TRUE, hover=TRUE)
  output$mot_tbl_wrapper <- renderUI({
    if (nrow(motifs()) > 0) {
      tableOutput("mot_tbl") %>% tagAppendAttributes(style = "font-size:0.9em; margin-top:8px;")
    } else NULL
  })

  observeEvent(input$clear_example_data, {
    example_csv(NULL); example_fasta_seq(NULL); example_fasta_flank_seq(NULL)
    domains(data.frame(Name=character(),Start=integer(),End=numeric(),stringsAsFactors=FALSE))
    motifs(data.frame(Name=character(),Pattern=character(),stringsAsFactors=FALSE))
    result_data(NULL); showNotification("Example data cleared.", type="message")
  })
  
  observeEvent(input$show_bib_bg, {
    showModal(modalDialog(
      title = "BibTeX",
      easyClose = TRUE,
      size = "m",
      footer = tagList(
        modalButton("Close"),
        tags$button(type="button", class="btn btn-primary copy-bib", "Copy")
      ),
      tags$pre(
        class = "bib-block",
        "@article{Paper,
  author  = {Schubert, T. and Tietzel, A. and Pottayil, H. and Caro, P. and Gilmore, R. B. and Franke, F. and Althammer, F. and Schaaf, C.P.},
  title   = {A blueprint for protein-centric genotype-phenotype investigations in single-exon disease genes applied to MAGEL2 and Schaaf-Yang syndrome},
  journal = {Journal},
  year    = {2025},
  volume  = {VOLUME},
  number  = {ISSUE},
  pages   = {PAGES},
  doi     = {10.xxxx/xxxxx}
}"
      ),
      tags$script(HTML("
      document.addEventListener('click', function(e){
        if (e.target && e.target.classList.contains('copy-bib')) {
          var modal = e.target.closest('.modal-content');
          var pre = modal && modal.querySelector('pre');
          if (pre) {
            navigator.clipboard.writeText(pre.innerText).then(function(){
              var old = e.target.textContent;
              e.target.textContent = 'Copied';
              setTimeout(function(){ e.target.textContent = old; }, 1200);
            });
          }
        }
      });
    "))
    ))
  })
  
  observeEvent(input$show_titer_bibtex_bg, {
    showModal(modalDialog(
      title = "TITER – BibTeX",
      easyClose = TRUE,
      size = "m",
      footer = tagList(
        modalButton("Close"),
        tags$button(type="button", class="btn btn-primary copy-bib", "Copy")
      ),
      tags$pre(
        class = "bib-block",
        "@article{TITER2017,
  author  = {Zhang, S. and Hu, H. and Jiang, T. and Zhang, L. and Zeng, J.},
  title   = {TITER: predicting translation initiation sites by deep learning},
  journal = {Bioinformatics},
  year    = {2017},
  volume  = {33},
  number  = {14},
  pages   = {i234--i242},
  doi     = {10.1093/bioinformatics/btx247}
}"
      ),
      tags$script(HTML("
      document.addEventListener('click', function(e){
        if (e.target && e.target.classList.contains('copy-bib')) {
          var modal = e.target.closest('.modal-content');
          var pre = modal && modal.querySelector('pre');
          if (pre) {
            navigator.clipboard.writeText(pre.innerText).then(function(){
              var old = e.target.textContent;
              e.target.textContent = 'Copied';
              setTimeout(function(){ e.target.textContent = old; }, 1200);
            });
          }
        }
      });
    "))
    ))
  })
  
  observeEvent(input$upload_help_trigger, {
    showModal(modalDialog(
      title = "Format & download help",
      easyClose = TRUE, size = "m",
      footer = modalButton("Close"),
      tags$div(
        tags$h6("Table requirements", class = "modal-subheading"),
        tags$p(class = "muted", style = "font-size:0.9rem; margin-bottom:1rem;",
               "Your CSV must have at least two columns: ",
               tags$b("DNA_variant"), " (HGVS, e.g. c.123delA) and ",
               tags$b("patient_ID"), " as the individual identifier."),
        tags$h6("Obtaining a FASTA", class = "modal-subheading"),
        tags$p(class = "muted", style = "font-size:0.9rem; margin-bottom:0.5rem;", "Example workflow:"),
        tags$ul(class = "muted", style = "font-size:0.9rem; margin-bottom:1rem; padding-left:1.2rem;",
                tags$li("Find the gene on ensembl.org"),
                tags$li("Click Download sequence"),
                tags$li("Included sequences: all; Format: FASTA; Flanks: 0"),
                tags$li("Download")),
        tags$h6("FASTA with flanks (for TITER)", class = "modal-subheading"),
        tags$p(class = "muted", style = "font-size:0.9rem; margin-bottom:0.5rem;",
               "Same workflow as above, but set Flanks to 100 instead of 0."),
        tags$h6("Why flanks for TITER?", class = "modal-subheading"),
        tags$p(class = "muted", style = "font-size:0.9rem; margin-bottom:0;",
               "TITER uses a 203 bp window around the central bases. A CDS file with ±100 bp flanks (5'/3' UTR) is needed to evaluate start sites near the 5' end.")
      )
    ))
  })

  observeEvent(input$show_bib, {
    showModal(modalDialog(
      title = "BibTeX — Schubert et al. (2025)",
      easyClose = TRUE, size = "m",
      footer = tagList(
        modalButton("Close"),
        tags$button(type="button", class="btn btn-primary", "Copy",
                    `data-copy-target` = "#bib_schubert_text")
      ),
      tags$pre(
        id   = "bib_schubert_text",
        class = "bib-block",
        "@article{Schubert2025,
  author  = {Schubert, T. and Tietzel, A. and Pottayil, H. and Caro, P. and Gilmore, R. B. and Franke, F. and Althammer, F. and Schaaf, C. P.},
  title   = {A blueprint for protein-centric genotype-phenotype investigations in single-exon disease genes applied to MAGEL2 and Schaaf-Yang syndrome},
  journal = {Journal},
  year    = {2025},
  volume  = {VOLUME},
  number  = {ISSUE},
  pages   = {PAGES},
  doi     = {10.xxxx/xxxxx}
}"
      )
    ))
  })
  
  observeEvent(input$show_titer_bibtex, {
    showModal(modalDialog(
      title = "BibTeX — Zhang et al. (2017)",
      easyClose = TRUE, size = "m",
      footer = tagList(
        modalButton("Close"),
        tags$button(type="button", class="btn btn-primary", "Copy",
                    `data-copy-target` = "#bib_zhang_text")
      ),
      tags$pre(
        id   = "bib_zhang_text",
        class = "bib-block",
        "@article{TITER2017,
  author  = {Zhang, S. and Hu, H. and Jiang, T. and Zhang, L. and Zeng, J.},
  title   = {TITER: predicting translation initiation sites by deep learning},
  journal = {Bioinformatics},
  year    = {2017},
  volume  = {33},
  number  = {14},
  pages   = {i234--i242},
  doi     = {10.1093/bioinformatics/btx247}
}"
      )
    ))
  })
  
  observeEvent(input$open_legal, {
    showModal(modalDialog(
      title = "Legal notice (Impressum)",
      easyClose = TRUE,
      size = "l",
      footer = modalButton("Close"),
      
      tags$div(
        style = "line-height:1.45;",
        
        tags$h5("Angaben gemäß § 5 TMG", style="margin-top:0;"),
        tags$p(
          "Institut für Humangenetik", tags$br(),
          "Tim Schubert", tags$br(),
          "Im Neuenheimer Feld 366, 5. OG", tags$br(),
          "69120 Heidelberg", tags$br(),
          "Germany"
        ),
        
        tags$h5("Kontakt"),
        tags$p(
          "E-Mail: tim dot schubert at med dot uni-heidelberg dot de"
        ),
        
        tags$h5("Verantwortlich für den Inhalt nach § 55 Abs. 2 RStV"),
        tags$p(
          "Institut für Humangenetik", tags$br(),
          "Tim Schubert", tags$br(),
          "Im Neuenheimer Feld 366, 5. OG", tags$br(),
          "69120 Heidelberg", tags$br(),
          "Germany"
        ),
        
        tags$h5("Verbraucherstreitbeilegung / Universalschlichtungsstelle"),
        tags$p(class="muted",
               "Wir sind nicht bereit oder verpflichtet, an Streitbeilegungsverfahren vor einer Verbraucherschlichtungsstelle teilzunehmen."
        ),
        
        tags$h5("Haftung für Inhalte"),
        tags$p(
          class="muted",
          "Als Diensteanbieter sind wir gemäß § 7 Abs.1 TMG für eigene Inhalte auf diesen Seiten nach den allgemeinen Gesetzen verantwortlich. Nach §§ 8 bis 10 TMG sind wir als Diensteanbieter jedoch nicht verpflichtet, übermittelte oder gespeicherte fremde Informationen zu überwachen oder nach Umständen zu forschen, die auf eine rechtswidrige Tätigkeit hinweisen.

Verpflichtungen zur Entfernung oder Sperrung der Nutzung von Informationen nach den allgemeinen Gesetzen bleiben hiervon unberührt. Eine diesbezügliche Haftung ist jedoch erst ab dem Zeitpunkt der Kenntnis einer konkreten Rechtsverletzung möglich. Bei Bekanntwerden von entsprechenden Rechtsverletzungen werden wir diese Inhalte umgehend entfernen."
        ),
        
        tags$h5("Haftung für Links"),
        tags$p(
          class="muted",
          "Unser Angebot enthält Links zu externen Websites Dritter, auf deren Inhalte wir keinen Einfluss haben. Deshalb können wir für diese fremden Inhalte auch keine Gewähr übernehmen. Für die Inhalte der verlinkten Seiten ist stets der jeweilige Anbieter oder Betreiber der Seiten verantwortlich. Die verlinkten Seiten wurden zum Zeitpunkt der Verlinkung auf mögliche Rechtsverstöße überprüft. Rechtswidrige Inhalte waren zum Zeitpunkt der Verlinkung nicht erkennbar.

Eine permanente inhaltliche Kontrolle der verlinkten Seiten ist jedoch ohne konkrete Anhaltspunkte einer Rechtsverletzung nicht zumutbar. Bei Bekanntwerden von Rechtsverletzungen werden wir derartige Links umgehend entfernen."
        ),
        
        tags$h5("Urheberrecht"),
        tags$p(
          class="muted",
          "Die durch die Seitenbetreiber erstellten Inhalte und Werke auf diesen Seiten unterliegen dem deutschen Urheberrecht. Die Vervielfältigung, Bearbeitung, Verbreitung und jede Art der Verwertung außerhalb der Grenzen des Urheberrechtes bedürfen der schriftlichen Zustimmung des jeweiligen Autors bzw. Erstellers. Downloads und Kopien dieser Seite sind nur für den privaten, nicht kommerziellen Gebrauch gestattet."
        ),
        tags$p(
          class="muted",
          "Soweit die Inhalte auf dieser Seite nicht vom Betreiber erstellt wurden, werden die Urheberrechte Dritter beachtet. Insbesondere werden Inhalte Dritter als solche gekennzeichnet. Sollten Sie trotzdem auf eine Urheberrechtsverletzung aufmerksam werden, bitten wir um einen entsprechenden Hinweis. Bei Bekanntwerden von Rechtsverletzungen werden wir derartige Inhalte umgehend entfernen. Quelle: e-recht24.de"
        )
      )
    ))
  })
  
  output$results_body <- renderUI({
    df <- result_data()
    
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
      tags$div(
        class="muted",
        style="padding:26px 6px; text-align:center;",
        tags$p(style="margin:0 0 6px 0; font-weight:600; font-size:1.05rem;",
               "No results to be displayed, yet."),
        tags$p(style="margin:0;",
               "Insert data on the left and click “Analyze”. Your results will appear here.")
        
        
      )
    } else {
      DTOutput("res_tbl")
    }
  })
  
  
  
  
  observeEvent(input$run, {
    shinyjs::disable("run"); shinyjs::hide("res_tbl")
    req(!is.null(input$csv) || !is.null(example_csv()))
    vd <- if (!is.null(example_csv())) {
      example_csv()
    } else if (!is.null(csv_df_norm())) {
      csv_df_norm()
    } else {
      read.csv(input$csv$datapath, stringsAsFactors=FALSE)
    }
    
    if (nrow(vd)==0) {
      showModal(modalDialog(title="Empty CSV","Your variant list file has zero rows.", easyClose=TRUE, footer=modalButton("OK")))
      result_data(NULL); shinyjs::enable("run"); return(NULL)
    }
    
    if (!all(c("patient_ID","DNA_variant") %in% colnames(vd))) {
      showModal(modalDialog(title="Missing column","Your CSV must contain columns `patient_ID` and `DNA_variant`.", easyClose=TRUE, footer=modalButton("OK")))
      result_data(NULL); shinyjs::enable("run"); return(NULL)
    }
    
    miss_pid <- sum(is.na(vd$patient_ID) | !nzchar(trimws(vd$patient_ID)))
    if (miss_pid > 0) {
      showModal(modalDialog(title="Missing values", paste0("patient_ID has ", miss_pid, " missing/empty values."), easyClose=TRUE, footer=modalButton("OK")))
      result_data(NULL); shinyjs::enable("run"); return(NULL)
    }
    
    miss_var <- sum(is.na(vd$DNA_variant) | !nzchar(trimws(vd$DNA_variant)))
    if (miss_var > 0) {
      showModal(modalDialog(title="Missing values", paste0("DNA_variant has ", miss_var, " missing/empty values."), easyClose=TRUE, footer=modalButton("OK")))
      result_data(NULL); shinyjs::enable("run"); return(NULL)
    }
    
    dup_pid <- sum(duplicated(vd$patient_ID))
    if (dup_pid > 0) {
      showModal(modalDialog(title="Duplicated IDs", paste0("patient_ID contains ", dup_pid, " duplicated value(s)."), easyClose=TRUE, footer=modalButton("OK")))
      result_data(NULL); shinyjs::enable("run"); return(NULL)
    }
    
    issues <- vapply(vd$DNA_variant, variant_issue, character(1))
    bad <- which(!is.na(issues) & issues != "")
    if (length(bad)) {
      i <- bad[1]
      showModal(modalDialog(
        title="Variant format issue",
        paste0("Example row ", i, " (", vd$DNA_variant[i], "): ", issues[i]),
        easyClose=TRUE,
        footer=modalButton("OK")
      ))
      result_data(NULL); shinyjs::enable("run"); return(NULL)
    }
    
    if (isTRUE(input$use_titer)) {
      if (is.null(input$fasta_flank) && is.null(example_fasta_flank_seq())) {
        showModal(modalDialog(title="Missing FASTA","Upload the FASTA with ±100 bp flanks for TITER.", easyClose=TRUE, footer=modalButton("OK")))
        result_data(NULL); shinyjs::enable("run"); return(NULL)
      }
    }
    
    
    lines <- if (!is.null(input$fasta)) readLines(input$fasta$datapath, warn=FALSE) else c(">example", example_fasta_seq())
    if (!any(startsWith(lines, ">"))) {showModal(modalDialog(title="Invalid FASTA","FASTA must have at least one header line beginning with `>`.", easyClose=TRUE, footer=modalButton("OK"))); result_data(NULL); return(NULL)}
    seq_nt <- toupper(paste(lines[!startsWith(lines, ">")], collapse=""))
    if (nchar(seq_nt) %% 3 != 0) {
      showModal(modalDialog(title="Sequence not multiple of 3","Reference length is not divisible by 3.", easyClose=TRUE, footer=modalButton("OK")))
      result_data(NULL); shinyjs::enable("run"); return(NULL)
    }
    
    probs <- validate_variants_against_ref(vd$DNA_variant, seq_nt)
    bad_sem <- which(nzchar(probs))
    if (length(bad_sem)) {
      ex <- head(bad_sem, 8)
      showModal(modalDialog(
        title = "Variant validation failed",
        easyClose = TRUE,
        size = "l",
        footer = modalButton("OK"),
        tags$p("Some variants are incompatible with the provided CDS FASTA."),
        tags$ul(lapply(ex, function(i) tags$li(paste0("Row ", i, " (", vd$DNA_variant[i], "): ", probs[i])))),
        if (length(bad_sem) > length(ex)) tags$p(class="muted", paste0("... and ", length(bad_sem) - length(ex), " more.")) else NULL
      ))
      result_data(NULL); shinyjs::enable("run"); return(NULL)
    }
    
    aa_len <- nchar(Biostrings::translate(Biostrings::DNAString(seq_nt), if.fuzzy.codon="X"))
    wt_aa_len(aa_len)
    
    
    withProgress(message="Running analysis...", value=0, {
      incProgress(0.1, detail="Parsing mutations")
      parsed <- tryCatch(
        lapply(vd$DNA_variant, parse_mutation_strict, wildtype_seq = seq_nt),
        error = function(e) e
      )
      if (inherits(parsed, "error")) {
        showModal(modalDialog(title="Variant parsing failed", parsed$message, easyClose=TRUE, footer=modalButton("OK")))
        result_data(NULL); shinyjs::enable("run"); return(NULL)
      }
      
      df <- vd %>% dplyr::mutate(
        Locus=sapply(parsed, `[[`,"Locus"),
        Mutation_Type=sapply(parsed, `[[`,"Mutation_Type"),
        Deleted_Bases=sapply(parsed, `[[`,"Deleted_Bases"),
        Duplicated_Bases=sapply(parsed, `[[`,"Duplicated_Bases"),
        Inserted_Bases=sapply(parsed, `[[`,"Inserted_Bases"),
        Mutated_Sequence=sapply(parsed, `[[`,"Mutated_Sequence")
      )
      
      incProgress(0.2, detail="Translating proteins")
      prot <- tryCatch(
        mapply(translate_prot_strict, df$Mutated_Sequence, df$DNA_variant, USE.NAMES = FALSE),
        error = function(e) e
      )
      if (inherits(prot, "error")) {
        showModal(modalDialog(title="Protein translation failed", prot$message, easyClose=TRUE, footer=modalButton("OK")))
        result_data(NULL); shinyjs::enable("run"); return(NULL)
      }
      df$Mutated_Protein <- prot
      df$Protein_Length_aa <- nchar(df$Mutated_Protein)
      
      
      incProgress(0.2, detail="Refining mutation types")
      df <- df %>% dplyr::rowwise() %>% dplyr::mutate(Mutation_Type = dplyr::case_when(
        Mutation_Type=="Point" ~ refine_point(Locus, Mutated_Sequence, seq_nt, GENETIC_CODE),
        Mutation_Type %in% c("Deletion","Duplication","Insertion","Indel") ~ refine_indel(Mutation_Type, Deleted_Bases, Duplicated_Bases, Inserted_Bases),
        TRUE ~ Mutation_Type
      )) %>% dplyr::ungroup()
      
      incProgress(0.2, detail="Annotating domains & motifs")
      doms <- domains(); mot <- motifs()
      df <- df %>% dplyr::mutate(
        AA_Position = floor((Locus - 1)/3) + 1,
        Domain_Location_Of_Variant = mapply(function(p){
          hits <- doms$Name[p >= doms$Start & p <= doms$End]; if (!length(hits)) "" else paste(hits, collapse='; ')
        }, AA_Position),
        Lost_Functional_Domains = mapply(function(p, mt){
          hits <- doms$Name[p >= doms$Start & p <= doms$End]; if (!length(hits)) return(NA_character_)
          mt <- tolower(mt)
          if (mt %in% c('nonsense','frameshifting indel')) {
            i <- which(doms$Name %in% hits); lost <- doms$Name[i:length(doms$Name)]
          } else if (mt %in% c('missense','in-frame indel')) { lost <- hits } else return(NA_character_)
          paste(lost, collapse='; ')
        }, AA_Position, Mutation_Type)
      )
      
      if (nrow(mot)>0) {
        for (i in seq_len(nrow(mot))) {
          nm <- mot$Name[i]; pat <- mot$Pattern[i]
          df[[nm]] <- ifelse(is.na(df$Mutated_Protein),"", ifelse(stringr::str_detect(df$Mutated_Protein, pat),"1","0"))
        }
        df <- df %>% dplyr::rowwise() %>% dplyr::mutate(`Lost Motifs`={
          vals <- dplyr::c_across(dplyr::all_of(mot$Name)); miss <- mot$Name[vals=="0"]; if (!length(miss)) "" else paste(miss, collapse='; ')
        }) %>% dplyr::ungroup()
      } else df$`Lost Motifs` <- ""
      
      # Non-canonical TIS via TITER
      if (isTRUE(input$use_titer)) {
        incProgress(0.05, detail = "Running non-canonical TIS analysis (TITER)")
        titer_status <- ensure_titer_python_ready()
        if (!isTRUE(titer_status$ok)) {
          showNotification(
            paste0(titer_status$message, " Skipping TITER analysis; TITER columns are set to NA."),
            type = "warning",
            duration = 12
          )
          df <- add_missing_titer_columns(df)
        } else {
          titer_ok <- tryCatch({
            # Resolve & copy bundled `titer/` to a writable temp workdir
            app_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
            titer_src <- normalizePath(file.path(app_dir, "titer"), winslash = "/", mustWork = FALSE)
            if (!dir.exists(titer_src)) {
              stop("The 'titer' folder is missing from the deployed bundle. ",
                   "Ensure it sits next to app.R and is NOT ignored by .gitignore/.Rbuildignore.")
            }
            temp_titer_dir <- file.path(tempdir(), "titer")
            if (dir.exists(temp_titer_dir)) unlink(temp_titer_dir, recursive = TRUE)
            dir.create(temp_titer_dir, recursive = TRUE, showWarnings = FALSE)
            ok <- file.copy(list.files(titer_src, full.names = TRUE, all.files = TRUE, no.. = TRUE),
                            temp_titer_dir, recursive = TRUE)
            if (!all(ok)) stop("Failed to copy TITER files into tempdir().")

            # Ensure data dir exists
            dir.create(file.path(temp_titer_dir, "data"), recursive = TRUE, showWarnings = FALSE)

            # Put the flanked FASTA into data/
            if (!is.null(input$fasta_flank)) {
              file.copy(input$fasta_flank$datapath,
                        file.path(temp_titer_dir, "data", basename(input$fasta_flank$name)),
                        overwrite = TRUE)
            } else {
              file.copy("www/example_reference_flank.fasta",
                        file.path(temp_titer_dir, "data", "example_reference_flank.fasta"),
                        overwrite = TRUE)
            }

            # Write the CSV TITER expects
            write.csv(df %>% dplyr::select(patient_ID, Mutated_Sequence, DNA_variant),
                      file.path(temp_titer_dir, "data", "variant_list_with_mutated_sequences.csv"),
                      row.names = FALSE)

            # Run Python
            script_path <- file.path(temp_titer_dir, "codes", "analyze_patients_for_variant_specific_additional_TIS.py")
            reticulate::source_python(script_path)
            summary_path <- file.path(temp_titer_dir, "data", "summary_patients_most_likely_additional_TIS.csv")
            if (!file.exists(summary_path))
              stop("TITER did not produce summary file. Check that the Python script completed successfully.")
            titer_summary <- read.csv(summary_path, stringsAsFactors = FALSE)
            missing_summary <- setdiff(titer_summary_columns, names(titer_summary))
            if (length(missing_summary) > 0) {
              stop("TITER output is missing required columns: ", paste(missing_summary, collapse = ", "))
            }
            list(summary = titer_summary)
          }, error = function(e) {
            showNotification(
              paste0("TITER analysis failed; TITER columns are set to NA. ", conditionMessage(e)),
              type = "warning",
              duration = 12
            )
            NULL
          })

          if (!is.null(titer_ok)) {
            titer_summary <- titer_ok$summary %>% dplyr::select(-dplyr::any_of("DNA_variant"))
            df <- dplyr::left_join(df, titer_summary, by = "patient_ID") %>%
              dplyr::mutate(
                Protein_from_most_likely_non_canonical_TIS          = sapply(RNA_sequence_most_likely_non_canonical_TIS, translate_prot),
                Protein_from_most_likely_non_canonical_in_frame_TIS = sapply(RNA_sequence_most_likely_non_canonical_in_frame_TIS, translate_prot),
                Protein_Length_from_most_likely_non_canonical_TIS          = nchar(Protein_from_most_likely_non_canonical_TIS),
                Protein_Length_from_most_likely_non_canonical_in_frame_TIS = nchar(Protein_from_most_likely_non_canonical_in_frame_TIS)
              ) %>% {
                df_inner <- .
                df_inner <- df_inner %>% dplyr::mutate(
                  Mutation_AA_Pos_Canonical = floor((Locus - 1) / 3) + 1,
                  Main_TIS_bp    = as.numeric(Most_likely_non_canonical_TIS_CDS_Position),
                  Main_TIS_Codon = floor((Main_TIS_bp - 1) / 3) + 1,
                  Mutation_AA_Pos_Main_TIS = dplyr::if_else(
                    Mutation_Type == "Frameshifting indel" & Mutation_AA_Pos_Canonical >= Main_TIS_Codon,
                    Mutation_AA_Pos_Canonical - (Main_TIS_Codon - 1), NA_integer_),
                  Unaltered_Raw_Main  = Mutation_AA_Pos_Main_TIS - 1,
                  Frameshift_Raw_Main = Protein_Length_from_most_likely_non_canonical_TIS - Unaltered_Raw_Main,
                  Unaltered_Length_from_most_likely_non_canonical_TIS = dplyr::if_else(
                    Mutation_Type == "Frameshifting indel" &
                      !is.na(Unaltered_Raw_Main) & Unaltered_Raw_Main >= 0 &
                      Unaltered_Raw_Main <= Protein_Length_from_most_likely_non_canonical_TIS,
                    Unaltered_Raw_Main, NA_integer_),
                  Frameshift_Length_from_most_likely_non_canonical_TIS = dplyr::if_else(
                    Mutation_Type == "Frameshifting indel" &
                      !is.na(Frameshift_Raw_Main) & Frameshift_Raw_Main >= 0,
                    Frameshift_Raw_Main, NA_integer_),
                  InFrame_TIS_bp    = as.numeric(Most_likely_non_canonical_in_frame_TIS_CDS_Position),
                  InFrame_TIS_Codon = floor((InFrame_TIS_bp - 1) / 3) + 1,
                  Mutation_AA_Pos_Alternative = dplyr::if_else(
                    Mutation_Type == "Frameshifting indel" & Mutation_AA_Pos_Canonical >= InFrame_TIS_Codon,
                    Mutation_AA_Pos_Canonical - (InFrame_TIS_Codon - 1), NA_integer_),
                  Unaltered_Raw_InFrame  = Mutation_AA_Pos_Alternative - 1,
                  Frameshift_Raw_InFrame = Protein_Length_from_most_likely_non_canonical_in_frame_TIS - Unaltered_Raw_InFrame,
                  Unaltered_Length_from_most_likely_non_canonical_in_frame_TIS = dplyr::if_else(
                    Mutation_Type == "Frameshifting indel" &
                      !is.na(Unaltered_Raw_InFrame) & Unaltered_Raw_InFrame >= 0 &
                      Unaltered_Raw_InFrame <= Protein_Length_from_most_likely_non_canonical_in_frame_TIS,
                    Unaltered_Raw_InFrame, NA_integer_),
                  Frameshift_Length_from_most_likely_non_canonical_in_frame_TIS = dplyr::if_else(
                    Mutation_Type == "Frameshifting indel" &
                      !is.na(Frameshift_Raw_InFrame) & Frameshift_Raw_InFrame >= 0,
                    Frameshift_Raw_InFrame, NA_integer_)
                ) %>% dplyr::select(
                  -Mutation_AA_Pos_Canonical, -Main_TIS_bp, -Main_TIS_Codon,
                  -Mutation_AA_Pos_Main_TIS, -Unaltered_Raw_Main, -Frameshift_Raw_Main,
                  -InFrame_TIS_bp, -InFrame_TIS_Codon, -Mutation_AA_Pos_Alternative,
                  -Unaltered_Raw_InFrame, -Frameshift_Raw_InFrame
                )
                df_inner
              }
            df <- add_missing_titer_columns(df)
          } else {
            df <- add_missing_titer_columns(df)
          }
        }
      }
      incProgress(0.1, detail="Calculating frameshift metrics")
      df <- df %>% dplyr::mutate(
        Unaltered_Length_aa  = dplyr::if_else(Mutation_Type=='Frameshifting indel', (floor((Locus-1)/3)+1)-1, NA_integer_),
        Frameshift_Length_aa = dplyr::if_else(Mutation_Type=='Frameshifting indel', Protein_Length_aa - ((floor((Locus-1)/3)+1)-1), NA_integer_)
      )
      
      colnames(df) <- gsub("_"," ", colnames(df))
      df <- df %>% dplyr::select(-c(`Mutated Sequence`,`Mutated Protein`), dplyr::everything(), `Mutated Sequence`,`Mutated Protein`)
      result_data(df); shinyjs::enable("run")
    })
  })
  
  output$res_tbl <- renderDT({
    req(result_data())
    
    datatable(
      result_data(),
      rownames = FALSE,
      filter   = "top",
      escape   = TRUE,
      options  = list(
        pageLength = 10,
        scrollX    = TRUE,
        autoWidth  = TRUE,
        dom        = "tip",
        createdRow = JS("
        function(row, data, dataIndex){
          $('td', row).each(function(i){
            var txt = $(this).text();
            $(this).attr('title', txt);
          });
        }
      "),
        initComplete = JS("
        function(){
          var api = this;
          var wrapper = $(api.table().container());
          var applyPaginateStyle = function(){
            if (!wrapper.length) return;
            var isLight = document.documentElement.getAttribute('data-mfap-theme') === 'light';
            var gray = isLight ? '#94a3b8' : '#495057';
            var blue = isLight ? '#2563eb' : '#4dabf7';
            var grayDisabled = isLight ? '#cbd5e1' : '#343a40';
            var grayColor = isLight ? '#64748b' : 'rgba(255,255,255,0.5)';
            var styleId = 'mfap-dt-pagination-style';
            var styleEl = document.getElementById(styleId);
            if (!styleEl) {
              styleEl = document.createElement('style');
              styleEl.id = styleId;
              document.head.appendChild(styleEl);
            }
            var sel = '#res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button';
            styleEl.textContent = sel + '{ background:' + gray + ' !important; border-color:' + gray + ' !important; color:#fff !important; }' +
              sel + '.current{ background:' + blue + ' !important; border-color:' + blue + ' !important; color:#fff !important; }' +
              sel + ':hover:not(.disabled){ background:' + blue + ' !important; border-color:' + blue + ' !important; color:#fff !important; }' +
              sel + '.disabled{ background:' + grayDisabled + ' !important; border-color:' + grayDisabled + ' !important; color:' + grayColor + ' !important; }';
            wrapper.find('.dataTables_paginate .paginate_button').each(function(){
              var el = $(this);
              var e = el[0];
              e.classList.remove('btn-success');
              if (el.hasClass('disabled')) {
                e.style.setProperty('background-color', grayDisabled, 'important');
                e.style.setProperty('border-color', grayDisabled, 'important');
                e.style.setProperty('color', grayColor, 'important');
              } else if (el.hasClass('current')) {
                e.style.setProperty('background-color', blue, 'important');
                e.style.setProperty('border-color', blue, 'important');
                e.style.setProperty('color', '#fff', 'important');
              } else {
                e.style.setProperty('background-color', gray, 'important');
                e.style.setProperty('border-color', gray, 'important');
                e.style.setProperty('color', '#fff', 'important');
              }
            });
            wrapper.find('.dataTables_paginate .paginate_button:not(.disabled)').off('mouseenter mouseleave').on('mouseenter', function(){
              this.style.setProperty('background-color', blue, 'important');
              this.style.setProperty('border-color', blue, 'important');
              this.style.setProperty('color', '#fff', 'important');
            }).on('mouseleave', function(){
              if (!$(this).hasClass('current')) {
                this.style.setProperty('background-color', gray, 'important');
                this.style.setProperty('border-color', gray, 'important');
                this.style.setProperty('color', '#fff', 'important');
              }
            });
          };
          applyPaginateStyle();
          setTimeout(applyPaginateStyle, 200);
        }
      "),
        drawCallback = JS("
        function(){
          var api = this;
          var wrapper = $(api.table().container());
          var applyPaginateStyle = function(){
            if (!wrapper.length) return;
            var isLight = document.documentElement.getAttribute('data-mfap-theme') === 'light';
            var gray = isLight ? '#94a3b8' : '#495057';
            var blue = isLight ? '#2563eb' : '#4dabf7';
            var grayDisabled = isLight ? '#cbd5e1' : '#343a40';
            var grayColor = isLight ? '#64748b' : 'rgba(255,255,255,0.5)';
            var styleId = 'mfap-dt-pagination-style';
            var styleEl = document.getElementById(styleId);
            if (!styleEl) {
              styleEl = document.createElement('style');
              styleEl.id = styleId;
              document.head.appendChild(styleEl);
            }
            var sel = '#res_tbl .dataTables_wrapper .dataTables_paginate .paginate_button';
            styleEl.textContent = sel + '{ background:' + gray + ' !important; border-color:' + gray + ' !important; color:#fff !important; }' +
              sel + '.current{ background:' + blue + ' !important; border-color:' + blue + ' !important; color:#fff !important; }' +
              sel + ':hover:not(.disabled){ background:' + blue + ' !important; border-color:' + blue + ' !important; color:#fff !important; }' +
              sel + '.disabled{ background:' + grayDisabled + ' !important; border-color:' + grayDisabled + ' !important; color:' + grayColor + ' !important; }';
            wrapper.find('.dataTables_paginate .paginate_button').each(function(){
              var el = $(this);
              var e = el[0];
              e.classList.remove('btn-success');
              if (el.hasClass('disabled')) {
                e.style.setProperty('background-color', grayDisabled, 'important');
                e.style.setProperty('border-color', grayDisabled, 'important');
                e.style.setProperty('color', grayColor, 'important');
              } else if (el.hasClass('current')) {
                e.style.setProperty('background-color', blue, 'important');
                e.style.setProperty('border-color', blue, 'important');
                e.style.setProperty('color', '#fff', 'important');
              } else {
                e.style.setProperty('background-color', gray, 'important');
                e.style.setProperty('border-color', gray, 'important');
                e.style.setProperty('color', '#fff', 'important');
              }
            });
            wrapper.find('.dataTables_paginate .paginate_button:not(.disabled)').off('mouseenter mouseleave').on('mouseenter', function(){
              this.style.setProperty('background-color', blue, 'important');
              this.style.setProperty('border-color', blue, 'important');
              this.style.setProperty('color', '#fff', 'important');
            }).on('mouseleave', function(){
              if (!$(this).hasClass('current')) {
                this.style.setProperty('background-color', gray, 'important');
                this.style.setProperty('border-color', gray, 'important');
                this.style.setProperty('color', '#fff', 'important');
              }
            });
          };
          applyPaginateStyle();
          setTimeout(applyPaginateStyle, 150);
        }
      ")
      )
    )
  })
  
  
  
  output$dl <- downloadHandler(
    filename = function() paste0('MAP_results_', Sys.Date(), '.csv'),
    content  = function(f) write.csv(result_data(), f, row.names=FALSE)
  )
  
  output$results_header <- renderUI({
    df <- result_data()
    has_res <- !is.null(df) && is.data.frame(df) && nrow(df) > 0
    
    tags$div(
      style = "display:flex; align-items:center; justify-content:space-between; gap:10px; width:100%;",
      tags$span("Results"),
      if (has_res) downloadButton("dl", "Download CSV", class = "btn btn-mfap-primary btn-sm") else NULL
    )
  })
  
  
  output$csv_status <- renderUI({
    msg <- csv_msg()
    if (is.null(msg) || msg == "") return(NULL)
    ok <- isTRUE(csv_ok())
    tags$div(
      class = if (ok) "text-success" else "text-warning",
      style = "font-size:0.88em; margin-top:6px;",
      msg
    )
  })
  
  output$fasta_status <- renderUI({
    msg <- fasta_msg()
    if (is.null(msg) || msg == "") return(NULL)
    ok <- isTRUE(fasta_ok())
    tags$div(
      class = if (ok) "text-success" else "text-warning",
      style = "font-size:0.88em; margin-top:6px;",
      msg
    )
  })
  
  output$fasta_flank_status <- renderUI({
    req(input$use_titer)
    msg <- flank_msg()
    if (is.null(msg) || msg == "") return(NULL)
    ok <- isTRUE(flank_ok())
    tags$div(
      class = if (ok) "text-success" else "text-warning",
      style = "font-size:0.88em; margin-top:6px;",
      msg
    )
  })
  
  output$csv_ui <- renderUI({
    if (is.null(example_csv())) {
      fileInput("csv","Variant table (.csv)", accept=c(".csv","text/csv"))
    } else tags$div(style="color:#ccc;font-size:0.9em;margin-bottom:10px;","Using example Variant table")
  })
  output$fasta_ui <- renderUI({
    if (is.null(example_fasta_seq())) {
      fileInput("fasta","Reference sequence (FASTA)", accept=c(".fa",".fasta","text/plain"))
    } else tags$div(style="color:#ccc;font-size:0.9em;margin-bottom:10px;","Using example Reference sequence")
  })
  output$fasta_flank_ui <- renderUI({
    req(input$use_titer)
    if (is.null(example_fasta_flank_seq())) {
      tagList(
        fileInput("fasta_flank", "Reference FASTA with 100 bp flanks", accept = c(".fa", ".fasta", "text/plain")),
        tags$div(class="muted", style="font-size:0.8em; margin-top:5px;", "Upload CDS ±100 bp.")
      )
    } else {
      tags$div(style="color:#ccc;font-size:0.8em;margin-bottom:10px;",
               "Using example FASTA with 100 bp flanks for TITER.")
    }
  })
}

shinyApp(ui, server)
