library(shiny)
library(shinyjs)
library(DT)
library(stringr)
library(dplyr)
library(Biostrings)
library(bslib)

# ———————————————— PYTHON VENV SETUP ————————————————
venv_name   <- "r-reticulate"
venv_path   <- file.path("~", ".virtualenvs", venv_name)
system_python <- "/usr/bin/python3"
if (!dir.exists(venv_path)) {
  message("Creating Python venv at ", venv_path)
  reticulate::virtualenv_create(envname = venv_name, python = system_python)
  reticulate::virtualenv_install(
    envname  = venv_name,
    packages = c("pip","wheel","setuptools","numpy","biopython","tensorflow"),
    ignore_installed = FALSE
  )
}
venv_python <- file.path(venv_path, "bin", "python")
Sys.setenv(RETICULATE_PYTHON = venv_python)
message("RETICULATE_PYTHON set to ", Sys.getenv("RETICULATE_PYTHON"))
library(reticulate)

use_python(venv_python, required = TRUE)
py_available(initialize = TRUE)

#------------------------------------------------------------------------------#
# 0. parse_mutation 
#------------------------------------------------------------------------------#
parse_mutation <- function(mutation, wildtype_seq) {
  info <- list(
    Locus=NA_integer_, Mutation_Type=NA_character_, New_Base=NA_character_,
    Original_Base=NA_character_, Deleted_Bases=NA_character_,
    Duplicated_Bases=NA_character_, Inserted_Bases=NA_character_,
    Mutated_Sequence=NA_character_
  )
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

#------------------------------------------------------------------------------#
# 1. translation
#------------------------------------------------------------------------------#
translate_prot <- function(dna_seq) {
  if (is.na(dna_seq) || nchar(dna_seq)<3) return(NA_character_)
  prot_full <- suppressWarnings(as.character(translate(DNAString(dna_seq), if.fuzzy.codon="X")))
  if (grepl("\\*", prot_full)) substr(prot_full, 1, regexpr("\\*", prot_full)[1]-1) else prot_full
}

#------------------------------------------------------------------------------#
# 2. Refinement helpers
#------------------------------------------------------------------------------#
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
  if (net %% 3 == 0) "In-Frame indel" else "Frameshifting indel"
}

#------------------------------------------------------------------------------#
# UI 
#------------------------------------------------------------------------------#
ui <- bslib::page_navbar(
  useShinyjs(),
  id = "topnav",
  title = tags$div(
    class = "brand",
    style = "display:flex;flex-direction:column;line-height:1.05;align-items:flex-start;padding-right:32px;",
    tags$span(class="brand-title",   "MfAP 🗺"),
    tags$span(class="brand-subtitle","Molecular feature Association Pipeline")
  ),
  theme = bslib::bs_theme(version = 5, bootswatch = "darkly"),
  window_title = "MfAP — Molecular feature Association Pipeline",
  header = tags$head(
    tags$style(HTML("
      .page-container{max-width:1180px;margin:0 auto;padding:0 14px;}
      .navbar .container-fluid{max-width:1180px;margin:0 auto;padding:0 14px;}
      .lit-item p{ margin:0 0 6px 0; }
      .navbar{background:#1d1f21 !important;border-bottom:1px solid rgba(255,255,255,.08)!important;box-shadow:none !important;}
      .navbar-brand{margin-right:auto;}
      .brand-title{
        font-weight: 900;
        font-size: 1.85rem;
        letter-spacing: .2px;
        color:#fff !important;
      }
      .brand-subtitle{
        font-size: 1.08rem;
        color:#e5e7eb !important;
        margin-top: 2px;
      }
      .navbar-nav{ margin-left:auto; gap: 30px; }
      /* responsive tweak */
      @media (max-width: 768px){
        .brand-title{   font-size: 1.55rem; }
        .brand-subtitle{font-size: 0.98rem; }
        .navbar-nav{ gap: 18px; }
      }
      .navbar-nav .nav-link{padding:.5rem .6rem;color:#e9ecef !important;border-radius:0 !important;border-bottom:2px solid transparent;}
      .navbar-nav .nav-link:hover{color:#4dabf7 !important;background:transparent !important;border-bottom-color:#4dabf7;}
      .navbar-nav .nav-link.active{color:#4dabf7 !important;background:transparent !important;border-bottom-color:#4dabf7;}

      .card{border:1px solid rgba(255,255,255,0.08)!important;border-radius:10px!important;background:#202326;}
      .card-header{background:#2a2e32!important;border-bottom:1px solid rgba(255,255,255,0.06)!important;color:#fff;padding:10px 14px!important;font-weight:600;}
      .card-body{padding:14px!important;}
      .card-footer{background:#202326;border-top:1px solid rgba(255,255,255,.06);padding:12px 14px;}

      .sidebar{position:sticky;top:16px;max-height:calc(100vh - 32px);overflow:auto;}
      .action-dock{
  position: sticky;
  bottom: 0;
  background: rgba(32,35,38,0.35);      /* more see-through */
  backdrop-filter: blur(8px) saturate(120%);
  -webkit-backdrop-filter: blur(8px) saturate(120%); /* Safari */
  border: 1px solid rgba(255,255,255,0.12);
  box-shadow: 0 6px 24px rgba(0,0,0,0.25);
  padding: 10px 12px;
  border-radius: 12px;
  margin-top: 10px;
  z-index: 3;
}
      .action-row{display:flex;gap:8px;}
      .action-dock .btn{width:100%;}

      details{background:#262a2e;border:1px solid rgba(255,255,255,0.08);border-radius:8px;padding:8px 10px;}
      details>summary{cursor:pointer;font-weight:600;color:#e9ecef;list-style:none;}
      details>summary::-webkit-details-marker{display:none;}
      details[open]{background:#2d3236;}
      .muted{color:#adb5bd;}
      .download-row{display:flex;justify-content:flex-end;margin:0;}

      .about-mfap .card-body{padding:12px 14px !important;}
      .about-mfap .lead-line{font-size:1.15rem;color:#cfd4da;margin-bottom:10px;}
      .about-mfap .notes{background:#23272b;border:1px solid rgba(255,255,255,.06);
                         border-radius:10px;padding:12px 14px;margin-top:8px;}
      .about-mfap .notes p{margin:0 0 6px 0;font-size:.98rem;line-height:1.35;}
      .about-mfap .cite-row{display:flex;align-items:center;gap:.5rem;margin-top:10px;}
      .about-mfap .bib-link{cursor:pointer;text-decoration:underline;}

      
      #bibtex_block{white-space:pre-wrap;font-size:.86rem;background:#1f1f1f;
      padding:10px;border-radius:8px;border:1px solid #2a2a2a;}

      footer, .site-footer { clear: both; }

      .site-footer{
        background:transparent;
        color:#fff;
        border-top:1px solid rgba(255,255,255,0.10);
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
      .site-footer__center{ flex:1 1 auto; text-align:center; color:#fff; }
      .site-footer__link{ color:#fff; text-decoration:underline; }

      @media (max-width:768px){
        .site-footer__logo{ height:42px; }
        .site-footer__center{ font-size:0.88rem; }
      }
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
  ),
  
  # ————— Analyze
  bslib::nav_panel(
    "Analysis",
    div(class="page-container",
        sidebarLayout(
          sidebarPanel(
            width = 3, class="sidebar",
            style = "max-height: calc(100vh - 220px); overflow-y: auto; overflow-x: hidden;",
            bslib::card(
              bslib::card_header("Upload"),
              bslib::card_body(
                uiOutput("csv_ui"),
                tags$details(
                  tags$summary("File requirements ▿"),
                  tags$div(class="muted", style="font-size:0.9em;",
                           HTML(paste(
                             "<b>Your table needs at least two columns with these names:</b>",
                             "<ul>",
                             "<li>'DNA_variant' in HGVS-annotation, e.g. 'c.123delA'",
                             "<li>'patient_ID' as the individual identifier</li>",
                             "</ul>", sep = ""
                           )))
                ),
                div(style="height:8px;"),
                uiOutput("fasta_ui"),
                tags$details(
                  tags$summary("Obtaining a FASTA ▿"),
                  tags$div(class="muted", style="font-size:0.9em;",
                           HTML(paste(
                             "Example workflow:",
                             "<ul style='margin-bottom:0;'>",
                             "<li>Find your gene on <i>ensembl.org</i></li>",
                             "<li>Click <b>Download sequence</b></li>",
                             "<li>Included sequences: select <b>all</b></li>",
                             "<li>Format: <b>FASTA</b></li>",
                             "<li>Flanks: set to <b>0</b></li>",
                             "<li>Click <b>Download</b></li>",
                             "</ul>", sep = ""
                           )))
                ),
                div(style="height:8px;"),
                checkboxInput("use_titer", "Include alternative translation initiation site prediction (TITER, developed by Zhang et al. 2017)", FALSE),
                conditionalPanel(
                  condition = "input.use_titer == true",
                  uiOutput("fasta_flank_ui"),
                  tags$details(
                    tags$summary("Obtaining a FASTA with flanks ▿"),
                    tags$div(class="muted", style="font-size:0.9em;",
                             HTML("Follow the same workflow as for the normal FASTA file. Simply set Flanks to 100 instead of 0."))
                  ),
                  tags$details(
                    tags$summary("Why do you need other FASTA files with flanks for TITER? ▿"),
                    tags$div(class="muted", style="font-size:0.9em;",
                             HTML("TITER evaluates sequence context within a 203 bp window around the central bases. To evaluate potential start sites close to the 5' end, we need a CDS file with ±100 bp flanks (part of 5' and 3' UTR)."))
                  )
                ),
                div(style="height:8px;"),
                actionButton("clear_uploads", "Clear uploads", class="btn btn-secondary")
              )
            ),
            bslib::card(
              bslib::card_header("Domains"),
              bslib::card_body(
                fluidRow(column(6, textInput("dom_name","Name")),
                         column(6, numericInput("dom_start","Start (AA)",1,min=1))),
                fluidRow(column(6, numericInput("dom_end","End (AA)",1,min=1)),
                         column(6, checkboxInput("dom_to_end","to end", FALSE))),
                fluidRow(column(6, actionButton("add_dom","Add", class="btn btn-info")),
                         column(6, actionButton("clr_dom","Clear", class="btn btn-secondary"))),
                tableOutput("dom_tbl") %>% tagAppendAttributes(style="font-size:0.9em; margin-top:8px;")
              )
            ),
            bslib::card(
              bslib::card_header("Motifs"),
              bslib::card_body(
                fluidRow(column(6, textInput("mot_name","Name")),
                         column(6, textInput("mot_pattern","Pattern (AA sequence)"))),
                fluidRow(column(6, actionButton("add_mot","Add", class="btn btn-info")),
                         column(6, actionButton("clr_mot","Clear", class="btn btn-secondary"))),
                tableOutput("mot_tbl") %>% tagAppendAttributes(style="font-size:0.9em; margin-top:8px;")
              )
            ),
            bslib::card(
              bslib::card_header("Demo"),
              bslib::card_body(
                tags$p(class="muted", style="font-size:0.9em;",
                       "Load a fictional cohort with example variants, reference sequences, and example domains/motifs."),
                fluidRow(
                  column(6, actionButton("load_example_data", "Load demo cohort", class="btn btn-dark")),
                  column(6, actionButton("clear_example_data", "Clear demo", class="btn btn-secondary"))
                )
              )
            ),
            div(class="action-dock",
                div(class="action-row",
                    actionButton("run", "Run analysis", icon=icon("play"), class="btn btn-success"),
                    actionButton("cancel", "Cancel", class="btn btn-danger")
                )
            )
          ),
          mainPanel(
            width = 9,
            style = "max-height: calc(100vh - 220px); overflow-y: auto; overflow-x: hidden;",
            bslib::card(
              class = "about-mfap",
              bslib::card_header("About MfAP"),
              bslib::card_body(
                tags$p(
                  class = "lead-line",
                  "Use MfAP to analyze protein-level consequences of cDNA variants in single-exon genes. ",
                  a("Background →", href = "#",
                    onclick = "
            var link = document.querySelector('.navbar .nav-link[data-value=\"Background\"]');
            if (link) { link.click(); }
            return false;")
                ),
                tags$div(
                  class = "notes",
                  tags$p(
                    strong("Important notes:"),
                    "You are solely responsible for compliance with privacy and data-protection laws and policies. ",
                    "This app is solely designed for research purposes."
                  ),
                  tags$div(
                    class = "cite-row",
                    tags$span(strong("Please cite:"), "Schubert et al. (2025) ·"),
                    actionLink("show_bib", label = "BibTeX", class = "bib-link")
                  ),
                  tags$div(
                    class = "cite-row",
                    tags$span("If using TITER, please additionally cite Zhang et al. (2017)"),
                    "·",
                    actionLink("show_titer_bibtex", label = "BibTeX", class = "bib-link")
                  )
                )
              )
            ),
            bslib::card(
              bslib::card_header("Results"),
              bslib::card_body(DTOutput("res_tbl")),
              bslib::card_footer(div(class="download-row", uiOutput("download_btn")))
            )
          )
        )
    )
  ),
  
  # ————— Background
  bslib::nav_panel(
    "Background",
    div(class="page-container",
        fluidRow(
          column(12,
                 bslib::card(
                   bslib::card_header("Overview"),
                   bslib::card_body(
                     p("MfAP is a framework developed by Tim Schubert and colleagues at the Institute of Human Genetics of Heidelberg University (Germany), which allows users to calculate and predict DNA and protein 
                       level consequences of genetic variation in single-exon genes. 
                       MfAP requires minimal user input: a list of variants in HGVS cDNA 
                       notation and the reference DNA sequence of the gene of interest.")
                   )
                 ),
                 bslib::card(
                   bslib::card_header("Rationale"),
                   bslib::card_body(
                     p("Protein-truncating variants (PTVs) canonically lead to nonsense-mediated mRNA decay. 
                       This process necessitates an exon-exon junction, which is lacking in single-exon genes. 
                       Therefore, PTVs in single-exon genes can result in a truncated protein that may have dominant-negative, gain-of-function, or neomorphic effects.
                       In our accompanying manuscript, we demonstrate that MfAP-predicted protein-level attributes of such
                       truncated proteins provide a stronger causal link between genotype and phenotypic severity."),
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
                     p("MfAP outputs calculations and predictions of the effects of variants on the DNA and protein level, including predictions of non-canonical translation initiation sites (TISs) and resulting proteins. For more information and inspiration for downstream analyses, please read the accompanying paper."),
                     tags$figure(
                       tags$img(
                         src = "https://lh3.googleusercontent.com/d/1h-h2yfvsKnxscA3ko8VyeFXkPaa9iAkE",
                         alt = "MfAP workflow schematic",
                         style = "display:block;margin:10px auto;max-width:60%;height:auto;border-radius:10px;border:1px solid rgba(255,255,255,.08);",
                         loading = "lazy"
                       ),
                       tags$figcaption(
                         "Figure: MfAP workflow. Figure created with Biorender.com",
                         style="text-align:center;color:#cfd4da;font-size:.9em;margin-top:6px;"
                       )
                     )
                   )
                 ),
                 bslib::card(
                   bslib::card_header("Literature"),
                   bslib::card_body(
                     tags$p(
                       style = "margin:0 0 8px 0;",
                       HTML("Schubert T, Tietzel A, Pottayil H, Caro P, Gilmore RB, Franke F, Althammer F, Schaaf CP. 2025. "),
                       em("A blueprint for protein-centric genotype-phenotype investigations in single-exon disease genes applied to MAGEL2 and Schaaf-Yang syndrome."),
                       HTML(" <i>Am J Hum Genet</i> <b>VOLUME</b>(ISSUE):PAGES · "),
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
                 )
          )
        )
    )
  ),
  
  # --- Footer spacer & footer ---
  tags$div(style="height:100px;"),
  tags$footer(
    class = "site-footer",
    style = "background:transparent; color:#fff;",
    tags$div(
      class = "site-footer__inner",
      tags$img(src="https://lh3.googleusercontent.com/d/12vCwvlmN8xrKnlU1GsDGJeDBgq_yl0sM", class="site-footer__logo"),
      tags$div(
        class="site-footer__center",
        tags$p(HTML("&copy; 2025 Tim Schubert"), style="margin:0;"),
        tags$p(
          tags$a(href="https://www.apache.org/licenses/LICENSE-2.0.html",
                 target="_blank", "Apache License 2.0", class="site-footer__link"),
          style="margin:0;"
        )
      ),
      tags$img(src="https://lh3.googleusercontent.com/d/1BPl641wAs67xCTAKEEETVO0zjzOYUrda", class="site-footer__logo")
    )
  )
)

#------------------------------------------------------------------------------#
# Server
#------------------------------------------------------------------------------#
server <- function(input, output, session) {
  example_csv <- reactiveVal(NULL)
  example_fasta_seq <- reactiveVal(NULL)
  example_fasta_flank_seq <- reactiveVal(NULL)
  example_domains <- data.frame(Name=c("DomainA","DomainB"), Start=c(5,50), End=c(20,100), stringsAsFactors=FALSE)
  example_motifs  <- data.frame(Name=c("MotifX","MotifY"),  Pattern=c("ACDE","WXYZ"), stringsAsFactors=FALSE)
  domains <- reactiveVal(data.frame(Name=character(),Start=integer(),End=numeric(),stringsAsFactors=FALSE))
  motifs  <- reactiveVal(data.frame(Name=character(),Pattern=character(),stringsAsFactors=FALSE))
  wt_aa_len <- reactiveVal(NULL)
  
  result_data <- reactiveVal(NULL)
  cancel_requested <- reactiveVal(FALSE)
  observeEvent(input$cancel, {
    cancel_requested(TRUE); result_data(NULL); shinyjs::hide("res_tbl")
  })
  
  observe({
    has_csv   <- !is.null(input$csv)   || !is.null(example_csv())
    has_fasta <- !is.null(input$fasta) || !is.null(example_fasta_seq())
    has_flank <- if (isTRUE(input$use_titer)) {!is.null(input$fasta_flank) || !is.null(example_fasta_flank_seq())} else TRUE
    shinyjs::toggleState("run", has_csv && has_fasta && has_flank)
  })
  
  observeEvent(input$clear_uploads, {
    shinyjs::reset("csv"); shinyjs::reset("fasta"); shinyjs::reset("fasta_flank")
    example_csv(NULL); example_fasta_seq(NULL); example_fasta_flank_seq(NULL)
    result_data(NULL); shinyjs::disable("run"); shinyjs::disable("cancel")
  })
  
  observeEvent(input$load_example_data, {
    example_csv(read.csv("www/example_variants.csv", stringsAsFactors=FALSE))
    lines <- readLines("www/example_reference.fasta", warn=FALSE)
    example_fasta_seq(toupper(paste(lines[!startsWith(lines, ">")], collapse="")))
    lines_flank <- readLines("www/example_reference_flank.fasta", warn=FALSE)
    example_fasta_flank_seq(toupper(paste(lines_flank[!startsWith(lines_flank, ">")], collapse="")))
    domains(example_domains); motifs(example_motifs)
    showNotification("Example data loaded.", type="message")
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
    req(input$dom_name, input$dom_start)
    
    max_aa <- wt_aa_len()
    
    # Decide end value robustly
    if (isTRUE(input$dom_to_end)) {
      if (is.null(max_aa) || is.na(max_aa) || !is.finite(max_aa)) {
        showNotification("Protein length unknown. Load FASTA or run analysis first (so 'to end' can resolve).", type = "error")
        return()
      }
      end_val <- as.integer(max_aa)
    } else {
      if (is.null(input$dom_end) || is.na(input$dom_end)) {
        showNotification("Please provide a domain end position.", type = "error")
        return()
      }
      end_val <- as.integer(input$dom_end)
    }
    
    start_val <- as.integer(input$dom_start)
    
    # Basic sanity checks
    if (!is.null(max_aa) && end_val > max_aa) {
      showNotification(paste0("Domain end (", end_val, ") exceeds protein length (", max_aa, ")."), type = "error")
      return()
    }
    if (!is.null(max_aa) && (start_val < 1 || end_val < 1 || start_val > end_val)) {
      showNotification("Domain coordinates must be within [1, end] and Start ≤ End.", type = "error")
      return()
    }
    
    df <- domains()
    if (nrow(df) > 0 && any(df$Name == input$dom_name & df$Start == start_val & df$End == end_val)) {
      showNotification("Exact duplicate domain exists.", type = "error")
      return()
    }
    
    # Ensure correct column types
    new_row <- data.frame(
      Name  = as.character(input$dom_name),
      Start = start_val,
      End   = end_val,
      stringsAsFactors = FALSE
    )
    
    domains(rbind(df, new_row))
  })
  observeEvent(input$clr_dom, domains(data.frame(Name=character(),Start=integer(),End=numeric(),stringsAsFactors=FALSE)))
  output$dom_tbl <- renderTable(domains(), striped=TRUE, hover=TRUE)
  
  observeEvent(input$add_mot, {
    req(input$mot_name, input$mot_pattern)
    if (!str_detect(input$mot_pattern, "^[ACDEFGHIKLMNPQRSTVWYBZX]+$")) {showNotification("Pattern must use only IUPAC amino acid letters A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,B,Z,X", type="error"); return()}
    df <- motifs()
    if (any(df$Name==input$mot_name & df$Pattern==input$mot_pattern)) {showNotification("Exact duplicate motif exists.", type="error"); return()}
    motifs(rbind(df, data.frame(Name=input$mot_name, Pattern=input$mot_pattern, stringsAsFactors=FALSE)))
  })
  observeEvent(input$clr_mot, motifs(data.frame(Name=character(),Pattern=character(),stringsAsFactors=FALSE)))
  output$mot_tbl <- renderTable(motifs(), striped=TRUE, hover=TRUE)
  
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
  
  observeEvent(input$run, {
    cancel_requested(FALSE); shinyjs::disable("run"); shinyjs::enable("cancel"); shinyjs::hide("res_tbl")
    req(!is.null(input$csv) || !is.null(example_csv()))
    vd <- if (!is.null(input$csv)) read.csv(input$csv$datapath, stringsAsFactors=FALSE) else example_csv()
    if (nrow(vd)==0) {showModal(modalDialog(title="Empty CSV","Your variant list file has zero rows.", easyClose=TRUE, footer=modalButton("OK"))); result_data(NULL); return(NULL)}
    if (!"DNA_variant" %in% colnames(vd)) {showModal(modalDialog(title="Missing column","Your CSV must contain a column named `DNA_variant`.", easyClose=TRUE, footer=modalButton("OK"))); result_data(NULL); return(NULL)}
    if (isTRUE(input$use_titer)) {
      if (!"patient_ID" %in% colnames(vd)) {showModal(modalDialog(title="Missing column","TITER analysis requires a `patient_ID` column.", easyClose=TRUE, footer=modalButton("OK"))); result_data(NULL); return(NULL)}
      if (is.null(input$fasta_flank) && is.null(example_fasta_flank_seq())) {showModal(modalDialog(title="Missing FASTA","Upload the FASTA with ±100 bp flanks for TITER.", easyClose=TRUE, footer=modalButton("OK"))); result_data(NULL); return(NULL)}
    }
    
    lines <- if (!is.null(input$fasta)) readLines(input$fasta$datapath, warn=FALSE) else c(">example", example_fasta_seq())
    if (!any(startsWith(lines, ">"))) {showModal(modalDialog(title="Invalid FASTA","FASTA must have at least one header line beginning with `>`.", easyClose=TRUE, footer=modalButton("OK"))); result_data(NULL); return(NULL)}
    seq_nt <- toupper(paste(lines[!startsWith(lines, ">")], collapse=""))
    if (nchar(seq_nt) %% 3 != 0) {showModal(modalDialog(title="Sequence not multiple of 3","Reference length is not divisible by 3.", easyClose=TRUE, footer=modalButton("OK"))); result_data(NULL); return(NULL)}
    aa_len <- nchar(Biostrings::translate(Biostrings::DNAString(seq_nt), if.fuzzy.codon="X")); wt_aa_len(aa_len)
    
    withProgress(message="Running analysis...", value=0, {
      incProgress(0.1, detail="Parsing mutations")
      if (cancel_requested()) {showNotification("Analysis cancelled.", type="warning"); shinyjs::disable("cancel"); result_data(NULL); return(NULL)}
      parsed <- lapply(vd$DNA_variant, parse_mutation, wildtype_seq=seq_nt)
      df <- vd %>% dplyr::mutate(
        Locus=sapply(parsed, `[[`,"Locus"),
        Mutation_Type=sapply(parsed, `[[`,"Mutation_Type"),
        Deleted_Bases=sapply(parsed, `[[`,"Deleted_Bases"),
        Duplicated_Bases=sapply(parsed, `[[`,"Duplicated_Bases"),
        Inserted_Bases=sapply(parsed, `[[`,"Inserted_Bases"),
        Mutated_Sequence=sapply(parsed, `[[`,"Mutated_Sequence")
      )
      
      incProgress(0.2, detail="Translating proteins")
      if (cancel_requested()) {showNotification("Analysis cancelled.", type="warning"); shinyjs::disable("cancel"); result_data(NULL); return(NULL)}
      df <- df %>% dplyr::mutate(Mutated_Protein=sapply(Mutated_Sequence, translate_prot), Protein_Length_aa=nchar(Mutated_Protein))
      
      incProgress(0.2, detail="Refining mutation types")
      if (cancel_requested()) {showNotification("Analysis cancelled.", type="warning"); shinyjs::disable("cancel"); result_data(NULL); return(NULL)}
      df <- df %>% dplyr::rowwise() %>% dplyr::mutate(Mutation_Type = dplyr::case_when(
        Mutation_Type=="Point" ~ refine_point(Locus, Mutated_Sequence, seq_nt, GENETIC_CODE),
        Mutation_Type %in% c("Deletion","Duplication","Insertion","Indel") ~ refine_indel(Mutation_Type, Deleted_Bases, Duplicated_Bases, Inserted_Bases),
        TRUE ~ Mutation_Type
      )) %>% dplyr::ungroup()
      
      incProgress(0.2, detail="Annotating domains & motifs")
      if (cancel_requested()) {showNotification("Analysis cancelled.", type="warning"); shinyjs::disable("cancel"); result_data(NULL); return(NULL)}
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
      
      # ----------------------- Non-canonical TIS via TITER -----------------------
      if (isTRUE(input$use_titer)) {
        incProgress(0.05, detail = "Running non-canonical TIS analysis (TITER)")
        if (cancel_requested()) {
          showNotification("Analysis cancelled.", type = "warning")
          shinyjs::disable("cancel"); result_data(NULL); return(NULL)
        }
        
        # Use the already-selected interpreter (RETICULATE_PYTHON) — do NOT switch envs here
        # reticulate::use_virtualenv("titer-venv", required = TRUE)  # <-- REMOVE
        
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
        if (!reticulate::py_available(initialize = TRUE)) {
          showNotification("Python not available; skipping TITER analysis.", type = "error")
        } else {
          # Optional sanity check for required modules (helpful on Connect)
          needed <- c("numpy", "tensorflow", "Bio")
          missing <- needed[!vapply(needed, reticulate::py_module_available, logical(1))]
          if (length(missing)) {
            stop("Missing Python modules for TITER: ", paste(missing, collapse = ", "),
                 ". Make sure your startup venv install step ran on Connect.")
          }
          
          reticulate::source_python(script_path)
          
          titer_summary <- read.csv(file.path(temp_titer_dir, "data",
                                              "summary_patients_most_likely_additional_TIS.csv"),
                                    stringsAsFactors = FALSE) %>% dplyr::select(-DNA_variant)
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
        }
      }
      # --------------------------------------------------------------------------
      
      incProgress(0.1, detail="Calculating frameshift metrics")
      if (cancel_requested()) {showNotification("Analysis cancelled.", type="warning"); shinyjs::disable("cancel"); result_data(NULL); return(NULL)}
      df <- df %>% dplyr::mutate(
        Unaltered_Length_aa  = dplyr::if_else(Mutation_Type=='Frameshifting indel', (floor((Locus-1)/3)+1)-1, NA_integer_),
        Frameshift_Length_aa = dplyr::if_else(Mutation_Type=='Frameshifting indel', Protein_Length_aa - ((floor((Locus-1)/3)+1)-1), NA_integer_)
      )
      
      colnames(df) <- gsub("_"," ", colnames(df))
      df <- df %>% dplyr::select(-c(`Mutated Sequence`,`Mutated Protein`), dplyr::everything(), `Mutated Sequence`,`Mutated Protein`)
      result_data(df); shinyjs::disable("cancel"); shinyjs::enable("run"); shinyjs::show("res_tbl")
    })
  })
  
  output$res_tbl <- renderDT({
    req(result_data()); datatable(result_data(), rownames=FALSE, filter='top',
                                  options=list(pageLength=10, scrollX=TRUE))
  })
  
  output$dl <- downloadHandler(
    filename = function() paste0('MAP_results_', Sys.Date(), '.csv'),
    content  = function(f) write.csv(result_data(), f, row.names=FALSE)
  )
  
  output$download_btn <- renderUI({
    df <- result_data()
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
    downloadButton("dl", "Download CSV", class="btn btn-success")
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
