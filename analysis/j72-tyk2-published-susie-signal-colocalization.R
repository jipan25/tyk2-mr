suppressPackageStartupMessages({library(data.table); library(coloc); library(ggplot2)})
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args)) args[[1]] else "data/upstream/clinical"
root <- file.path(data_dir, "molecular_qtl_followup")
out_dir <- file.path(root, "results_signal_coloc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
meta <- fread(file.path(root, "resource_metadata", "selected_finemapping_downloads.tsv"))
meta <- meta[dataset_id %in% c("QTD000024", "QTD000034", "QTD000329", "QTD000359")]
stopifnot(nrow(meta) == 4L)
required_bf <- file.path(root, "derived_finemap", paste0(meta$dataset_id, ".TYK2_KEAP1.lbf_variable.tsv.gz"))
if (any(!file.exists(required_bf))) stop("A prespecified molecular BF input is still missing: ", paste(meta$dataset_id[!file.exists(required_bf)], collapse = ", "))
primary_event <- "ENSG00000105397.grp_2.contained.ENST00000525976"
normalize_variant <- function(x) gsub("_", ":", sub("^chr", "", x))
logsumexp <- function(x) {m <- max(x); m + log(sum(exp(x - m)))}
posterior <- function(x) exp(x - logsumexp(x))

gwas <- list()
gwas_cs <- list()
gwas_audit <- list()
for (endpoint in c("M13_RHEUMA", "RHEUMA_SEROPOS_STRICT")) {
  prefix <- paste0("finngen_R13_", endpoint, ".SUSIE")
  x <- fread(file.path(root, "derived_finemap", paste0(prefix, ".TYK2_region.tsv.gz")))
  x <- x[region == "chr19:8852442-11852442"]
  x[, key := normalize_variant(v)]
  stopifnot(!anyDuplicated(x$key), all(x$trait == endpoint), nrow(x) > 1000)
  cs <- fread(file.path(root, "raw_finemap", paste0(prefix, ".cred.bgz")))
  cs <- cs[region == "chr19:8852442-11852442" & !low_purity & cs_min_r2 >= 0.5 & cs_log10bf >= 0.9]
  for (signal in cs$cs) {
    lbf <- x[[paste0("lbf_variable", signal)]]
    alpha <- x[[paste0("alpha", signal)]]
    reconstructed_alpha <- posterior(lbf)
    stopifnot(max(abs(reconstructed_alpha - alpha)) < 1e-8, abs(sum(alpha) - 1) < 1e-8)
    gwas_audit[[length(gwas_audit) + 1L]] <- data.table(
      outcome_id = endpoint, signal = signal, n_variants = nrow(x),
      cs_min_r2 = cs[cs == signal, cs_min_r2],
      max_alpha_reconstruction_difference = max(abs(reconstructed_alpha - alpha)),
      lead_variant = x$key[which.max(alpha)], lead_alpha = max(alpha)
    )
  }
  gwas[[endpoint]] <- x
  gwas_cs[[endpoint]] <- cs
}

coloc_rows <- list()
overlap_rows <- list()
qtl_audit <- list()
for (i in seq_len(nrow(meta))) {
  m <- meta[i]
  bf_path <- file.path(root, "derived_finemap", paste0(m$dataset_id, ".TYK2_KEAP1.lbf_variable.tsv.gz"))
  qtl <- fread(bf_path)
  stopifnot(all(c("molecular_trait_id", "variant", "lbf_variable1") %in% names(qtl)))
  qtl[, key := normalize_variant(variant)]
  cs <- fread(file.path(root, "raw_finemap", paste0(m$dataset_id, ".TYK2_KEAP1.credible_sets.tsv")))
  cs[, signal_index := as.integer(sub(".*_L", "", cs_id))]
  signals <- unique(cs[, .(molecular_trait_id, cs_id, signal_index, cs_min_r2, cs_size)])
  signals[, qualified_primary := cs_min_r2 >= 0.5]
  signals[, `:=`(dataset_id = m$dataset_id, context = m$sample_group,
                 prespecified_event = molecular_trait_id == primary_event)]
  qtl_audit[[i]] <- signals
  signals <- signals[qualified_primary == TRUE]
  for (j in seq_len(nrow(signals))) {
    s <- signals[j]
    q <- qtl[molecular_trait_id == s$molecular_trait_id]
    stopifnot(nrow(q) > 1000, !anyDuplicated(q$key))
    bf_qtl <- q[[paste0("lbf_variable", s$signal_index)]]
    alpha_qtl <- posterior(bf_qtl)
    for (endpoint in names(gwas)) {
      g <- gwas[[endpoint]]
      common <- intersect(g$key, q$key)
      g_idx <- match(common, g$key)
      q_idx <- match(common, q$key)
      for (gs in gwas_cs[[endpoint]]$cs) {
        bf_gwas <- g[[paste0("lbf_variable", gs)]]
        alpha_gwas <- posterior(bf_gwas)
        mass_gwas <- sum(alpha_gwas[g_idx])
        mass_qtl <- sum(alpha_qtl[q_idx])
        clpp <- sum(alpha_gwas[g_idx] * alpha_qtl[q_idx])
        overlap_rows[[length(overlap_rows) + 1L]] <- data.table(
          outcome_id = endpoint, dataset_id = m$dataset_id, study = m$study_label,
          context = m$sample_group, molecular_trait_id = s$molecular_trait_id,
          qtl_signal = s$signal_index, gwas_signal = gs,
          n_qtl_variants = nrow(q), n_gwas_variants = nrow(g), n_shared_variants = length(common),
          gwas_posterior_mass_shared = mass_gwas, qtl_posterior_mass_shared = mass_qtl,
          CLPP = clpp, qualified_overlap = mass_gwas >= 0.9 & mass_qtl >= 0.9
        )
        if (mass_gwas < 0.9 || mass_qtl < 0.9) next
        b1 <- matrix(bf_gwas, nrow = 1, dimnames = list(paste0("RA_L", gs), g$key))
        b2 <- matrix(bf_qtl, nrow = 1, dimnames = list(s$cs_id, q$key))
        for (prior12 in c(1e-6, 5e-6, 1e-5)) {
          fit <- coloc.bf_bf(b1, b2, p1 = 1e-4, p2 = 1e-4, p12 = prior12,
                            overlap.min = 0.9, trim_by_posterior = TRUE)
          result <- as.data.table(fit$summary)
          if (!nrow(result)) stop("A qualified signal pair unexpectedly produced no coloc summary")
          result[, `:=`(
            outcome_id = endpoint, dataset_id = m$dataset_id, study = m$study_label,
            context = m$sample_group, quant_method = m$quant_method,
            molecular_trait_id = s$molecular_trait_id,
            qtl_signal = s$signal_index, gwas_signal = gs,
            qtl_cs_min_r2 = s$cs_min_r2,
            gwas_cs_min_r2 = gwas_cs[[endpoint]][cs == gs, cs_min_r2],
            p1 = 1e-4, p2 = 1e-4, p12 = prior12, CLPP = clpp,
            n_shared_variants = length(common),
            gwas_posterior_mass_shared = mass_gwas, qtl_posterior_mass_shared = mass_qtl,
            primary_event = s$molecular_trait_id == primary_event,
            source_method = "Published in-sample SuSiE single-effect BF; coloc::coloc.bf_bf"
          )]
          coloc_rows[[length(coloc_rows) + 1L]] <- result
        }
      }
    }
  }
}

if (!length(coloc_rows)) stop("No qualified colocalization comparisons; input or overlap gate remains unresolved")
results <- rbindlist(coloc_rows, fill = TRUE)
overlap <- rbindlist(overlap_rows)
fwrite(results, file.path(out_dir, "TYK2_Published_SuSiE_Signal_Colocalization_Prior_Sensitivity.csv"))
fwrite(overlap, file.path(out_dir, "TYK2_Published_SuSiE_Signal_Overlap_Audit.csv"))
fwrite(rbindlist(gwas_audit), file.path(out_dir, "FinnGen_SuSiE_BF_Alpha_Reconstruction_Audit.csv"))
fwrite(rbindlist(qtl_audit, fill = TRUE), file.path(out_dir, "QTL_SuSiE_Credible_Set_Quality_Audit.csv"))
primary <- results[p12 == 1e-5 & primary_event == TRUE]
fwrite(primary, file.path(out_dir, "TYK2_Main_Event_Signal_Colocalization_Primary.csv"))
print(primary)
writeLines(c("# Published-SuSiE single-effect colocalization", "",
  "Published variant-level single-effect log Bayes factors were used directly. No external LD matrix or positive square root of r2 was constructed, and no new fine-mapping was claimed.",
  "Primary gates: credible-set minimum r2>=0.5; posterior mass >=0.9 in both datasets within the shared variant set; p1=p2=1e-4, p12=1e-5. Prior sensitivity: p12=1e-6 and 5e-6.", "",
  paste(capture.output(print(primary)), collapse = "\n"), "",
  "## Interpretation",
  "The prespecified TYK2 txrevise phenotype measures relative event usage. Gene-level QTL coverage and lack of qualified total-expression credible sets are reported separately.",
  "GTEx tissues and BLUEPRINT cell types are correlated contexts within their respective studies; they are not four independent cohorts.",
  "Shared genetic signals do not establish that event usage mediates RA risk or exclude a coding-protein mechanism at rs34536443.",
  "This is a new colocalization calculation using previously published in-sample fine-mapping outputs, not a new individual-level or externally-LD-refitted SuSiE analysis."),
  file.path(out_dir, "TYK2_Published_SuSiE_Signal_Colocalization_Summary.md"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Signal colocalization results: ", out_dir)
