suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(MRPRESSO)
  library(TwoSampleMR)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else "data/upstream/clinical"
project_dir <- if (length(args) >= 2) args[[2]] else "data/upstream/project"
raw_dir <- file.path(data_dir, "raw_complete_summary_stats")
out_dir <- file.path(data_dir, "results_complete_summary_stats")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

eqtl <- fread(file.path(data_dir, "derived", "eqtlgen_TYK2_all_cis.tsv.gz"))
stopifnot(nrow(eqtl) == 6995L, all(eqtl$GeneSymbol == "TYK2"), !anyDuplicated(eqtl$SNP))
af_cache <- file.path(data_dir, "derived", "eqtlgen_TYK2_allele_frequency.tsv.gz")
af_file <- if (length(args) >= 3) args[[3]] else af_cache
if (!file.exists(af_file)) stop("Supply the eQTLGen cis-EQTL-AF.txt path as argument 3 for the first run")
af <- fread(af_file)
af <- af[SNP %in% eqtl$SNP]
if (normalizePath(af_file) != normalizePath(af_cache, mustWork = FALSE)) {
  fwrite(af, af_cache, sep = "\t")
}
eqtl <- merge(eqtl, af[SNP %in% eqtl$SNP], by = "SNP")
eqtl[, frequency_alleles_valid :=
  (AssessedAllele == AlleleA & OtherAllele == AlleleB) |
    (AssessedAllele == AlleleB & OtherAllele == AlleleA)]
if (!all(eqtl$frequency_alleles_valid)) stop("eQTL allele-frequency reference mismatch")
eqtl[, assessed_af := fifelse(AssessedAllele == AlleleB, as.numeric(AlleleB_all), 1 - as.numeric(AlleleB_all))]
eqtl[, maf := pmin(assessed_af, 1 - assessed_af)]
eqtl[, se := 1 / sqrt(2 * maf * (1 - maf) * (NrSamples + Zscore^2))]
eqtl[, beta := Zscore * se]
eqtl <- eqtl[is.finite(beta) & is.finite(se) & se > 0 & maf > 0 & maf <= 0.5]

instruments <- as.data.table(readRDS(file.path(project_dir, "import", "tyk2_exp_clumped.rds")))
exposure <- as.data.frame(instruments)
exposure$exposure <- "TYK2 eQTLGen whole blood"
exposure$id.exposure <- "TYK2"
manifest <- fread(file.path(raw_dir, "finngen_R13_manifest.tsv"))
metadata <- manifest[phenocode %in% c("M13_RHEUMA", "RHEUMA_SEROPOS_STRICT")]
stopifnot(nrow(metadata) == 2L)

prior_grid <- CJ(p1 = c(1e-5, 1e-4), p2 = c(1e-5, 1e-4), p12 = c(1e-6, 1e-5, 1e-4))
prior_grid <- prior_grid[p12 <= pmin(p1, p2)]
all_coloc <- list()
all_coverage <- list()
all_mr <- list()
all_het <- list()
all_egger <- list()
all_single <- list()
all_loo <- list()
all_presso <- list()
all_exclusions <- list()
all_harmonised <- list()
all_top <- list()
all_comparison <- list()
all_aligned <- list()
all_regional_gwas <- list()

for (i in seq_len(nrow(metadata))) {
  meta <- metadata[i]
  endpoint <- meta$phenocode
  fg <- fread(file.path(raw_dir, paste0("finngen_R13_", endpoint, "_TYK2_complete_region.tsv.gz")))
  required <- c("#chrom", "pos", "ref", "alt", "rsids", "beta", "sebeta", "pval", "af_alt")
  stopifnot(all(required %in% names(fg)), all(fg[["#chrom"]] == 19L))
  fg[, varid := paste(`#chrom`, pos, ref, alt, sep = ":")]
  stopifnot(!anyDuplicated(fg$varid))
  fg[, qc_keep := is.finite(beta) & is.finite(sebeta) & sebeta > 0 &
       is.finite(pval) & pval > 0 & pval <= 1 & is.finite(af_alt) & af_alt > 0 & af_alt < 1]
  all_regional_gwas[[i]] <- fg[qc_keep == TRUE, .(
    outcome_id = endpoint, position_hg38 = pos,
    layer = "Clinically defined RA GWAS", neglog10p = -log10(pval)
  )]

  # Resolve rsID aliases using the actual allele pair; no constant coordinate offset is used for matching.
  aliases <- fg[, .(SNP = unlist(strsplit(as.character(rsids), "[,;]"))), by = varid]
  aliases <- unique(aliases[grepl("^rs[0-9]+$", SNP)])
  fg_rs <- merge(aliases, fg[qc_keep == TRUE], by = "varid")
  aligned <- merge(eqtl[, .(
    SNP, position_hg19 = SNPPos, beta_eqtl = beta, se_eqtl = se,
    p_eqtl = Pvalue, maf_eqtl = maf, n_eqtl = NrSamples,
    ea_eqtl = AssessedAllele, oa_eqtl = OtherAllele
  )], fg_rs[, .(
    SNP, varid, position_hg38 = pos, ref, alt,
    beta_ra = beta, se_ra = sebeta, p_ra = pval,
    maf_ra = pmin(af_alt, 1 - af_alt), eaf_ra = af_alt
  )], by = "SNP")
  aligned[, allele_orientation := fifelse(
    ea_eqtl == alt & oa_eqtl == ref, 1L,
    fifelse(ea_eqtl == ref & oa_eqtl == alt, -1L, 0L)
  )]
  n_rsid_overlap <- uniqueN(aligned$SNP)
  n_nonmatching_allele_records <- nrow(aligned[allele_orientation == 0L])
  n_allele_mismatch <- length(setdiff(aligned$SNP, aligned[allele_orientation != 0L, SNP]))
  excluded <- eqtl[!SNP %in% aligned[allele_orientation != 0L, SNP], .(
    SNP, position_hg19 = SNPPos, exposure_effect_allele = AssessedAllele,
    exposure_other_allele = OtherAllele, p_eqtl = Pvalue
  )]
  excluded[, `:=`(
    outcome_id = endpoint,
    exclusion_reason = fifelse(SNP %in% aligned$SNP, "rsID present but no matching allele pair", "rsID absent from QC-passing official regional statistics")
  )]
  all_exclusions[[i]] <- excluded
  aligned <- aligned[allele_orientation != 0L]
  if (anyDuplicated(aligned$SNP) || anyDuplicated(aligned$varid)) stop("Ambiguous rsID/allele mapping for ", endpoint)
  aligned[, `:=`(beta_ra_aligned = beta_ra * allele_orientation, outcome_id = endpoint)]
  if (nrow(aligned) < 5000L) stop("Insufficient complete eQTL-GWAS overlap: ", endpoint)
  all_aligned[[i]] <- aligned

  old_json <- fromJSON(file.path(data_dir, "raw", paste0("finngen_r13_", endpoint, "_TYK2_region.json")))$data
  old <- as.data.table(old_json)
  all_coverage[[i]] <- data.table(
    outcome_id = endpoint,
    ncase = meta$num_cases, ncontrol = meta$num_controls,
    n_official_regional_variants = nrow(fg),
    n_qc_variants = sum(fg$qc_keep),
    n_regional_plotting_variants = nrow(old),
    n_eqtlgen_variants = nrow(eqtl),
    n_rsid_overlap = n_rsid_overlap,
    n_allele_mismatches = n_allele_mismatch,
    n_nonmatching_alias_allele_records = n_nonmatching_allele_records,
    n_final_coloc_variants = nrow(aligned),
    n_original_instruments_recovered = uniqueN(fg_rs[SNP %in% instruments$SNP, SNP]),
    region_grch38 = "19:9365604-11365604",
    matching_rule = "rsID and exact unordered allele pair; no significance filtering",
    source_url = meta$path_https
  )

  dataset_eqtl <- list(
    beta = aligned$beta_eqtl, varbeta = aligned$se_eqtl^2,
    snp = aligned$SNP, position = aligned$position_hg38,
    MAF = aligned$maf_eqtl, N = round(median(aligned$n_eqtl)), sdY = 1, type = "quant"
  )
  n_outcome <- meta$num_cases + meta$num_controls
  dataset_ra <- list(
    beta = aligned$beta_ra_aligned, varbeta = aligned$se_ra^2,
    snp = aligned$SNP, position = aligned$position_hg38,
    MAF = aligned$maf_ra, N = n_outcome, s = meta$num_cases / n_outcome, type = "cc"
  )
  rows <- list()
  for (j in seq_len(nrow(prior_grid))) {
    prior <- prior_grid[j]
    fit <- coloc.abf(dataset_eqtl, dataset_ra, p1 = prior$p1, p2 = prior$p2, p12 = prior$p12)
    rows[[j]] <- data.table(
      outcome_id = endpoint, n_snps = nrow(aligned),
      p1 = prior$p1, p2 = prior$p2, p12 = prior$p12,
      PP.H0 = unname(fit$summary[["PP.H0.abf"]]),
      PP.H1 = unname(fit$summary[["PP.H1.abf"]]),
      PP.H2 = unname(fit$summary[["PP.H2.abf"]]),
      PP.H3 = unname(fit$summary[["PP.H3.abf"]]),
      PP.H4 = unname(fit$summary[["PP.H4.abf"]])
    )
    if (prior$p1 == 1e-4 && prior$p2 == 1e-4 && prior$p12 == 1e-5) {
      top <- as.data.table(fit$results)
      setorder(top, -SNP.PP.H4)
      top[, outcome_id := endpoint]
      all_top[[i]] <- top[1:min(.N, 100L)]
    }
  }
  all_coloc[[i]] <- rbindlist(rows)

  outcome_inst <- fg_rs[SNP %in% instruments$SNP]
  if (nrow(outcome_inst) != 5L) stop("Expected exactly five direct instruments for ", endpoint)
  outcome_dat <- data.frame(
    SNP = outcome_inst$SNP, beta.outcome = outcome_inst$beta,
    se.outcome = outcome_inst$sebeta,
    effect_allele.outcome = outcome_inst$alt, other_allele.outcome = outcome_inst$ref,
    eaf.outcome = outcome_inst$af_alt, pval.outcome = outcome_inst$pval,
    samplesize.outcome = n_outcome, ncase.outcome = meta$num_cases,
    ncontrol.outcome = meta$num_controls,
    outcome = if (endpoint == "M13_RHEUMA") "FinnGen R13 rheumatoid arthritis" else "FinnGen R13 seropositive RA, strict",
    id.outcome = endpoint, stringsAsFactors = FALSE
  )
  h <- harmonise_data(exposure, outcome_dat, action = 2)
  h <- h[h$mr_keep, , drop = FALSE]
  h <- merge(h, outcome_inst[, .(
    SNP, varid_hg38 = varid,
    data_route = "Official FinnGen R13 BGZF summary statistics, Tabix regional extraction",
    se_source = "reported in official summary statistics"
  )], by = "SNP")
  stopifnot(nrow(h) == 5L)
  set.seed(20261002 + i)
  all_harmonised[[i]] <- as.data.table(h)
  all_mr[[i]] <- as.data.table(mr(h))
  all_het[[i]] <- as.data.table(mr_heterogeneity(h))
  all_egger[[i]] <- as.data.table(mr_pleiotropy_test(h))
  all_single[[i]] <- as.data.table(mr_singlesnp(h))
  all_loo[[i]] <- as.data.table(mr_leaveoneout(h))
  all_presso[[endpoint]] <- mr_presso(
    BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
    SdOutcome = "se.outcome", SdExposure = "se.exposure",
    OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = h,
    NbDistribution = 10000, SignifThreshold = 0.05, seed = 20260911 + i
  )

  old_h <- fread(file.path(data_dir, "results_finngen_r13_five_instruments", "tyk2_finngen_r13_five_instrument_harmonised.csv"))
  cmp <- merge(
    as.data.table(h)[, .(SNP, official_beta = beta.outcome, official_se = se.outcome, official_eaf = eaf.outcome)],
    old_h[id.outcome == endpoint, .(SNP, previous_beta = beta.outcome, previous_se = se.outcome,
                                  previous_eaf = eaf.outcome, previous_se_source = se_source)],
    by = "SNP"
  )
  cmp[, `:=`(
    outcome_id = endpoint,
    beta_difference = official_beta - previous_beta,
    se_relative_difference = abs(official_se - previous_se) / official_se,
    eaf_difference = official_eaf - previous_eaf
  )]
  all_comparison[[i]] <- cmp
}

coverage <- rbindlist(all_coverage)
coloc_results <- rbindlist(all_coloc)
mr_results <- rbindlist(all_mr, fill = TRUE)
fwrite(coverage, file.path(out_dir, "TYK2_R13_Complete_Region_Coverage.csv"))
fwrite(coloc_results, file.path(out_dir, "TYK2_R13_Complete_Region_Coloc_Prior_Sensitivity.csv"))
fwrite(rbindlist(all_top, fill = TRUE), file.path(out_dir, "TYK2_R13_Complete_Region_Coloc_Top_SNPs.csv"))
fwrite(rbindlist(all_aligned), file.path(out_dir, "TYK2_R13_Complete_Region_Aligned.csv.gz"))
fwrite(mr_results, file.path(out_dir, "TYK2_R13_Official_SE_MR_Methods.csv"))
fwrite(rbindlist(all_het), file.path(out_dir, "TYK2_R13_Official_SE_Heterogeneity.csv"))
fwrite(rbindlist(all_egger), file.path(out_dir, "TYK2_R13_Official_SE_Egger_Intercept.csv"))
fwrite(rbindlist(all_single), file.path(out_dir, "TYK2_R13_Official_SE_Single_SNP.csv"))
fwrite(rbindlist(all_loo), file.path(out_dir, "TYK2_R13_Official_SE_Leave_One_Out.csv"))
fwrite(rbindlist(all_harmonised), file.path(out_dir, "TYK2_R13_Official_SE_Harmonised.csv"))
fwrite(rbindlist(all_comparison), file.path(out_dir, "TYK2_R13_Official_vs_Previous_API_Instrument_Audit.csv"))
fwrite(rbindlist(all_exclusions), file.path(out_dir, "TYK2_R13_Complete_Region_Excluded_eQTLs.csv"))
saveRDS(all_presso, file.path(out_dir, "TYK2_R13_Official_SE_MR_PRESSO.rds"))
capture.output(print(all_presso), file = file.path(out_dir, "TYK2_R13_Official_SE_MR_PRESSO.txt"))
presso_summary <- lapply(seq_len(nrow(metadata)), function(i) {
  endpoint <- metadata$phenocode[[i]]
  result <- all_presso[[endpoint]][["MR-PRESSO results"]]
  distortion <- result[["Distortion Test"]]
  indices <- distortion[["Outliers Indices"]]
  data.table(
    outcome_id = endpoint,
    global_p = as.numeric(result[["Global Test"]][["Pvalue"]]),
    outlier_snps = paste(all_harmonised[[i]]$SNP[as.integer(indices)], collapse = ";"),
    n_outliers = length(indices),
    distortion_p = as.numeric(distortion[["Pvalue"]]),
    interpretation = "Sparse-instrument sensitivity diagnostic; corrected estimate is not definitive"
  )
})
fwrite(rbindlist(presso_summary), file.path(out_dir, "TYK2_R13_Official_SE_MR_PRESSO_Summary.csv"))

plot_data <- rbindlist(list(
  rbindlist(all_aligned)[, .(outcome_id, position_hg38, layer = "Whole-blood TYK2 eQTL", neglog10p = -log10(p_eqtl))],
  rbindlist(all_regional_gwas)
))
p <- ggplot(plot_data, aes(position_hg38 / 1e6, neglog10p, color = layer)) +
  geom_point(size = 0.55, alpha = 0.55) +
  geom_point(data = rbindlist(all_regional_gwas)[position_hg38 == 10352442L],
             color = "#C74343", shape = 18, size = 2.2) +
  geom_vline(xintercept = 10.352442, linetype = 2, color = "#C74343", linewidth = 0.6) +
  facet_grid(layer ~ outcome_id, scales = "free_y") +
  scale_color_manual(values = c("Whole-blood TYK2 eQTL" = "#2878A8", "Clinically defined RA GWAS" = "#168B6B"), guide = "none") +
  labs(title = "Complete-summary-statistic TYK2 regional association patterns",
       subtitle = "GWAS: all 15,625 regional variants; eQTL: 6,291 allele-matched variants",
       caption = "Dashed line: main FinnGen RA credible-set variant rs34536443, absent from the eQTLGen TYK2 cis records",
       x = "Chromosome 19 position (GRCh38, Mb)", y = "-log10(P)") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), strip.background = element_blank(),
        plot.caption = element_text(hjust = 0))
ggsave(file.path(out_dir, "Figure_TYK2_R13_Complete_Regional_Associations.pdf"), p, width = 11, height = 7)
ggsave(file.path(out_dir, "Figure_TYK2_R13_Complete_Regional_Associations.png"), p, width = 11, height = 7, dpi = 240)

primary <- coloc_results[p1 == 1e-4 & p2 == 1e-4 & p12 == 1e-5]
writeLines(c(
  "# Complete FinnGen R13 TYK2-region summary-statistic analysis",
  "",
  "Input: the official unfiltered BGZF summary-statistic files queried using their Tabix indexes. Region: GRCh38 chr19:9,365,604-11,365,604. Exposure: all 6,995 TYK2 cis-eQTL records in the official unthresholded eQTLGen release. Matching: rsID plus exact allele pair.",
  "",
  "## Coverage",
  paste(capture.output(print(coverage)), collapse = "\n"),
  "",
  "## Primary single-causal-variant ABF",
  paste(capture.output(print(primary)), collapse = "\n"),
  "",
  "## Official-SE MR",
  paste(capture.output(print(mr_results[method == "Inverse variance weighted"])), collapse = "\n"),
  "",
  "## Prior sensitivity",
  paste(capture.output(print(coloc_results)), collapse = "\n"),
  "",
  "## Interpretation boundary",
  "Complete regional inputs resolve the plotting-API ascertainment limitation. The ABF still assumes at most one causal variant per trait. FinnGen's in-sample fine-mapping identifies multiple regional signals, so a low H4 cannot exclude a shared secondary signal.",
  "Separate molecular-QTL and GWAS LD matrices with signed allele-aligned correlations are required for defensible coloc.susie. The public SISu42 LD API exposes only r2 and D-prime, so it cannot substitute for a signed matrix."
), file.path(out_dir, "TYK2_R13_Complete_Summary_Analysis.md"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Complete summary-statistic analysis saved to ", out_dir)
