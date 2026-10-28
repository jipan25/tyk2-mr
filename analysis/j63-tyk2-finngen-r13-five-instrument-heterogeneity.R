suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else "data/upstream/clinical"
project_dir <- if (length(args) >= 2) args[[2]] else getwd()
input_dir <- if (length(args) >= 3) args[[3]] else file.path(data_dir, "results_finngen_r13_five_instruments")
raw_dir <- file.path(data_dir, "raw")
out_dir <- if (length(args) >= 4) args[[4]] else file.path(data_dir, "results_finngen_r13_five_instrument_heterogeneity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

complete_input <- file.exists(file.path(input_dir, "TYK2_R13_Official_SE_MR_Methods.csv"))
input_path <- function(legacy_name, official_name) {
  file.path(input_dir, if (complete_input) official_name else legacy_name)
}

single <- fread(input_path("tyk2_finngen_r13_five_instrument_single_snp.csv", "TYK2_R13_Official_SE_Single_SNP.csv"))
single <- single[grepl("^rs[0-9]+$", SNP)]
loo <- fread(input_path("tyk2_finngen_r13_five_instrument_leave_one_out.csv", "TYK2_R13_Official_SE_Leave_One_Out.csv"))
loo <- loo[SNP != "All"]
heterogeneity <- fread(input_path("tyk2_finngen_r13_five_instrument_heterogeneity.csv", "TYK2_R13_Official_SE_Heterogeneity.csv"))
heterogeneity <- heterogeneity[method == "Inverse variance weighted"]
mr <- fread(input_path("tyk2_finngen_r13_five_instrument_mr_methods.csv", "TYK2_R13_Official_SE_MR_Methods.csv"))
mr_ivw <- mr[method == "Inverse variance weighted"]
harmonised <- fread(input_path("tyk2_finngen_r13_five_instrument_harmonised.csv", "TYK2_R13_Official_SE_Harmonised.csv"))
recovered <- fread(file.path(raw_dir, "finngen_r13_exact_variant_api_recovered_instruments.csv"))

outcome_meta <- data.table(
  id.outcome = c("M13_RHEUMA", "RHEUMA_SEROPOS_STRICT"),
  outcome_short = c("FinnGen R13 broad RA", "FinnGen R13 seropositive RA, strict"),
  json_file = file.path(raw_dir, c(
    "finngen_r13_M13_RHEUMA_TYK2_region.json",
    "finngen_r13_RHEUMA_SEROPOS_STRICT_TYK2_region.json"
  ))
)

q_rows <- list()
loo_rows <- list()
annotation_rows <- list()
summary_rows <- list()

for (i in seq_len(nrow(outcome_meta))) {
  meta <- outcome_meta[i]
  x <- single[id.outcome == meta$id.outcome]
  full <- mr_ivw[id.outcome == meta$id.outcome]
  reported_q <- heterogeneity[id.outcome == meta$id.outcome]
  if (nrow(x) != 5L || nrow(full) != 1L || nrow(reported_q) != 1L) {
    stop("Unexpected input dimensions for ", meta$id.outcome)
  }

  x[, weight := 1 / se^2]
  fixed_beta <- x[, sum(weight * b) / sum(weight)]
  x[, q_contribution := weight * (b - fixed_beta)^2]
  x[, q_contribution_percent := 100 * q_contribution / sum(q_contribution)]
  x[, `:=`(
    fixed_ivw_beta = fixed_beta,
    reported_random_ivw_beta = full$b,
    reported_random_ivw_se = full$se,
    reported_random_ivw_p = full$pval,
    reported_total_q = reported_q$Q,
    calculated_total_q = sum(q_contribution),
    outcome_short = meta$outcome_short
  )]
  if (abs(sum(x$q_contribution) - reported_q$Q) > 1e-6) {
    stop("Q contribution sum does not match reported Q for ", meta$id.outcome)
  }

  loo_x <- loo[id.outcome == meta$id.outcome]
  for (snp in x$SNP) {
    kept <- x[SNP != snp]
    kept_beta <- kept[, sum(weight * b) / sum(weight)]
    residual_q <- kept[, sum(weight * (b - kept_beta)^2)]
    residual_df <- nrow(kept) - 1L
    current_loo <- loo_x[SNP == snp]
    current_q <- x[SNP == snp]
    loo_rows[[length(loo_rows) + 1L]] <- data.table(
      id.outcome = meta$id.outcome,
      outcome_short = meta$outcome_short,
      omitted_snp = snp,
      full_random_ivw_beta = full$b,
      full_random_ivw_se = full$se,
      full_random_ivw_p = full$pval,
      loo_beta = current_loo$b,
      loo_se = current_loo$se,
      loo_p = current_loo$p,
      beta_change = current_loo$b - full$b,
      absolute_relative_beta_change = abs(current_loo$b - full$b) / abs(full$b),
      direction_preserved = sign(current_loo$b) == sign(full$b),
      omitted_q_contribution = current_q$q_contribution,
      omitted_q_percent = current_q$q_contribution_percent,
      residual_q = residual_q,
      residual_q_df = residual_df,
      residual_q_p = pchisq(residual_q, df = residual_df, lower.tail = FALSE),
      major_q_driver = current_q$q_contribution_percent >= 50,
      material_beta_influence = abs(current_loo$b - full$b) / abs(full$b) >= 0.25,
      nominal_significance_preserved = current_loo$p < 0.05
    )
  }

  h <- harmonised[id.outcome == meta$id.outcome]
  regional <- as.data.table(fromJSON(meta$json_file, simplifyDataFrame = TRUE)$data)
  regional <- regional[rsid %in% h$SNP, .(
    SNP = rsid,
    most_severe_consequence = most_severe,
    nearest_gene = nearest_genes
  )]
  rescued_ann <- unique(recovered[outcome_id == meta$id.outcome, .(
    SNP = snp,
    most_severe_consequence = most_severe,
    nearest_gene
  )])
  ann <- merge(
    h[, .(
      SNP, varid_hg38, data_route, se_source,
      exposure_effect_allele = effect_allele.exposure,
      exposure_other_allele = other_allele.exposure,
      exposure_eaf = eaf.exposure,
      exposure_beta = beta.exposure,
      exposure_se = se.exposure,
      exposure_p = pval.exposure,
      outcome_effect_allele = effect_allele.outcome,
      outcome_other_allele = other_allele.outcome,
      outcome_eaf = eaf.outcome,
      outcome_beta_harmonized = beta.outcome,
      outcome_se = se.outcome,
      outcome_p = pval.outcome
    )],
    rbindlist(list(regional, rescued_ann), use.names = TRUE),
    by = "SNP",
    all.x = TRUE
  )
  ann[, calculated_wald_ratio := outcome_beta_harmonized / exposure_beta]
  ann <- merge(
    ann,
    x[, .(SNP, reported_wald_ratio = b, reported_wald_se = se,
           reported_wald_p = p, q_contribution, q_contribution_percent)],
    by = "SNP"
  )
  ann[, `:=`(id.outcome = meta$id.outcome, outcome_short = meta$outcome_short)]
  annotation_rows[[i]] <- ann

  q_rows[[i]] <- x
  summary_rows[[i]] <- data.table(
    id.outcome = meta$id.outcome,
    outcome_short = meta$outcome_short,
    n_snps = nrow(x),
    random_ivw_beta = full$b,
    random_ivw_se = full$se,
    random_ivw_p = full$pval,
    total_q = reported_q$Q,
    q_df = reported_q$Q_df,
    q_p = reported_q$Q_pval,
    top_q_snp = x$SNP[which.max(x$q_contribution)],
    top_q_percent = max(x$q_contribution_percent),
    n_inverse_single_snp = sum(x$b < 0),
    all_leave_one_out_directions_preserved = all(loo_x$b < 0),
    leave_one_out_nominal_p_count = sum(loo_x$p < 0.05)
  )
}

q_decomposition <- rbindlist(q_rows, fill = TRUE)
loo_influence <- rbindlist(loo_rows, fill = TRUE)
annotations <- rbindlist(annotation_rows, fill = TRUE)
summary <- rbindlist(summary_rows, fill = TRUE)

fwrite(q_decomposition, file.path(out_dir, "TYK2_R13_Five_Instrument_SNP_Q_Contributions.csv"))
fwrite(loo_influence, file.path(out_dir, "TYK2_R13_Five_Instrument_Leave_One_Out_Influence.csv"))
fwrite(annotations, file.path(out_dir, "TYK2_R13_Five_Instrument_Annotation_and_Direction.csv"))
fwrite(summary, file.path(out_dir, "TYK2_R13_Five_Instrument_Heterogeneity_Summary.csv"))

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(list(
    Summary = as.data.frame(summary),
    Q_Contributions = as.data.frame(q_decomposition),
    Leave_One_Out = as.data.frame(loo_influence),
    SNP_Annotation = as.data.frame(annotations)
  ), file.path(out_dir, "TYK2_R13_Five_Instrument_Heterogeneity.xlsx"))
}

q_plot <- copy(q_decomposition)
q_plot[, SNP := factor(SNP, levels = rev(unique(SNP)))]
snp_colors <- c(
  "rs1058154" = "#D96CE2",
  "rs11882133" = "#0EA5D9",
  "rs12720250" = "#08B879",
  "rs144309607" = "#AAA600",
  "rs34725611" = "#F56E64"
)
p_q <- ggplot(q_plot, aes(q_contribution_percent, SNP, color = SNP)) +
  geom_segment(aes(x = 0, xend = q_contribution_percent, yend = SNP), linewidth = 0.9, alpha = 0.7) +
  geom_point(size = 3.2) +
  geom_text(aes(label = sprintf("%.1f%%", q_contribution_percent)), nudge_x = 2, hjust = 0, size = 3.1) +
  facet_wrap(~outcome_short, ncol = 1) +
  coord_cartesian(xlim = c(0, max(q_plot$q_contribution_percent) + 12)) +
  scale_color_manual(values = snp_colors) +
  labs(title = "A. SNP contributions to Cochran heterogeneity", x = "Contribution to total Q (%)", y = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"), strip.background = element_blank(), strip.text = element_text(face = "bold"))

p_influence <- ggplot(loo_influence, aes(100 * absolute_relative_beta_change, omitted_snp, color = omitted_snp)) +
  geom_vline(xintercept = 25, linetype = 2, color = "#8A969C") +
  geom_segment(aes(x = 0, xend = 100 * absolute_relative_beta_change, yend = omitted_snp), linewidth = 0.9, alpha = 0.7) +
  geom_point(size = 3.2) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * absolute_relative_beta_change)), nudge_x = 1.5, hjust = 0, size = 3.1) +
  facet_wrap(~outcome_short, ncol = 1) +
  coord_cartesian(xlim = c(0, max(100 * loo_influence$absolute_relative_beta_change) + 8)) +
  scale_color_manual(values = snp_colors) +
  labs(title = "B. Influence on pooled random-effects IVW beta", x = "Absolute change in IVW beta (%)", y = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"), strip.background = element_blank(), strip.text = element_text(face = "bold"))

draw_figure <- function() {
  grid.newpage()
  grid.text("TYK2-FinnGen R13 five-instrument heterogeneity audit", x = 0.03, y = 0.975,
            just = "left", gp = gpar(fontsize = 17, fontface = "bold", col = "#17242B"))
  print(p_q, newpage = FALSE, vp = viewport(x = 0.27, y = 0.49, width = 0.50, height = 0.88))
  print(p_influence, newpage = FALSE, vp = viewport(x = 0.77, y = 0.49, width = 0.46, height = 0.88))
}

pdf(file.path(out_dir, "Supplementary_Figure_TYK2_R13_Five_Instrument_Heterogeneity.pdf"), width = 12, height = 7)
draw_figure()
dev.off()
png(file.path(out_dir, "Supplementary_Figure_TYK2_R13_Five_Instrument_Heterogeneity.png"), width = 3600, height = 2100, res = 300)
draw_figure()
dev.off()

loo_plot <- copy(loo_influence)
loo_plot[, `:=`(
  lower = loo_beta - 1.96 * loo_se,
  upper = loo_beta + 1.96 * loo_se,
  label = paste0("Omit ", omitted_snp)
)]
full_rows <- summary[, .(
  id.outcome,
  outcome_short,
  omitted_snp = "Full model",
  loo_beta = random_ivw_beta,
  loo_se = random_ivw_se,
  loo_p = random_ivw_p,
  lower = random_ivw_beta - 1.96 * random_ivw_se,
  upper = random_ivw_beta + 1.96 * random_ivw_se,
  label = "Full random-effects IVW"
)]
forest <- rbindlist(list(
  loo_plot[, .(id.outcome, outcome_short, omitted_snp, loo_beta, loo_se, loo_p, lower, upper, label)],
  full_rows
), fill = TRUE)
label_order <- c("Full random-effects IVW", paste0("Omit ", unique(q_decomposition$SNP)))
forest[, label := factor(label, levels = rev(label_order))]

p_loo <- ggplot(forest, aes(loo_beta, label, color = omitted_snp == "Full model")) +
  geom_vline(xintercept = 0, linetype = 2, color = "#8A969C") +
  geom_segment(aes(x = lower, xend = upper, yend = label), linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("P=%.3g", loo_p)), nudge_y = 0.22, size = 3) +
  facet_wrap(~outcome_short, ncol = 2, scales = "free_x") +
  scale_color_manual(values = c("TRUE" = "#168B6B", "FALSE" = "#2878A8"), guide = "none") +
  labs(
    title = "TYK2-FinnGen R13 five-instrument leave-one-out analysis",
    subtitle = "All pooled estimates preserve the inverse direction; each leave-one-out estimate retains four instruments",
    x = "MR beta (95% CI)", y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), strip.background = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(out_dir, "Supplementary_Figure_TYK2_R13_Five_Instrument_Leave_One_Out.pdf"), p_loo, width = 11, height = 6.5)
ggsave(file.path(out_dir, "Supplementary_Figure_TYK2_R13_Five_Instrument_Leave_One_Out.png"), p_loo, width = 11, height = 6.5, dpi = 300)

writeLines(c(
  "# TYK2 FinnGen R13 five-instrument heterogeneity audit",
  "",
  paste(capture.output(print(summary)), collapse = "\n"),
  "",
  "## SNP-level Q contributions",
  paste(capture.output(print(q_decomposition[, .(id.outcome, SNP, b, se, p, q_contribution, q_contribution_percent)])), collapse = "\n"),
  "",
  "## Leave-one-out influence",
  paste(capture.output(print(loo_influence)), collapse = "\n"),
  "",
  "## Interpretation",
  if (complete_input) {
    "All five direct instruments use beta, standard error, and alternate-allele frequency reported in the official FinnGen R13 summary statistics. Historical plotting-API and exact-API estimates remain archived."
  } else {
    "The exact-variant API recovery changes the analysis from three to five direct instruments. Results must therefore distinguish the original regional-API analysis from this exact-API-augmented analysis."
  },
  "The five-instrument analysis should not be described as uniformly concordant at the single-SNP level if any ratio estimate changes sign. Leave-one-out estimates provide the more relevant direction-stability diagnostic.",
  "rs11882133 is annotated by FinnGen within KEAP1 despite being a strong cis-eQTL for TYK2; this is a specific horizontal-pleiotropy concern that requires explicit leave-one-out and variant-level interpretation rather than post hoc exclusion."
), file.path(out_dir, "TYK2_R13_Five_Instrument_Heterogeneity_Summary.md"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Completed five-instrument heterogeneity audit in: ", out_dir)
