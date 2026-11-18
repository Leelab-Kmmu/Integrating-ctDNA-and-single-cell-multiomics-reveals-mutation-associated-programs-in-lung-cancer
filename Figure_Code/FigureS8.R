#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(survminer)
  library(ggplot2)
})

set.seed(20260910)
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args)) args[[1L]] else "E:/Lung/EGA/PDF/gene_list/TCGA"
input_file <- file.path(workdir, "results", "TCGA-LUAD_19genes_OS_data.csv")
outdir <- file.path(workdir, "TCGA_LUAD_LASSO_Figure6")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

genes <- c("E2F3", "AKR1C2", "AKR1C3", "CDKN2A", "COL4A4", "CPS1", "DIP2C",
           "GSTCD", "IGF2BP3", "IGSF1", "KIAA0319", "NPAS2", "NPNT", "PCSK6",
           "PPARGC1A", "PTPRH", "RGS17", "VCAN", "WEE1")

dat <- read.csv(input_file, check.names = FALSE)
genes_use <- intersect(genes, names(dat))
dat <- dat[complete.cases(dat[, c("OS_time_days", "OS", genes_use)]), , drop = FALSE]
dat <- dat[dat$OS_time_days > 0 & dat$OS %in% c(0, 1), , drop = FALSE]
x_raw <- as.matrix(dat[, genes_use, drop = FALSE])
x <- scale(x_raw)
y <- Surv(dat$OS_time_days, dat$OS)

cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10,
                   type.measure = "deviance", standardize = FALSE)
fullfit <- glmnet(x, y, family = "cox", alpha = 1, standardize = FALSE)

extract_coef <- function(cvfit, s, label) {
  b <- as.matrix(coef(cvfit, s = s))[, 1]
  data.frame(gene = names(b), coefficient = as.numeric(b),
             selected = as.numeric(b) != 0, lambda_rule = label,
             stringsAsFactors = FALSE)
}
coef_min_all <- extract_coef(cvfit, "lambda.min", "lambda.min")
coef_1se_all <- extract_coef(cvfit, "lambda.1se", "lambda.1se")
coef_all <- rbind(coef_min_all, coef_1se_all)
write.csv(coef_all, file.path(outdir, "TCGA_LUAD_LASSO_all_coefficients.csv"), row.names = FALSE)
write.csv(coef_min_all[coef_min_all$selected, ],
          file.path(outdir, "TCGA_LUAD_LASSO_lambda_min_nonzero_coefficients.csv"), row.names = FALSE)
write.csv(coef_1se_all[coef_1se_all$selected, ],
          file.path(outdir, "TCGA_LUAD_LASSO_lambda_1se_nonzero_coefficients.csv"), row.names = FALSE)

risk_coef <- setNames(coef_min_all$coefficient, coef_min_all$gene)
risk_coef <- risk_coef[risk_coef != 0]
if (!length(risk_coef)) stop("lambda.min selected no genes")
risk_score <- as.numeric(x[, names(risk_coef), drop = FALSE] %*% risk_coef)
cutoff <- median(risk_score, na.rm = TRUE)
risk_group <- factor(ifelse(risk_score > cutoff, "High", "Low"), levels = c("Low", "High"))

risk_df <- data.frame(
  patient_id = dat$patient_id, OS_time_days = dat$OS_time_days, OS_event = dat$OS,
  time = dat$OS_time_days / 30.4375, time_unit = "months",
  risk_score = risk_score, risk_group = risk_group, stringsAsFactors = FALSE
)
write.csv(risk_df, file.path(outdir, "TCGA_LUAD_LASSO_risk_score_samples.csv"), row.names = FALSE)

fit_km <- survfit(Surv(time, OS_event) ~ risk_group, data = risk_df)
fit_group <- coxph(Surv(time, OS_event) ~ risk_group, data = risk_df)
fit_score <- coxph(Surv(time, OS_event) ~ as.numeric(scale(risk_score)), data = risk_df)
sg <- summary(fit_group)
ss <- summary(fit_score)
lr <- survdiff(Surv(time, OS_event) ~ risk_group, data = risk_df)
logrank_p <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)

stats <- data.frame(
  dataset = "TCGA-LUAD", endpoint = "OS", n = nrow(risk_df), events = sum(risk_df$OS_event == 1),
  selected_genes = paste(names(risk_coef), collapse = ";"), selected_gene_n = length(risk_coef),
  lambda_min = cvfit$lambda.min, lambda_1se = cvfit$lambda.1se, median_cutoff = cutoff,
  HR_high_vs_low = sg$conf.int[1, "exp(coef)"], lower95 = sg$conf.int[1, "lower .95"],
  upper95 = sg$conf.int[1, "upper .95"], cox_group_p = sg$coefficients[1, "Pr(>|z|)"],
  logrank_p = logrank_p, HR_per_SD_score = ss$conf.int[1, "exp(coef)"],
  score_lower95 = ss$conf.int[1, "lower .95"], score_upper95 = ss$conf.int[1, "upper .95"],
  cox_score_p = ss$coefficients[1, "Pr(>|z|)"], stringsAsFactors = FALSE
)
write.csv(stats, file.path(outdir, "TCGA_LUAD_LASSO_risk_score_statistics.csv"), row.names = FALSE)

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
hr_text <- paste0(
  "Ref: Low risk\n",
  "n=", stats$n, ", events=", stats$events,
  "\nLow=", sum(risk_group == "Low"), ", High=", sum(risk_group == "High"),
  "\nHigh risk: HR=", sprintf("%.2f", stats$HR_high_vs_low),
  " (95% CI ", sprintf("%.2f", stats$lower95), "-", sprintf("%.2f", stats$upper95), ")",
  "\nCox P=", fmt_p(stats$cox_group_p), "\nLog-rank P=", fmt_p(stats$logrank_p)
)

p_km <- ggsurvplot(
  fit_km, data = risk_df, pval = TRUE, conf.int = FALSE, risk.table = TRUE,
  risk.table.height = 0.28, palette = c("#4DBBD5", "#E64B35"),
  xlab = "Time (months)", ylab = "Overall survival probability",
  legend.title = "Risk group", legend.labs = c("Low risk", "High risk"),
  title = paste0("TCGA-LUAD OS by ", length(risk_coef), "-gene LASSO risk score")
)
p_km$plot <- p_km$plot + theme_bw() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5, size = 11),
        axis.text = element_text(size = 9, color = "black"), axis.title = element_text(size = 10),
        legend.title = element_text(size = 9), legend.text = element_text(size = 8)) +
  annotate("text", x = max(risk_df$time) * 0.04, y = 0.08, label = hr_text,
           hjust = 0, vjust = 0, size = 2.0)
p_km$table <- p_km$table + theme_bw() +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9))

pdf(file.path(outdir, "TCGA_LUAD_LASSO_risk_score_KM_Figure6_style.pdf"), width = 5.8, height = 7.0)
print(p_km)
dev.off()
png(file.path(outdir, "TCGA_LUAD_LASSO_risk_score_KM_Figure6_style.png"), width = 1740, height = 2100, res = 300)
print(p_km)
dev.off()

pdf(file.path(outdir, "TCGA_LUAD_LASSO_CV_curve.pdf"), width = 6, height = 5.5)
plot(cvfit); title("TCGA-LUAD 10-fold CV LASSO-Cox")
dev.off()
pdf(file.path(outdir, "TCGA_LUAD_LASSO_coefficient_path.pdf"), width = 7, height = 5.5)
plot(fullfit, xvar = "lambda", label = TRUE); title("TCGA-LUAD LASSO-Cox coefficient paths")
dev.off()

coef_plot_df <- coef_min_all[coef_min_all$selected, , drop = FALSE]
coef_plot_df$gene <- factor(coef_plot_df$gene, levels = coef_plot_df$gene[order(coef_plot_df$coefficient)])
p_coef <- ggplot(coef_plot_df, aes(coefficient, gene, fill = coefficient > 0)) +
  geom_col(width = 0.72) + geom_vline(xintercept = 0, color = "grey35") +
  scale_fill_manual(values = c(`TRUE` = "#E64B35", `FALSE` = "#4DBBD5"), guide = "none") +
  labs(title = "TCGA-LUAD LASSO-Cox coefficients", subtitle = "Selected at lambda.min",
       x = "Coefficient", y = NULL) + theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5), axis.text.y = element_text(color = "black"))
ggsave(file.path(outdir, "TCGA_LUAD_LASSO_lambda_min_coefficients.pdf"), p_coef, width = 6.5, height = 5.5)
ggsave(file.path(outdir, "TCGA_LUAD_LASSO_lambda_min_coefficients.png"), p_coef, width = 6.5, height = 5.5, dpi = 300)

writeLines(c(
  "TCGA-LUAD LASSO-Cox risk model",
  "Training cohort: TCGA-LUAD only; endpoint: overall survival.",
  "Predictors: z-scored expression of the 19 candidate genes.",
  "Selection: 10-fold cross-validation; primary coefficients at lambda.min.",
  "Risk score = sum(lambda.min coefficient * z-scored gene expression).",
  "Median training-cohort risk score defines Low and High groups.",
  "This training-set KM estimate is optimistic; GEO should be used as independent validation."
), file.path(outdir, "README_TCGA_LUAD_LASSO.txt"))

message("Selected genes: ", paste(names(risk_coef), collapse = ", "))
message("Done: ", normalizePath(outdir, winslash = "/"))
