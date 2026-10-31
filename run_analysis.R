#!/usr/bin/env Rscript

# TCGA-LUAD Somatic Mutation & Survival Analysis
# ================================================
# Author: Sohel Ahmed, Ph.D.
# Date: October 2026

# ---- Configuration ----
PROJECT_DIR <- path.expand("~/tcga-luad-pipeline")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
DATA_DIR    <- file.path(RESULTS_DIR, "data")
MUT_DIR     <- file.path(RESULTS_DIR, "mutation")
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
SURV_DIR    <- file.path(RESULTS_DIR, "survival")

for (d in c(DATA_DIR, MUT_DIR, FIGURES_DIR, SURV_DIR)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("All output will be saved under:", RESULTS_DIR, "\n\n")

# ---- Step 1: Download TCGA-LUAD Mutation Data ----
cat("=== Step 1: Downloading TCGA-LUAD mutation data ===\n")

library(TCGAbiolinks)
library(maftools)

cached_maf <- file.path(DATA_DIR, "tcga_luad_maf.rds")

if (file.exists(cached_maf)) {
    cat("Using cached MAF data:", cached_maf, "\n")
    maf <- readRDS(cached_maf)
} else {
    query_mut <- GDCquery(
        project = "TCGA-LUAD",
        data.category = "Simple Nucleotide Variation",
        data.type = "Masked Somatic Mutation",
        access = "open",
        workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
    )

    GDCdownload(query_mut, method = "api", directory = file.path(DATA_DIR, "GDCdata"))
    maf_data <- GDCprepare(query_mut, directory = file.path(DATA_DIR, "GDCdata"))

    maf <- read.maf(maf = maf_data)
    saveRDS(maf, cached_maf)
    cat("MAF data saved to:", cached_maf, "\n")
}

cat("Samples:", length(getSampleSummary(maf)$Tumor_Sample_Barcode), "\n")
cat("Genes mutated:", length(getGeneSummary(maf)$Hugo_Symbol), "\n\n")

# ---- Step 2: Mutation Summary ----
cat("=== Step 2: Generating Mutation Summary ===\n")

gene_summary <- getGeneSummary(maf)
write.csv(gene_summary, file.path(MUT_DIR, "gene_mutation_summary.csv"), row.names = FALSE)

sample_summary <- getSampleSummary(maf)
write.csv(sample_summary, file.path(MUT_DIR, "sample_mutation_summary.csv"), row.names = FALSE)

cat("Top 10 mutated genes:\n")
print(head(gene_summary[, c("Hugo_Symbol", "MutatedSamples", "AlteredSamples")], 10))
cat("\n")

# ---- Step 3: Oncoplot ----
cat("=== Step 3: Generating Oncoplot ===\n")

png(file.path(FIGURES_DIR, "oncoplot_top20.png"), width = 14, height = 10, units = "in", res = 300)
oncoplot(maf = maf, top = 20, fontSize = 0.6,
         titleText = "TCGA-LUAD: Top 20 Mutated Genes")
dev.off()
cat("Saved: oncoplot_top20.png\n")

# ---- Step 4: Variant Classification Plot ----
cat("=== Step 4: Variant Classification Summary ===\n")

png(file.path(FIGURES_DIR, "variant_classification.png"), width = 12, height = 8, units = "in", res = 300)
plotmafSummary(maf = maf, rmOutlier = TRUE,
               addStat = "median", dashboard = TRUE,
               titvRaw = FALSE)
dev.off()
cat("Saved: variant_classification.png\n")

# ---- Step 5: Lollipop Plots for Key Genes ----
cat("=== Step 5: Lollipop Plots ===\n")

key_genes <- c("TP53", "KRAS", "EGFR", "STK11", "KEAP1")

for (gene in key_genes) {
    tryCatch({
        png(file.path(FIGURES_DIR, paste0("lollipop_", gene, ".png")),
            width = 10, height = 4, units = "in", res = 300)
        lollipopPlot(maf = maf, gene = gene, AACol = "HGVSp_Short",
                     showMutationRate = TRUE,
                     labelPos = "all", labPosSize = 0.8)
        dev.off()
        cat("Saved: lollipop_", gene, ".png\n")
    }, error = function(e) {
        cat("Skipped lollipop for", gene, ":", e$message, "\n")
        try(dev.off(), silent = TRUE)
    })
}

# ---- Step 6: Mutual Exclusivity / Co-occurrence ----
cat("=== Step 6: Somatic Interactions ===\n")

png(file.path(FIGURES_DIR, "somatic_interactions.png"), width = 10, height = 10, units = "in", res = 300)
interactions <- somaticInteractions(maf = maf, top = 15, pvalue = c(0.05, 0.01))
dev.off()
cat("Saved: somatic_interactions.png\n")

write.csv(interactions, file.path(MUT_DIR, "somatic_interactions.csv"), row.names = FALSE)

# ---- Step 7: Tumor Mutation Burden ----
cat("=== Step 7: Computing TMB ===\n")

library(data.table)

tmb <- sample_summary[, .(Tumor_Sample_Barcode, total)]
tmb$TMB <- tmb$total / 38
tmb$TMB_class <- ifelse(tmb$TMB >= 10, "TMB-High", "TMB-Low")

write.csv(tmb, file.path(MUT_DIR, "tmb_scores.csv"), row.names = FALSE)

cat("TMB distribution:\n")
cat("  TMB-High (>=10 mut/Mb):", sum(tmb$TMB_class == "TMB-High"), "\n")
cat("  TMB-Low  (<10 mut/Mb):", sum(tmb$TMB_class == "TMB-Low"), "\n\n")

library(ggplot2)

png(file.path(FIGURES_DIR, "tmb_distribution.png"), width = 8, height = 6, units = "in", res = 300)
print(
    ggplot(tmb, aes(x = TMB)) +
        geom_histogram(aes(fill = TMB_class), bins = 50, alpha = 0.8) +
        geom_vline(xintercept = 10, linetype = "dashed", color = "red", linewidth = 1) +
        scale_fill_manual(values = c("TMB-High" = "#B2182B", "TMB-Low" = "#2166AC")) +
        labs(title = "Tumor Mutation Burden Distribution - TCGA-LUAD",
             x = "TMB (mutations/Mb)", y = "Number of Patients",
             fill = "TMB Class") +
        annotate("text", x = 12, y = Inf, vjust = 2, label = "10 mut/Mb threshold",
                 color = "red", fontface = "italic", size = 4) +
        theme_minimal(base_size = 14)
)
dev.off()
cat("Saved: tmb_distribution.png\n")

# ---- Step 8: Survival Analysis ----
cat("=== Step 8: Survival Analysis ===\n")

library(survival)
library(survminer)

cached_clin <- file.path(DATA_DIR, "tcga_luad_clinical.rds")

if (file.exists(cached_clin)) {
    clinical <- readRDS(cached_clin)
} else {
    clinical <- GDCquery_clinic("TCGA-LUAD", type = "clinical")
    saveRDS(clinical, cached_clin)
}

clin <- clinical[, c("submitter_id", "vital_status", "days_to_death", "days_to_last_follow_up",
                      "ajcc_pathologic_stage")]
clin$time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death, clin$days_to_last_follow_up)
clin$status <- ifelse(clin$vital_status == "Dead", 1, 0)
clin <- clin[!is.na(clin$time) & clin$time > 0, ]

cat("Patients with survival data:", nrow(clin), "\n")

# 8a: Overall Survival
fit_os <- survfit(Surv(time, status) ~ 1, data = clin)

png(file.path(FIGURES_DIR, "km_overall.png"), width = 8, height = 6, units = "in", res = 300)
print(ggsurvplot(fit_os, data = clin,
           title = "Overall Survival - TCGA-LUAD",
           xlab = "Time (days)", ylab = "Survival Probability",
           risk.table = TRUE, conf.int = TRUE,
           palette = "#2166AC", ggtheme = theme_minimal(base_size = 14)))
dev.off()
cat("Saved: km_overall.png\n")

# 8b: Survival by TMB
tmb$patient_id <- substr(tmb$Tumor_Sample_Barcode, 1, 12)
clin_tmb <- merge(clin, tmb[, c("patient_id", "TMB_class")],
                  by.x = "submitter_id", by.y = "patient_id")

fit_tmb <- survfit(Surv(time, status) ~ TMB_class, data = clin_tmb)

png(file.path(FIGURES_DIR, "km_tmb.png"), width = 8, height = 6, units = "in", res = 300)
print(ggsurvplot(fit_tmb, data = clin_tmb, pval = TRUE,
           title = "Survival by Tumor Mutation Burden - TCGA-LUAD",
           xlab = "Time (days)", ylab = "Survival Probability",
           risk.table = TRUE, conf.int = TRUE,
           palette = c("#2166AC", "#B2182B"),
           legend.labs = c("TMB-High", "TMB-Low"),
           ggtheme = theme_minimal(base_size = 14)))
dev.off()
cat("Saved: km_tmb.png\n")

# 8c: Survival by Key Gene Mutations
key_mut_genes <- c("TP53", "KRAS", "EGFR", "STK11")
cox_results <- data.frame()

for (gene in key_mut_genes) {
    mutated_samples <- unique(substr(
        subsetMaf(maf, genes = gene)@data$Tumor_Sample_Barcode, 1, 12))

    clin[[gene]] <- ifelse(clin$submitter_id %in% mutated_samples, "Mutant", "Wild-Type")

    fit_gene <- survfit(as.formula(paste0("Surv(time, status) ~ ", gene)), data = clin)

    png(file.path(SURV_DIR, paste0("km_", gene, ".png")), width = 8, height = 6, units = "in", res = 300)
    print(ggsurvplot(fit_gene, data = clin, pval = TRUE,
               title = paste0("Survival by ", gene, " Mutation - TCGA-LUAD"),
               xlab = "Time (days)", ylab = "Survival Probability",
               risk.table = TRUE, conf.int = TRUE,
               palette = c("#B2182B", "#2166AC"),
               legend.labs = c("Mutant", "Wild-Type"),
               ggtheme = theme_minimal(base_size = 14)))
    dev.off()
    cat("Saved: km_", gene, ".png\n")

    cox_fit <- coxph(as.formula(paste0("Surv(time, status) ~ ", gene)), data = clin)
    s <- summary(cox_fit)
    cox_results <- rbind(cox_results, data.frame(
        Gene = gene,
        HR = round(s$conf.int[1], 3),
        HR_lower = round(s$conf.int[3], 3),
        HR_upper = round(s$conf.int[4], 3),
        p_value = round(s$coefficients[5], 4)
    ))
}

write.csv(cox_results, file.path(SURV_DIR, "cox_regression_results.csv"), row.names = FALSE)
cat("\nCox Regression Results:\n")
print(cox_results)

# Forest plot
png(file.path(FIGURES_DIR, "forest_plot_hr.png"), width = 8, height = 5, units = "in", res = 300)
print(
    ggplot(cox_results, aes(x = HR, y = reorder(Gene, HR))) +
        geom_point(size = 3, color = "#B2182B") +
        geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2, color = "#B2182B") +
        geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
        labs(title = "Hazard Ratios - Key Gene Mutations in TCGA-LUAD",
             x = "Hazard Ratio (95% CI)", y = "") +
        theme_minimal(base_size = 14)
)
dev.off()
cat("Saved: forest_plot_hr.png\n")

cat("\n============================================\n")
cat("Analysis complete! All outputs in:\n")
cat(RESULTS_DIR, "\n")
cat("============================================\n")
