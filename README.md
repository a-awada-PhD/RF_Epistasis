# Custom Epistasis Pipeline

# Preprocess and filter files using PLINK to QC standards you require in bash

# Rscript

library(data.table)
library(ranger)

geno <- fread("EUR_rf.raw")

pheno <- fread("EUR_pheno.txt")
colnames(pheno) <- c("FID", "IID", "VL")

sex <- fread("EUR_sex.txt", fill = Inf)
colnames(sex) <- c("FID", "IID", "Sex")

pcs <- fread("EUR_pca.eigenvec")
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:5))

# Merge
df <- merge(geno, pheno, by = c("FID", "IID"))
df <- merge(df, sex, by = c("FID", "IID"))
df <- merge(
  df,
  pcs[, .(FID, IID, PC1, PC2, PC3, PC4, PC5)],
  by = c("FID", "IID")
)

fwrite(df[, .(FID)], "EUR_FIDs.txt", col.names = FALSE)

# Check who still has VL but is not in genotype
geno_ids <- unique(geno[, .(FID, IID)])
pheno_vl_ids <- unique(pheno[!is.na(VL), .(FID, IID)])
pheno_with_vl_not_in_geno <- fsetdiff(pheno_vl_ids, geno_ids)
cat("VL people not in genotype after ID fix:", nrow(pheno_with_vl_not_in_geno), "\n")
head(pheno_with_vl_not_in_geno)

# Remove unnecessary PLINK columns if present
drop_cols <- intersect(c("PAT", "MAT", "SEX", "PHENOTYPE"), names(df))
if (length(drop_cols) > 0) {
  df[, (drop_cols) := NULL]
}

df[, VL := as.numeric(VL)]
df[, Sex := as.numeric(Sex)]

# Keep only complete cases for outcome and sex
df <- df[!is.na(VL) & !is.na(Sex)]

# Exclude non-SNP columns
snp_cols <- setdiff(names(df), c("FID", "IID", "VL", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5"))
rm(geno)

table(df$Sex)

# Keep only SNPs with no missing values
keep_snps <- snp_cols[colSums(is.na(df[, ..snp_cols])) == 0]

df <- df[, c("FID", "IID", "VL", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5", keep_snps), with = FALSE]

# Remove SNPs with no variation
snp_cols <- setdiff(names(df), c("FID", "IID", "VL", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5"))
variable_snps <- snp_cols[sapply(df[, ..snp_cols], function(x) length(unique(x)) > 1)]
df <- df[, c("FID", "IID", "VL", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5", variable_snps), with = FALSE]

# Random forest with averaged impurity importance
rf_df <- copy(df[, !c("FID", "IID"), with = FALSE])

n_runs <- 50
importance_sum <- NULL

for (i in seq_len(n_runs)) {
  
  cat("Impurity RF run", i, "\n")
  gc()
  
  rf <- ranger(
    dependent.variable.name = "VL",
    data = rf_df,
    num.trees = 50,
    importance = "impurity",
    seed = i,
    num.threads = 1,
    write.forest = FALSE,
    save.memory = TRUE,
    mtry = 100,
    min.node.size = 10
  )
  
  if (is.null(importance_sum)) {
    importance_sum <- rf$variable.importance
  } else {
    importance_sum <- importance_sum + rf$variable.importance
  }
  
  rm(rf)
  gc()
}

# average impurity importance
importance_avg <- importance_sum / n_runs

imp <- data.table(
  SNP = names(importance_avg),
  Importance = as.numeric(importance_avg)
)

imp <- imp[order(-Importance)]

fwrite(
  imp,
  "EUR_rf_full_5PCs.txt",
  sep = "\t"
)

imp <- fread("EUR_rf_full_5PCs.txt")
imp <- imp[!SNP %in% c("Sex", "PC1", "PC2", "PC3", "PC4", "PC5")]

# Keep top SNPs for permutation step
top_n <- 1000
top_snps <- imp$SNP[1:min(top_n, nrow(imp))]

df_perm <- df[, c("FID", "IID", "VL", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5", top_snps), with = FALSE]
rm(df)
#--------------------------------------------------
# Random forest with averaged permutation importance
perm_sum <- NULL
n_runs <- 10

for (i in 1:n_runs) {
  
  cat("Permutation RF run", i, "\n")
  gc()
  
  rf <- ranger(
    VL ~ . - FID - IID,
    data = df_perm,
    num.trees = 500,
    importance = "permutation",
    seed = i,
    num.threads = 1,
    write.forest = FALSE,
    save.memory = TRUE,
    mtry = 100,
    min.node.size = 10
  )
  
  perm_vec <- rf$variable.importance
  
  if (is.null(perm_sum)) {
    perm_sum <- perm_vec
  } else {
    perm_sum <- perm_sum + perm_vec
  }
  
  rm(rf)
  gc()
}

# average permutation importance
perm_avg <- perm_sum / n_runs

imp_final <- data.table(
  SNP = names(perm_avg),
  Importance = as.numeric(perm_avg)
)

imp_final <- imp_final[order(-Importance)]

fwrite(
  imp_final,
  "EUR_rf_importance_0.1_5PCs.txt",
  sep = "\t"
)
imp_final <- imp_final[!SNP %in% c("Sex", "PC1", "PC2", "PC3", "PC4", "PC5")]

# Select top SNPs
n_top <- min(500, nrow(imp_final))
top_snps <- imp_final$SNP[1:n_top]

# -----------------------------
# THIS EXAMPLE USES HLA-B*57 DOSAGE AS A COVARIATE, YOU MAY CHOOSE TO EXCLUDE IN PAIRWISE TESTING

# LOAD B57 DOSAGE DATA, CALCULATE IT BEFOREHAND USING SNP2HLA 
# -----------------------------

B57_dose <- fread("B57_dosage.txt")

# Fix IDs in B57 dosage file
B57_dose[, FID := fix_double_id(FID)]
B57_dose[, IID := fix_double_id(IID)]

# Rename the dosage column to B57_dose for clarity
# Adjust the column name if yours differs from "HLA_B*57_T"
setnames(B57_dose, "HLA_B*57_T", "B57_dose", skip_absent = TRUE)

# Merge B57 dosage into df_perm
df_perm <- merge(df_perm, B57_dose[, .(FID, IID, B57_dose)], by = c("FID", "IID"), all.x = TRUE)

# Check merge result
cat("Samples with B57 dosage:", sum(!is.na(df_perm$B57_dose)), "/", nrow(df_perm), "\n")

# -----------------------------
# Pairwise epistasis testing with Sex + PCs + B57_dose in final model
# -----------------------------
results <- list()
k <- 1

if (length(top_snps) >= 2) {
  for (i in 1:(length(top_snps) - 1)) {
    for (j in (i + 1):length(top_snps)) {
      
      snp1 <- top_snps[i]
      snp2 <- top_snps[j]
      
      # Skip if either SNP is missing from df for any reason
      if (!(snp1 %in% names(df_perm)) || !(snp2 %in% names(df_perm))) next
      
      # Skip invariant SNPs just to be safe
      if (length(unique(df_perm[[snp1]])) < 2 || length(unique(df_perm[[snp2]])) < 2) next
      
      # Include Sex, PCs, AND B57_dose in the interaction model
      formula <- as.formula(
        paste("VL ~ Sex + PC1 + PC2 + PC3 + PC4 + PC5 + B57_dose +", snp1, "*", snp2)
      )
      
      model <- tryCatch(
        lm(formula, data = df_perm),
        error = function(e) NULL
      )
      
      pval <- NA_real_
      
      if (!is.null(model)) {
        coef_tab <- summary(model)$coefficients
        
        interaction1 <- paste0(snp1, ":", snp2)
        interaction2 <- paste0(snp2, ":", snp1)
        
        if (interaction1 %in% rownames(coef_tab)) {
          pval <- coef_tab[interaction1, "Pr(>|t|)"]
        } else if (interaction2 %in% rownames(coef_tab)) {
          pval <- coef_tab[interaction2, "Pr(>|t|)"]
        }
      }
      
      results[[k]] <- data.table(
        SNP1 = snp1,
        SNP2 = snp2,
        Interaction_P = pval
      )
      
      k <- k + 1
    }
  }
}

results <- rbindlist(results, fill = TRUE)

# Save results
fwrite(
  results,
  "EUR_epistasis_B57_conditioned.txt",
  sep = "\t"
)

results[, FDR := p.adjust(Interaction_P, method = "fdr")]

# Sort by raw p-value
results <- results[order(Interaction_P)]

bim <- fread("EUR_final.bim", header = FALSE)
colnames(bim) <- c("CHR", "SNP_ID", "CM", "BP", "A1", "A2")

# Remove allele suffix from PLINK .raw SNP names
results[, SNP1_base := sub("_[^_]+$", "", SNP1)]
results[, SNP2_base := sub("_[^_]+$", "", SNP2)]

# Match SNPs to BIM file
results[, CHR1 := bim$CHR[match(SNP1_base, bim$SNP_ID)]]
results[, CHR2 := bim$CHR[match(SNP2_base, bim$SNP_ID)]]

# Replace unmatched values
results[is.na(CHR1), CHR1 := "N/A"]
results[is.na(CHR2), CHR2 := "N/A"]

# Reorder columns
setcolorder(
  results,
  c("SNP1", "SNP2", "Interaction_P", "FDR", "CHR1", "CHR2",
    setdiff(names(results), c("SNP1", "SNP2", "Interaction_P", "FDR", "CHR1", "CHR2")))
)

# Remove helper columns
results[, c("SNP1_base", "SNP2_base") := NULL]

# Save results
fwrite(
  results,
  "EUR_epistasis_5PCs_B57_Conditioned.txt",
  sep = "\t"
)




