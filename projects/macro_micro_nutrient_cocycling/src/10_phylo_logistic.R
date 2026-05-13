library(ape)
library(phylolm)

args <- commandArgs(trailingOnly = TRUE)
tree_path   <- args[1]  # pruned tree
trait_path  <- args[2]  # species_gene_families.csv
pairs_path  <- args[3]  # CSV defining which pairs to test (gene_x, gene_y, col_x, col_y)
output_path <- args[4]  # output CSV

cat("Loading pruned tree...\n")
tree <- read.tree(tree_path)
cat(sprintf("Tree has %d tips\n", Ntip(tree)))

cat("Loading trait data...\n")
traits <- read.csv(trait_path, stringsAsFactors = FALSE)

acc_from_id <- function(sid) {
  parts <- strsplit(sid, "--", fixed = TRUE)
  sapply(parts, function(x) x[length(x)])
}
traits$accession <- acc_from_id(traits$gtdb_species_clade_id)

rownames(traits) <- traits$accession
traits_matched <- traits[tree$tip.label, ]
stopifnot(nrow(traits_matched) == Ntip(tree))

cat("Loading pair definitions...\n")
pairs <- read.csv(pairs_path, stringsAsFactors = FALSE)
cat(sprintf("Testing %d pairs\n", nrow(pairs)))

results <- data.frame(
  pair_label = character(),
  gene_x = character(),
  gene_y = character(),
  n_species = integer(),
  n_both = integer(),
  uncorrected_logOR = numeric(),
  uncorrected_p = numeric(),
  phylo_alpha = numeric(),
  phylo_logOR = numeric(),
  phylo_se = numeric(),
  phylo_z = numeric(),
  phylo_p = numeric(),
  phylo_ci_low = numeric(),
  phylo_ci_high = numeric(),
  converged = logical(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(pairs)) {
  label <- pairs$label[i]
  gx    <- pairs$gene_x[i]
  gy    <- pairs$gene_y[i]
  cx    <- pairs$col_x[i]
  cy    <- pairs$col_y[i]

  x <- as.numeric(traits_matched[[cx]])
  y <- as.numeric(traits_matched[[cy]])
  names(x) <- tree$tip.label
  names(y) <- tree$tip.label

  n_species <- sum(!is.na(x) & !is.na(y))
  n_both    <- sum(x == 1 & y == 1, na.rm = TRUE)

  tbl <- table(factor(x, levels = 0:1), factor(y, levels = 0:1))
  ft <- fisher.test(tbl)
  uncorr_logOR <- log(ft$estimate)
  uncorr_p     <- ft$p.value

  cat(sprintf("  [%d/%d] %s: n=%d, n_both=%d, uncorr logOR=%.3f ... ",
              i, nrow(pairs), label, n_species, n_both, uncorr_logOR))

  df <- data.frame(y_var = y, x_var = x)
  rownames(df) <- tree$tip.label

  tryCatch({
    fit <- phyloglm(y_var ~ x_var, data = df, phy = tree,
                    method = "logistic_MPLE",
                    btol = 50, log.alpha.bound = 8)

    coef_x <- summary(fit)$coefficients["x_var", ]
    alpha   <- fit$alpha
    logOR   <- coef_x["Estimate"]
    se      <- coef_x["StdErr"]
    z_val   <- coef_x["z.value"]
    p_val   <- coef_x["p.value"]
    ci_lo   <- logOR - 1.96 * se
    ci_hi   <- logOR + 1.96 * se
    conv    <- fit$convergence == 0

    cat(sprintf("phylo logOR=%.3f (p=%.2e), alpha=%.4f\n", logOR, p_val, alpha))

    results <- rbind(results, data.frame(
      pair_label = label, gene_x = gx, gene_y = gy,
      n_species = n_species, n_both = n_both,
      uncorrected_logOR = uncorr_logOR, uncorrected_p = uncorr_p,
      phylo_alpha = alpha, phylo_logOR = logOR, phylo_se = se,
      phylo_z = z_val, phylo_p = p_val,
      phylo_ci_low = ci_lo, phylo_ci_high = ci_hi,
      converged = conv
    ))
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    results <<- rbind(results, data.frame(
      pair_label = label, gene_x = gx, gene_y = gy,
      n_species = n_species, n_both = n_both,
      uncorrected_logOR = uncorr_logOR, uncorrected_p = uncorr_p,
      phylo_alpha = NA, phylo_logOR = NA, phylo_se = NA,
      phylo_z = NA, phylo_p = NA,
      phylo_ci_low = NA, phylo_ci_high = NA,
      converged = FALSE
    ))
  })
}

write.csv(results, output_path, row.names = FALSE)
cat(sprintf("\nResults saved to %s\n", output_path))
cat("\nSummary:\n")
print(results[, c("pair_label", "uncorrected_logOR", "phylo_logOR", "phylo_p", "converged")])
