library(ape)
library(phytools)

args <- commandArgs(trailingOnly = TRUE)
tree_path   <- args[1]  # pruned tree (from prior run)
trait_path  <- args[2]  # species_gene_families.csv
tax_path    <- args[3]  # species_taxonomy.csv
output_path <- args[4]  # output CSV
n_subsample <- as.integer(ifelse(length(args) >= 5, args[5], 3000))

stopifnot(
  "Tree file not found"  = file.exists(tree_path),
  "Trait file not found"  = file.exists(trait_path),
  "Taxonomy file not found" = file.exists(tax_path)
)

set.seed(42)

cat("Loading pruned tree...\n")
tree <- read.tree(tree_path)
cat(sprintf("Tree has %d tips\n", Ntip(tree)))

cat("Loading trait data...\n")
traits <- read.csv(trait_path, stringsAsFactors = FALSE)

env_map <- read.csv(file.path(dirname(trait_path), "env_species_mapping.csv"), stringsAsFactors = FALSE)
soil_plant_ids <- env_map$gtdb_species_clade_id[
  env_map$primary_env %in% c("soil/rhizosphere", "plant-associated")]
traits <- traits[traits$gtdb_species_clade_id %in% soil_plant_ids, ]
cat(sprintf("v4 soil+plant filter: %d species\n", nrow(traits)))

acc_from_id <- function(sid) {
  parts <- strsplit(sid, "--", fixed = TRUE)
  sapply(parts, function(x) x[length(x)])
}
traits$accession <- acc_from_id(traits$gtdb_species_clade_id)
rownames(traits) <- traits$accession

common <- intersect(tree$tip.label, traits$accession)
traits_matched <- traits[common, ]

cat("Loading taxonomy for stratified sampling...\n")
tax <- read.csv(tax_path, stringsAsFactors = FALSE)
tax$accession <- acc_from_id(tax$gtdb_species_clade_id)
rownames(tax) <- tax$accession
tax_matched <- tax[common, ]

phyla <- tax_matched$phylum
phylum_tab <- table(phyla)
cat(sprintf("Species across %d phyla\n", length(phylum_tab)))

n_target <- min(n_subsample, length(common))
cat(sprintf("Subsampling %d species (stratified by phylum)...\n", n_target))

sampled <- c()
for (p in names(phylum_tab)) {
  p_tips <- common[phyla == p]
  n_take <- max(1, round(length(p_tips) / length(common) * n_target))
  n_take <- min(n_take, length(p_tips))
  sampled <- c(sampled, sample(p_tips, n_take))
}
if (length(sampled) > n_target) {
  sampled <- sample(sampled, n_target)
}
cat(sprintf("Sampled %d species\n", length(sampled)))

sub_tree <- keep.tip(tree, sampled)
cat(sprintf("Subsampled tree has %d tips\n", Ntip(sub_tree)))

sub_path <- sub("\\.tree$", sprintf("_sub%d.tree", length(sampled)), tree_path)
write.tree(sub_tree, sub_path)

trait_cols <- c(
  "has_phoA", "has_pstA", "has_pstB", "has_pstC", "has_pstS",
  "has_phnC", "has_phnD", "has_phnE",
  "has_nifH", "has_nifD",
  "has_copA", "has_corA",
  "has_phzF", "has_phzA", "has_phzB", "has_phzD", "has_phzG", "has_phzS", "has_phzM",
  "has_phoD_pfam", "has_feoB_pfam", "has_HMA_pfam"
)

gene_labels <- c(
  "phoA", "pstA", "pstB", "pstC", "pstS",
  "phnC", "phnD", "phnE",
  "nifH", "nifD",
  "copA", "corA",
  "phzF", "phzA", "phzB", "phzD", "phzG", "phzS", "phzM",
  "phoD", "feoB", "HMA"
)

group_labels <- c(
  rep("P-acquisition", 8),
  rep("N-fixation", 2),
  rep("Metal-handling", 2),
  rep("Phenazine", 7),
  "P-acquisition", rep("Metal-handling", 2)
)

traits_sub <- traits[sub_tree$tip.label, ]

results <- data.frame(
  gene = character(), group = character(), prevalence = numeric(),
  lambda = numeric(), lambda_p = numeric(), lambda_loglik = numeric(),
  n_subsample = integer(),
  stringsAsFactors = FALSE
)

cat("\nComputing Pagel's lambda for 22 gene families...\n")
for (i in seq_along(trait_cols)) {
  col  <- trait_cols[i]
  gene <- gene_labels[i]
  grp  <- group_labels[i]

  x <- as.numeric(traits_sub[[col]])
  names(x) <- sub_tree$tip.label
  prev <- mean(x, na.rm = TRUE)

  if (prev < 0.005 || prev > 0.995) {
    cat(sprintf("  %s: prevalence %.4f — skipping (near-fixed)\n", gene, prev))
    results <- rbind(results, data.frame(
      gene = gene, group = grp, prevalence = prev,
      lambda = NA, lambda_p = NA, lambda_loglik = NA,
      n_subsample = length(sampled)
    ))
    next
  }

  cat(sprintf("  %s (prevalence=%.3f)... ", gene, prev))

  tryCatch({
    lam_fit <- phylosig(sub_tree, x, method = "lambda", test = TRUE)
    lambda_val <- lam_fit$lambda
    lambda_p   <- lam_fit$P
    lambda_ll  <- lam_fit$logL

    cat(sprintf("lambda=%.4f (p=%.2e)\n", lambda_val, lambda_p))
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    lambda_val <<- NA
    lambda_p   <<- NA
    lambda_ll  <<- NA
  })

  results <- rbind(results, data.frame(
    gene = gene, group = grp, prevalence = prev,
    lambda = lambda_val, lambda_p = lambda_p, lambda_loglik = lambda_ll,
    n_subsample = length(sampled)
  ))
}

write.csv(results, output_path, row.names = FALSE)
cat(sprintf("\nResults saved to %s\n", output_path))
cat("\nSummary:\n")
print(results[, c("gene", "group", "prevalence", "lambda", "lambda_p")])
