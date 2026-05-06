# =============================================================================
# poisson_regression_model.R
# Poisson regression model for mutation bias analysis across 14 bacterial species
#
# Model formula: n_mut ~ mut_class + species * is_igr + offset(log(opp))
# where:
#   n_mut     = observed mutation count per class per category per species
#   mut_class = substitution class (C>A, C>G, C>T, T>A, T>C, T>G)
#   species   = bacterial species (14 levels)
#   is_igr    = TRUE if intergenic, FALSE if coding (synonymous or nonsynonymous)
#   opp       = mutational opportunity (number of available sites on both strands)
#
# The offset(log(opp)) term normalises for differences in opportunity between
# species and mutation classes, converting raw counts to rates.
# The species * is_igr interaction tests whether the coding/intergenic
# rate difference varies across species.
#
# glm() usage: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/glm.html
# Poisson family: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/family.html
# BH FDR correction: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
# =============================================================================

BASE <- "/home/jovyan/shared-team/2025-masters-project/people/alison"

# ── Load data ─────────────────────────────────────────────────────────────────
data <- read.csv(file.path(BASE, "mutation_data_14sp.csv"))

# Derive binary coding/intergenic variable from category column
data$is_igr <- data$category == "intergenic"

# Sense check
head(data)
nrow(data)

# ── Fit Poisson regression model ──────────────────────────────────────────────
model <- glm(n_mut ~ mut_class + species * is_igr + offset(log(opp)),
             data   = data,
             family = poisson)

summary(model)

# ── FDR correction ────────────────────────────────────────────────────────────
# Benjamini-Hochberg correction applied across all model terms.
p_vals <- summary(model)$coefficients[, 4]
p_adj  <- p.adjust(p_vals, method = "BH")

# ── Export results to CSV ─────────────────────────────────────────────────────
# Combines raw coefficients, standard errors, z-values, raw p-values,
# and BH-adjusted p-values into a single results table.
results <- as.data.frame(summary(model)$coefficients)
colnames(results) <- c("estimate", "std_error", "z_value", "p_value")
results$p_adj_BH  <- p_adj[rownames(results)]
results$term      <- rownames(results)
results           <- results[, c("term", "estimate", "std_error", "z_value",
                                 "p_value", "p_adj_BH")]

write.csv(results,
          file      = file.path(BASE, "model_results_14sp.csv"),
          row.names = FALSE)

cat("Model results written to model_results_14sp.csv\n")