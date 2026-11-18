library(data.table)
library(ggplot2)

# ============================================================
# 1. Data
# ============================================================

plot_dt <- final_results[
  !is.na(N_targets) &
    N_targets > 0
]

# Threshold
cutoff <- 20

# Proportion with N_targets >= 20
pct_ge20 <- mean(
  plot_dt$N_targets >= cutoff
) * 100

pct_ge20