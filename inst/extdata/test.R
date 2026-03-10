library(patchwork)
library(tidyverse)
# Create parameter grid with 3 values for each parameter
param_grid <- expand.grid(
    edk50 = c(0.3, 0.7),
    kde = c(0.05, 0.15),
    ks = c(0.03, 0.07),
    gamma = c(0.7, 1.3)
)

# Generate data for each parameter combination
results_list <- lapply(seq_len(nrow(param_grid)), function(i) {
    params <- param_grid[i, ]
    model <- kpd_mod(
        edk50 = params$edk50,
        kde = params$kde,
        ks = params$ks,
        gamma = params$gamma
    )
    res <- simulate_steph_curve(model)
    res$edk50 <- params$edk50
    res$kde <- params$kde
    res$ks <- params$ks
    res$gamma <- params$gamma
    res
})

# Combine all results into one data frame
combined_data <- bind_rows(results_list)

# Create facet grid plot
combined_data |>
    mutate(
        edk50 = factor(edk50),
        kde = factor(kde),
        ks = factor(ks),
        gamma = factor(gamma)
    ) |>
    plot_pH_time() +
    facet_grid(edk50 + kde ~ ks + gamma, labeller = "label_both")

### 2nd parameterization

# Generate data for each parameter combination
results_list <- lapply(seq_len(nrow(param_grid)), function(i) {
    params <- param_grid[i, ]
    model <- kpd_mod2(
        edk50 = params$edk50,
        kde = params$kde,
        ks = params$ks,
        gamma = params$gamma
    )
    res <- simulate_steph_curve(model)
    res$edk50 <- params$edk50
    res$kde <- params$kde
    res$ks <- params$ks
    res$gamma <- params$gamma
    res
})

# Combine all results into one data frame
combined_data <- bind_rows(results_list)

# Create facet grid plot
combined_data |>
    mutate(
        edk50 = factor(edk50),
        kde = factor(kde),
        ks = factor(ks),
        gamma = factor(gamma)
    ) |>
    plot_pH_time() +
    facet_grid(edk50 + kde ~ ks + gamma, labeller = "label_both")
