# ENRTsensitivity

[![](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
[![](https://img.shields.io/badge/status-dev-orange.svg)]()

R package for **S**ensitivity **A**nalysis for **C**ontamination in **E**gocentric-**N**etwork **R**andomized **T**rials.

[cite_start]This package implements the sensitivity analysis (SA) framework developed in Weinstein & Nevo (2025)[cite: 30, 31].

## Overview

[cite_start]Egocentric-network randomized trials (ENRTs) are a powerful design for estimating causal effects in the presence of network interference[cite: 18, 20]. [cite_start]However, the egocentric sampling process can lead to **network contamination**, where latent ego-ego or alter-ego edges are missing from the observed data[cite: 22, 23, 118].

[cite_start]This contamination can bias the standard estimators for the **Indirect Effect (IE)** and **Direct Effect (DE)**[cite: 29, 131, 135, 140].

The `ENRTsensitivity` package provides a framework to assess the robustness of causal estimates to such contamination. [cite_start]It implements the bias-corrected estimators from the paper [cite: 162, 190] and provides two methods for analysis:

1.  [cite_start]**Grid Sensitivity Analysis (GSA):** Systematically compute estimates over a grid of sensitivity parameters[cite: 268].
2.  [cite_start]**Probabilistic Bias Analysis (PBA):** Compute the distribution of bias-corrected estimates based on prior distributions for the sensitivity parameters[cite: 270].

## Installation

You can install the development version of `ENRTsensitivity` from [GitHub](httpsCongratulations on getting this far! Here is a complete README.md file for your package.

This file is written in Markdown (that's what the .md means). You should:

In your RStudio project, go to File > New File > Text File.
Copy and paste the entire contents of the block below into that new file.
Save the file in the root directory of your package (the ENRT_sensitivity/ folder, not the R/ folder) with the exact name README.md.

After you push this file to GitHub, it will appear as the homepage for your repository.

---
ACTION: Copy all the text below this line.

# ENRTsensitivity

[![](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md) [![](https://img.shields.io/badge/status-dev-orange.svg)]()

R package for **S**ensitivity **A**nalysis for **C**ontamination in **E**gocentric-**N**etwork **R**andomized **T**rials.

[cite_start]This package implements the sensitivity analysis (SA) framework developed in Weinstein & Nevo (2025)[cite: 30, 31].

## Overview

[cite_start]Egocentric-network randomized trials (ENRTs) are a powerful design for estimating causal effects in the presence of network interference[cite: 18, 20]. [cite_start]However, the egocentric sampling process can lead to **network contamination**, where latent ego-ego or alter-ego edges are missing from the observed data[cite: 22, 23, 118].

[cite_start]This contamination can bias the standard estimators for the **Indirect Effect (IE)** and **Direct Effect (DE)**[cite: 29, 131, 135, 140].

The `ENRTsensitivity` package provides a framework to assess the robustness of causal estimates to such contamination. [cite_start]It implements the bias-corrected estimators from the paper [cite: 162, 190] and provides two methods for analysis:

1.  [cite_start]**Grid Sensitivity Analysis (GSA):** Systematically compute estimates over a grid of sensitivity parameters[cite: 268].
2.  [cite_start]**Probabilistic Bias Analysis (PBA):** Compute the distribution of bias-corrected estimates based on prior distributions for the sensitivity parameters[cite: 270].

## Installation

You can install the development version of `ENRTsensitivity` from [GitHub](httpshttps://github.com/YOUR-USERNAME/ENRT_sensitivity) with:

```r
# install.packages("devtools")
devtools::install_github("YOUR-USERNAME/ENRT_sensitivity")
```
**Note:** Remember to replace `YOUR-USERNAME` with your actual GitHub username.

## Quick Start: Grid Sensitivity Analysis (GSA)

Here is a minimal example of performing a GSA using the `enrt_sa()` function.

```r
library(ENRTsensitivity)

# 1. Simulate minimal data
set.seed(123)
n_e <- 50
n_a <- 100
pz <- 0.5

Z_e <- rbinom(n_e, 1, pz)
Y_e <- rbinom(n_e, 1, 0.2 + 0.1 * Z_e)

ego_id_a <- sample(1:n_e, n_a, replace = TRUE)
F_a <- Z_e[ego_id_a] # Observed exposure
Y_a <- rbinom(n_a, 1, 0.3 - 0.1 * F_a)

# 2. Define Sensitivity Parameter Grids
# [cite_start]We'll use the homogeneous model (Example 2 from the paper) [cite: 232]

# IE: Grid of expected missing alter-ego edges (m_a)
m_a_grid <- c(0, 50, 100)
pi_list_ae_homo <- pi_homo(m_vec = m_a_grid, n_e = n_e, n_a = n_a,
                           type = "alter", pz = pz)

# DE: Grid of expected missing ego-ego edges (m_e)
m_e_grid <- c(0, 25, 50)
pi_list_ee_homo <- pi_homo(m_vec = m_e_grid, n_e = n_e,
                           type = "ego", pz = pz)

# [cite_start]DE: Grid for kappa (treatment-exposure interaction) [cite: 180]
kappa_grid <- c(1, 1.5, 2)

# 3. Format Parameter Lists for enrt_sa
# The names "Homo" and "Hetero" will be used in the plots
pi_lists_ae <- list(Homo = pi_list_ae_homo)
pi_lists_ee <- list(Homo = pi_list_ee_homo)

# 4. Run the Grid Sensitivity Analysis (GSA)
# We use the simple, non-augmented estimators (no outcome models)
sa_results <- enrt_sa(
  Y_e = Y_e, Y_a = Y_a,
  Z_e = Z_e, F_a = F_a, ego_id_a = ego_id_a,
  pi_lists_ego_ego = pi_lists_ee,
  pi_lists_alter_ego = pi_lists_ae,
  kappa_vec = kappa_grid,
  pz = pz,
  bootstrap = FALSE # Use fast empirical variance
)

# 5. View Results
# The naive estimate (no contamination)
print(sa_results$null_results$IE)

# The full grid of bias-corrected estimates
print(sa_results$sa_results$IE)

# Show the plot for Indirect Effects
# [cite_start]This recreates a plot similar to Figure 2 in the paper [cite: 379]
print(sa_results$ie_rd_plot)

# Show the plot for Direct Effects
# [cite_start]This recreates a plot similar to Figure 3 in the paper [cite: 406]
print(sa_results$de_rd_plot)
```

## Probabilistic Bias Analysis (PBA)

You can also perform a PBA using the `enrt_pba()` function by defining prior distributions for the sensitivity parameters.

```r
# 1. Define Prior Sampling Functions

# Prior for IE: Poisson(50) for m_a
prior_ie_func <- function() {
  round(rpois(1, 50))
}

# Prior for DE: Poisson(25) for m_e and Uniform(1, 3) for kappa
prior_de_func <- function() {
  list(
    pi_param = round(rpois(1, 25)), # Sampled m_e
    kappa = runif(1, 1, 3)           # Sampled kappa
  )
}

# 2. Define Base Arguments for pi_homo()
pi_args_ie_list <- list(n_e = n_e, n_a = n_a, type = "alter", pz = pz)
pi_args_de_list <- list(n_e = n_e, type = "ego", pz = pz)

# 3. Run the PBA
pba_results <- enrt_pba(
  Y_e = Y_e, Y_a = Y_a,
  Z_e = Z_e, F_a = F_a, ego_id_a = ego_id_a,
  bootstrap = FALSE, # Use Normal Approximation
  B = 500,           # Number of Monte Carlo samples
  pz = pz,
  prior_func_ie = prior_ie_func,
  pi_func_ie = pi_homo,
  pi_args_ie = pi_args_ie_list,
  pi_param_name_ie = "m_vec", # pi_func_ie's arg to sample
  prior_func_de = prior_de_func,
  pi_func_de = pi_homo,
  pi_args_de = pi_args_de_list,
  pi_param_name_de = "m_vec" # pi_func_de's arg to sample
)

# 4. View Results
# Shows 95% intervals for bias-only and total uncertainty
print(pba_results$IE_results)
print(pba_results$DE_results)
```

## Citation

If you use this package in your research, please cite the accompanying paper:

> Weinstein, B. and Nevo, D. (2025). Sensitivity analysis for contamination in egocentric-network randomized trials with interference. *[Journal Name or "Working Paper, Tel Aviv University"]*.

```bibtex
@article{weinstein2025enrt,
  title   = {Sensitivity analysis for contamination in egocentric-network randomized trials with interference},
  author  = {Weinstein, Bar and Nevo, Daniel},
  year    = {2025},
  journal = {Working Paper}
}
```

## License

MIT
) with:

```r
# install.packages("devtools")
devtools::install_github("YOUR-USERNAME/ENRT_sensitivity")
```
**Note:** Remember to replace `YOUR-USERNAME` with your actual GitHub username.

## Quick Start: Grid Sensitivity Analysis (GSA)

Here is a minimal example of performing a GSA using the `enrt_sa()` function.

```r
library(ENRTsensitivity)

# 1. Simulate minimal data
set.seed(123)
n_e <- 50
n_a <- 100
pz <- 0.5

Z_e <- rbinom(n_e, 1, pz)
Y_e <- rbinom(n_e, 1, 0.2 + 0.1 * Z_e)

ego_id_a <- sample(1:n_e, n_a, replace = TRUE)
F_a <- Z_e[ego_id_a] # Observed exposure
Y_a <- rbinom(n_a, 1, 0.3 - 0.1 * F_a)

# 2. Define Sensitivity Parameter Grids
# [cite_start]We'll use the homogeneous model (Example 2 from the paper) [cite: 232]

# IE: Grid of expected missing alter-ego edges (m_a)
m_a_grid <- c(0, 50, 100)
pi_list_ae_homo <- pi_homo(m_vec = m_a_grid, n_e = n_e, n_a = n_a,
                           type = "alter", pz = pz)

# DE: Grid of expected missing ego-ego edges (m_e)
m_e_grid <- c(0, 25, 50)
pi_list_ee_homo <- pi_homo(m_vec = m_e_grid, n_e = n_e,
                           type = "ego", pz = pz)

# [cite_start]DE: Grid for kappa (treatment-exposure interaction) [cite: 180]
kappa_grid <- c(1, 1.5, 2)

# 3. Format Parameter Lists for enrt_sa
# The names "Homo" and "Hetero" will be used in the plots
pi_lists_ae <- list(Homo = pi_list_ae_homo)
pi_lists_ee <- list(Homo = pi_list_ee_homo)

# 4. Run the Grid Sensitivity Analysis (GSA)
# We use the simple, non-augmented estimators (no outcome models)
sa_results <- enrt_sa(
  Y_e = Y_e, Y_a = Y_a,
  Z_e = Z_e, F_a = F_a, ego_id_a = ego_id_a,
  pi_lists_ego_ego = pi_lists_ee,
  pi_lists_alter_ego = pi_lists_ae,
  kappa_vec = kappa_grid,
  pz = pz,
  bootstrap = FALSE # Use fast empirical variance
)

# 5. View Results
# The naive estimate (no contamination)
print(sa_results$null_results$IE)

# The full grid of bias-corrected estimates
print(sa_results$sa_results$IE)

# Show the plot for Indirect Effects
# [cite_start]This recreates a plot similar to Figure 2 in the paper [cite: 379]
print(sa_results$ie_rd_plot)

# Show the plot for Direct Effects
# [cite_start]This recreates a plot similar to Figure 3 in the paper [cite: 406]
print(sa_results$de_rd_plot)
```

## Probabilistic Bias Analysis (PBA)

You can also perform a PBA using the `enrt_pba()` function by defining prior distributions for the sensitivity parameters.

```r
# 1. Define Prior Sampling Functions

# Prior for IE: Poisson(50) for m_a
prior_ie_func <- function() {
  round(rpois(1, 50))
}

# Prior for DE: Poisson(25) for m_e and Uniform(1, 3) for kappa
prior_de_func <- function() {
  list(
    pi_param = round(rpois(1, 25)), # Sampled m_e
    kappa = runif(1, 1, 3)           # Sampled kappa
  )
}

# 2. Define Base Arguments for pi_homo()
pi_args_ie_list <- list(n_e = n_e, n_a = n_a, type = "alter", pz = pz)
pi_args_de_list <- list(n_e = n_e, type = "ego", pz = pz)

# 3. Run the PBA
pba_results <- enrt_pba(
  Y_e = Y_e, Y_a = Y_a,
  Z_e = Z_e, F_a = F_a, ego_id_a = ego_id_a,
  bootstrap = FALSE, # Use Normal Approximation
  B = 500,           # Number of Monte Carlo samples
  pz = pz,
  prior_func_ie = prior_ie_func,
  pi_func_ie = pi_homo,
  pi_args_ie = pi_args_ie_list,
  pi_param_name_ie = "m_vec", # pi_func_ie's arg to sample
  prior_func_de = prior_de_func,
  pi_func_de = pi_homo,
  pi_args_de = pi_args_de_list,
  pi_param_name_de = "m_vec" # pi_func_de's arg to sample
)

# 4. View Results
# Shows 95% intervals for bias-only and total uncertainty
print(pba_results$IE_results)
print(pba_results$DE_results)
```

## Citation

If you use this package in your research, please cite the accompanying paper:

> Weinstein, B. and Nevo, D. (2025). Sensitivity analysis for contamination in egocentric-network randomized trials with interference. *[Journal Name or "Working Paper, Tel Aviv University"]*.

```bibtex
@article{weinstein2025enrt,
  title   = {Sensitivity analysis for contamination in egocentric-network randomized trials with interference},
  author  = {Weinstein, Bar and Nevo, Daniel},
  year    = {2025},
  journal = {Working Paper}
}
```

## License

MIT
