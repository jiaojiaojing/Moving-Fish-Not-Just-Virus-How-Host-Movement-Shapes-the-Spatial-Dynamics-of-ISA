# APPENDIX FIGURE 8 — HOST MOVEMENT AND STOCKING INTENSITY
# Extracted from Experiment 22 of the historical 2026-10-04 analysis script.
# Related main manuscript: Figure 3.
# Historical plot exports are separate outcome heatmaps; assembly into the
# manuscript three-panel layout is NOT performed here.
# This script is not yet runtime-validated against the final PDF.

# ==============================================================================
# ISA two-patch mechanistic analysis plan
# Human-mediated host movement (M^H) vs waterborne pathogen movement (M^P)
#
# Recreated and updated version
#
# Model:
#   dS_i/dt = Lambda_i - beta_i S_i P_i - (mu_i + h_i) S_i
#              - sum_j M^H_ij S_i + sum_j M^H_ji S_j
#
#   dI_i/dt = beta_i S_i P_i - (alpha_i + mu_i + h_i) I_i
#              - sum_j M^H_ij I_i + sum_j M^H_ji I_j
#
#   dP_i/dt = chi_i I_i - delta_i P_i
#              - sum_j M^P_ij P_i + sum_j M^P_ji P_j
#
# Two-patch convention:
#   MH12 = human-mediated fish movement rate from patch 1 -> patch 2
#   MH21 = human-mediated fish movement rate from patch 2 -> patch 1
#   MP12 = waterborne pathogen movement rate from patch 1 -> patch 2
#   MP21 = waterborne pathogen movement rate from patch 2 -> patch 1
#
# IMPORTANT:
# - Lambda is interpreted as external stocking / host input, not natural birth.
# - h is routine harvest/background management removal.
# - alpha is disease-related mortality/removal.
# - MH acts non-selectively on both S and I.
# - Internal movement cancels at the total-host level, but patch depletion can
#   still occur. This updated version explicitly tracks N1 and N2.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Packages and global settings
# ------------------------------------------------------------------------------

required_packages <- c("deSolve", "ggplot2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(deSolve)
library(ggplot2)

options(stringsAsFactors = FALSE)

set.seed(123)


# ------------------------------------------------------------------------------
# 1. Reference parameterization
# ------------------------------------------------------------------------------

# These values are placeholders / working reference values from the earlier plan.
# Replace with final empirically justified values before publication.
cycle_months <- 24
Lambda_ref <- mean(c(493233, 1108106))/ cycle_months  # external stocking / host input
mu_ref     <- 0.0127                     # natural mortality rate
h_ref      <- 1/cycle_months                          # routine harvest; placeholder until estimated
beta_ref   <- mean(c(1e-8, 1e-6))        # environment-to-fish transmission
alpha_ref  <- 1.43  
V_effect <- 20000 # disease-related removal; VERIFY provenance/units
chi_ref    <- mean(c(1e4, 1e6))/V_effect        # pathogen shedding
delta_ref  <- mean(c(1.5, 10.5))         # environmental pathogen decay/loss

MP12_ref   <- 0.05
MP21_ref   <- 0.01

MH_ref     <- 0.05

# Simulation horizon
t_end      <- 120
dt         <- 0.1
times      <- seq(0, t_end, by = dt)

# Small positive value for numerical checks
eps <- 1e-12


# ------------------------------------------------------------------------------
# 2. Parameter constructor
# ------------------------------------------------------------------------------

make_pars <- function(
  Lambda1 = Lambda_ref,
  Lambda2 = Lambda_ref,
  beta1   = beta_ref,
  beta2   = beta_ref,
  mu1     = mu_ref,
  mu2     = mu_ref,
  h1      = h_ref,
  h2      = h_ref,
  alpha1  = alpha_ref,
  alpha2  = alpha_ref,
  chi1    = chi_ref,
  chi2    = chi_ref,
  delta1  = delta_ref,
  delta2  = delta_ref,
  MH12    = MH_ref,
  MH21    = 0,
  MP12    = MP12_ref,
  MP21    = MP21_ref
) {
  c(
    Lambda1 = Lambda1, Lambda2 = Lambda2,
    beta1 = beta1, beta2 = beta2,
    mu1 = mu1, mu2 = mu2,
    h1 = h1, h2 = h2,
    alpha1 = alpha1, alpha2 = alpha2,
    chi1 = chi1, chi2 = chi2,
    delta1 = delta1, delta2 = delta2,
    MH12 = MH12, MH21 = MH21,
    MP12 = MP12, MP21 = MP21
  )
}


# ------------------------------------------------------------------------------
# 3. ODE system
# ------------------------------------------------------------------------------

isa_two_patch <- function(t, state, pars) {
  with(as.list(c(state, pars)), {

    # Infection incidence
    inc1 <- beta1 * S1 * P1
    inc2 <- beta2 * S2 * P2

    # Host equations
    dS1 <- Lambda1 - inc1 - (mu1 + h1) * S1 - MH12 * S1 + MH21 * S2
    dS2 <- Lambda2 - inc2 - (mu2 + h2) * S2 - MH21 * S2 + MH12 * S1

    dI1 <- inc1 - (alpha1 + mu1 + h1) * I1 - MH12 * I1 + MH21 * I2
    dI2 <- inc2 - (alpha2 + mu2 + h2) * I2 - MH21 * I2 + MH12 * I1

    # Environmental pathogen equations
    dP1 <- chi1 * I1 - delta1 * P1 - MP12 * P1 + MP21 * P2
    dP2 <- chi2 * I2 - delta2 * P2 - MP21 * P2 + MP12 * P1

    list(
      c(dS1, dI1, dP1, dS2, dI2, dP2),
      incidence1 = inc1,
      incidence2 = inc2
    )
  })
}


# ------------------------------------------------------------------------------
# 4. Initial conditions
# ------------------------------------------------------------------------------

make_init <- function(
  S1 = 1e6,
  I1 = 100,
  P1 = 1e5,
  S2 = 1e6,
  I2 = 10,
  P2 = 1e4
) {
  c(S1 = S1, I1 = I1, P1 = P1, S2 = S2, I2 = I2, P2 = P2)
}

init_ref <- make_init()


# ------------------------------------------------------------------------------
# 5. Simulation and summary helpers
# ------------------------------------------------------------------------------

simulate_model <- function(pars,
                           init = init_ref,
                           times_vec = times,
                           method = "lsoda") {

  out <- ode(
    y = init,
    times = times_vec,
    func = isa_two_patch,
    parms = pars,
    method = method
  )

  out <- as.data.frame(out)

  # Explicit host totals and prevalence
  out$N1 <- out$S1 + out$I1
  out$N2 <- out$S2 + out$I2
  out$NT <- out$N1 + out$N2

  out$prev1 <- ifelse(out$N1 > 0, out$I1 / out$N1, NA_real_)
  out$prev2 <- ifelse(out$N2 > 0, out$I2 / out$N2, NA_real_)
  out$prev_total <- ifelse(out$NT > 0,
                           (out$I1 + out$I2) / out$NT,
                           NA_real_)

  out
}


tail_mean <- function(x, fraction = 0.1) {
  n <- length(x)
  k <- max(1, ceiling(n * fraction))
  mean(tail(x, k), na.rm = TRUE)
}


summarize_run <- function(out, init = NULL, tail_fraction = 0.1) {

  if (is.null(init)) {
    init <- c(
      S1 = out$S1[1], I1 = out$I1[1], P1 = out$P1[1],
      S2 = out$S2[1], I2 = out$I2[1], P2 = out$P2[1]
    )
  }

  N1_0 <- init["S1"] + init["I1"]
  N2_0 <- init["S2"] + init["I2"]

  data.frame(
    I1_final = tail(out$I1, 1),
    I2_final = tail(out$I2, 1),
    P1_final = tail(out$P1, 1),
    P2_final = tail(out$P2, 1),
    N1_final = tail(out$N1, 1),
    N2_final = tail(out$N2, 1),
    NT_final = tail(out$NT, 1),
    prev1_final = tail(out$prev1, 1),
    prev2_final = tail(out$prev2, 1),
    prev_total_final = tail(out$prev_total, 1),

    I1_tail_mean = tail_mean(out$I1, tail_fraction),
    I2_tail_mean = tail_mean(out$I2, tail_fraction),
    prev1_tail_mean = tail_mean(out$prev1, tail_fraction),
    prev2_tail_mean = tail_mean(out$prev2, tail_fraction),
    prev_total_tail_mean = tail_mean(out$prev_total, tail_fraction),

    peak_I1 = max(out$I1, na.rm = TRUE),
    peak_I2 = max(out$I2, na.rm = TRUE),
    peak_total_I = max(out$I1 + out$I2, na.rm = TRUE),

    min_N1 = min(out$N1, na.rm = TRUE),
    min_N2 = min(out$N2, na.rm = TRUE),
    N1_retention = min(out$N1, na.rm = TRUE) / max(N1_0, eps),
    N2_retention = min(out$N2, na.rm = TRUE) / max(N2_0, eps)
  )
}


# Patch-depletion flag:
# A prevalence decrease is potentially misleading if a patch loses most hosts.
flag_depletion <- function(summary_df, threshold = 0.25) {
  summary_df$depleted_patch1 <- summary_df$N1_retention < threshold
  summary_df$depleted_patch2 <- summary_df$N2_retention < threshold
  summary_df$any_depletion <- summary_df$depleted_patch1 |
                              summary_df$depleted_patch2
  summary_df
}


# ------------------------------------------------------------------------------
# 6. Disease-free equilibrium (DFE) and R0
# ------------------------------------------------------------------------------

# At the DFE, I1 = I2 = P1 = P2 = 0.
# S1 and S2 satisfy the linear host-balance system:
#
# (mu1+h1+MH12) S1 - MH21 S2 = Lambda1
# -MH12 S1 + (mu2+h2+MH21) S2 = Lambda2

dfe_susceptibles <- function(pars) {

  p <- as.list(pars)

  A <- matrix(
    c(
      p$mu1 + p$h1 + p$MH12, -p$MH21,
      -p$MH12, p$mu2 + p$h2 + p$MH21
    ),
    nrow = 2,
    byrow = TRUE
  )

  b <- c(p$Lambda1, p$Lambda2)

  Sstar <- solve(A, b)

  names(Sstar) <- c("S1_DFE", "S2_DFE")
  Sstar
}


# Next-generation matrix for infected-state vector:
# x = (I1, I2, P1, P2)
#
# F contains only newly infected fish:
#   new I1 = beta1 * S1_DFE * P1
#   new I2 = beta2 * S2_DFE * P2
#
# V contains all other transfers/removals among infected compartments.

calc_R0 <- function(pars) {

  p <- as.list(pars)
  Sstar <- dfe_susceptibles(pars)

  S1s <- Sstar["S1_DFE"]
  S2s <- Sstar["S2_DFE"]

  gamma1 <- p$alpha1 + p$mu1 + p$h1
  gamma2 <- p$alpha2 + p$mu2 + p$h2

  F <- matrix(0, nrow = 4, ncol = 4)

  # x order = I1, I2, P1, P2
  F[1, 3] <- p$beta1 * S1s
  F[2, 4] <- p$beta2 * S2s

  V <- matrix(
    c(
      gamma1 + p$MH12, -p$MH21,                0,              0,
      -p$MH12,            gamma2 + p$MH21,     0,              0,
      -p$chi1,             0,                   p$delta1 + p$MP12, -p$MP21,
       0,                 -p$chi2,             -p$MP12,       p$delta2 + p$MP21
    ),
    nrow = 4,
    byrow = TRUE
  )

  K <- F %*% solve(V)

  eig <- eigen(K, only.values = TRUE)$values
  R0 <- max(Mod(eig))

  as.numeric(R0)

# ==============================================================================
# EXPERIMENT 22
# Host movement x stocking intensity
#
# Question:
# Does stocking intensity modify the epidemiological effect of
# host movement from the source patch (Patch 1) to the receiving patch
# (Patch 2)?
#
# Design:
#   x-axis: MH12
#   y-axis: r_Lambda = Lambda / Lambda_ref
#
# Stocking intensity is varied equally in both patches:
#
#   Lambda1 = Lambda2 = Lambda_ref * r_Lambda
#
# All other biological parameters remain at their reference values.
# ==============================================================================


# ------------------------------------------------------------------------------
# 22A. Parameter ranges
# ------------------------------------------------------------------------------

# Host movement from source patch (1) to receiving patch (2)
MH12_Lambda_grid <- seq(
  from = 0,
  to = 0.10,
  length.out = 51
)

# Stocking-intensity multiplier
#
# r_Lambda = 1   : reference stocking intensity
# r_Lambda < 1   : lower stocking intensity
# r_Lambda > 1   : higher stocking intensity

r_Lambda_grid <- seq(
  from = 0.25,
  to = 2.0,
  length.out = 51
)


# ------------------------------------------------------------------------------
# 22B. Construct the two-dimensional parameter grid
# ------------------------------------------------------------------------------

MHxLambda_grid <- expand.grid(
  MH12 = MH12_Lambda_grid,
  r_Lambda = r_Lambda_grid
)

MHxLambda_grid$Lambda1 <-
  Lambda_ref * MHxLambda_grid$r_Lambda

MHxLambda_grid$Lambda2 <-
  Lambda_ref * MHxLambda_grid$r_Lambda


# ------------------------------------------------------------------------------
# 22C. Run the full nonlinear model
# ------------------------------------------------------------------------------

MHxLambda_results <- do.call(
  rbind,
  lapply(seq_len(nrow(MHxLambda_grid)), function(k) {
    
    # Current parameter combination
    mh <- MHxLambda_grid$MH12[k]
    rL <- MHxLambda_grid$r_Lambda[k]
    
    Lambda_current <- Lambda_ref * rL
    
    
    # ------------------------------------------------------------
    # Model parameters
    # ------------------------------------------------------------
    
    pars <- make_pars(
      Lambda1 = Lambda_current,
      Lambda2 = Lambda_current,
      
      MH12 = mh,
      MH21 = 0,
      
      beta1 = beta_ref,
      beta2 = beta_ref,
      
      MP12 = MP12_ref,
      MP21 = MP21_ref
    )
    
    
    # ------------------------------------------------------------
    # Simulate nonlinear model
    # ------------------------------------------------------------
    
    out <- simulate_model(
      pars = pars,
      init = init_ref,
      times_vec = times
    )
    
    
    # ------------------------------------------------------------
    # Existing summary statistics
    # ------------------------------------------------------------
    
    s <- summarize_run(
      out = out,
      init = init_ref
    )
    
    s <- flag_depletion(s)
    
    
    # ------------------------------------------------------------
    # Long-term fish abundance
    # ------------------------------------------------------------
    
    NT_tail_mean <- tail_mean(
      out$NT,
      fraction = 0.1
    )
    
    N1_tail_mean <- tail_mean(
      out$N1,
      fraction = 0.1
    )
    
    N2_tail_mean <- tail_mean(
      out$N2,
      fraction = 0.1
    )
    
    
    # ------------------------------------------------------------
    # System-level R0
    # ------------------------------------------------------------
    
    R0_value <- calc_R0(pars)
    
    
    # ------------------------------------------------------------
    # Return results
    # ------------------------------------------------------------
    
    cbind(
      MH12 = mh,
      MH21 = 0,
      
      r_Lambda = rL,
      
      Lambda1 = Lambda_current,
      Lambda2 = Lambda_current,
      
      R0 = R0_value,
      
      NT_tail_mean = NT_tail_mean,
      N1_tail_mean = N1_tail_mean,
      N2_tail_mean = N2_tail_mean,
      
      s
    )
  })
)

rownames(MHxLambda_results) <- NULL


# ------------------------------------------------------------------------------
# 22D. Basic checks
# ------------------------------------------------------------------------------

head(MHxLambda_results)

summary(MHxLambda_results)

cat(
  "\nDuplicated column names:",
  anyDuplicated(names(MHxLambda_results)),
  "\n"
)


# ==============================================================================
# SAVE Experiment 22: MH x stocking intensity
# ==============================================================================

output_dir <- "ISA_two_patch_outputs"

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

saveRDS(
  MHxLambda_results,
  file = file.path(
    output_dir,
    "MHxLambda_results.rds"
  )
)

write.csv(
  MHxLambda_results,
  file = file.path(
    output_dir,
    "MHxLambda_results.csv"
  ),
  row.names = FALSE
)

cat(
  "\nExperiment 22 saved successfully.\n"
)

print(
  file.exists(
    file.path(
      output_dir,
      "MHxLambda_results.rds"
    )
  )
)


# ==============================================================================
# READ Experiment 22
#
# Future sessions can start here.
# There is NO need to rerun the simulation once the RDS file has been saved.
# ==============================================================================

MHxLambda_results <- readRDS(
  file.path(
    output_dir,
    "MHxLambda_results.rds"
  )
)

cat(
  "\nSaved Experiment 22 results loaded successfully.\n",
  "Number of parameter combinations:",
  nrow(MHxLambda_results),
  "\n"
)


# ------------------------------------------------------------------------------
# 22E. Reference parameter combination
#
# MH12    = 0.05
# rLambda = 1
# ------------------------------------------------------------------------------

reference_MHxLambda <- data.frame(
  MH12 = MH_ref,
  r_Lambda = 1
)


# ==============================================================================
# HEATMAP 1
# System-level R0
# ==============================================================================

p_MHxLambda_R0 <- ggplot(
  MHxLambda_results,
  aes(
    x = MH12,
    y = r_Lambda,
    fill = R0
  )
) +
  
  geom_tile() +
  
  scale_fill_gradientn(
    colours = c(
      "navy",
      "dodgerblue",
      "yellow",
      "orange",
      "red"
    )
  ) +
  
  geom_contour(
    aes(z = R0),
    breaks = 1,
    color = "black",
    linewidth = 0.7
  ) +
  
  geom_point(
    data = reference_MHxLambda,
    aes(
      x = MH12,
      y = r_Lambda
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 3.5,
    stroke = 1
  ) +
  
  labs(
    x = expression(M[12]^H),
    
    y = expression(
      r[Lambda] == Lambda / Lambda[ref]
    ),
    
    fill = expression(R[0]),
    
    title = "Host movement and stocking intensity: invasion potential"
  ) +
  
  theme_classic(base_size = 13)

print(p_MHxLambda_R0)


# ==============================================================================
# HEATMAP 2
# Long-term total prevalence
# ==============================================================================

p_MHxLambda_prevT <- ggplot(
  MHxLambda_results,
  aes(
    x = MH12,
    y = r_Lambda,
    fill = prev_total_tail_mean
  )
) +
  
  geom_tile() +
  
  scale_fill_gradientn(
    colours = c(
      "navy",
      "dodgerblue",
      "yellow",
      "orange",
      "red"
    )
  ) +
  
  geom_point(
    data = reference_MHxLambda,
    aes(
      x = MH12,
      y = r_Lambda
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 3.5,
    stroke = 1
  ) +
  
  labs(
    x = expression(M[12]^H),
    
    y = expression(
      r[Lambda] == Lambda / Lambda[ref]
    ),
    
    fill = "Total\nprevalence",
    
    title = "Long-term total disease prevalence"
  ) +
  
  theme_classic(base_size = 13)

print(p_MHxLambda_prevT)


# ==============================================================================
# HEATMAP 3
# Long-term total fish abundance
# ==============================================================================

p_MHxLambda_NT <- ggplot(
  MHxLambda_results,
  aes(
    x = MH12,
    y = r_Lambda,
    fill = NT_tail_mean
  )
) +
  
  geom_tile() +
  
  scale_fill_gradientn(
    colours = c(
      "navy",
      "dodgerblue",
      "yellow",
      "orange",
      "red"
    )
  ) +
  
  geom_point(
    data = reference_MHxLambda,
    aes(
      x = MH12,
      y = r_Lambda
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 3.5,
    stroke = 1
  ) +
  
  labs(
    x = expression(M[12]^H),
    
    y = expression(
      r[Lambda] == Lambda / Lambda[ref]
    ),
    
    fill = "Total fish\nabundance",
    
    title = "Long-term total fish abundance"
  ) +
  
  theme_classic(base_size = 13)

print(p_MHxLambda_NT)


# ==============================================================================
# HEATMAP 4
# Source-patch long-term prevalence
# ==============================================================================

p_MHxLambda_prev1 <- ggplot(
  MHxLambda_results,
  aes(
    x = MH12,
    y = r_Lambda,
    fill = prev1_tail_mean
  )
) +
  
  geom_tile() +
  
  scale_fill_gradientn(
    colours = c(
      "navy",
      "dodgerblue",
      "yellow",
      "orange",
      "red"
    )
  ) +
  
  geom_point(
    data = reference_MHxLambda,
    aes(
      x = MH12,
      y = r_Lambda
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 3.5,
    stroke = 1
  ) +
  
  labs(
    x = expression(M[12]^H),
    
    y = expression(
      r[Lambda] == Lambda / Lambda[ref]
    ),
    
    fill = "Source-patch\nprevalence",
    
    title = "Long-term prevalence in the source patch"
  ) +
  
  theme_classic(base_size = 13)

print(p_MHxLambda_prev1)

##############save and read
# ==============================================================================
# EXPERIMENT 22 DIAGNOSTIC
# Source-patch prevalence at selected stocking intensities
#
# Panel A: r_Lambda = 0.5
# Panel B: r_Lambda = 1.0
#
# Uses the existing MHxLambda_results object.
# NO model rerun is required.
# ==============================================================================


# ------------------------------------------------------------------------------
# Select the nearest simulated stocking-intensity values
# ------------------------------------------------------------------------------

rLambda_A <- MHxLambda_results$r_Lambda[
  which.min(
    abs(MHxLambda_results$r_Lambda - 0.5)
  )
]

rLambda_B <- MHxLambda_results$r_Lambda[
  which.min(
    abs(MHxLambda_results$r_Lambda - 1.0)
  )
]

cat(
  "\nStocking levels used for diagnostic:\n",
  "Panel A =", rLambda_A, "\n",
  "Panel B =", rLambda_B, "\n"
)


# ------------------------------------------------------------------------------
# Extract the two slices
# ------------------------------------------------------------------------------

source_prev_A <- MHxLambda_results[
  abs(MHxLambda_results$r_Lambda - rLambda_A) < 1e-10,
  c(
    "MH12",
    "r_Lambda",
    "prev1_tail_mean"
  )
]

source_prev_B <- MHxLambda_results[
  abs(MHxLambda_results$r_Lambda - rLambda_B) < 1e-10,
  c(
    "MH12",
    "r_Lambda",
    "prev1_tail_mean"
  )
]


# Make sure host movement is ordered correctly
source_prev_A <- source_prev_A[
  order(source_prev_A$MH12),
]

source_prev_B <- source_prev_B[
  order(source_prev_B$MH12),
]


# ------------------------------------------------------------------------------
# Panel A: r_Lambda = 0.5
# ------------------------------------------------------------------------------

p_source_prev_rLambda_05 <- ggplot(
  source_prev_A,
  aes(
    x = MH12,
    y = prev1_tail_mean
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  labs(
    x = expression(M[12]^H),
    y = "Source-patch long-term prevalence",
    title = expression(
      bold("(A)") ~ r[Lambda] == 0.5
    )
  ) +
  
  theme_classic(base_size = 13)


# ------------------------------------------------------------------------------
# Panel B: r_Lambda = 1.0
# ------------------------------------------------------------------------------

p_source_prev_rLambda_10 <- ggplot(
  source_prev_B,
  aes(
    x = MH12,
    y = prev1_tail_mean
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  labs(
    x = expression(M[12]^H),
    y = "Source-patch long-term prevalence",
    title = expression(
      bold("(B)") ~ r[Lambda] == 1.0
    )
  ) +
  
  theme_classic(base_size = 13)


# ------------------------------------------------------------------------------
# Combine panels A and B
# ------------------------------------------------------------------------------

p_source_prev_stocking_AB <-
  p_source_prev_rLambda_05 +
  p_source_prev_rLambda_10 +
  patchwork::plot_layout(
    ncol = 2
  )

print(
  p_source_prev_stocking_AB
)


# ==============================================================================
# SAVE FIGURE
# ==============================================================================

output_dir <- "ISA_two_patch_outputs"

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

ggsave(
  filename = file.path(
    output_dir,
    "MHxLambda_source_patch_prevalence_stocking050_100.png"
  ),
  plot = p_source_prev_stocking_AB,
  width = 10,
  height = 4.5,
  dpi = 300
)

# ==============================================================================
# HEATMAP 5
# Receiving-patch long-term prevalence
# ==============================================================================

p_MHxLambda_prev2 <- ggplot(
  MHxLambda_results,
  aes(
    x = MH12,
    y = r_Lambda,
    fill = prev2_tail_mean
  )
) +
  
  geom_tile() +
  
  scale_fill_gradientn(
    colours = c(
      "navy",
      "dodgerblue",
      "yellow",
      "orange",
      "red"
    )
  ) +
  
  geom_point(
    data = reference_MHxLambda,
    aes(
      x = MH12,
      y = r_Lambda
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 3.5,
    stroke = 1
  ) +
  
  labs(
    x = expression(M[12]^H),
    
    y = expression(
      r[Lambda] == Lambda / Lambda[ref]
    ),
    
    fill = "Receiving-patch\nprevalence",
    
    title = "Long-term prevalence in the receiving patch"
  ) +
  
  theme_classic(base_size = 13)

print(p_MHxLambda_prev2)

######save all output
# ==============================================================================
# SAVE ALL EXPERIMENT 22 FIGURES
# Host movement x stocking intensity
# ==============================================================================

output_dir <- "ISA_two_patch_outputs"

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------------------------
# 1. System-level R0 heatmap
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "MHxLambda_R0.png"
  ),
  plot = p_MHxLambda_R0,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ------------------------------------------------------------------------------
# 2. Long-term total disease prevalence heatmap
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "MHxLambda_total_prevalence.png"
  ),
  plot = p_MHxLambda_prevT,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ------------------------------------------------------------------------------
# 3. Long-term total fish abundance heatmap
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "MHxLambda_total_fish_abundance.png"
  ),
  plot = p_MHxLambda_NT,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ------------------------------------------------------------------------------
# 4. Host-exporting-patch long-term prevalence heatmap
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "MHxLambda_exporting_patch_prevalence.png"
  ),
  plot = p_MHxLambda_prev1,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ------------------------------------------------------------------------------
# 5. Host-receiving-patch long-term prevalence heatmap
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "MHxLambda_receiving_patch_prevalence.png"
  ),
  plot = p_MHxLambda_prev2,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ------------------------------------------------------------------------------
# 6. Diagnostic: exporting-patch prevalence at selected stocking intensities
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "MHxLambda_exporting_patch_prevalence_diagnostic.png"
  ),
  plot = p_source_prev_stocking_AB,
  width = 10,
  height = 4.5,
  dpi = 300,
  bg = "white"
)


# ------------------------------------------------------------------------------
# Confirm saved files
# ------------------------------------------------------------------------------

saved_figures <- c(
  "MHxLambda_R0.png",
  "MHxLambda_total_prevalence.png",
  "MHxLambda_total_fish_abundance.png",
  "MHxLambda_exporting_patch_prevalence.png",
  "MHxLambda_receiving_patch_prevalence.png",
  "MHxLambda_exporting_patch_prevalence_diagnostic.png"
)

cat(
  "\nExperiment 22 figures saved to:\n",
  normalizePath(output_dir),
  "\n\n"
)

print(
  data.frame(
    figure = saved_figures,
    saved = file.exists(
      file.path(
        output_dir,
        saved_figures
      )
    )
  )
)
# ============================================================
