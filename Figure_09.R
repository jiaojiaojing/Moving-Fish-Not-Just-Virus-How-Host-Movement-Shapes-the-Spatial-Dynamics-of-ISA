# APPENDIX FIGURE 9 — HOST MOVEMENT AND PATCH-1 TRANSMISSION HETEROGENEITY
# Extracted from Experiment 21 of the historical 2026-10-04 analysis script.
# Related main manuscript: Figure 4.
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
# ==============================================================================
# 21. Combined effects of host movement and patch-1 transmission intensity
# ==============================================================================
#
# Scientific question:
# Does the epidemiological effect of host movement from patch 1 to patch 2
# depend on the local environmental transmission intensity in patch 1?
#
# Mechanistic interpretation:
#
# Susceptible fish in patch 1 can either:
#
#   (1) become infected locally:
#
#         S1 -- beta1 * P1 --> I1
#
#   or
#
#   (2) leave patch 1 through host movement:
#
#         S1 -- MH12 --> S2
#
# Therefore, this experiment examines whether the effect of host export
# depends on the relative strength of local infection pressure in patch 1.
#
# x-axis:
#   MH12 = host movement from patch 1 -> patch 2
#
# y-axis:
#   r_beta1 = beta1 / beta_ref
#
# Patch 2 transmission is fixed:
#   beta2 = beta_ref
#
# Waterborne pathogen movement is fixed at the GLOBAL reference values:
#   MP12 = MP12_ref = 0.05
#   MP21 = MP21_ref = 0.01
#
# Reverse host movement remains zero:
#   MH21 = 0
#
# Main outcomes:
#   1. R0
#   2. Long-term total prevalence
#   3. Long-term total fish abundance
#   4. Patch-1 long-term prevalence
#   5. Patch-2 long-term prevalence
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 21A. Parameter space
# ------------------------------------------------------------------------------

# Host movement from patch 1 -> patch 2.
# Keep the same focused range used in Experiment 20.

MH12_beta1_grid <- seq(
  from = 0,
  to = 0.10,
  length.out = 51
)

# Reverse host movement remains zero.
MH21_beta1_fixed <- 0


# Patch-1 transmission multiplier:
#
#             beta1
# r_beta1 = ----------
#            beta_ref
#
# r_beta1 = 1   : reference transmission rate in patch 1
# r_beta1 < 1   : weaker local transmission in patch 1
# r_beta1 > 1   : stronger local transmission in patch 1
#
# Use a log-spaced range because beta is a rate parameter and the
# mechanistic question concerns relative changes around beta_ref.

r_beta1_grid <- exp(
  seq(
    log(0.1),
    log(10),
    length.out = 51
  )
)


# Patch 2 transmission remains fixed at the reference value.

beta2_beta1_fixed <- beta_ref


# Waterborne pathogen movement remains at the GLOBAL reference values.

MP12_beta1_fixed <- MP12_ref
MP21_beta1_fixed <- MP21_ref


# ------------------------------------------------------------------------------
# 21B. Construct the two-dimensional parameter grid
# ------------------------------------------------------------------------------

MHxBeta1_grid <- expand.grid(
  MH12 = MH12_beta1_grid,
  r_beta1 = r_beta1_grid
)


# Convert the multiplier into the actual beta1 value used by the ODE model.

MHxBeta1_grid$beta1 <-
  beta_ref * MHxBeta1_grid$r_beta1

MHxBeta1_grid$beta2 <-
  beta2_beta1_fixed


# ------------------------------------------------------------------------------
# 21C. Run the model across the full parameter space
# ------------------------------------------------------------------------------

MHxBeta1_results <- do.call(
  rbind,
  lapply(seq_len(nrow(MHxBeta1_grid)), function(k) {
    
    # Current parameter combination
    mh <- MHxBeta1_grid$MH12[k]
    rb <- MHxBeta1_grid$r_beta1[k]
    
    beta1_k <- beta_ref * rb
    beta2_k <- beta2_beta1_fixed
    
    
    # All other biological parameters remain at their reference values.
    #
    # IMPORTANT:
    # beta1 changes.
    # beta2 remains fixed at beta_ref.
    #
    # Waterborne transport uses the GLOBAL reference values,
    # not the exploratory values used in Experiment 20.
    
    pars <- make_pars(
      beta1 = beta1_k,
      beta2 = beta2_k,
      
      MH12 = mh,
      MH21 = MH21_beta1_fixed,
      
      MP12 = MP12_beta1_fixed,
      MP21 = MP21_beta1_fixed
    )
    
    
    # Simulate the full nonlinear two-patch model.
    
    out <- simulate_model(
      pars = pars,
      init = init_ref,
      times_vec = times
    )
    
    
    # Existing summary function:
    # final states, tail-mean prevalence, peak infection,
    # and patch-retention diagnostics.
    
    s <- summarize_run(
      out = out,
      init = init_ref
    )
    
    s <- flag_depletion(s)
    
    
    # ----------------------------------------------------------
    # Long-term state variables
    # ----------------------------------------------------------
    
    S1_tail_mean <- tail_mean(
      out$S1,
      fraction = 0.1
    )
    
    
    P1_tail_mean <- tail_mean(
      out$P1,
      fraction = 0.1
    )
    
    S2_tail_mean <- tail_mean(
      out$S2,
      fraction = 0.1
    )
   
    
    P2_tail_mean <- tail_mean(
      out$P2,
      fraction = 0.1
    )
    
    
    # ----------------------------------------------------------
    # Long-term fish abundance
    # ----------------------------------------------------------
    
    N1_tail_mean <- tail_mean(
      out$N1,
      fraction = 0.1
    )
    
    N2_tail_mean <- tail_mean(
      out$N2,
      fraction = 0.1
    )
    
    NT_tail_mean <- tail_mean(
      out$NT,
      fraction = 0.1
    )
    
    
    # ----------------------------------------------------------
    # Long-term patch-1 incidence
    #
    # incidence1 is already returned by isa_two_patch()
    # through the ODE output.
    # ----------------------------------------------------------
    
    incidence1_tail_mean <- tail_mean(
      out$incidence1,
      fraction = 0.1
    )
    
    
    # ----------------------------------------------------------
    # System-level R0
    # ----------------------------------------------------------
    
    R0_value <- calc_R0(pars)
    
    
    # ----------------------------------------------------------
    # Return all relevant outcomes
    # ----------------------------------------------------------
    
    cbind(
      MH12 = mh,
      MH21 = MH21_beta1_fixed,
      
      r_beta1 = rb,
      
      beta1 = beta1_k,
      beta2 = beta2_k,
      
      MP12 = MP12_beta1_fixed,
      MP21 = MP21_beta1_fixed,
      
      R0 = R0_value,
      
      S1_tail_mean = S1_tail_mean,
      
      P1_tail_mean = P1_tail_mean,
      
      S2_tail_mean = S2_tail_mean,
      
      P2_tail_mean = P2_tail_mean,
      
      N1_tail_mean = N1_tail_mean,
      N2_tail_mean = N2_tail_mean,
      NT_tail_mean = NT_tail_mean,
      
      incidence1_tail_mean = incidence1_tail_mean,
      
      s
    )
  })
)

rownames(MHxBeta1_results) <- NULL


# ------------------------------------------------------------------------------
# 21D. Basic checks
# ------------------------------------------------------------------------------

head(MHxBeta1_results)

summary(MHxBeta1_results)


# Check actual beta1 range

cat(
  "\nPatch-1 beta range:",
  range(MHxBeta1_results$beta1),
  "\n"
)


# Check that beta2 stayed fixed

cat(
  "\nUnique beta2 values:",
  unique(MHxBeta1_results$beta2),
  "\n"
)


# Check that waterborne transport stayed at the global reference values

cat(
  "\nMP12:",
  unique(MHxBeta1_results$MP12),
  "\nMP21:",
  unique(MHxBeta1_results$MP21),
  "\n"
)


# ------------------------------------------------------------------------------
# 21E. Reference condition
# ------------------------------------------------------------------------------

# Reference condition:
#
# MH12   = MH_ref = 0.05
# r_beta1 = 1
# beta1  = beta_ref
# beta2  = beta_ref

reference_MHxBeta1 <- data.frame(
  MH12 = MH_ref,
  r_beta1 = 1
)


reference_MHxBeta1_result <- MHxBeta1_results[
  abs(MHxBeta1_results$MH12 - MH_ref) < 1e-10 &
    abs(MHxBeta1_results$r_beta1 - 1) < 1e-10,
]

print(reference_MHxBeta1_result)


# ==============================================================================
# HEATMAP 1
# System-level R0
# ==============================================================================

p_MHxBeta1_R0 <- ggplot(
  MHxBeta1_results,
  aes(
    x = MH12,
    y = r_beta1,
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
  
  # R0 = 1 invasion threshold
  geom_contour(
    aes(z = R0),
    breaks = 1,
    color = "black",
    linewidth = 0.7
  ) +
  
  # Reference parameter combination
  geom_point(
    data = reference_MHxBeta1,
    aes(
      x = MH12,
      y = r_beta1
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 3.5,
    stroke = 1
  ) +
  
  # Log scale is appropriate because r_beta1 spans 0.1--10
  scale_y_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10)
  ) +
  
  labs(
    x = expression(
      "Host movement " *
        M[12]^H *
        " (month"^{-1} * ")"
    ),
    
    y = expression(
      "Patch 1 transmission multiplier " *
        r[beta*1] == beta[1] / beta[ref]
    ),
    
    fill = expression(R[0]),
    
    title = expression(
      "Combined effects of host movement and patch-1 transmission on " *
        R[0]
    )
  ) +
  
  theme_classic(base_size = 13)

print(p_MHxBeta1_R0)


# ==============================================================================
# HEATMAP 2
# Long-term total disease prevalence
# ==============================================================================

p_MHxBeta1_prev <- ggplot(
  MHxBeta1_results,
  aes(
    x = MH12,
    y = r_beta1,
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
    data = reference_MHxBeta1,
    aes(
      x = MH12,
      y = r_beta1
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 3.5,
    stroke = 1
  ) +
  
  scale_y_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10)
  ) +
  
  labs(
    x = expression(M[12]^H),
    
    y = expression(
      r[beta*1] == beta[1] / beta[ref]
    ),
    
    fill = "Long-term\nprevalence",
    
    title = "Long-term total disease prevalence"
  ) +
  
  theme_classic(base_size = 13)

print(p_MHxBeta1_prev)


# ==============================================================================
# HEATMAP 3
# Long-term total fish abundance
# ==============================================================================

p_MHxBeta1_NT <- ggplot(
  MHxBeta1_results,
  aes(
    x = MH12,
    y = r_beta1,
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
    data = reference_MHxBeta1,
    aes(
      x = MH12,
      y = r_beta1
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 3.5,
    stroke = 1
  ) +
  
  scale_y_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10)
  ) +
  
  labs(
    x = expression(M[12]^H),
    
    y = expression(
      r[beta*1] == beta[1] / beta[ref]
    ),
    
    fill = "Long-term\nfish abundance",
    
    title = "Long-term total fish abundance"
  ) +
  
  theme_classic(base_size = 13)

print(p_MHxBeta1_NT)


# ==============================================================================
# HEATMAP 4
# Patch-1 long-term prevalence
# ==============================================================================

p_MHxBeta1_prev1 <- ggplot(
  MHxBeta1_results,
  aes(
    x = MH12,
    y = r_beta1,
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
    data = reference_MHxBeta1,
    aes(
      x = MH12,
      y = r_beta1
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 3.5,
    stroke = 1
  ) +
  
  scale_y_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10)
  ) +
  
  labs(
    x = expression(M[12]^H),
    
    y = expression(
      r[beta*1] == beta[1] / beta[ref]
    ),
    
    fill = "Patch 1\nprevalence",
    
    title = "Long-term prevalence in patch 1"
  ) +
  
  theme_classic(base_size = 13)

print(p_MHxBeta1_prev1)

###################save and read the values of heatmap

# ==============================================================================
# SAVE AND RELOAD EXPERIMENT 21 RESULTS
# Host movement x Patch-1 transmission intensity
#
# FIRST RUN:
#   Run the SAVE section once after the full simulation has finished.
#
# FUTURE RUNS:
#   Skip the SAVE section.
#   Start directly from the READ section below.
# ==============================================================================


# ==============================================================================
# A. SAVE RESULTS
# RUN THIS SECTION ONLY AFTER A NEW FULL SIMULATION
# ==============================================================================

output_dir <- "ISA_two_patch_outputs"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save the complete R object for future use
saveRDS(
  MHxBeta1_results,
  file = file.path(
    output_dir,
    "MHxBeta1_results.rds"
  )
)

# Also save a CSV copy for inspection / external use
write.csv(
  MHxBeta1_results,
  file = file.path(
    output_dir,
    "MHxBeta1_results.csv"
  ),
  row.names = FALSE
)

cat(
  "\nExperiment 21 results saved successfully.\n",
  "RDS file: ",
  file.path(output_dir, "MHxBeta1_results.rds"),
  "\nCSV file: ",
  file.path(output_dir, "MHxBeta1_results.csv"),
  "\n"
)


# ==============================================================================
# B. READ SAVED RESULTS
#
# FOR ALL FUTURE ANALYSES, START HERE.
# You do NOT need to rerun the simulation or SAVE section above.
# ==============================================================================

output_dir <- "ISA_two_patch_outputs"

MHxBeta1_results <- readRDS(
  file.path(
    output_dir,
    "MHxBeta1_results.rds"
  )
)

cat(
  "\nSaved Experiment 21 results loaded successfully.\n",
  "Number of parameter combinations:",
  nrow(MHxBeta1_results),
  "\n"
)

# Check that the expected columns are present
print(names(MHxBeta1_results))

# Check that there are no duplicated column names
cat(
  "Duplicated column names:",
  anyDuplicated(names(MHxBeta1_results)),
  "\n"
)

########followup test when beta1 is relatively smaller at 0.2 and 0.5
# ==============================================================================
# PATCH-1 PREVALENCE AT LOW beta1
#
# Panel A: r_beta1 approximately 0.2
# Panel B: r_beta1 approximately 0.5
#
# x-axis: MH12
# y-axis: Patch-1 long-term prevalence
#
# Uses saved MHxBeta1_results only.
# NO simulation is rerun.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Find the simulated r_beta1 values closest to 0.2 and 0.5
# ------------------------------------------------------------------------------

available_r_beta1 <- sort(
  unique(MHxBeta1_results$r_beta1)
)

r_beta1_A <- available_r_beta1[
  which.min(abs(available_r_beta1 - 0.2))
]

r_beta1_B <- available_r_beta1[
  which.min(abs(available_r_beta1 - 0.5))
]

cat(
  "\nPanel A: target r_beta1 = 0.2",
  "\nActual simulated value =", r_beta1_A,
  "\n\nPanel B: target r_beta1 = 0.5",
  "\nActual simulated value =", r_beta1_B,
  "\n"
)


# ------------------------------------------------------------------------------
# 2. Extract the two slices
# ------------------------------------------------------------------------------

prev1_beta_A <- MHxBeta1_results[
  abs(MHxBeta1_results$r_beta1 - r_beta1_A) < 1e-12,
]

prev1_beta_B <- MHxBeta1_results[
  abs(MHxBeta1_results$r_beta1 - r_beta1_B) < 1e-12,
]

# Sort by host movement
prev1_beta_A <- prev1_beta_A[
  order(prev1_beta_A$MH12),
]

prev1_beta_B <- prev1_beta_B[
  order(prev1_beta_B$MH12),
]


# ------------------------------------------------------------------------------
# 3. Check the prevalence ranges
# ------------------------------------------------------------------------------

cat(
  "\nPatch-1 prevalence range -- Panel A:\n"
)

print(
  range(
    prev1_beta_A$prev1_tail_mean,
    na.rm = TRUE
  )
)

cat(
  "\nPatch-1 prevalence range -- Panel B:\n"
)

print(
  range(
    prev1_beta_B$prev1_tail_mean,
    na.rm = TRUE
  )
)


# ------------------------------------------------------------------------------
# 4. Panel A: r_beta1 approximately 0.2
# ------------------------------------------------------------------------------

p_prev1_beta_A <- ggplot(
  prev1_beta_A,
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
    y = "Patch-1 long-term prevalence",
    
    title = expression(
      bold("(A)") ~ r[beta*1] == 0.2
    )
  ) +
  
  theme_classic(base_size = 13)


# ------------------------------------------------------------------------------
# 5. Panel B: r_beta1 approximately 0.5
# ------------------------------------------------------------------------------

p_prev1_beta_B <- ggplot(
  prev1_beta_B,
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
    y = "Patch-1 long-term prevalence",
    
    title = expression(
      bold("(B)") ~ r[beta*1] == 0.5
    )
  ) +
  
  theme_classic(base_size = 13)


# ------------------------------------------------------------------------------
# 6. Display A and B together
# ------------------------------------------------------------------------------

library(patchwork)

p_prev1_low_beta_AB <-
  p_prev1_beta_A +
  p_prev1_beta_B +
  plot_layout(ncol = 2)

print(p_prev1_low_beta_AB)
# ==============================================================================
# HEATMAP 5
# Patch-2 long-term prevalence
# ==============================================================================

p_MHxBeta1_prev2 <- ggplot(
  MHxBeta1_results,
  aes(
    x = MH12,
    y = r_beta1,
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
    data = reference_MHxBeta1,
    aes(
      x = MH12,
      y = r_beta1
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 3.5,
    stroke = 1
  ) +
  
  scale_y_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10)
  ) +
  
  labs(
    x = expression(M[12]^H),
    
    y = expression(
      r[beta*1] == beta[1] / beta[ref]
    ),
    
    fill = "Patch 2\nprevalence",
    
    title = "Long-term prevalence in patch 2"
  ) +
  
  theme_classic(base_size = 13)

print(p_MHxBeta1_prev2)


# ==============================================================================
# OPTIONAL MECHANISTIC HEATMAP
# Patch-1 long-term incidence
# ==============================================================================

p_MHxBeta1_inc1 <- ggplot(
  MHxBeta1_results,
  aes(
    x = MH12,
    y = r_beta1,
    fill = incidence1_tail_mean
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
    data = reference_MHxBeta1,
    aes(
      x = MH12,
      y = r_beta1
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 3.5,
    stroke = 1
  ) +
  
  scale_y_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10)
  ) +
  
  labs(
    x = expression(M[12]^H),
    
    y = expression(
      r[beta*1] == beta[1] / beta[ref]
    ),
    
    fill = "Patch 1\nincidence",
    
    title = expression(
      "Long-term patch-1 incidence: " *
        beta[1] * S[1] * P[1]
    )
  ) +
  
  theme_classic(base_size = 13)

print(p_MHxBeta1_inc1)


# ==============================================================================
# 21F. Depletion diagnostic
# ==============================================================================

cat(
  "\nNumber of parameter combinations with patch depletion:",
  sum(MHxBeta1_results$any_depletion, na.rm = TRUE),
  "out of",
  nrow(MHxBeta1_results),
  "\n"
)


# ==============================================================================
# 21G. Export results and figures
# ==============================================================================

# Reuse the output directory already defined earlier in the script.

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


write.csv(
  MHxBeta1_results,
  file.path(
    output_dir,
    "MHx_patch1_transmission_results.csv"
  ),
  row.names = FALSE
)


ggsave(
  file.path(
    output_dir,
    "MHx_patch1_transmission_R0.png"
  ),
  p_MHxBeta1_R0,
  width = 7,
  height = 5.5,
  dpi = 300
)


ggsave(
  file.path(
    output_dir,
    "MHx_patch1_transmission_total_prevalence.png"
  ),
  p_MHxBeta1_prev,
  width = 7,
  height = 5.5,
  dpi = 300
)


ggsave(
  file.path(
    output_dir,
    "MHx_patch1_transmission_total_fish_abundance.png"
  ),
  p_MHxBeta1_NT,
  width = 7,
  height = 5.5,
  dpi = 300
)


ggsave(
  file.path(
    output_dir,
    "MHx_patch1_transmission_patch1_prevalence.png"
  ),
  p_MHxBeta1_prev1,
  width = 7,
  height = 5.5,
  dpi = 300
)


ggsave(
  file.path(
    output_dir,
    "MHx_patch1_transmission_patch2_prevalence.png"
  ),
  p_MHxBeta1_prev2,
  width = 7,
  height = 5.5,
  dpi = 300
)


ggsave(
  file.path(
    output_dir,
    "MHx_patch1_transmission_patch1_incidence.png"
  ),
  p_MHxBeta1_inc1,
  width = 7,
  height = 5.5,
  dpi = 300
)


# ==============================================================================
# End of Experiment 21
# ==============================================================================

