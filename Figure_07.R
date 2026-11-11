# APPENDIX FIGURE 7 — HOST MOVEMENT AND WATERBORNE TRANSPORT
# Extracted from Experiment 20 of the historical 2026-10-04 analysis script.
# Related main manuscript: Figure 2.
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
# EXPERIMENT 20
# HOST MOVEMENT x WATERBORNE PATHOGEN MOVEMENT
# ==============================================================================
#
# Scientific question:
#
# How does the epidemiological effect of directional human-mediated
# host movement interact with asymmetric waterborne pathogen transport?
#
#
# Experimental dimensions:
#
#   x-axis:
#
#       M12^H
#
#       Human-mediated host movement from the host-exporting patch
#       (Patch 1) to the host-receiving patch (Patch 2).
#
#
#   y-axis:
#
#                    M12^P
#       r_P = ----------------
#                    M21^P
#
#
# Interpretation:
#
#   r_P < 1 :
#       stronger waterborne pathogen transport from Patch 2 -> Patch 1
#
#   r_P = 1 :
#       symmetric waterborne pathogen transport
#
#   r_P > 1 :
#       stronger waterborne pathogen transport from Patch 1 -> Patch 2
#
#
# Human-mediated host movement:
#
#   M12^H = 0 to 0.10 month^-1
#   M21^H = 0
#
#
# Waterborne pathogen movement:
#
#   M21^P is fixed locally at 0.05 month^-1
#
#   M12^P = r_P * M21^P
#
#   r_P = 0 to 100
#
# Therefore:
#
#   M12^P = 0 to 5 month^-1
#
#
# IMPORTANT:
#
# This experiment is a deliberately broad stress-test of waterborne
# pathogen redistribution.
#
# The fixed M21^P = 0.05 used here is NOT the global reference value
# MP21_ref = 0.01.
#
# Consequently, the global waterborne reference combination
#
#       MP12_ref = 0.05
#       MP21_ref = 0.01
#
# does NOT lie on this two-dimensional experimental slice.
#
#
# Reference marker shown on the heatmaps:
#
#       M12^H = MH_ref = 0.05 month^-1
#       r_P   = 1
#
# Under the Experiment-20 fixed value M21^P = 0.05,
# this corresponds to:
#
#       M12^P = 0.05
#       M21^P = 0.05
#
# Thus, the "x" represents:
#
#       reference host movement
#       +
#       symmetric waterborne pathogen transport
#
# It should NOT be interpreted as the full global model baseline.
#
#
# Main outcomes:
#
#   1. System-level R0
#   2. Long-term total disease prevalence
#   3. Long-term total fish abundance
#   4. Host-exporting-patch long-term prevalence
#   5. Host-receiving-patch long-term prevalence
#
# ==============================================================================



# ==============================================================================
# 20A. PARAMETER SPACE
# ==============================================================================


# ------------------------------------------------------------------------------
# Human-mediated host movement
#
# Patch 1 -> Patch 2
# ------------------------------------------------------------------------------

MH12_HP_grid <- seq(
  from = 0,
  to = 0.10,
  length.out = 51
)


# Reverse human-mediated host movement remains zero.

MH21_HP_fixed <- 0



# ------------------------------------------------------------------------------
# Waterborne pathogen movement
# ------------------------------------------------------------------------------

# IMPORTANT:
#
# This is the LOCAL Experiment-20 value.
#
# It deliberately differs from:
#
# MP21_ref = 0.01
#
# used in the global reference parameterization.

MP21_HP_fixed <- 0.05



# ------------------------------------------------------------------------------
# Waterborne pathogen transport ratio
#
#                 M12^P
#       r_P = -----------
#                 M21^P
#
# ------------------------------------------------------------------------------

rP_grid <- seq(
  from = 0,
  to = 100,
  length.out = 51
)



# ==============================================================================
# 20B. CONSTRUCT THE TWO-DIMENSIONAL PARAMETER GRID
# ==============================================================================

MHxRP_grid <- expand.grid(
  
  MH12 = MH12_HP_grid,
  
  rP = rP_grid,
  
  KEEP.OUT.ATTRS = FALSE,
  
  stringsAsFactors = FALSE
)



# ------------------------------------------------------------------------------
# Convert r_P into actual waterborne pathogen movement rates
# ------------------------------------------------------------------------------

MHxRP_grid$MP21 <- MP21_HP_fixed


MHxRP_grid$MP12 <-
  MHxRP_grid$rP *
  MHxRP_grid$MP21



# ==============================================================================
# 20C. CHECK THE PARAMETER SPACE BEFORE RUNNING THE MODEL
# ==============================================================================

cat(
  "\n============================================================\n",
  "EXPERIMENT 20 PARAMETER SPACE\n",
  "============================================================\n"
)


cat(
  "\nHuman-mediated host movement M12^H:\n"
)

print(
  range(
    MHxRP_grid$MH12
  )
)


cat(
  "\nWaterborne movement ratio r_P:\n"
)

print(
  range(
    MHxRP_grid$rP
  )
)


cat(
  "\nFixed M21^P:\n"
)

print(
  unique(
    MHxRP_grid$MP21
  )
)


cat(
  "\nResulting M12^P range:\n"
)

print(
  range(
    MHxRP_grid$MP12
  )
)


cat(
  "\nTotal parameter combinations:\n"
)

print(
  nrow(
    MHxRP_grid
  )
)



# Expected:
#
# 51 x 51 = 2601

stopifnot(
  nrow(MHxRP_grid) == 2601
)



# ==============================================================================
# 20D. RUN THE FULL MODEL
# ==============================================================================

MHxRP_results <- do.call(
  
  rbind,
  
  lapply(
    
    seq_len(
      nrow(MHxRP_grid)
    ),
    
    function(k) {
      
      
      # ------------------------------------------------------------------------
      # Current parameter combination
      # ------------------------------------------------------------------------
      
      mh <- MHxRP_grid$MH12[k]
      
      rp <- MHxRP_grid$rP[k]
      
      
      mp21 <- MP21_HP_fixed
      
      mp12 <- rp * mp21
      
      
      
      # ------------------------------------------------------------------------
      # Construct model parameters
      #
      # All biological parameters not explicitly changed here remain
      # at their global reference values through make_pars().
      # ------------------------------------------------------------------------
      
      pars_k <- make_pars(
        
        MH12 = mh,
        
        MH21 = MH21_HP_fixed,
        
        
        MP12 = mp12,
        
        MP21 = mp21
      )
      
      
      
      # ------------------------------------------------------------------------
      # Run nonlinear ODE model
      # ------------------------------------------------------------------------
      
      out_k <- simulate_model(
        
        pars = pars_k,
        
        init = init_ref,
        
        times_vec = times
      )
      
      
      
      # ------------------------------------------------------------------------
      # Standard model summary
      #
      # Includes:
      #
      #   final states
      #   long-term prevalence
      #   peak infection
      #   retention/depletion diagnostics
      #
      # ------------------------------------------------------------------------
      
      s_k <- summarize_run(
        
        out = out_k,
        
        init = init_ref
      )
      
      
      s_k <- flag_depletion(
        s_k
      )
      
      
      
      # ------------------------------------------------------------------------
      # Long-term fish abundance
      #
      # Mean over final 10% of simulation.
      # ------------------------------------------------------------------------
      
      N1_tail_mean <- tail_mean(
        
        out_k$N1,
        
        fraction = 0.1
      )
      
      
      N2_tail_mean <- tail_mean(
        
        out_k$N2,
        
        fraction = 0.1
      )
      
      
      NT_tail_mean <- tail_mean(
        
        out_k$NT,
        
        fraction = 0.1
      )
      
      
      
      # ------------------------------------------------------------------------
      # Long-term environmental pathogen concentration
      #
      # Retained as a mechanistic diagnostic.
      # ------------------------------------------------------------------------
      
      P1_tail_mean <- tail_mean(
        
        out_k$P1,
        
        fraction = 0.1
      )
      
      
      P2_tail_mean <- tail_mean(
        
        out_k$P2,
        
        fraction = 0.1
      )
      
      
      
      # ------------------------------------------------------------------------
      # Long-term susceptible abundance
      #
      # Retained as a mechanistic diagnostic.
      # ------------------------------------------------------------------------
      
      S1_tail_mean <- tail_mean(
        
        out_k$S1,
        
        fraction = 0.1
      )
      
      
      S2_tail_mean <- tail_mean(
        
        out_k$S2,
        
        fraction = 0.1
      )
      
      
      
      # ------------------------------------------------------------------------
      # System-level R0
      #
      # IMPORTANT:
      #
      # R0 is calculated analytically at the disease-free equilibrium
      # using the next-generation matrix.
      #
      # It is NOT calculated from the 120-month ODE trajectory and does
      # NOT depend on the initial I or P values.
      # ------------------------------------------------------------------------
      
      R0_value <- calc_R0(
        pars_k
      )
      
      
      
      # ------------------------------------------------------------------------
      # Return results
      # ------------------------------------------------------------------------
      
      data.frame(
        
        MH12 = mh,
        
        MH21 = MH21_HP_fixed,
        
        
        rP = rp,
        
        
        MP12 = mp12,
        
        MP21 = mp21,
        
        
        R0 = R0_value,
        
        
        S1_tail_mean = S1_tail_mean,
        
        S2_tail_mean = S2_tail_mean,
        
        
        P1_tail_mean = P1_tail_mean,
        
        P2_tail_mean = P2_tail_mean,
        
        
        N1_tail_mean = N1_tail_mean,
        
        N2_tail_mean = N2_tail_mean,
        
        NT_tail_mean = NT_tail_mean,
        
        
        s_k,
        
        check.names = FALSE
      )
    }
  )
)


rownames(
  MHxRP_results
) <- NULL



# ==============================================================================
# 20E. POST-SIMULATION SANITY CHECKS
# ==============================================================================

stopifnot(
  nrow(MHxRP_results) == 2601
)



cat(
  "\n============================================================\n",
  "EXPERIMENT 20 MODEL CHECK\n",
  "============================================================\n"
)



cat(
  "\nMH12 range:\n"
)

print(
  range(
    MHxRP_results$MH12
  )
)



cat(
  "\nr_P range:\n"
)

print(
  range(
    MHxRP_results$rP
  )
)



cat(
  "\nM12^P range:\n"
)

print(
  range(
    MHxRP_results$MP12
  )
)



cat(
  "\nUnique M21^P:\n"
)

print(
  unique(
    MHxRP_results$MP21
  )
)



cat(
  "\nUnique M21^H:\n"
)

print(
  unique(
    MHxRP_results$MH21
  )
)



# ==============================================================================
# 20F. DEFINE THE EXPERIMENT-20 REFERENCE MARKER
# ==============================================================================

# ------------------------------------------------------------------------------
# IMPORTANT:
#
# The marker used in Experiment 20 is NOT the full global baseline.
#
# It represents:
#
#       reference human-mediated host movement
#
#                   MH12 = 0.05
#
# plus:
#
#       symmetric waterborne pathogen transport
#
#                   r_P = 1
#
# Because M21^P is fixed at 0.05 in Experiment 20:
#
#                   M12^P = 0.05
#                   M21^P = 0.05
#
# ------------------------------------------------------------------------------


reference_HP_target <- data.frame(
  
  MH12 = MH_ref,
  
  rP = 1
)



# ------------------------------------------------------------------------------
# Find the actual simulated grid point nearest the intended marker.
#
# This avoids relying on exact floating-point equality.
# ------------------------------------------------------------------------------

reference_HP_index <- which.min(
  
  (
    MHxRP_results$MH12 -
      MH_ref
  )^2 +
    
    (
      MHxRP_results$rP -
        1
    )^2
)



reference_HP_result <-
  
  MHxRP_results[
    
    reference_HP_index,
    
    ,
    
    drop = FALSE
  ]



# ------------------------------------------------------------------------------
# Use the ACTUAL simulated coordinates for plotting the cross.
# ------------------------------------------------------------------------------

reference_HP <- data.frame(
  
  MH12 =
    reference_HP_result$MH12,
  
  rP =
    reference_HP_result$rP
)



# ------------------------------------------------------------------------------
# Print exactly what the cross represents.
# ------------------------------------------------------------------------------

cat(
  "\n============================================================\n",
  "EXPERIMENT 20 REFERENCE MARKER\n",
  "============================================================\n"
)


cat(
  "\nThe heatmap cross represents:\n"
)


print(
  
  reference_HP_result[
    ,
    c(
      "MH12",
      "rP",
      "MP12",
      "MP21",
      "R0"
    )
  ]
)



cat(
  "\nIMPORTANT:\n",
  "This marker is the reference host-movement setting combined\n",
  "with symmetric waterborne pathogen transport within the\n",
  "Experiment-20 parameter slice.\n",
  "It is NOT the full global baseline parameter combination.\n"
)



# ==============================================================================
# 20G. OUTPUT DIRECTORY
# ==============================================================================

output_dir <- "ISA_two_patch_outputs"


if (
  !dir.exists(
    output_dir
  )
) {
  
  dir.create(
    
    output_dir,
    
    recursive = TRUE
  )
}



# ==============================================================================
# 20H. SAVE COMPLETE NUMERICAL RESULTS
# ==============================================================================

saveRDS(
  
  MHxRP_results,
  
  file = file.path(
    
    output_dir,
    
    "MHx_waterborne_ratio_results.rds"
  )
)



write.csv(
  
  MHxRP_results,
  
  file = file.path(
    
    output_dir,
    
    "MHx_waterborne_ratio_results.csv"
  ),
  
  row.names = FALSE
)



# ==============================================================================
# 20I. COMMON HEATMAP SETTINGS
# ==============================================================================


# ------------------------------------------------------------------------------
# Common colour palette
# ------------------------------------------------------------------------------

heat_colours <- c(
  
  "navy",
  
  "dodgerblue",
  
  "yellow",
  
  "orange",
  
  "red"
)



# ------------------------------------------------------------------------------
# Common x-axis
# ------------------------------------------------------------------------------

MH_HP_x_scale <- scale_x_continuous(
  
  breaks = seq(
    0,
    0.10,
    by = 0.02
  ),
  
  limits = c(
    0,
    0.10
  ),
  
  expand = c(
    0,
    0
  )
)



# ------------------------------------------------------------------------------
# Common y-axis
#
# r_P is deliberately shown on a linear scale because the grid
# includes r_P = 0.
# ------------------------------------------------------------------------------

rP_y_scale <- scale_y_continuous(
  
  breaks = c(
    0,
    20,
    40,
    60,
    80,
    100
  ),
  
  limits = c(
    0,
    100
  ),
  
  expand = c(
    0,
    0
  )
)



# ------------------------------------------------------------------------------
# Axis labels
# ------------------------------------------------------------------------------

x_label_HP <- expression(
  
  "Human-mediated host movement, " ~
    
    M[12]^H ~
    
    " (month"^{-1}*")"
)



y_label_rP <- expression(
  
  "Waterborne pathogen movement ratio, " ~
    
    r[P] == M[12]^P / M[21]^P
)



# ------------------------------------------------------------------------------
# Common theme
# ------------------------------------------------------------------------------

heat_theme_HP <- theme_classic(
  base_size = 13
) +
  
  theme(
    
    panel.grid.major = element_blank(),
    
    panel.grid.minor = element_blank(),
    
    
    axis.title = element_text(
      size = 12
    ),
    
    
    axis.text = element_text(
      size = 10
    ),
    
    
    axis.line = element_line(
      
      linewidth = 0.65,
      
      colour = "black"
    ),
    
    
    axis.ticks = element_line(
      
      linewidth = 0.55,
      
      colour = "black"
    ),
    
    
    axis.ticks.length = unit(
      0.10,
      "cm"
    ),
    
    
    legend.title = element_text(
      size = 11
    ),
    
    
    legend.text = element_text(
      size = 10
    ),
    
    
    plot.title = element_blank()
  )



# ==============================================================================
# 20J. HEATMAP 1
# SYSTEM-LEVEL R0
# ==============================================================================

p_MHxRP_R0 <- ggplot(
  
  MHxRP_results,
  
  aes(
    
    x = MH12,
    
    y = rP,
    
    fill = R0
  )
) +
  
  
  geom_tile() +
  
  
  scale_fill_gradientn(
    
    colours = heat_colours,
    
    name = expression(
      R[0]
    )
  ) +
  
  
  # ---------------------------------------------------------------------------
# R0 = 1 invasion threshold
# ---------------------------------------------------------------------------

geom_contour(
  
  aes(
    z = R0
  ),
  
  breaks = 1,
  
  colour = "black",
  
  linewidth = 0.8
) +
  
  
  # ---------------------------------------------------------------------------
# Experiment-20 reference marker
# ---------------------------------------------------------------------------

geom_point(
  
  data = reference_HP,
  
  aes(
    
    x = MH12,
    
    y = rP
  ),
  
  inherit.aes = FALSE,
  
  shape = 4,
  
  size = 4,
  
  stroke = 1.2,
  
  colour = "black"
) +
  
  
  MH_HP_x_scale +
  
  rP_y_scale +
  
  
  labs(
    
    x = x_label_HP,
    
    y = y_label_rP
  ) +
  
  
  heat_theme_HP



print(
  p_MHxRP_R0
)



# ==============================================================================
# 20K. HEATMAP 2
# LONG-TERM TOTAL DISEASE PREVALENCE
# ==============================================================================

p_MHxRP_prev <- ggplot(
  
  MHxRP_results,
  
  aes(
    
    x = MH12,
    
    y = rP,
    
    fill = prev_total_tail_mean * 100
  )
) +
  
  
  geom_tile() +
  
  
  scale_fill_gradientn(
    
    colours = heat_colours,
    
    name =
      "Long-term total\nprevalence (%)"
  ) +
  
  
  # Experiment-20 reference marker
  
  geom_point(
    
    data = reference_HP,
    
    aes(
      
      x = MH12,
      
      y = rP
    ),
    
    inherit.aes = FALSE,
    
    shape = 4,
    
    size = 4,
    
    stroke = 1.2,
    
    colour = "black"
  ) +
  
  
  MH_HP_x_scale +
  
  rP_y_scale +
  
  
  labs(
    
    x = x_label_HP,
    
    y = y_label_rP
  ) +
  
  
  heat_theme_HP



print(
  p_MHxRP_prev
)



# ==============================================================================
# 20L. HEATMAP 3
# LONG-TERM TOTAL FISH ABUNDANCE
# ==============================================================================

p_MHxRP_NT <- ggplot(
  
  MHxRP_results,
  
  aes(
    
    x = MH12,
    
    y = rP,
    
    fill = NT_tail_mean / 1e6
  )
) +
  
  
  geom_tile() +
  
  
  scale_fill_gradientn(
    
    colours = heat_colours,
    
    name = expression(
      "Long-term total fish abundance (" %*% 10^6 * ")"
    )
  ) +
  
  
  # Experiment-20 reference marker
  
  geom_point(
    
    data = reference_HP,
    
    aes(
      
      x = MH12,
      
      y = rP
    ),
    
    inherit.aes = FALSE,
    
    shape = 4,
    
    size = 4,
    
    stroke = 1.2,
    
    colour = "black"
  ) +
  
  
  MH_HP_x_scale +
  
  rP_y_scale +
  
  
  labs(
    
    x = x_label_HP,
    
    y = y_label_rP
  ) +
  
  
  heat_theme_HP



print(
  p_MHxRP_NT
)



# ==============================================================================
# 20M. HEATMAP 4
# HOST-EXPORTING-PATCH LONG-TERM PREVALENCE
# ==============================================================================

p_MHxRP_prev1 <- ggplot(
  
  MHxRP_results,
  
  aes(
    
    x = MH12,
    
    y = rP,
    
    fill = prev1_tail_mean * 100
  )
) +
  
  
  geom_tile() +
  
  
  scale_fill_gradientn(
    
    colours = heat_colours,
    
    name =
      "Host-exporting-patch\nprevalence (%)"
  ) +
  
  
  # Experiment-20 reference marker
  
  geom_point(
    
    data = reference_HP,
    
    aes(
      
      x = MH12,
      
      y = rP
    ),
    
    inherit.aes = FALSE,
    
    shape = 4,
    
    size = 4,
    
    stroke = 1.2,
    
    colour = "black"
  ) +
  
  
  MH_HP_x_scale +
  
  rP_y_scale +
  
  
  labs(
    
    x = x_label_HP,
    
    y = y_label_rP
  ) +
  
  
  heat_theme_HP



print(
  p_MHxRP_prev1
)



# ==============================================================================
# 20N. HEATMAP 5
# HOST-RECEIVING-PATCH LONG-TERM PREVALENCE
# ==============================================================================

p_MHxRP_prev2 <- ggplot(
  
  MHxRP_results,
  
  aes(
    
    x = MH12,
    
    y = rP,
    
    fill = prev2_tail_mean * 100
  )
) +
  
  
  geom_tile() +
  
  
  scale_fill_gradientn(
    
    colours = heat_colours,
    
    name =
      "Host-receiving-patch\nprevalence (%)"
  ) +
  
  
  # Experiment-20 reference marker
  
  geom_point(
    
    data = reference_HP,
    
    aes(
      
      x = MH12,
      
      y = rP
    ),
    
    inherit.aes = FALSE,
    
    shape = 4,
    
    size = 4,
    
    stroke = 1.2,
    
    colour = "black"
  ) +
  
  
  MH_HP_x_scale +
  
  rP_y_scale +
  
  
  labs(
    
    x = x_label_HP,
    
    y = y_label_rP
  ) +
  
  
  heat_theme_HP



print(
  p_MHxRP_prev2
)



# ==============================================================================
# 20O. OPTIONAL MECHANISTIC DIAGNOSTICS
# ==============================================================================


# ------------------------------------------------------------------------------
# Depletion diagnostic
# ------------------------------------------------------------------------------

cat(
  "\n============================================================\n",
  "DEPLETION DIAGNOSTIC\n",
  "============================================================\n"
)


cat(
  
  "\nNumber of parameter combinations with patch depletion:",
  
  sum(
    MHxRP_results$any_depletion,
    na.rm = TRUE
  ),
  
  "out of",
  
  nrow(
    MHxRP_results
  ),
  
  "\n"
)



# ------------------------------------------------------------------------------
# R0 range
# ------------------------------------------------------------------------------

cat(
  "\nR0 range:\n"
)

print(
  
  range(
    
    MHxRP_results$R0,
    
    na.rm = TRUE
  )
)



# ------------------------------------------------------------------------------
# Total prevalence range
# ------------------------------------------------------------------------------

cat(
  "\nLong-term total prevalence range:\n"
)

print(
  
  range(
    
    MHxRP_results$prev_total_tail_mean,
    
    na.rm = TRUE
  )
)



# ------------------------------------------------------------------------------
# Total fish abundance range
# ------------------------------------------------------------------------------

cat(
  "\nLong-term total fish abundance range:\n"
)

print(
  
  range(
    
    MHxRP_results$NT_tail_mean,
    
    na.rm = TRUE
  )
)



# ------------------------------------------------------------------------------
# Host-exporting-patch prevalence range
# ------------------------------------------------------------------------------

cat(
  "\nHost-exporting-patch long-term prevalence range:\n"
)

print(
  
  range(
    
    MHxRP_results$prev1_tail_mean,
    
    na.rm = TRUE
  )
)



# ------------------------------------------------------------------------------
# Host-receiving-patch prevalence range
# ------------------------------------------------------------------------------

cat(
  "\nHost-receiving-patch long-term prevalence range:\n"
)

print(
  
  range(
    
    MHxRP_results$prev2_tail_mean,
    
    na.rm = TRUE
  )
)



# ==============================================================================
# 20P. SAVE ALL FIGURES
# ==============================================================================

ggsave(
  
  filename = file.path(
    
    output_dir,
    
    "MHx_waterborne_ratio_R0.png"
  ),
  
  plot = p_MHxRP_R0,
  
  width = 7.5,
  
  height = 5.8,
  
  dpi = 400
)



ggsave(
  
  filename = file.path(
    
    output_dir,
    
    "MHx_waterborne_ratio_total_prevalence.png"
  ),
  
  plot = p_MHxRP_prev,
  
  width = 7.5,
  
  height = 5.8,
  
  dpi = 400
)



ggsave(
  
  filename = file.path(
    
    output_dir,
    
    "MHx_waterborne_ratio_total_fish_abundance.png"
  ),
  
  plot = p_MHxRP_NT,
  
  width = 7.5,
  
  height = 5.8,
  
  dpi = 400
)



ggsave(
  
  filename = file.path(
    
    output_dir,
    
    "MHx_waterborne_ratio_exporting_patch_prevalence.png"
  ),
  
  plot = p_MHxRP_prev1,
  
  width = 7.5,
  
  height = 5.8,
  
  dpi = 400
)



ggsave(
  
  filename = file.path(
    
    output_dir,
    
    "MHx_waterborne_ratio_receiving_patch_prevalence.png"
  ),
  
  plot = p_MHxRP_prev2,
  
  width = 7.5,
  
  height = 5.8,
  
  dpi = 400
)



# ==============================================================================
# 20Q. FINAL SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n",
  "EXPERIMENT 20 COMPLETE\n",
  "============================================================\n"
)


cat(
  
  "\nParameter combinations:",
  
  nrow(
    MHxRP_results
  ),
  
  "\n"
)



cat(
  "\nExperimental waterborne setting:\n"
)


cat(
  
  "Fixed M21^P =",
  
  MP21_HP_fixed,
  
  "month^-1\n"
)


cat(
  
  "r_P range =",
  
  min(
    rP_grid
  ),
  
  "to",
  
  max(
    rP_grid
  ),
  
  "\n"
)


cat(
  
  "M12^P range =",
  
  min(
    MHxRP_results$MP12
  ),
  
  "to",
  
  max(
    MHxRP_results$MP12
  ),
  
  "month^-1\n"
)



cat(
  "\nHeatmap reference marker:\n"
)


cat(
  
  "M12^H =",
  
  reference_HP$MH12,
  
  "month^-1\n"
)


cat(
  
  "r_P =",
  
  reference_HP$rP,
  
  "\n"
)


cat(
  
  "M12^P =",
  
  reference_HP_result$MP12,
  
  "month^-1\n"
)


cat(
  
  "M21^P =",
  
  reference_HP_result$MP21,
  
  "month^-1\n"
)


cat(
  
  "R0 at marker =",
  
  reference_HP_result$R0,
  
  "\n"
)



cat(
  "\nNOTE:\n",
  "The heatmap marker is NOT the full global model baseline.\n",
  "It represents reference host movement plus symmetric\n",
  "waterborne pathogen transport within Experiment 20.\n"
)



cat(
  
  "\nResults saved to:",
  
  output_dir,
  
  "\n"
)



# ==============================================================================
# END OF EXPERIMENT 20
