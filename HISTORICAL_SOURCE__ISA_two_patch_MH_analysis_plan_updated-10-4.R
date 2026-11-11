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
}


######realistic Re
calc_Re <- function(pars, S1_current, S2_current) {
  
  p <- as.list(pars)
  
  gamma1 <- p$alpha1 + p$mu1 + p$h1
  gamma2 <- p$alpha2 + p$mu2 + p$h2
  
  # F evaluated at the CURRENT susceptible abundances
  F <- matrix(0, nrow = 4, ncol = 4)
  
  # infected-state order: I1, I2, P1, P2
  F[1, 3] <- p$beta1 * S1_current
  F[2, 4] <- p$beta2 * S2_current
  
  # Same transition matrix used in calc_R0()
  # NOTE: chi1 and chi2 already include division by V_effect
  V <- matrix(
    c(
      gamma1 + p$MH12, -p$MH21, 0, 0,
      -p$MH12, gamma2 + p$MH21, 0, 0,
      -p$chi1, 0, p$delta1 + p$MP12, -p$MP21,
      0, -p$chi2, -p$MP12, p$delta2 + p$MP21
    ),
    nrow = 4,
    byrow = TRUE
  )
  
  K <- F %*% solve(V)
  
  eig <- eigen(K, only.values = TRUE)$values
  
  Re <- max(Mod(eig))
  
  as.numeric(Re)
}

# ------------------------------------------------------------------------------
# 7. Reference checks: R0 and time-varying Re(t)
# ------------------------------------------------------------------------------

# Four reference connectivity scenarios
pars_no_spatial <- make_pars(
  MH12 = 0, MH21 = 0,
  MP12 = 0, MP21 = 0
)

pars_water_only <- make_pars(
  MH12 = 0, MH21 = 0,
  MP12 = MP12_ref, MP21 = MP21_ref
)

pars_host_only <- make_pars(
  MH12 = MH_ref, MH21 = 0,
  MP12 = 0, MP21 = 0
)

pars_combined <- make_pars(
  MH12 = MH_ref, MH21 = 0,
  MP12 = MP12_ref, MP21 = MP21_ref
)


# ------------------------------------------------------------------------------
# 7A. Basic reproduction number R0
# ------------------------------------------------------------------------------

reference_checks <- data.frame(
  scenario = c(
    "no_spatial",
    "waterborne_only",
    "host_movement_only",
    "combined"
  ),
  R0 = c(
    calc_R0(pars_no_spatial),
    calc_R0(pars_water_only),
    calc_R0(pars_host_only),
    calc_R0(pars_combined)
  )
)

print(reference_checks)


# ------------------------------------------------------------------------------
# 7B. Simulate each reference scenario and calculate Re(t)
# ------------------------------------------------------------------------------

scenario_pars <- list(
  no_spatial         = pars_no_spatial,
  waterborne_only    = pars_water_only,
  host_movement_only = pars_host_only,
  combined           = pars_combined
)

Re_runs <- do.call(
  rbind,
  lapply(names(scenario_pars), function(sc) {
    
    pars <- scenario_pars[[sc]]
    
    # Run the dynamic model
    out <- simulate_model(
      pars = pars,
      init = init_ref,
      times_vec = times
    )
    
    # Calculate Re(t) using current S1(t) and S2(t)
    out$Re <- mapply(
      function(S1_now, S2_now) {
        calc_Re(
          pars = pars,
          S1_current = S1_now,
          S2_current = S2_now
        )
      },
      out$S1,
      out$S2
    )
    
    # R0 is constant within each scenario
    out$R0 <- calc_R0(pars)
    
    out$scenario <- sc
    
    out
  })
)

rownames(Re_runs) <- NULL


# ------------------------------------------------------------------------------
# 7C. Selected time-point table
# ------------------------------------------------------------------------------

Re_check_times <- c(0, 1, 3, 6, 12, 18, 24, 48, 72, 96, 120)

Re_table <- do.call(
  rbind,
  lapply(names(scenario_pars), function(sc) {
    
    tmp <- Re_runs[Re_runs$scenario == sc, ]
    
    do.call(
      rbind,
      lapply(Re_check_times, function(tt) {
        
        idx <- which.min(abs(tmp$time - tt))
        
        data.frame(
          scenario = sc,
          time = tmp$time[idx],
          
          R0 = tmp$R0[idx],
          Re = tmp$Re[idx],
          
          S1 = tmp$S1[idx],
          S2 = tmp$S2[idx],
          
          I1 = tmp$I1[idx],
          I2 = tmp$I2[idx],
          
          P1 = tmp$P1[idx],
          P2 = tmp$P2[idx],
          
          prev1 = tmp$prev1[idx],
          prev2 = tmp$prev2[idx],
          prev_total = tmp$prev_total[idx]
        )
      })
    )
  })
)

rownames(Re_table) <- NULL

print(Re_table)


# ------------------------------------------------------------------------------
# 7D. Plot Re(t)
# ------------------------------------------------------------------------------

p_Re <- ggplot(
  Re_runs,
  aes(
    x = time,
    y = Re,
    color = scenario,
    linetype = scenario
  )
) +
  geom_line(linewidth = 1.1) +
  
  # Epidemic threshold
  geom_hline(
    yintercept = 1,
    color = "grey40",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  scale_color_manual(
    values = c(
      "no_spatial"         = "#0072B2",
      "waterborne_only"    = "#56B4E9",
      "host_movement_only" = "#D55E00",
      "combined"           = "#CC79A7"
    ),
    labels = c(
      "no_spatial"         = "No spatial movement",
      "waterborne_only"    = "Waterborne only",
      "host_movement_only" = "Host movement only",
      "combined"           = "Combined"
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "no_spatial"         = "solid",
      "waterborne_only"    = "dashed",
      "host_movement_only" = "dotdash",
      "combined"           = "longdash"
    ),
    labels = c(
      "no_spatial"         = "No spatial movement",
      "waterborne_only"    = "Waterborne only",
      "host_movement_only" = "Host movement only",
      "combined"           = "Combined"
    )
  ) +
  
  labs(
    x = "Time (months)",
    y = expression(R[e](t)),
    color = "Scenario",
    linetype = "Scenario",
    title = expression(
      "Time-varying effective reproduction number " * R[e](t)
    )
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

print(p_Re)

# ------------------------------------------------------------------------------
# 8. Main analysis: one-way human-mediated movement sweep, patch 1 -> patch 2
# ------------------------------------------------------------------------------

MH_grid <- seq(0, 0.5, length.out = 51)

main_sweep <- do.call(
  rbind,
  lapply(MH_grid, function(mh) {

    pars <- make_pars(MH12 = mh, MH21 = 0,
                      MP12 = MP12_ref, MP21 = MP21_ref)

    out <- simulate_model(pars, init = init_ref)
    s <- summarize_run(out, init_ref)
    s <- flag_depletion(s)

    cbind(
      MH12 = mh,
      MH21 = 0,
      R0 = calc_R0(pars),
      s
    )
  })
)

# Example plot: equilibrium/tail prevalence vs MH
p_main_prev <- ggplot(main_sweep, aes(x = MH12)) +
  geom_line(aes(y = prev1_tail_mean, linetype = "Patch 1")) +
  geom_line(aes(y = prev2_tail_mean, linetype = "Patch 2")) +
  geom_line(aes(y = prev_total_tail_mean, linetype = "Total")) +
  labs(
    x = expression(M[12]^H),
    y = "Tail-mean prevalence",
    linetype = NULL,
    title = "Human-mediated fish movement and long-term prevalence"
  ) +
  theme_classic()

print(p_main_prev)

# Host abundance diagnostic
p_main_N <- ggplot(main_sweep, aes(x = MH12)) +
  geom_line(aes(y = N1_final, linetype = "Patch 1")) +
  geom_line(aes(y = N2_final, linetype = "Patch 2")) +
  labs(
    x = expression(M[12]^H),
    y = "Final host abundance",
    linetype = NULL,
    title = "Patch abundance diagnostic during movement sweep"
  ) +
  theme_classic()

print(p_main_N)


# ------------------------------------------------------------------------------
# 9. Source-state dependence and transient dilution/amplification
# ------------------------------------------------------------------------------

# Scenario A: source patch 1 is more infected than destination patch 2
init_source_worse <- make_init(
  S1 = 1e6, I1 = 1000, P1 = 1e6,
  S2 = 1e6, I2 = 10,   P2 = 1e4
)

# Scenario B: destination patch 2 is more infected than source patch 1
init_dest_worse <- make_init(
  S1 = 1e6, I1 = 10,   P1 = 1e4,
  S2 = 1e6, I2 = 1000, P2 = 1e6
)

source_state_sweep <- do.call(
  rbind,
  lapply(c("source_worse", "destination_worse"), function(sc) {

    init <- if (sc == "source_worse") init_source_worse else init_dest_worse

    do.call(
      rbind,
      lapply(MH_grid, function(mh) {

        pars <- make_pars(MH12 = mh, MH21 = 0,
                          MP12 = MP12_ref, MP21 = MP21_ref)

        out <- simulate_model(pars, init = init)
        s <- summarize_run(out, init)
        s <- flag_depletion(s)

        cbind(
          scenario = sc,
          MH12 = mh,
          R0 = calc_R0(pars),
          s
        )
      })
    )
  })
)

# Short-term transient example
transient_mh_values <- c(0, 0.05, 0.15, 0.30)

transient_runs <- do.call(
  rbind,
  lapply(transient_mh_values, function(mh) {

    pars <- make_pars(MH12 = mh, MH21 = 0,
                      MP12 = MP12_ref, MP21 = MP21_ref)

    out <- simulate_model(pars, init = init_dest_worse)
    out$MH12 <- mh
    out
  })
)

p_transient <- ggplot(transient_runs,
                      aes(x = time, y = prev2,
                          group = factor(MH12),
                          linetype = factor(MH12))) +
  geom_line() +
  labs(
    x = "Time",
    y = "Patch 2 prevalence",
    linetype = expression(M[12]^H),
    title = "Possible short-term dilution followed by longer-term amplification"
  ) +
  theme_classic()

print(p_transient)


# ------------------------------------------------------------------------------
# 10. Directionality: human movement aligned with vs opposed to waterborne flow
# ------------------------------------------------------------------------------

# Stronger waterborne direction 1 -> 2:
MP_strong <- 0.10
MP_weak   <- 0.01

directionality <- do.call(
  rbind,
  lapply(c("aligned", "opposed"), function(dir) {

    do.call(
      rbind,
      lapply(MH_grid, function(mh) {

        if (dir == "aligned") {
          # both dominant directions 1 -> 2
          pars <- make_pars(
            MH12 = mh, MH21 = 0,
            MP12 = MP_strong, MP21 = MP_weak
          )
        } else {
          # human movement 1 -> 2, but stronger waterborne flow 2 -> 1
          pars <- make_pars(
            MH12 = mh, MH21 = 0,
            MP12 = MP_weak, MP21 = MP_strong
          )
        }

        out <- simulate_model(pars, init = init_ref)
        s <- summarize_run(out, init_ref)
        s <- flag_depletion(s)

        cbind(
          direction = dir,
          MH12 = mh,
          R0 = calc_R0(pars),
          s
        )
      })
    )
  })
)

p_direction <- ggplot(directionality,
                      aes(x = MH12, y = prev_total_tail_mean,
                          linetype = direction)) +
  geom_line() +
  labs(
    x = expression(M[12]^H),
    y = "Total tail-mean prevalence",
    linetype = "Connectivity",
    title = "Human movement aligned with vs opposed to waterborne transport"
  ) +
  theme_classic()

print(p_direction)


# ------------------------------------------------------------------------------
# 11. R0 threshold analysis
# ------------------------------------------------------------------------------

r0_sweep <- data.frame(
  MH12 = MH_grid,
  R0 = sapply(MH_grid, function(mh) {
    pars <- make_pars(MH12 = mh, MH21 = 0,
                      MP12 = MP12_ref, MP21 = MP21_ref)
    calc_R0(pars)
  })
)

p_r0 <- ggplot(r0_sweep, aes(x = MH12, y = R0)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = 2) +
  labs(
    x = expression(M[12]^H),
    y = expression(R[0]),
    title = expression(R[0]~"threshold across human-mediated movement")
  ) +
  theme_classic()

print(p_r0)


# ------------------------------------------------------------------------------
# 12. Optional bidirectional human-movement grid
# ------------------------------------------------------------------------------

MH_bi_grid <- seq(0, 0.5, length.out = 26)

bidirectional_grid <- expand.grid(
  MH12 = MH_bi_grid,
  MH21 = MH_bi_grid
)

bidirectional_results <- do.call(
  rbind,
  lapply(seq_len(nrow(bidirectional_grid)), function(k) {

    mh12 <- bidirectional_grid$MH12[k]
    mh21 <- bidirectional_grid$MH21[k]

    pars <- make_pars(
      MH12 = mh12, MH21 = mh21,
      MP12 = MP12_ref, MP21 = MP21_ref
    )

    out <- simulate_model(pars, init = init_ref)
    s <- summarize_run(out, init_ref)
    s <- flag_depletion(s)

    cbind(
      MH12 = mh12,
      MH21 = mh21,
      R0 = calc_R0(pars),
      s
    )
  })
)


# ------------------------------------------------------------------------------
# 13. Hypothesis-driven robustness analyses
# ------------------------------------------------------------------------------

# 13A. MH x MP
# Question:
# Does human-mediated movement matter most when natural waterborne connectivity
# is relatively weak?

MP_grid <- seq(0, 0.15, length.out = 25)
MH_robust_grid <- seq(0, 0.5, length.out = 31)

MHxMP_grid <- expand.grid(
  MH12 = MH_robust_grid,
  MP12 = MP_grid
)

MHxMP_results <- do.call(
  rbind,
  lapply(seq_len(nrow(MHxMP_grid)), function(k) {

    mh <- MHxMP_grid$MH12[k]
    mp <- MHxMP_grid$MP12[k]

    pars <- make_pars(
      MH12 = mh,
      MH21 = 0,
      MP12 = mp,
      MP21 = MP21_ref
    )

    out <- simulate_model(pars, init = init_ref)
    s <- summarize_run(out, init_ref)
    s <- flag_depletion(s)

    cbind(
      MH12 = mh,
      MP12 = mp,
      R0 = calc_R0(pars),
      s
    )
  })
)


# 13B. MH x beta
# Question:
# Is any dilution-to-amplification transition controlled by how efficiently
# environmental pathogen exposure infects incoming susceptible fish?

beta_mult_grid <- exp(seq(log(0.1), log(10), length.out = 25))

MHxbeta_grid <- expand.grid(
  MH12 = MH_robust_grid,
  beta_multiplier = beta_mult_grid
)

MHxbeta_results <- do.call(
  rbind,
  lapply(seq_len(nrow(MHxbeta_grid)), function(k) {

    mh <- MHxbeta_grid$MH12[k]
    bm <- MHxbeta_grid$beta_multiplier[k]

    pars <- make_pars(
      MH12 = mh,
      MH21 = 0,
      beta1 = beta_ref * bm,
      beta2 = beta_ref * bm,
      MP12 = MP12_ref,
      MP21 = MP21_ref
    )

    out <- simulate_model(pars, init = init_ref)
    s <- summarize_run(out, init_ref)
    s <- flag_depletion(s)

    cbind(
      MH12 = mh,
      beta_multiplier = bm,
      beta = beta_ref * bm,
      R0 = calc_R0(pars),
      s
    )
  })
)


# 13C. Optional MH x delta
# Use only if environmental persistence becomes mechanistically important.

delta_mult_grid <- exp(seq(log(0.25), log(4), length.out = 21))

MHxdelta_grid <- expand.grid(
  MH12 = MH_robust_grid,
  delta_multiplier = delta_mult_grid
)

MHxdelta_results <- do.call(
  rbind,
  lapply(seq_len(nrow(MHxdelta_grid)), function(k) {

    mh <- MHxdelta_grid$MH12[k]
    dm <- MHxdelta_grid$delta_multiplier[k]

    pars <- make_pars(
      MH12 = mh,
      MH21 = 0,
      delta1 = delta_ref * dm,
      delta2 = delta_ref * dm,
      MP12 = MP12_ref,
      MP21 = MP21_ref
    )

    out <- simulate_model(pars, init = init_ref)
    s <- summarize_run(out, init_ref)
    s <- flag_depletion(s)

    cbind(
      MH12 = mh,
      delta_multiplier = dm,
      delta = delta_ref * dm,
      R0 = calc_R0(pars),
      s
    )
  })
)


# ------------------------------------------------------------------------------
# 14. Management interaction: MH x routine harvest h
# ------------------------------------------------------------------------------

# This directly addresses a management question:
# Can faster routine harvest/removal offset epidemiological risk from stronger
# human-mediated fish movement?

h_grid <- seq(0, 0.5, length.out = 31)

MHxh_grid <- expand.grid(
  MH12 = MH_robust_grid,
  h = h_grid
)

MHxh_results <- do.call(
  rbind,
  lapply(seq_len(nrow(MHxh_grid)), function(k) {

    mh <- MHxh_grid$MH12[k]
    h  <- MHxh_grid$h[k]

    pars <- make_pars(
      MH12 = mh,
      MH21 = 0,
      h1 = h,
      h2 = h,
      MP12 = MP12_ref,
      MP21 = MP21_ref
    )

    out <- simulate_model(pars, init = init_ref)
    s <- summarize_run(out, init_ref)
    s <- flag_depletion(s)

    cbind(
      MH12 = mh,
      h = h,
      R0 = calc_R0(pars),
      s
    )
  })
)

# R0 heatmap; contour at R0 = 1 gives a candidate management threshold
p_MHxh_R0 <- ggplot(MHxh_results, aes(x = MH12, y = h, fill = R0)) +
  geom_tile() +
  geom_contour(aes(z = R0), breaks = 1, inherit.aes = TRUE) +
  labs(
    x = expression(M[12]^H),
    y = "Routine harvest/removal rate h",
    fill = expression(R[0]),
    title = expression("Management interaction: " * M^H * " x h")
  ) +
  theme_classic()

print(p_MHxh_R0)


# ------------------------------------------------------------------------------
# 15. External stocking / openness robustness
# ------------------------------------------------------------------------------

# Question:
# Are qualitative MH effects dependent on the magnitude of external host input?

Lambda_multipliers <- c(0.5, 1, 2)

Lambda_robustness <- do.call(
  rbind,
  lapply(Lambda_multipliers, function(lm) {

    do.call(
      rbind,
      lapply(MH_grid, function(mh) {

        pars <- make_pars(
          Lambda1 = Lambda_ref * lm,
          Lambda2 = Lambda_ref * lm,
          MH12 = mh,
          MH21 = 0,
          MP12 = MP12_ref,
          MP21 = MP21_ref
        )

        out <- simulate_model(pars, init = init_ref)
        s <- summarize_run(out, init_ref)
        s <- flag_depletion(s)

        cbind(
          Lambda_multiplier = lm,
          MH12 = mh,
          R0 = calc_R0(pars),
          s
        )
      })
    )
  })
)


# ------------------------------------------------------------------------------
# 16. Diagnostics: distinguish disease control from population depletion
# ------------------------------------------------------------------------------

diagnose_result <- function(result_df,
                            depletion_threshold = 0.25,
                            prevalence_col = "prev_total_tail_mean") {

  result_df$depletion_warning <-
    result_df$N1_retention < depletion_threshold |
    result_df$N2_retention < depletion_threshold

  result_df
}

main_sweep <- diagnose_result(main_sweep)
source_state_sweep <- diagnose_result(source_state_sweep)
directionality <- diagnose_result(directionality)
bidirectional_results <- diagnose_result(bidirectional_results)
MHxMP_results <- diagnose_result(MHxMP_results)
MHxbeta_results <- diagnose_result(MHxbeta_results)
MHxdelta_results <- diagnose_result(MHxdelta_results)
MHxh_results <- diagnose_result(MHxh_results)
Lambda_robustness <- diagnose_result(Lambda_robustness)

# Quick count of parameter combinations with strong patch depletion
depletion_counts <- data.frame(
  analysis = c(
    "main_sweep",
    "source_state",
    "directionality",
    "bidirectional",
    "MHxMP",
    "MHxbeta",
    "MHxdelta",
    "MHxh",
    "Lambda_robustness"
  ),
  n_depletion_warnings = c(
    sum(main_sweep$depletion_warning),
    sum(source_state_sweep$depletion_warning),
    sum(directionality$depletion_warning),
    sum(bidirectional_results$depletion_warning),
    sum(MHxMP_results$depletion_warning),
    sum(MHxbeta_results$depletion_warning),
    sum(MHxdelta_results$depletion_warning),
    sum(MHxh_results$depletion_warning),
    sum(Lambda_robustness$depletion_warning)
  )
)

print(depletion_counts)


# ------------------------------------------------------------------------------
# 17. Optional mechanistic balance checks
# ------------------------------------------------------------------------------

# Total host abundance satisfies:
#
# d(N1+N2)/dt =
# Lambda1 + Lambda2
# - (mu1+h1)S1 - (alpha1+mu1+h1)I1
# - (mu2+h2)S2 - (alpha2+mu2+h2)I2
#
# Internal human-mediated movement cancels exactly.
#
# This helper numerically checks the total-host mass balance along a simulation.

check_total_host_balance <- function(out, pars) {

  p <- as.list(pars)

  expected_dNT <- with(out,
    p$Lambda1 + p$Lambda2 -
      (p$mu1 + p$h1) * S1 -
      (p$alpha1 + p$mu1 + p$h1) * I1 -
      (p$mu2 + p$h2) * S2 -
      (p$alpha2 + p$mu2 + p$h2) * I2
  )

  numerical_dNT <- c(NA, diff(out$NT) / diff(out$time))

  data.frame(
    time = out$time,
    expected_dNT = expected_dNT,
    numerical_dNT = numerical_dNT
  )
}

example_balance <- check_total_host_balance(
  simulate_model(pars_combined, init_ref),
  pars_combined
)


# ------------------------------------------------------------------------------
# 18. Export results
# ------------------------------------------------------------------------------

output_dir <- "ISA_two_patch_outputs"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

write.csv(reference_checks,
          file.path(output_dir, "reference_checks.csv"),
          row.names = FALSE)

write.csv(main_sweep,
          file.path(output_dir, "main_MH_sweep.csv"),
          row.names = FALSE)

write.csv(source_state_sweep,
          file.path(output_dir, "source_state_dependence.csv"),
          row.names = FALSE)

write.csv(directionality,
          file.path(output_dir, "directionality.csv"),
          row.names = FALSE)

write.csv(r0_sweep,
          file.path(output_dir, "R0_MH_sweep.csv"),
          row.names = FALSE)

write.csv(bidirectional_results,
          file.path(output_dir, "bidirectional_MH_grid.csv"),
          row.names = FALSE)

write.csv(MHxMP_results,
          file.path(output_dir, "MHxMP.csv"),
          row.names = FALSE)

write.csv(MHxbeta_results,
          file.path(output_dir, "MHxbeta.csv"),
          row.names = FALSE)

write.csv(MHxdelta_results,
          file.path(output_dir, "MHxdelta.csv"),
          row.names = FALSE)

write.csv(MHxh_results,
          file.path(output_dir, "MHxh.csv"),
          row.names = FALSE)

write.csv(Lambda_robustness,
          file.path(output_dir, "Lambda_robustness.csv"),
          row.names = FALSE)

write.csv(depletion_counts,
          file.path(output_dir, "depletion_diagnostics.csv"),
          row.names = FALSE)


# Save plots
ggsave(file.path(output_dir, "main_MH_prevalence.png"),
       p_main_prev, width = 7, height = 5, dpi = 300)

ggsave(file.path(output_dir, "main_MH_host_abundance.png"),
       p_main_N, width = 7, height = 5, dpi = 300)

ggsave(file.path(output_dir, "transient_patch2_prevalence.png"),
       p_transient, width = 7, height = 5, dpi = 300)

ggsave(file.path(output_dir, "directionality.png"),
       p_direction, width = 7, height = 5, dpi = 300)

ggsave(file.path(output_dir, "R0_threshold.png"),
       p_r0, width = 7, height = 5, dpi = 300)

ggsave(file.path(output_dir, "MHxharvest_R0.png"),
       p_MHxh_R0, width = 7, height = 5, dpi = 300)


# ------------------------------------------------------------------------------
# 19. Minimal recommended workflow for the paper
# ------------------------------------------------------------------------------

cat("
Recommended interpretation workflow
-----------------------------------
1. Run reference checks only as coding/sanity checks.
2. Use the one-way MH sweep as the central mechanistic result.
3. Compare source-worse vs destination-worse initial states.
4. Inspect transients for short-term dilution vs long-term amplification.
5. Compare aligned vs opposed MH and MP directions.
6. Examine R0 crossing 1.
7. Use only hypothesis-driven robustness analyses in the main manuscript.
8. Always inspect N1 and N2 alongside prevalence to detect patch depletion.
9. Treat the bidirectional grid as optional/supplementary unless it reveals a
   qualitatively new result.
10. Verify empirical values and units for Lambda, alpha, h, beta, chi, and delta
    before interpreting numerical thresholds.
")

# End of script


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
# Experiment 23
# Focused higher-order interaction analysis
#
# Host movement x source-patch transmission x stocking
# x directional waterborne pathogen movement
#
# Figure:
#   x-axis        = MH12
#   y-axis        = R0
#   facet rows    = r_beta1 = 1, 5
#   facet columns = r_Lambda = 1, 5
#   line type     = rP = 1, 15, 50
#
# Definitions:
#   r_beta1  = beta1 / beta_ref
#   r_Lambda = Lambda / Lambda_ref
#   rP       = MP12 / MP21
#
# beta2 is fixed at beta_ref.
# Lambda1 and Lambda2 vary together.
# MH21 = 0.
# MP21 is fixed locally at 0.05 / month.
# ============================================================
# EXPERIMENT 23
# Focused combined interaction analysis
#
# Human-mediated host movement x
# transmission heterogeneity x
# stocking intensity x
# waterborne pathogen movement
# ============================================================


# ------------------------------------------------------------
# 23A. Define parameter levels
# ------------------------------------------------------------

# Human-mediated host movement:
# host-exporting patch (Patch 1) ->
# host-receiving patch (Patch 2)

MH12_multi_grid <- seq(
  from = 0,
  to = 0.10,
  length.out = 51
)


# Transmission intensity in Patch 1
#
# r_beta1 = beta1 / beta_ref
#
# Patch 2 remains fixed at beta_ref.

r_beta1_levels <- c(
  1,
  5
)


# Stocking intensity in both patches
#
# r_Lambda = Lambda / Lambda_ref

r_Lambda_levels <- c(
  1,
  5
)


# Directional waterborne pathogen movement
#
# rP = MP12 / MP21

rP_levels <- c(
  1,
  15,
  50
)


# Reverse pathogen movement is fixed locally

MP21_multi_fixed <- 0.05


# Therefore:
#
# rP = 1
# MP12 = 1 x 0.05 = 0.05
#
# rP = 15
# MP12 = 15 x 0.05 = 0.75
#
# rP = 50
# MP12 = 50 x 0.05 = 2.50



# ------------------------------------------------------------
# 23B. Create complete parameter grid
# ------------------------------------------------------------

multi_factor_grid <- expand.grid(
  MH12 = MH12_multi_grid,
  r_beta1 = r_beta1_levels,
  r_Lambda = r_Lambda_levels,
  rP = rP_levels,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)


# Calculate forward pathogen movement
#
# rP = MP12 / MP21
#
# therefore:
#
# MP12 = rP x MP21

multi_factor_grid$MP12 <-
  multi_factor_grid$rP * MP21_multi_fixed



# ------------------------------------------------------------
# 23C. Calculate R0 for every parameter combination
# ------------------------------------------------------------

multi_factor_results <- do.call(
  rbind,
  lapply(
    seq_len(nrow(multi_factor_grid)),
    function(i) {
      
      g <- multi_factor_grid[i, ]
      
      
      pars_i <- make_pars(
        
        # Stocking intensity:
        # varied equally in both patches
        
        Lambda1 = Lambda_ref * g$r_Lambda,
        Lambda2 = Lambda_ref * g$r_Lambda,
        
        
        # Transmission:
        # Patch 1 varies;
        # Patch 2 remains at reference
        
        beta1 = beta_ref * g$r_beta1,
        beta2 = beta_ref,
        
        
        # Human-mediated host movement:
        # Patch 1 -> Patch 2
        
        MH12 = g$MH12,
        MH21 = 0,
        
        
        # Waterborne pathogen movement
        
        MP12 = g$MP12,
        MP21 = MP21_multi_fixed
      )
      
      
      # System-level basic reproduction number
      
      R0_i <- calc_R0(pars_i)
      
      
      data.frame(
        
        MH12 = g$MH12,
        
        r_beta1 = g$r_beta1,
        beta1 = beta_ref * g$r_beta1,
        beta2 = beta_ref,
        
        r_Lambda = g$r_Lambda,
        Lambda1 = Lambda_ref * g$r_Lambda,
        Lambda2 = Lambda_ref * g$r_Lambda,
        
        rP = g$rP,
        MP12 = g$MP12,
        MP21 = MP21_multi_fixed,
        
        R0 = R0_i
      )
    }
  )
)



# ------------------------------------------------------------
# 23D. Sanity checks
# ------------------------------------------------------------

# Expected number of parameter combinations:
#
# 51 MH12 values
# x 2 beta levels
# x 2 stocking levels
# x 3 rP levels
# = 612

stopifnot(
  nrow(multi_factor_results) == 612
)


stopifnot(
  all(
    sort(unique(multi_factor_results$rP)) ==
      c(1, 15, 50)
  )
)


stopifnot(
  all(
    sort(unique(multi_factor_results$r_beta1)) ==
      c(1, 5)
  )
)


stopifnot(
  all(
    sort(unique(multi_factor_results$r_Lambda)) ==
      c(1, 5)
  )
)



# ------------------------------------------------------------
# 23E. Define the four panels
# ------------------------------------------------------------

# A: r_beta1 = 1, r_Lambda = 1
# B: r_beta1 = 1, r_Lambda = 5
# C: r_beta1 = 5, r_Lambda = 1
# D: r_beta1 = 5, r_Lambda = 5

multi_factor_results$panel_f <- with(
  multi_factor_results,
  ifelse(
    r_beta1 == 1 & r_Lambda == 1,
    "A",
    ifelse(
      r_beta1 == 1 & r_Lambda == 5,
      "B",
      ifelse(
        r_beta1 == 5 & r_Lambda == 1,
        "C",
        "D"
      )
    )
  )
)


multi_factor_results$panel_f <- factor(
  multi_factor_results$panel_f,
  levels = c(
    "A",
    "B",
    "C",
    "D"
  )
)


# Check that only four panels exist

stopifnot(
  identical(
    levels(multi_factor_results$panel_f),
    c("A", "B", "C", "D")
  )
)


# Check number of observations in each panel
#
# Each panel:
# 51 MH values x 3 rP levels = 153

print(
  table(multi_factor_results$panel_f)
)



# ------------------------------------------------------------
# 23F. Save numerical results
# ------------------------------------------------------------

dir.create(
  "ISA_two_patch_outputs",
  showWarnings = FALSE
)


saveRDS(
  multi_factor_results,
  file = file.path(
    "ISA_two_patch_outputs",
    paste0(
      "focused_interaction_R0_",
      "beta1_1_5_",
      "Lambda_1_5_",
      "rP_1_15_50.rds"
    )
  )
)


write.csv(
  multi_factor_results,
  file = file.path(
    "ISA_two_patch_outputs",
    paste0(
      "focused_interaction_R0_",
      "beta1_1_5_",
      "Lambda_1_5_",
      "rP_1_15_50.csv"
    )
  ),
  row.names = FALSE
)



# ------------------------------------------------------------
# 23G. Mathematical panel labels
# ------------------------------------------------------------

panel_labeller <- as_labeller(
  c(
    A = "bold(A)~':'~r[beta*1]==1*','~r[Lambda]==1",
    B = "bold(B)~':'~r[beta*1]==1*','~r[Lambda]==5",
    C = "bold(C)~':'~r[beta*1]==5*','~r[Lambda]==1",
    D = "bold(D)~':'~r[beta*1]==5*','~r[Lambda]==5"
  ),
  label_parsed
)



# ------------------------------------------------------------
# 23H. Custom rP key in Panel A ONLY
# ------------------------------------------------------------

# Mathematical definition:
#
# r_P = M^P_12 / M^P_21
#
# Written on ONE line rather than as a stacked fraction.
#
# The P symbols are italic.


key_title <- data.frame(
  
  panel_f = factor(
    "A",
    levels = c(
      "A",
      "B",
      "C",
      "D"
    )
  ),
  
  x = 0.010,
  y = 19.0,
  
  label =
    "italic(r)[italic(P)] == italic(M)[12]^italic(P) / italic(M)[21]^italic(P)*':'"
)



# Solid line: rP = 1

key_line_1 <- data.frame(
  
  panel_f = factor(
    "A",
    levels = c(
      "A",
      "B",
      "C",
      "D"
    )
  ),
  
  x = 0.012,
  xend = 0.040,
  
  y = 14.5,
  yend = 14.5
)



# Dashed line: rP = 15

key_line_15 <- data.frame(
  
  panel_f = factor(
    "A",
    levels = c(
      "A",
      "B",
      "C",
      "D"
    )
  ),
  
  x = 0.012,
  xend = 0.040,
  
  y = 12.2,
  yend = 12.2
)



# Dotted line: rP = 50

key_line_50 <- data.frame(
  
  panel_f = factor(
    "A",
    levels = c(
      "A",
      "B",
      "C",
      "D"
    )
  ),
  
  x = 0.012,
  xend = 0.040,
  
  y = 9.9,
  yend = 9.9
)



# Numerical labels for the custom key

key_numbers <- data.frame(
  
  panel_f = factor(
    rep("A", 3),
    levels = c(
      "A",
      "B",
      "C",
      "D"
    )
  ),
  
  x = rep(
    0.045,
    3
  ),
  
  y = c(
    14.5,
    12.2,
    9.9
  ),
  
  label = c(
    "1",
    "15",
    "50"
  )
)



# ------------------------------------------------------------
# 23I. Manual panel-axis elements
# ------------------------------------------------------------

# ============================================================
# Upper-row x-axes: Panels A and B
# ============================================================
#
# Panels A and B do not repeat the numerical x-axis labels,
# but retain the horizontal x-axis baseline and tick marks.


top_axis_lines <- data.frame(
  panel_f = factor(
    c("A", "B"),
    levels = c("A", "B", "C", "D")
  ),
  
  x = c(0, 0),
  xend = c(0.10, 0.10),
  
  y = c(0, 0),
  yend = c(0, 0)
)


# X-axis tick positions for Panels A and B

top_tick_positions <- seq(
  from = 0,
  to = 0.10,
  by = 0.02
)


top_axis_ticks <- expand.grid(
  panel_f = c("A", "B"),
  x = top_tick_positions,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)


top_axis_ticks$panel_f <- factor(
  top_axis_ticks$panel_f,
  levels = c("A", "B", "C", "D")
)


top_axis_ticks$xend <- top_axis_ticks$x

top_axis_ticks$y <- 0

# Short ticks extending upward into the upper panels
top_axis_ticks$yend <- 0.38



# ============================================================
# Y-axes: Panels B and D
# ============================================================
#
# Panels B and D use the same y-axis scale as A and C.
#
# B and D retain:
#   - the same vertical y-axis line
#   - outward-facing tick marks
#
# but DO NOT repeat the numerical y-axis labels.


BD_y_axis_lines <- data.frame(
  panel_f = factor(
    c("B", "D"),
    levels = c("A", "B", "C", "D")
  ),
  
  x = c(0, 0),
  xend = c(0, 0),
  
  y = c(0, 0),
  yend = c(23, 23)
)


# Same y-axis tick positions as A and C

BD_y_tick_positions <- seq(
  from = 0,
  to = 20,
  by = 5
)


BD_y_axis_ticks <- expand.grid(
  panel_f = c("B", "D"),
  y = BD_y_tick_positions,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)


BD_y_axis_ticks$panel_f <- factor(
  BD_y_axis_ticks$panel_f,
  levels = c("A", "B", "C", "D")
)


# Outward-facing y-axis ticks.
#
# x = 0 is the y-axis.
# Negative xend sends the tick LEFT, outside the panel.
#
# Keep these short so that they visually match the
# standard ggplot ticks used for Panels A and C.

BD_y_axis_ticks$x <- 0
BD_y_axis_ticks$xend <- -0.0010

BD_y_axis_ticks$yend <- BD_y_axis_ticks$y



# ------------------------------------------------------------
# 23J. Create final four-panel R0 figure
# ------------------------------------------------------------

p_multi_R0 <- ggplot(
  multi_factor_results,
  aes(
    x = MH12,
    y = R0,
    linetype = factor(rP),
    group = factor(rP)
  )
) +
  
  
  # ==========================================================
# Main R0 curves
# ==========================================================

geom_line(
  colour = "black",
  linewidth = 1.05
) +
  
  
  # ==========================================================
# R0 = 1 invasion threshold
# ==========================================================

geom_hline(
  yintercept = 1,
  colour = "red",
  linetype = "solid",
  linewidth = 0.75
) +
  
  
  # ==========================================================
# Four panels ONLY: A, B, C, D
# ==========================================================

facet_wrap(
  ~ panel_f,
  nrow = 2,
  ncol = 2,
  labeller = panel_labeller
) +
  
  
  # ==========================================================
# rP line types
# ==========================================================

scale_linetype_manual(
  values = c(
    "1"  = "solid",
    "15" = "dashed",
    "50" = "dotted"
  ),
  guide = "none"
) +
  
  
  # ==========================================================
# X-AXIS BASELINES FOR PANELS A AND B
# ==========================================================

geom_segment(
  data = top_axis_lines,
  aes(
    x = x,
    xend = xend,
    y = y,
    yend = yend
  ),
  inherit.aes = FALSE,
  colour = "black",
  linewidth = 0.65,
  lineend = "butt"
) +
  
  
  # X-axis ticks for Panels A and B
  
  geom_segment(
    data = top_axis_ticks,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.55,
    lineend = "butt"
  ) +
  
  
  # ==========================================================
# Y-AXIS BASELINES FOR PANELS B AND D
# ==========================================================
#
# Same line thickness as the standard A/C y-axis.

geom_segment(
  data = BD_y_axis_lines,
  aes(
    x = x,
    xend = xend,
    y = y,
    yend = yend
  ),
  inherit.aes = FALSE,
  colour = "black",
  linewidth = 0.65,
  lineend = "butt"
) +
  
  
  # ==========================================================
# Y-AXIS TICKS FOR PANELS B AND D
# ==========================================================
#
# Same thickness as A/C ticks.
# Ticks point OUTWARD to the left.
# No numerical labels are added.

geom_segment(
  data = BD_y_axis_ticks,
  aes(
    x = x,
    xend = xend,
    y = y,
    yend = yend
  ),
  inherit.aes = FALSE,
  colour = "black",
  linewidth = 0.55,
  lineend = "butt"
) +
  
  
  # ==========================================================
# CUSTOM rP KEY IN PANEL A ONLY
# ==========================================================

# r_P = M^P_12 / M^P_21 :

geom_text(
  data = key_title,
  aes(
    x = x,
    y = y,
    label = label
  ),
  inherit.aes = FALSE,
  parse = TRUE,
  hjust = 0,
  size = 4.5
) +
  
  
  # Solid line -> rP = 1
  
  geom_segment(
    data = key_line_1,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.15,
    linetype = "solid",
    lineend = "butt"
  ) +
  
  
  # Dashed line -> rP = 15
  
  geom_segment(
    data = key_line_15,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.15,
    linetype = "dashed",
    lineend = "butt"
  ) +
  
  
  # Dotted line -> rP = 50
  
  geom_segment(
    data = key_line_50,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.15,
    linetype = "dotted",
    lineend = "butt"
  ) +
  
  
  # Numerical labels: 1 / 15 / 50
  
  geom_text(
    data = key_numbers,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 4.2
  ) +
  
  
  # ==========================================================
# X-AXIS
# ==========================================================

scale_x_continuous(
  breaks = seq(
    from = 0,
    to = 0.10,
    by = 0.02
  ),
  
  limits = c(
    0,
    0.10
  ),
  
  expand = expansion(
    mult = c(
      0,
      0.02
    )
  )
) +
  
  
  # ==========================================================
# Y-AXIS
# ==========================================================
#
# Same y scale for all four panels.

scale_y_continuous(
  limits = c(
    0,
    23
  ),
  
  breaks = seq(
    from = 0,
    to = 20,
    by = 5
  ),
  
  expand = expansion(
    mult = c(
      0,
      0.02
    )
  )
) +
  
  
  # ==========================================================
# Allow B/D outward y-axis ticks to extend outside panel
# ==========================================================

coord_cartesian(
  clip = "off"
) +
  
  
  # ==========================================================
# AXIS TITLES
# ==========================================================

labs(
  x = expression(
    "Human-mediated host movement, " ~
      M[12]^H ~
      " (month"^{-1}*")"
  ),
  
  y = expression(
    R[0]
  )
) +
  
  
  # ==========================================================
# CLEAN MANUSCRIPT-STYLE THEME
# ==========================================================

theme_classic(
  base_size = 13
) +
  
  
  theme(
    
    # --------------------------------------------------------
    # No background grid
    # --------------------------------------------------------
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    
    # --------------------------------------------------------
    # No grey facet boxes
    # --------------------------------------------------------
    
    strip.background = element_blank(),
    
    
    # --------------------------------------------------------
    # Panel headings
    # --------------------------------------------------------
    
    strip.text = element_text(
      size = 12,
      margin = margin(
        b = 7
      )
    ),
    
    
    # --------------------------------------------------------
    # No standard legend
    # --------------------------------------------------------
    
    legend.position = "none",
    
    
    # --------------------------------------------------------
    # Axis titles
    # --------------------------------------------------------
    
    axis.title = element_text(
      size = 12
    ),
    
    
    # --------------------------------------------------------
    # Axis tick labels
    # --------------------------------------------------------
    
    axis.text = element_text(
      size = 10
    ),
    
    
    # --------------------------------------------------------
    # Standard axis lines
    #
    # A/C use these directly.
    # B/D manual y-axis lines use the SAME linewidth = 0.65.
    # --------------------------------------------------------
    
    axis.line = element_line(
      linewidth = 0.65,
      colour = "black"
    ),
    
    
    # --------------------------------------------------------
    # Standard axis ticks
    #
    # A/C use these directly.
    # B/D manual y ticks use the SAME linewidth = 0.55.
    # --------------------------------------------------------
    
    axis.ticks = element_line(
      linewidth = 0.55,
      colour = "black"
    ),
    
    
    # Fix standard A/C tick length explicitly
    
    axis.ticks.length = unit(
      0.10,
      "cm"
    ),
    
    
    # --------------------------------------------------------
    # Space between panels
    # --------------------------------------------------------
    
    panel.spacing = unit(
      1.0,
      "lines"
    ),
    
    
    # --------------------------------------------------------
    # Outer margins
    # --------------------------------------------------------
    
    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )



# ------------------------------------------------------------
# 23K. Display figure
# ------------------------------------------------------------

print(
  p_multi_R0
)



# ------------------------------------------------------------
# 23L. Save figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    "ISA_two_patch_outputs",
    "focused_interaction_R0_2x2_final.png"
  ),
  
  plot = p_multi_R0,
  
  width = 9,
  height = 7,
  dpi = 400
)

# ============================================================
# EXPERIMENT 23
# FINAL PLOTTING CODE FOR DYNAMIC OUTCOMES
#
# Assumes that multi_dynamic_results has already been created
# by the Experiment 23 ODE simulation block.
#
# Figures:
#   Figure 2 = Long-term total fish abundance
#   Figure 3 = Long-term total disease prevalence
#   Figure 4 = Patch-specific long-term disease prevalence
# ============================================================


# ============================================================
# COMMON SETTINGS
# ============================================================

# ------------------------------------------------------------
# Mathematical panel labels
# ------------------------------------------------------------

panel_labeller <- as_labeller(
  c(
    A = "bold(A)~':'~r[beta*1]==1*','~r[Lambda]==1",
    B = "bold(B)~':'~r[beta*1]==1*','~r[Lambda]==5",
    C = "bold(C)~':'~r[beta*1]==5*','~r[Lambda]==1",
    D = "bold(D)~':'~r[beta*1]==5*','~r[Lambda]==5"
  ),
  label_parsed
)


# ------------------------------------------------------------
# Common x-axis positions
# ------------------------------------------------------------

x_breaks_common <- seq(
  from = 0,
  to = 0.10,
  by = 0.02
)



# ============================================================
# FIGURE 2
# LONG-TERM TOTAL FISH ABUNDANCE
# ============================================================


# ------------------------------------------------------------
# Figure 2A. Determine common y-axis range
#
# Fish abundance is plotted in millions of fish.
# For example:
#
# 1,000,000 -> 1.0
# 1,500,000 -> 1.5
#
# The y-axis multiplier is shown separately as x 10^6.
# ------------------------------------------------------------

fish_plot_data <- multi_dynamic_results

fish_plot_data$NT_million <-
  fish_plot_data$NT_tail_mean / 1e6


fish_y_max <- max(
  fish_plot_data$NT_million,
  na.rm = TRUE
)


# Give some room above the curves for the Panel C legend.

fish_y_limit <- max(
  2.0,
  ceiling(
    fish_y_max * 4
  ) / 4
)


# ------------------------------------------------------------
# Figure 2B. Manual axes
# ------------------------------------------------------------

# Upper-row x-axis baselines: A and B

fish_top_x_lines <- data.frame(
  panel_f = factor(
    c("A", "B"),
    levels = c("A", "B", "C", "D")
  ),
  x = 0,
  xend = 0.10,
  y = 0,
  yend = 0
)


fish_top_x_ticks <- expand.grid(
  panel_f = c("A", "B"),
  x = x_breaks_common,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fish_top_x_ticks$panel_f <- factor(
  fish_top_x_ticks$panel_f,
  levels = c("A", "B", "C", "D")
)

fish_top_x_ticks$xend <- fish_top_x_ticks$x
fish_top_x_ticks$y <- 0

# Short upward ticks for upper panels
fish_top_x_ticks$yend <- fish_y_limit * 0.018


# B/D y-axis baselines

fish_BD_y_lines <- data.frame(
  panel_f = factor(
    c("B", "D"),
    levels = c("A", "B", "C", "D")
  ),
  x = 0,
  xend = 0,
  y = 0,
  yend = fish_y_limit
)


# Use simple common breaks in millions

fish_y_breaks <- pretty(
  c(
    0,
    fish_y_limit
  ),
  n = 4
)

fish_y_breaks <- fish_y_breaks[
  fish_y_breaks >= 0 &
    fish_y_breaks <= fish_y_limit
]


fish_BD_y_ticks <- expand.grid(
  panel_f = c("B", "D"),
  y = fish_y_breaks,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fish_BD_y_ticks$panel_f <- factor(
  fish_BD_y_ticks$panel_f,
  levels = c("A", "B", "C", "D")
)

fish_BD_y_ticks$x <- 0

# Outward = left
fish_BD_y_ticks$xend <- -0.0010

fish_BD_y_ticks$yend <-
  fish_BD_y_ticks$y



# ------------------------------------------------------------
# Figure 2C. rP legend in PANEL C
# ------------------------------------------------------------

fish_key_title <- data.frame(
  panel_f = factor(
    "C",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.010,
  y = fish_y_limit * 0.91,
  label =
    "italic(r)[italic(P)] == italic(M)[12]^italic(P) / italic(M)[21]^italic(P)*':'"
)


fish_key_1 <- data.frame(
  panel_f = factor(
    "C",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.012,
  xend = 0.040,
  y = fish_y_limit * 0.79,
  yend = fish_y_limit * 0.79
)


fish_key_15 <- data.frame(
  panel_f = factor(
    "C",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.012,
  xend = 0.040,
  y = fish_y_limit * 0.69,
  yend = fish_y_limit * 0.69
)


fish_key_50 <- data.frame(
  panel_f = factor(
    "C",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.012,
  xend = 0.040,
  y = fish_y_limit * 0.59,
  yend = fish_y_limit * 0.59
)


fish_key_numbers <- data.frame(
  panel_f = factor(
    rep("C", 3),
    levels = c("A", "B", "C", "D")
  ),
  x = rep(
    0.045,
    3
  ),
  y = fish_y_limit * c(
    0.79,
    0.69,
    0.59
  ),
  label = c(
    "1",
    "15",
    "50"
  )
)


# x 10^6 multiplier
#
# Put only in Panel A because all panels share the same scale.

fish_multiplier <- data.frame(
  panel_f = factor(
    "A",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.001,
  y = fish_y_limit * 0.985,
  label = "x~10^6"
)



# ------------------------------------------------------------
# Figure 2D. Plot
# ------------------------------------------------------------

p_total_fish <- ggplot(
  fish_plot_data,
  aes(
    x = MH12,
    y = NT_million,
    linetype = factor(rP),
    group = factor(rP)
  )
) +
  
  
  # Main curves
  
  geom_line(
    colour = "black",
    linewidth = 1.05
  ) +
  
  
  # Four panels
  
  facet_wrap(
    ~ panel_f,
    nrow = 2,
    ncol = 2,
    labeller = panel_labeller
  ) +
  
  
  # rP line types
  
  scale_linetype_manual(
    values = c(
      "1" = "solid",
      "15" = "dashed",
      "50" = "dotted"
    ),
    guide = "none"
  ) +
  
  
  # ----------------------------------------------------------
# Upper x-axes: A/B
# ----------------------------------------------------------

geom_segment(
  data = fish_top_x_lines,
  aes(
    x = x,
    xend = xend,
    y = y,
    yend = yend
  ),
  inherit.aes = FALSE,
  linewidth = 0.65,
  colour = "black",
  lineend = "butt"
) +
  
  
  geom_segment(
    data = fish_top_x_ticks,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    colour = "black",
    lineend = "butt"
  ) +
  
  
  # ----------------------------------------------------------
# B/D y-axes
# ----------------------------------------------------------

geom_segment(
  data = fish_BD_y_lines,
  aes(
    x = x,
    xend = xend,
    y = y,
    yend = yend
  ),
  inherit.aes = FALSE,
  linewidth = 0.65,
  colour = "black",
  lineend = "butt"
) +
  
  
  geom_segment(
    data = fish_BD_y_ticks,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    colour = "black",
    lineend = "butt"
  ) +
  
  
  # ----------------------------------------------------------
# Scientific multiplier x 10^6
# ----------------------------------------------------------

geom_text(
  data = fish_multiplier,
  aes(
    x = x,
    y = y,
    label = label
  ),
  inherit.aes = FALSE,
  parse = TRUE,
  hjust = 0,
  vjust = 1,
  size = 4.0
) +
  
  
  # ----------------------------------------------------------
# rP key: Panel C
# ----------------------------------------------------------

geom_text(
  data = fish_key_title,
  aes(
    x = x,
    y = y,
    label = label
  ),
  inherit.aes = FALSE,
  parse = TRUE,
  hjust = 0,
  size = 4.3
) +
  
  
  geom_segment(
    data = fish_key_1,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.15,
    linetype = "solid"
  ) +
  
  
  geom_segment(
    data = fish_key_15,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.15,
    linetype = "dashed"
  ) +
  
  
  geom_segment(
    data = fish_key_50,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.15,
    linetype = "dotted"
  ) +
  
  
  geom_text(
    data = fish_key_numbers,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 4.1
  ) +
  
  
  # ----------------------------------------------------------
# Scales
# ----------------------------------------------------------

scale_x_continuous(
  breaks = x_breaks_common,
  limits = c(
    0,
    0.10
  ),
  expand = expansion(
    mult = c(
      0,
      0.02
    )
  )
) +
  
  
  scale_y_continuous(
    limits = c(
      0,
      fish_y_limit
    ),
    breaks = fish_y_breaks,
    labels = function(x) {
      sprintf(
        "%.1f",
        x
      )
    },
    expand = expansion(
      mult = c(
        0,
        0.01
      )
    )
  ) +
  
  
  coord_cartesian(
    clip = "off"
  ) +
  
  
  labs(
    x = expression(
      "Human-mediated host movement, " ~
        M[12]^H ~
        " (month"^{-1}*")"
    ),
    y = "Long-term total fish abundance"
  ) +
  
  
  theme_classic(
    base_size = 13
  ) +
  
  
  theme(
    panel.grid = element_blank(),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 12,
      margin = margin(
        b = 7
      )
    ),
    
    legend.position = "none",
    
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
    
    panel.spacing = unit(
      1.0,
      "lines"
    ),
    
    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )


print(
  p_total_fish
)


ggsave(
  filename = file.path(
    "ISA_two_patch_outputs",
    "focused_interaction_total_fish_2x2_final.png"
  ),
  plot = p_total_fish,
  width = 9,
  height = 7,
  dpi = 400
)



# ============================================================
# FIGURE 3
# LONG-TERM TOTAL DISEASE PREVALENCE
# ============================================================


# ------------------------------------------------------------
# Figure 3A. Prepare plotting data
# ------------------------------------------------------------

total_prev_plot_data <- multi_dynamic_results

total_prev_plot_data$prev_percent <-
  total_prev_plot_data$prev_total_tail_mean * 100


prev_y_max <- max(
  total_prev_plot_data$prev_percent,
  na.rm = TRUE
)


# Round upward and provide some empty space.

prev_y_limit <- max(
  25,
  ceiling(
    prev_y_max / 5
  ) * 5
)


prev_y_breaks <- seq(
  from = 0,
  to = prev_y_limit,
  by = 5
)



# ------------------------------------------------------------
# Figure 3B. Manual axes
# ------------------------------------------------------------

prev_top_x_lines <- data.frame(
  panel_f = factor(
    c("A", "B"),
    levels = c("A", "B", "C", "D")
  ),
  x = 0,
  xend = 0.10,
  y = 0,
  yend = 0
)


prev_top_x_ticks <- expand.grid(
  panel_f = c("A", "B"),
  x = x_breaks_common,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

prev_top_x_ticks$panel_f <- factor(
  prev_top_x_ticks$panel_f,
  levels = c("A", "B", "C", "D")
)

prev_top_x_ticks$xend <- prev_top_x_ticks$x
prev_top_x_ticks$y <- 0
prev_top_x_ticks$yend <- prev_y_limit * 0.018


prev_BD_y_lines <- data.frame(
  panel_f = factor(
    c("B", "D"),
    levels = c("A", "B", "C", "D")
  ),
  x = 0,
  xend = 0,
  y = 0,
  yend = prev_y_limit
)


prev_BD_y_ticks <- expand.grid(
  panel_f = c("B", "D"),
  y = prev_y_breaks,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

prev_BD_y_ticks$panel_f <- factor(
  prev_BD_y_ticks$panel_f,
  levels = c("A", "B", "C", "D")
)

prev_BD_y_ticks$x <- 0
prev_BD_y_ticks$xend <- -0.0010
prev_BD_y_ticks$yend <- prev_BD_y_ticks$y



# ------------------------------------------------------------
# Figure 3C. rP key in PANEL A
# ------------------------------------------------------------

prev_key_title <- data.frame(
  panel_f = factor(
    "A",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.010,
  y = prev_y_limit * 0.88,
  label =
    "italic(r)[italic(P)] == italic(M)[12]^italic(P) / italic(M)[21]^italic(P)*':'"
)


prev_key_1 <- data.frame(
  panel_f = factor(
    "A",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.012,
  xend = 0.040,
  y = prev_y_limit * 0.72,
  yend = prev_y_limit * 0.72
)


prev_key_15 <- data.frame(
  panel_f = factor(
    "A",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.012,
  xend = 0.040,
  y = prev_y_limit * 0.61,
  yend = prev_y_limit * 0.61
)


prev_key_50 <- data.frame(
  panel_f = factor(
    "A",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.012,
  xend = 0.040,
  y = prev_y_limit * 0.50,
  yend = prev_y_limit * 0.50
)


prev_key_numbers <- data.frame(
  panel_f = factor(
    rep("A", 3),
    levels = c("A", "B", "C", "D")
  ),
  x = rep(
    0.045,
    3
  ),
  y = prev_y_limit * c(
    0.72,
    0.61,
    0.50
  ),
  label = c(
    "1",
    "15",
    "50"
  )
)



# ------------------------------------------------------------
# Figure 3D. Plot
# ------------------------------------------------------------

p_total_prev <- ggplot(
  total_prev_plot_data,
  aes(
    x = MH12,
    y = prev_percent,
    linetype = factor(rP),
    group = factor(rP)
  )
) +
  
  
  geom_line(
    colour = "black",
    linewidth = 1.05
  ) +
  
  
  facet_wrap(
    ~ panel_f,
    nrow = 2,
    ncol = 2,
    labeller = panel_labeller
  ) +
  
  
  scale_linetype_manual(
    values = c(
      "1" = "solid",
      "15" = "dashed",
      "50" = "dotted"
    ),
    guide = "none"
  ) +
  
  
  # Upper x-axes
  
  geom_segment(
    data = prev_top_x_lines,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 0.65,
    colour = "black"
  ) +
  
  
  geom_segment(
    data = prev_top_x_ticks,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    colour = "black"
  ) +
  
  
  # B/D y-axes
  
  geom_segment(
    data = prev_BD_y_lines,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 0.65,
    colour = "black"
  ) +
  
  
  geom_segment(
    data = prev_BD_y_ticks,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    colour = "black"
  ) +
  
  
  # rP title
  
  geom_text(
    data = prev_key_title,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = 0,
    size = 4.3
  ) +
  
  
  # rP lines
  
  geom_segment(
    data = prev_key_1,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 1.15,
    colour = "black",
    linetype = "solid"
  ) +
  
  
  geom_segment(
    data = prev_key_15,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 1.15,
    colour = "black",
    linetype = "dashed"
  ) +
  
  
  geom_segment(
    data = prev_key_50,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 1.15,
    colour = "black",
    linetype = "dotted"
  ) +
  
  
  geom_text(
    data = prev_key_numbers,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4.1
  ) +
  
  
  scale_x_continuous(
    breaks = x_breaks_common,
    limits = c(
      0,
      0.10
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  
  
  scale_y_continuous(
    breaks = prev_y_breaks,
    limits = c(
      0,
      prev_y_limit
    ),
    expand = expansion(
      mult = c(
        0,
        0.01
      )
    )
  ) +
  
  
  coord_cartesian(
    clip = "off"
  ) +
  
  
  labs(
    x = expression(
      "Human-mediated host movement, " ~
        M[12]^H ~
        " (month"^{-1}*")"
    ),
    y = "Long-term total disease prevalence (%)"
  ) +
  
  
  theme_classic(
    base_size = 13
  ) +
  
  
  theme(
    panel.grid = element_blank(),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 12,
      margin = margin(
        b = 7
      )
    ),
    
    legend.position = "none",
    
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
    
    panel.spacing = unit(
      1.0,
      "lines"
    ),
    
    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )


print(
  p_total_prev
)


ggsave(
  filename = file.path(
    "ISA_two_patch_outputs",
    "focused_interaction_total_prevalence_2x2_final.png"
  ),
  plot = p_total_prev,
  width = 9,
  height = 7,
  dpi = 400
)



# ============================================================
# FIGURE 4
# PATCH-SPECIFIC LONG-TERM DISEASE PREVALENCE
# ============================================================


# ------------------------------------------------------------
# Figure 4A. Convert to long format
# ------------------------------------------------------------

patch_prev_long <- rbind(
  
  data.frame(
    multi_dynamic_results[
      ,
      c(
        "MH12",
        "r_beta1",
        "r_Lambda",
        "rP",
        "panel_f"
      )
    ],
    
    patch = "Host-exporting patch",
    
    prevalence =
      multi_dynamic_results$prev1_tail_mean
  ),
  
  
  data.frame(
    multi_dynamic_results[
      ,
      c(
        "MH12",
        "r_beta1",
        "r_Lambda",
        "rP",
        "panel_f"
      )
    ],
    
    patch = "Host-receiving patch",
    
    prevalence =
      multi_dynamic_results$prev2_tail_mean
  )
)


patch_prev_long$patch <- factor(
  patch_prev_long$patch,
  levels = c(
    "Host-exporting patch",
    "Host-receiving patch"
  )
)


patch_prev_long$prev_percent <-
  patch_prev_long$prevalence * 100



# ------------------------------------------------------------
# Figure 4B. Common y scale
# ------------------------------------------------------------

patch_y_max <- max(
  patch_prev_long$prev_percent,
  na.rm = TRUE
)


patch_y_limit <- max(
  50,
  ceiling(
    patch_y_max / 10
  ) * 10
)


patch_y_breaks <- seq(
  from = 0,
  to = patch_y_limit,
  by = 10
)



# ------------------------------------------------------------
# Figure 4C. Manual axes
# ------------------------------------------------------------

patch_top_x_lines <- data.frame(
  panel_f = factor(
    c("A", "B"),
    levels = c("A", "B", "C", "D")
  ),
  x = 0,
  xend = 0.10,
  y = 0,
  yend = 0
)


patch_top_x_ticks <- expand.grid(
  panel_f = c("A", "B"),
  x = x_breaks_common,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

patch_top_x_ticks$panel_f <- factor(
  patch_top_x_ticks$panel_f,
  levels = c("A", "B", "C", "D")
)

patch_top_x_ticks$xend <- patch_top_x_ticks$x
patch_top_x_ticks$y <- 0
patch_top_x_ticks$yend <- patch_y_limit * 0.018


patch_BD_y_lines <- data.frame(
  panel_f = factor(
    c("B", "D"),
    levels = c("A", "B", "C", "D")
  ),
  x = 0,
  xend = 0,
  y = 0,
  yend = patch_y_limit
)


patch_BD_y_ticks <- expand.grid(
  panel_f = c("B", "D"),
  y = patch_y_breaks,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

patch_BD_y_ticks$panel_f <- factor(
  patch_BD_y_ticks$panel_f,
  levels = c("A", "B", "C", "D")
)

patch_BD_y_ticks$x <- 0
patch_BD_y_ticks$xend <- -0.0010
patch_BD_y_ticks$yend <- patch_BD_y_ticks$y



# ============================================================
# Figure 4D.
# rP LEGEND IN PANEL A
# ============================================================

patch_rP_title <- data.frame(
  panel_f = factor(
    "A",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.010,
  y = patch_y_limit * 0.88,
  label =
    "italic(r)[italic(P)] == italic(M)[12]^italic(P) / italic(M)[21]^italic(P)*':'"
)


patch_rP_1 <- data.frame(
  panel_f = factor(
    "A",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.012,
  xend = 0.040,
  y = patch_y_limit * 0.73,
  yend = patch_y_limit * 0.73
)


patch_rP_15 <- data.frame(
  panel_f = factor(
    "A",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.012,
  xend = 0.040,
  y = patch_y_limit * 0.62,
  yend = patch_y_limit * 0.62
)


patch_rP_50 <- data.frame(
  panel_f = factor(
    "A",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.012,
  xend = 0.040,
  y = patch_y_limit * 0.51,
  yend = patch_y_limit * 0.51
)


patch_rP_numbers <- data.frame(
  panel_f = factor(
    rep("A", 3),
    levels = c("A", "B", "C", "D")
  ),
  x = rep(
    0.045,
    3
  ),
  y = patch_y_limit * c(
    0.73,
    0.62,
    0.51
  ),
  label = c(
    "1",
    "15",
    "50"
  )
)



# ============================================================
# Figure 4E.
# PATCH COLOUR LEGEND IN PANEL C
# ============================================================

patch_colour_title <- data.frame(
  panel_f = factor(
    "C",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.010,
  y = patch_y_limit * 0.88,
  label = "Patch:"
)


# Host-exporting patch = RED

patch_export_key_line <- data.frame(
  panel_f = factor(
    "C",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.012,
  xend = 0.035,
  y = patch_y_limit * 0.73,
  yend = patch_y_limit * 0.73
)


patch_export_key_text <- data.frame(
  panel_f = factor(
    "C",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.040,
  y = patch_y_limit * 0.73,
  label = "Host-exporting patch"
)


# Host-receiving patch = BLUE

patch_receive_key_line <- data.frame(
  panel_f = factor(
    "C",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.012,
  xend = 0.035,
  y = patch_y_limit * 0.61,
  yend = patch_y_limit * 0.61
)


patch_receive_key_text <- data.frame(
  panel_f = factor(
    "C",
    levels = c("A", "B", "C", "D")
  ),
  x = 0.040,
  y = patch_y_limit * 0.61,
  label = "Host-receiving patch"
)



# ------------------------------------------------------------
# Figure 4F. Plot
# ------------------------------------------------------------

p_patch_prev <- ggplot(
  patch_prev_long,
  aes(
    x = MH12,
    y = prev_percent,
    
    # Patch identity = colour
    colour = patch,
    
    # Waterborne movement ratio = line type
    linetype = factor(rP),
    
    group = interaction(
      patch,
      rP
    )
  )
) +
  
  
  # ----------------------------------------------------------
# Main curves
# ----------------------------------------------------------

geom_line(
  linewidth = 1.05
) +
  
  
  # ----------------------------------------------------------
# Four panels
# ----------------------------------------------------------

facet_wrap(
  ~ panel_f,
  nrow = 2,
  ncol = 2,
  labeller = panel_labeller
) +
  
  
  # ----------------------------------------------------------
# Patch colours
#
# Patch 1 / host-exporting = RED
# Patch 2 / host-receiving = BLUE
# ----------------------------------------------------------

scale_colour_manual(
  values = c(
    "Host-exporting patch" = "red",
    "Host-receiving patch" = "blue"
  ),
  guide = "none"
) +
  
  
  # ----------------------------------------------------------
# rP line types
# ----------------------------------------------------------

scale_linetype_manual(
  values = c(
    "1" = "solid",
    "15" = "dashed",
    "50" = "dotted"
  ),
  guide = "none"
) +
  
  
  # ==========================================================
# Upper x-axes: A/B
# ==========================================================

geom_segment(
  data = patch_top_x_lines,
  aes(
    x = x,
    xend = xend,
    y = y,
    yend = yend
  ),
  inherit.aes = FALSE,
  linewidth = 0.65,
  colour = "black"
) +
  
  
  geom_segment(
    data = patch_top_x_ticks,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    colour = "black"
  ) +
  
  
  # ==========================================================
# B/D y-axes
# ==========================================================

geom_segment(
  data = patch_BD_y_lines,
  aes(
    x = x,
    xend = xend,
    y = y,
    yend = yend
  ),
  inherit.aes = FALSE,
  linewidth = 0.65,
  colour = "black"
) +
  
  
  geom_segment(
    data = patch_BD_y_ticks,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    colour = "black"
  ) +
  
  
  # ==========================================================
# rP LEGEND — PANEL A
# ==========================================================

geom_text(
  data = patch_rP_title,
  aes(
    x = x,
    y = y,
    label = label
  ),
  inherit.aes = FALSE,
  parse = TRUE,
  hjust = 0,
  size = 4.3,
  colour = "black"
) +
  
  
  # rP = 1
  
  geom_segment(
    data = patch_rP_1,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.15,
    linetype = "solid"
  ) +
  
  
  # rP = 15
  
  geom_segment(
    data = patch_rP_15,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.15,
    linetype = "dashed"
  ) +
  
  
  # rP = 50
  
  geom_segment(
    data = patch_rP_50,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.15,
    linetype = "dotted"
  ) +
  
  
  geom_text(
    data = patch_rP_numbers,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4.1,
    colour = "black"
  ) +
  
  
  # ==========================================================
# PATCH COLOUR LEGEND — PANEL C
# ==========================================================

geom_text(
  data = patch_colour_title,
  aes(
    x = x,
    y = y,
    label = label
  ),
  inherit.aes = FALSE,
  hjust = 0,
  fontface = "bold",
  size = 4.3,
  colour = "black"
) +
  
  
  # Red line = host-exporting patch
  
  geom_segment(
    data = patch_export_key_line,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "red",
    linewidth = 1.3,
    linetype = "solid"
  ) +
  
  
  geom_text(
    data = patch_export_key_text,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.8,
    colour = "black"
  ) +
  
  
  # Blue line = host-receiving patch
  
  geom_segment(
    data = patch_receive_key_line,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "blue",
    linewidth = 1.3,
    linetype = "solid"
  ) +
  
  
  geom_text(
    data = patch_receive_key_text,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.8,
    colour = "black"
  ) +
  
  
  # ==========================================================
# Scales
# ==========================================================

scale_x_continuous(
  breaks = x_breaks_common,
  limits = c(
    0,
    0.10
  ),
  expand = expansion(
    mult = c(
      0,
      0.02
    )
  )
) +
  
  
  scale_y_continuous(
    breaks = patch_y_breaks,
    limits = c(
      0,
      patch_y_limit
    ),
    expand = expansion(
      mult = c(
        0,
        0.01
      )
    )
  ) +
  
  
  coord_cartesian(
    clip = "off"
  ) +
  
  
  # ----------------------------------------------------------
# Axis labels
# ----------------------------------------------------------

labs(
  x = expression(
    "Human-mediated host movement, " ~
      M[12]^H ~
      " (month"^{-1}*")"
  ),
  
  y = "Long-term disease prevalence (%)"
) +
  
  
  # ----------------------------------------------------------
# Theme
# ----------------------------------------------------------

theme_classic(
  base_size = 13
) +
  
  
  theme(
    panel.grid = element_blank(),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 12,
      margin = margin(
        b = 7
      )
    ),
    
    # Absolutely no automatic legend at bottom
    legend.position = "none",
    
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
    
    panel.spacing = unit(
      1.0,
      "lines"
    ),
    
    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )



# ------------------------------------------------------------
# Display Figure 4
# ------------------------------------------------------------

print(
  p_patch_prev
)



# ------------------------------------------------------------
# Save Figure 4
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    "ISA_two_patch_outputs",
    "focused_interaction_patch_prevalence_2x2_final.png"
  ),
  plot = p_patch_prev,
  width = 9,
  height = 7,
  dpi = 400
)