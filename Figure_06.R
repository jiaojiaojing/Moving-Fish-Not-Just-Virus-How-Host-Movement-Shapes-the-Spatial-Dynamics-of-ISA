# Appendix Figure 6: mechanistic responses of the host-exporting patch
# MATCHES original R script: 120 months, dt=0.1, final 10% tail means.
# Reproducible two-patch ISA model, using deSolve and ggplot2.
# Install once if needed: install.packages(c('deSolve','ggplot2','dplyr','tidyr'))

library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

cycle_months <- 24
Lambda_ref <- mean(c(493233, 1108106)) / cycle_months
mu_ref <- 0.0127
h_ref <- 1 / cycle_months
beta_ref <- mean(c(1e-8, 1e-6))
alpha_ref <- 1.43
V_effect <- 20000
chi_ref <- mean(c(1e4, 1e6)) / V_effect
delta_ref <- mean(c(1.5, 10.5))
MP12_ref <- 0.05
MP21_ref <- 0.01
MH21_ref <- 0

# Full system. Rates are per month; P_i is environmental concentration.
isa_two_patch <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    dS1 <- Lambda1 - beta1*S1*P1 - (mu1+h1)*S1 - MH12*S1 + MH21*S2
    dI1 <- beta1*S1*P1 - (alpha1+mu1+h1)*I1 - MH12*I1 + MH21*I2
    dP1 <- chi1*I1 - delta1*P1 - MP12*P1 + MP21*P2
    dS2 <- Lambda2 - beta2*S2*P2 - (mu2+h2)*S2 + MH12*S1 - MH21*S2
    dI2 <- beta2*S2*P2 - (alpha2+mu2+h2)*I2 + MH12*I1 - MH21*I2
    dP2 <- chi2*I2 - delta2*P2 + MP12*P1 - MP21*P2
    list(c(dS1,dI1,dP1,dS2,dI2,dP2))
  })
}

initial_state <- c(S1=1e6, I1=100, P1=1e5,
                   S2=1e6, I2=10, P2=1e4)

make_parameters <- function(stock_ratio, host_movement) {
  c(Lambda1=Lambda_ref*stock_ratio, Lambda2=Lambda_ref*stock_ratio,
    mu1=mu_ref, mu2=mu_ref, h1=h_ref, h2=h_ref,
    beta1=beta_ref, beta2=beta_ref,
    alpha1=alpha_ref, alpha2=alpha_ref,
    chi1=chi_ref, chi2=chi_ref,
    delta1=delta_ref, delta2=delta_ref,
    MP12=MP12_ref, MP21=MP21_ref,
    MH12=host_movement, MH21=MH21_ref)
}

# Match the original manuscript: 120 months, output every 0.1 month,
# and the mean of the final 10% of saved observations.
movement_grid <- seq(0, 0.10, length.out = 51)
stock_ratios <- c(1, 2)
t_end <- 120
dt <- 0.1
time_points <- seq(0, t_end, by = dt)

tail_mean_original <- function(x, fraction = 0.1) {
  k <- max(1, ceiling(length(x) * fraction))
  mean(tail(x, k), na.rm = TRUE)
}

run_case <- function(stock_ratio, host_movement) {
  pars <- make_parameters(stock_ratio, host_movement)
  sol <- as.data.frame(ode(y = initial_state, times = time_points,
                           func = isa_two_patch, parms = pars,
                           method = 'lsoda'))
  sol$prev1 <- sol$I1 / (sol$S1 + sol$I1)
  sol$prev2 <- sol$I2 / (sol$S2 + sol$I2)
  sol$prev_total <- (sol$I1 + sol$I2) /
                    (sol$S1 + sol$I1 + sol$S2 + sol$I2)
  data.frame(stock_ratio = stock_ratio, MH12 = host_movement,
             S1 = tail_mean_original(sol$S1),
             I1 = tail_mean_original(sol$I1),
             P1 = tail_mean_original(sol$P1),
             prev1 = tail_mean_original(sol$prev1),
             prev2 = tail_mean_original(sol$prev2),
             prev_total = tail_mean_original(sol$prev_total),
             S1_final = tail(sol$S1, 1),
             I1_final = tail(sol$I1, 1),
             P1_final = tail(sol$P1, 1),
             prev1_final = tail(sol$prev1, 1))
}

cases <- expand.grid(stock_ratio = stock_ratios, MH12 = movement_grid)
results <- bind_rows(lapply(seq_len(nrow(cases)), function(i) {
  run_case(cases$stock_ratio[i], cases$MH12[i])
}))
write.csv(results, 'appendix_figure6_values.csv', row.names = FALSE)

cat('Selected movement rates and 120-month tail-mean prevalence:\n')
print(results %>% filter(abs(MH12 - 0.05) < 1e-8 |
                         abs(MH12 - 0.10) < 1e-8) %>%
      select(stock_ratio, MH12, S1, I1, P1, prev1, prev1_final),
      row.names = FALSE)

# Verify exact reference inputs against the user's original script.
stopifnot(abs(Lambda_ref - mean(c(493233, 1108106))/24) < 1e-10,
          MP12_ref == 0.05, MP21_ref == 0.01,
          all(stock_ratios == c(1,2)),
          length(movement_grid) == 51,
          t_end == 120, dt == 0.1)

plot_data <- results %>%
  mutate(Stocking=factor(stock_ratio, levels=c(1,2),
                         labels=c('Reference stocking (1x)',
                                  'Twice-reference stocking (2x)'))) %>%
  pivot_longer(cols=c(S1,I1,P1,prev1),
               names_to='variable',values_to='value') %>%
  mutate(variable=factor(variable,levels=c('S1','I1','P1','prev1'),
                         labels=c('Susceptible fish (S1)',
                                  'Infected fish (I1)',
                                  'Environmental pathogen (P1)',
                                  'Patch 1 prevalence')))

p <- ggplot(plot_data,aes(x=MH12,y=value)) +
  geom_line(linewidth=0.65,color='black') +
  facet_grid(variable~Stocking,scales='free_y') +
  scale_x_continuous(breaks=seq(0,0.1,0.025)) +
  labs(x=expression('Host movement rate '*M[12]^H*' (month'^-1*')'),
       y=NULL,
       title='Host-exporting patch: 120-month tail-mean responses') +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(),
        strip.background=element_rect(fill='grey94'),
        plot.title=element_text(face='bold',size=12))

ggsave('appendix_patch1_mechanism.png',p,width=10,height=10,dpi=350)
ggsave('appendix_patch1_mechanism.pdf',p,width=10,height=10)
cat('\nFiles created: appendix_patch1_mechanism.png, .pdf, and appendix_figure6_values.csv\n')
