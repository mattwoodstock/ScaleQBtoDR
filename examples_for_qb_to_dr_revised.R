# R Script to Reproduce Q/B vs. Daily Ration Scaling Discrepancy
#
# This script accompanies the manuscript:
# "From population consumption (Q/B) to individual daily ration:
#  a workflow for allometric parameterization in fisheries models"
#
# Woodstock & Bigman
#
# It reproduces the four simulation scenarios, generates the key metrics
# presented in Table 1, and creates the manuscript figures using a
# colorblind-friendly palette (Okabe-Ito).
#
# Colorblind palette reference:
# Okabe, M. & Ito, K. (2008). Color Universal Design (CUD).
# https://jfly.uni-koeln.de/color/

# --- 0. Load Libraries ---
# install.packages(c("ggplot2", "tidyr", "dplyr", "scales"))

library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# --- Colorblind-Friendly Palette (Okabe-Ito) ---
# Safe for deuteranopia, protanopia, and tritanopia
cb_black      <- "#000000"
cb_orange     <- "#E69F00"
cb_sky_blue   <- "#56B4E9"
cb_green      <- "#009E73"
cb_yellow     <- "#F0E442"
cb_blue       <- "#0072B2"
cb_vermillion <- "#D55E00"
cb_pink       <- "#CC79A7"

# --- 1. Main Simulation Function ---

#' Calculates population-level scaling metrics from individual parameters.
#'
#' Simulates a stable, age-structured population and calculates the true
#' population Q/B, the "Naive Daily Ration" (R_naive), and related metrics.
#'
#' @param M Annual instantaneous mortality rate.
#' @param k von Bertalanffy growth coefficient (year^-1).
#' @param W_inf Asymptotic weight (g).
#' @param t_0 Theoretical age at zero weight (years).
#' @param a Allometric consumption coefficient (ration for a 1-g fish; g/g/day).
#' @param b Allometric consumption exponent (dimensionless; typically negative).
#' @param N_1 Recruitment (abundance at Age-1).
#' @param max_age Maximum age class in the simulation.
#' @param return_pop_details If TRUE, returns full age-structured data frame.
#'
#' @return A summary data.frame row OR full age-structured data.frame.
run_population_simulation <- function(M, k, W_inf, t_0, a, b,
                                      N_1 = 1e6, max_age = 10,
                                      return_pop_details = FALSE) {

  ages <- 1:max_age
  pop  <- data.frame(age = ages)

  # 1. Population structure (Eq. 7)
  pop$N_i <- N_1 * exp(-M * (pop$age - 1))

  # 2. Growth (Eq. 8)
  pop$W_i <- W_inf * (1 - exp(-k * (pop$age - t_0)))^3

  # 3. Individual specific daily ration (Eq. 1)
  pop$R_spec_i <- a * pop$W_i^b

  # 4. Aggregate calculations
  pop$Biomass_i          <- pop$N_i * pop$W_i
  total_B                <- sum(pop$Biomass_i)
  pop$Consumption_daily_i <- pop$Biomass_i * pop$R_spec_i
  total_Q_daily          <- sum(pop$Consumption_daily_i)
  total_Q_annual         <- total_Q_daily * 365

  # 5. Key scaling metrics
  pop_QB  <- total_Q_annual / total_B          # Population Q/B (yr^-1) (Eq. 4)
  R_naive <- total_Q_daily  / total_B          # Naive daily ration (Eq. 5)
  total_N <- sum(pop$N_i)
  R_avg   <- sum(pop$N_i * pop$R_spec_i) / total_N  # True numerical avg (Eq. 6)

  if (return_pop_details) {
    pop <- pop %>%
      mutate(prop_N = N_i / total_N,
             prop_B = Biomass_i / total_B)
    return(pop)
  } else {
    R_spec_1 <- pop$R_spec_i[pop$age == 1]
    R_spec_5 <- pop$R_spec_i[pop$age == 5]
    W_1      <- pop$W_i[pop$age == 1]
    W_5      <- pop$W_i[pop$age == 5]

    results <- data.frame(
      Scenario    = NA,
      Key_Param   = paste0("M=", M, ", k=", k, ", b=", b),
      Pop_QB_yr   = pop_QB,
      R_naive_pct = R_naive * 100,
      R_spec_1_pct = R_spec_1 * 100,
      R_spec_5_pct = R_spec_5 * 100,
      R_avg_pct   = R_avg * 100,
      W_1_g       = W_1,
      W_5_g       = W_5
    )
    return(results)
  }
}

# --- 2. Define Simulation Scenarios ---

base_params <- list(M = 0.4, k = 0.3, W_inf = 1000, t_0 = -0.1,
                    a = 0.05, b = -0.27)

scenarios <- list(
  Base           = base_params,
  High_Mortality = modifyList(base_params, list(M   = 0.8)),
  Fast_Growth    = modifyList(base_params, list(k   = 0.6)),
  Strong_Allometry = modifyList(base_params, list(b = -0.40))
)

# Scenario display labels
scenario_labels <- c(
  Base             = "Base (M=0.4, k=0.3, b=-0.27)",
  High_Mortality   = "High Mortality (M=0.8)",
  Fast_Growth      = "Fast Growth (k=0.6)",
  Strong_Allometry = "Strong Allometry (b=-0.40)"
)

# --- 3. Run Simulations (Table 1) ---

results_list <- lapply(names(scenarios), function(name) {
  params    <- scenarios[[name]]
  res       <- do.call(run_population_simulation, params)
  res$Scenario <- name
  res       <- res[, c("Scenario", "Key_Param", "Pop_QB_yr",
                       "R_naive_pct", "R_spec_1_pct",
                       "R_spec_5_pct", "R_avg_pct", "W_1_g", "W_5_g")]
  return(res)
})

final_results              <- do.call(rbind, results_list)
rownames(final_results)    <- NULL
final_results$Scenario     <- factor(final_results$Scenario,
                                     levels = names(scenarios))
final_results$Scenario_label <- scenario_labels[as.character(final_results$Scenario)]
final_results$Scenario_label <- factor(final_results$Scenario_label,
                                       levels = scenario_labels)

# --- 4. Print Table 1 ---

cat("--- Simulation Results (Table 1) ---\n\n")
table_print <- final_results[, 1:7]
colnames(table_print) <- c("Scenario", "Key Parameter", "Pop. Q/B (yr-1)",
                           "R_naive (% BW/d)", "R_spec,1 (% BW/d)",
                           "R_spec,5 (% BW/d)", "R_avg (% BW/d)")
table_print[, 3:7] <- round(table_print[, 3:7], 2)
print(table_print, row.names = FALSE)

# --- 5. Generate Plot Data ---

# Figure 1: Theoretical allometric scaling curves
plot1_data <- data.frame(Weight = seq(10, 1000, by = 5)) %>%
  mutate(
    `Base (b = -0.27)`           = 100 * base_params$a * Weight^base_params$b,
    `Strong Allometry (b = -0.40)` = 100 * base_params$a * Weight^(-0.40)
  ) %>%
  pivot_longer(cols = -Weight, names_to = "Scenario", values_to = "Ration_pct") %>%
  mutate(Scenario = factor(Scenario,
                           levels = c("Base (b = -0.27)",
                                      "Strong Allometry (b = -0.40)")))

# Figure 2: Population structure
pop_base   <- do.call(run_population_simulation,
                      c(scenarios$Base,           return_pop_details = TRUE))
pop_high_m <- do.call(run_population_simulation,
                      c(scenarios$High_Mortality, return_pop_details = TRUE))

plot2_data <- bind_rows(
  mutate(pop_base,   Scenario = "Base (M = 0.4)"),
  mutate(pop_high_m, Scenario = "High Mortality (M = 0.8)")
) %>%
  pivot_longer(cols = c(prop_N, prop_B),
               names_to  = "Metric",
               values_to = "Proportion") %>%
  mutate(Metric = ifelse(Metric == "prop_N",
                         "Abundance (N)",
                         "Biomass (B)"),
         Metric = factor(Metric, levels = c("Abundance (N)", "Biomass (B)")))

# --- 6. Build Figures ---

# Shared theme
ms_theme <- theme_bw(base_size = 13) +
  theme(
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold"),
    strip.background   = element_rect(fill = "grey92", colour = "grey60"),
    strip.text         = element_text(face = "bold")
  )

# ----- Figure 1: Allometric Scaling -----
fig1 <- ggplot(plot1_data, aes(x = Weight, y = Ration_pct,
                               linetype = Scenario, colour = Scenario)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c(cb_blue, cb_vermillion)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(
    x        = "Fish Weight (g)",
    y        = "Specific Daily Ration (% BW day\u207B\u00B9)",
    colour   = "Allometric Scenario",
    linetype = "Allometric Scenario"
  ) +
  ms_theme

# ----- Figure 2: Population Structure -----
fig2 <- ggplot(plot2_data,
               aes(x = age, y = Proportion, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge",
           colour = "white", linewidth = 0.3) +
  facet_wrap(~Scenario) +
  scale_fill_manual(values = c("Abundance (N)" = cb_sky_blue,
                               "Biomass (B)"   = cb_blue)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x    = "Age Class (years)",
    y    = "Proportion of Total",
    fill = "Population Metric"
  ) +
  ms_theme

# ----- Figure 3: Scaling Discrepancy (Dumbbell) -----
fig3 <- ggplot(final_results, aes(y = Scenario_label)) +

  # Range bar: true allometric range (adult to juvenile)
  geom_segment(aes(x = R_spec_5_pct, xend = R_spec_1_pct,
                   yend = Scenario_label),
               colour   = "grey55",
               linewidth = 2) +

  # Juvenile (Age-1) ration
  geom_point(aes(x = R_spec_1_pct, colour = "Age 1",
                 shape = "Age 1"), size = 4) +

  # Adult (Age-5) ration
  geom_point(aes(x = R_spec_5_pct, colour = "Age 5",
                 shape = "Age 5"), size = 4) +

  # True numerical average
  geom_point(aes(x = R_avg_pct, colour = "True Avg.",
                 shape = "True Avg."), size = 4.5) +

  # Naive ration
  geom_point(aes(x = R_naive_pct, colour = "Naive Ration",
                 shape = "Naive Ration"), size = 4.5, stroke = 1.5) +

  scale_colour_manual(
    name   = "Metric",
    values = c(
      "Age 1"       = cb_sky_blue,
      "Age 5"          = cb_blue,
      "True Avg."      = cb_green,
      "Naive Ration" = cb_vermillion
    )
  ) +
  scale_shape_manual(
    name   = "Metric",
    values = c(
      "Age 1"       = 16,   # filled circle
      "Age 5"          = 16,   # filled circle (darker colour)
      "True Avg."      = 18,   # filled diamond
      "Naive Ration" = 4     # X
    )
  ) +
  labs(
    x = "Specific Daily Ration (% Body Weight day\u207B\u00B9)",
    y = NULL
  ) +
  ms_theme +
  theme(legend.key.width = unit(1.5, "lines"))

# --- 7. Save Figures ---

ggsave("figure1.png", fig1, width = 6, height = 5, dpi = 300)
ggsave("figure2.png", fig2, width = 7, height = 5, dpi = 300)
ggsave("figure3.png", fig3, width = 7, height = 5, dpi = 300)

message("Figures saved: figure1.png, figure2.png, figure3.png")
