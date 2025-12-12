# R Script to Reproduce Q/B vs. Daily Ration Scaling Discrepancy
#
# This script accompanies the manuscript "A methodological clarification: 
# The scaling discrepancy between population consumption (Q/B) and individual daily ration"
# It reproduces the four simulation scenarios, generates the key metrics
# presented in Table 1, and creates the manuscript's key figures.

# --- 0. Load Libraries ---
# Install these packages if you don't have them:
# install.packages(c("ggplot2", "tidyr", "dplyr"))

library(ggplot2)
library(tidyr)
library(dplyr)


# --- 1. Main Simulation Function ---

#' Calculates population-level scaling metrics from individual parameters.
#'
#' This function simulates a stable, age-structured population and calculates
#' the true population Q/B, the "Naive Daily Ration" (R_naive), and other
#' metrics as defined in the manuscript.
#'
#' @param M Annual instantaneous mortality rate.
#' @param k von Bertalanffy growth coefficient (k).
#' @param W_inf Asymptotic weight (g).
#' @param t_0 Theoretical age at zero weight (years).
#' @param a Allometric consumption coefficient (ration for 1-g fish).
#' @param b Allometric consumption exponent.
#' @param N_1 Recruitment (abundance at Age-1).
#' @param max_age Maximum age class in the simulation.
#' @param return_pop_details (Logical) If TRUE, returns the full age-structured
#'   data frame instead of the summary row.
#'
#' @return A data.frame row with all key parameters and calculated metrics
#'   OR a full data.frame with population details by age.
run_population_simulation <- function(M, k, W_inf, t_0, a, b, 
                                      N_1 = 1e6, max_age = 10,
                                      return_pop_details = FALSE) {
  
  # Create a data frame for the population, structured by age
  ages <- 1:max_age
  pop <- data.frame(age = ages)
  
  # --- 1. POPULATION STRUCTURE (Eq 7) ---
  # Abundance at age (N_i)
  pop$N_i <- N_1 * exp(-M * (pop$age - 1))
  
  # --- 2. GROWTH (Eq 8) ---
  # Weight at age (W_i)
  pop$W_i <- W_inf * (1 - exp(-k * (pop$age - t_0)))^3
  
  # --- 3. INDIVIDUAL CONSUMPTION (Eq 1) ---
  # Specific daily ration (g food / g fish / day) for each age class
  pop$R_spec_i <- a * pop$W_i^b
  
  # --- 4. AGGREGATE CALCULATIONS ---
  
  # Total population biomass (B) (Eq 2)
  # Biomass for each age class
  pop$Biomass_i <- pop$N_i * pop$W_i
  # Sum for total B
  total_B <- sum(pop$Biomass_i)
  
  # Total annual population consumption (Q) (Eq 3)
  # Daily consumption (g/day) for each age class
  pop$Consumption_daily_i <- pop$Biomass_i * pop$R_spec_i
  # Sum for total daily Q
  total_Q_daily <- sum(pop$Consumption_daily_i)
  # Annualize
  total_Q_annual <- total_Q_daily * 365
  
  # --- 5. CALCULATE KEY SCALING METRICS ---
  
  # Population Q/B (annual) (Eq 4)
  pop_QB <- total_Q_annual / total_B
  
  # "Naive Daily Ration" (R_naive) (Eq 5)
  # This is the biomass-weighted average ration
  # Note: (pop_QB / 365) is identical to (total_Q_daily / total_B)
  R_naive <- total_Q_daily / total_B
  
  # "True Average Daily Ration" (R_avg) (Eq 6)
  # This is the numerical-average ration across all individuals
  total_N <- sum(pop$N_i)
  R_avg <- sum(pop$N_i * pop$R_spec_i) / total_N
  
  # --- 6. COMPILE RESULTS ---
  
  if (return_pop_details) {
    # Add proportional metrics for plotting
    pop <- pop %>%
      mutate(
        prop_N = N_i / total_N,
        prop_B = Biomass_i / total_B
      )
    return(pop)
    
  } else {
    # Get specific rations for Age-1 and Age-5 fish for comparison
    R_spec_1 <- pop$R_spec_i[pop$age == 1]
    R_spec_5 <- pop$R_spec_i[pop$age == 5]
    
    # Get weights for context
    W_1 <- pop$W_i[pop$age == 1]
    W_5 <- pop$W_i[pop$age == 5]
    
    # Create a summary results data.frame
    results <- data.frame(
      Scenario = NA, # Placeholder, to be filled in later
      Key_Param = paste0("M=", M, ", k=", k, ", b=", b),
      Pop_QB_yr = pop_QB,
      R_naive_pct = R_naive * 100,
      R_spec_1_pct = R_spec_1 * 100,
      R_spec_5_pct = R_spec_5 * 100,
      R_avg_pct = R_avg * 100,
      W_1_g = W_1,
      W_5_g = W_5
    )
    return(results)
  }
}

# --- 2. Define Simulation Scenarios ---

# Base parameters
base_params <- list(
  M = 0.4,
  k = 0.3,
  W_inf = 1000,
  t_0 = -0.1,
  a = 0.05,
  b = -0.27
)

# Define the list of four scenarios from the manuscript
scenarios <- list(
  Base = base_params,
  
  High_Mortality = modifyList(base_params, list(M = 0.8)),
  
  Fast_Growth = modifyList(base_params, list(k = 0.6)),
  
  Strong_Allometry = modifyList(base_params, list(b = -0.40))
)

# --- 3. Run Simulations (for Table 1) ---

# Use lapply to run the simulation for each scenario in the list
results_list <- lapply(names(scenarios), function(name) {
  params <- scenarios[[name]]
  
  # Call the main function using do.call to pass the list of parameters
  res <- do.call(run_population_simulation, params)
  
  # Add the scenario name
  res$Scenario <- name
  
  # Re-order columns for clarity
  res <- res[, c("Scenario", "Key_Param", "Pop_QB_yr", "R_naive_pct", "R_spec_1_pct", "R_spec_5_pct", "R_avg_pct", "W_1_g", "W_5_g")]
  return(res)
})

# Combine the list of data.frame rows into a single data.frame
final_results <- do.call(rbind, results_list)
rownames(final_results) <- NULL # Clean up row names

# --- 4. Print Final Results (Table 1) ---

cat("--- Simulation Results (Reproduces Table 1) ---\n\n")

# Format for printing
table_1_data <- final_results[, 1:7]
colnames(table_1_data) <- c("Scenario", "Key Parameter", "Pop. Q/B (year-1)", "R_naive (% BW)", "R_spec,1 (Juv) (% BW)", "R_spec,5 (Ad) (% BW)", "R_avg (True Avg) (% BW)")

# Create a copy for formatting, round only numeric columns
table_1_formatted <- table_1_data
numeric_cols_table_1 <- 3:7 # Columns 3 through 7 are numeric
table_1_formatted[, numeric_cols_table_1] <- round(table_1_formatted[, numeric_cols_table_1], 2)

# Print the formatted table, without row numbers
print(table_1_formatted, row.names = FALSE)

cat("\n--- Weight-at-Age Context ---\n")

# Create a copy for formatting, round only numeric columns
context_data <- final_results[, c("Scenario", "W_1_g", "W_5_g")]
numeric_cols_context <- 2:3 # Columns 2 and 3 are numeric
context_data[, numeric_cols_context] <- round(context_data[, numeric_cols_context], 1)

# Print the formatted context data, without row numbers
print(context_data, row.names = FALSE)

cat("\nNote: R_naive (biomass-weighted avg) is consistently lower than R_avg (numerical avg) \n")
cat("and provides a poor estimate of the true allometric rations (R_spec,1, R_spec,5).\n")


# --- 5. Generate Data for Plots ---

# Plot 1 Data (Theoretical Allometry)
# Create a sequence of weights and calculate rations for two scenarios
plot1_data <- data.frame(Weight = seq(25, 1000, by = 10)) %>%
  mutate(
    Base = 100 * base_params$a * Weight^base_params$b,
    Strong_Allometry = 100 * base_params$a * Weight^scenarios$Strong_Allometry$b
  ) %>%
  pivot_longer(cols = c(Base, Strong_Allometry), names_to = "Scenario", values_to = "Ration_pct")

# Plot 2 Data (Population Structure)
# Rerun Base and High_Mortality scenarios to get detailed population data
pop_base <- do.call(run_population_simulation, c(scenarios$Base, return_pop_details = TRUE))
pop_high_m <- do.call(run_population_simulation, c(scenarios$High_Mortality, return_pop_details = TRUE))

plot2_data <- bind_rows(
  mutate(pop_base, Scenario = "Base (M=0.4)"),
  mutate(pop_high_m, Scenario = "High Mortality (M=0.8)")
) %>%
  pivot_longer(cols = c(prop_N, prop_B), names_to = "Metric", values_to = "Proportion") %>%
  mutate(Metric = ifelse(Metric == "prop_N", "Abundance (N)", "Biomass (B)"))


# Plot 3 Data (Main Results)
# The `final_results` data frame is already in the correct format to be plotted.
# We will just re-order the scenarios for a logical plot
final_results$Scenario <- factor(final_results$Scenario, 
                                 levels = c("Base", "High_Mortality", "Fast_Growth", "Strong_Allometry"))

# --- 6. Generate and Print Plots ---

# Plot 1: Theoretical Allometric Scaling
plot1 <- ggplot(plot1_data, aes(x = Weight, y = Ration_pct, linetype = Scenario)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Fish Weight (g)",
    y = "Specific Daily Ration (% BW/day)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

print(plot1)


# Plot 2: Population Structure (Abundance vs. Biomass)
plot2 <- ggplot(plot2_data, aes(x = age, y = Proportion, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Scenario) +
  labs(
    x = "Age Class (years)",
    y = "Proportion of Total"
  ) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw(base_size = 14)

print(plot2)


# Plot 3: The Scaling Discrepancy (Dumbbell Plot)
plot3 <- ggplot(final_results, aes(y = Scenario)) +
  
  # Dumbbell segment showing the "true" allometric range
  geom_segment(aes(x = R_spec_5_pct, xend = R_spec_1_pct, yend = Scenario),
               color = "grey", linewidth = 1.5) +
  
  # Point for true Juvenile ration
  geom_point(aes(x = R_spec_1_pct, color = "Juvenile (Age-1)"), size = 4) +
  
  # Point for true Adult ration
  geom_point(aes(x = R_spec_5_pct, color = "Adult (Age-5)"), size = 4) +
  
  # Point for the true numerical average
  geom_point(aes(x = R_avg_pct, color = "True Avg (R_avg)"), shape = 18, size = 5) +
  
  # The "X" of ERROR: The Naive Ration
  geom_point(aes(x = R_naive_pct, color = "Naive Ration (R_naive)"), 
             shape = 4, size = 5, stroke = 1.5) +
  
  scale_color_manual(
    name = "Metric",
    values = c(
      "Juvenile (Age-1)" = "skyblue",
      "Adult (Age-5)" = "navyblue",
      "True Avg (R_avg)" = "darkgreen",
      "Naive Ration (R_naive)" = "red"
    )
  ) +
  labs(
    x = "Specific Daily Ration\n(% Body Weight / day)",
    y = "Simulation Scenario"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",legend.title = element_blank())+
  guides(color = guide_legend(ncol = 2))

print(plot3)

ggsave("figure1.png",plot1,height = unit(5,"in"),width=unit(6,"in"),dpi=300)
ggsave("figure2.png",plot2,height = unit(5,"in"),width=unit(6,"in"),dpi=300)
ggsave("figure3.png",plot3,height = unit(5,"in"),width=unit(6,"in"),dpi=300)

