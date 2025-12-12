# --- 0. Load Libraries ---
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(writexl)

# --- 1. Define Main Simulation Function ---
run_simulation <- function(M, k, W_inf, b, t_0 = -0.1, a = 0.05, N_1 = 1e6, max_age = 10) {
  ages <- 1:max_age
  
  # Population Structure & Growth
  N_i <- N_1 * exp(-M * (ages - 1))
  W_i <- W_inf * (1 - exp(-k * (ages - t_0)))^3
  
  # Consumption (Allometric)
  R_spec_i <- a * W_i^b
  
  # Aggregate Metrics
  Biomass_i <- N_i * W_i
  Consumption_daily_i <- Biomass_i * R_spec_i
  
  total_B <- sum(Biomass_i)
  total_Q_daily <- sum(Consumption_daily_i)
  total_Q_annual <- total_Q_daily * 365
  total_N <- sum(N_i)
  
  # Calculated Rations (% Body Weight)
  R_naive <- (total_Q_daily / total_B) * 100        # Biomass-weighted (The Error)
  R_avg <- (sum(N_i * R_spec_i) / total_N) * 100    # Numerical Average (The Truth)
  
  return(data.frame(
    Pop_QB_Annual = total_Q_annual / total_B,
    R_naive_pct = R_naive,
    R_avg_pct = R_avg,
    Discrepancy_Ratio = R_avg / R_naive
  ))
}

# --- 2. Reproduce Manuscript Table 1 ---
scenarios <- list(
  Base = list(M=0.4, k=0.3, W_inf=1000, b=-0.27),
  High_Mortality = list(M=0.8, k=0.3, W_inf=1000, b=-0.27),
  Fast_Growth = list(M=0.4, k=0.6, W_inf=1000, b=-0.27),
  Strong_Allometry = list(M=0.4, k=0.3, W_inf=1000, b=-0.40)
)

table1_results <- bind_rows(lapply(names(scenarios), function(name) {
  params <- scenarios[[name]]
  res <- do.call(run_simulation, params)
  cbind(Scenario = name, res)
}))

print("--- Manuscript Table 1 ---")
print(table1_results)

# --- 3. Sensitivity Analysis on 'b' ---
# We vary b from -0.5 (strong allometry) to 0.0 (no allometry)
b_values <- seq(-0.5, 0.0, by = 0.01)
sensitivity_results <- bind_rows(lapply(b_values, function(val) {
  # Run simulation using Base parameters but changing 'b'
  res <- run_simulation(M=0.4, k=0.3, W_inf=1000, b=val)
  res$b_value <- val
  return(res)
}))

# Plot the Sensitivity
p <- ggplot(sensitivity_results, aes(x = b_value)) +
  geom_line(aes(y = R_avg_pct, color = "True Average (R_avg)"), size = 1.2) +
  geom_line(aes(y = R_naive_pct, color = "Naive Ration (R_naive)"), size = 1.2, linetype = "dashed") +
  scale_color_manual(values = c("True Average (R_avg)" = "darkgreen", "Naive Ration (R_naive)" = "red")) +
  labs(
    title = "Sensitivity of Daily Ration Estimates to Allometric Exponent (b)",
    subtitle = "As |b| increases, the naive estimate diverges further from the true average",
    x = "Allometric Exponent (b)",
    y = "Daily Ration (% Body Weight / Day)",
    color = "Metric"
  ) +
  theme_bw()

print(p)
ggsave("Sensitivity_Analysis_b.png", p, width=7, height=5)

# --- 4. Process FishBase Data ---
# Load provided CSV
fishbase_data <- read_csv("FishBase_QB_table.csv", show_col_types = FALSE)

fishbase_processed <- fishbase_data %>%
  filter(!is.na(PopQB)) %>%
  mutate(
    Naive_Daily_Ration_pct = (PopQB / 365) * 100,
    Label = paste(Genus, Species)
  ) %>%
  select(Label, PopQB, Naive_Daily_Ration_pct, Winf, K, Mortality) %>%
  arrange(desc(PopQB))

# --- 5. Export All to Excel ---
sheets <- list(
  "Manuscript_Table_1" = table1_results,
  "Sensitivity_Data" = sensitivity_results,
  "FishBase_Processed" = fishbase_processed
)
write_xlsx(sheets, "FishBase_and_Simulation_Analysis.xlsx")