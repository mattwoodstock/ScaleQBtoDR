# --- 0. Setup and Libraries ---
# install.packages(c("dplyr", "writexl", "readr"))
library(dplyr)
library(writexl)
library(readr)

# --- PART 1: SIMULATION TABLE (From Literature/Manuscript) ---
# This section reproduces Table 1 using the logic provided in your "examples for qb to dr.R"

run_simulation_row <- function(M, k, W_inf, b) {
  # Fixed parameters from the manuscript
  t_0 <- -0.1
  a <- 0.05
  N_1 <- 1e6
  max_age <- 10
  
  ages <- 1:max_age
  
  # 1. Population Structure
  N_i <- N_1 * exp(-M * (ages - 1))
  
  # 2. Growth
  W_i <- W_inf * (1 - exp(-k * (ages - t_0)))^3
  
  # 3. Individual Consumption (True Allometry)
  R_spec_i <- a * W_i^b
  
  # 4. Aggregate Calculations
  Biomass_i <- N_i * W_i
  Consumption_daily_i <- Biomass_i * R_spec_i
  
  total_B <- sum(Biomass_i)
  total_Q_daily <- sum(Consumption_daily_i)
  total_Q_annual <- total_Q_daily * 365
  total_N <- sum(N_i)
  
  # 5. Metrics
  pop_QB <- total_Q_annual / total_B
  R_naive <- (total_Q_daily / total_B) * 100 # %BW
  R_avg <- (sum(N_i * R_spec_i) / total_N) * 100 # %BW
  R_spec_1 <- (a * W_i[1]^b) * 100
  R_spec_5 <- (a * W_i[5]^b) * 100
  
  return(data.frame(
    M = M, k = k, b = b,
    Pop_QB_Annual = pop_QB,
    Naive_Ration_pct = R_naive,
    True_Avg_Ration_pct = R_avg,
    Juv_Ration_Age1_pct = R_spec_1,
    Adult_Ration_Age5_pct = R_spec_5
  ))
}

# Define Scenarios
scenarios <- list(
  Base = c(M=0.4, k=0.3, W_inf=1000, b=-0.27),
  High_Mortality = c(M=0.8, k=0.3, W_inf=1000, b=-0.27),
  Fast_Growth = c(M=0.4, k=0.6, W_inf=1000, b=-0.27),
  Strong_Allometry = c(M=0.4, k=0.3, W_inf=1000, b=-0.40)
)

# Run simulations
sim_rows <- lapply(names(scenarios), function(name) {
  p <- scenarios[[name]]
  res <- run_simulation_row(p['M'], p['k'], p['W_inf'], p['b'])
  cbind(Scenario = name, res)
})
df_sim_table <- do.call(rbind, sim_rows)


# --- PART 2: FISHBASE TABLE DEVELOPMENT ---
# This processes your provided FishBase CSV to generate daily ration estimates

# 1. Load the FishBase Data
# Ensure "FishBase_QB_table (2).csv" is in your working directory
# If not, set working directory: setwd("path/to/your/files")
fishbase_data <- read_csv("FishBase_QB_table.csv", show_col_types = FALSE)

# 2. Process and Calculate
# We calculate "Naive Daily Ration" because we lack the allometric exponent 'b' 
# for every species to calculate the "True" average.
# Naive Ration = (Population Q/B per year / 365) * 100
df_fishbase_processed <- fishbase_data %>%
  filter(!is.na(PopQB), !is.na(Winf), !is.na(Mortality)) %>% # Clean empty rows
  mutate(
    # Basic identifier
    Scientific_Name = paste(Genus, Species),
    
    # The calculation: Convert Annual Q/B to Daily % Body Weight
    Naive_Daily_Ration_pct = (PopQB / 365) * 100
  ) %>%
  # Select and reorder columns for the final clean table
  select(
    Scientific_Name,
    Common_Name = sciname, # Assuming sciname might contain common text or just verify column
    W_inf_g = Winf,
    K_coeff = K,
    Mortality_M = Mortality,
    Pop_QB_Annual = PopQB,
    Naive_Daily_Ration_pct
  ) %>%
  arrange(desc(Pop_QB_Annual)) # Sort by highest consumption

# --- PART 3: EXPORT TO EXCEL ---

# Create a list of sheets for the Excel file
sheets <- list(
  "Manuscript_Simulation_Table" = df_sim_table,
  "FishBase_Empirical_Data" = df_fishbase_processed
)

# Write to file
write_xlsx(sheets, path = "FishBase_and_Simulation_Rations.xlsx")

cat("Success! 'FishBase_and_Simulation_Rations.xlsx' has been created.\n")
cat("It contains two tabs:\n")
cat("1. The simulated scenarios from the manuscript (Table 1).\n")
cat("2. The empirical FishBase data with calculated Daily Rations.\n")