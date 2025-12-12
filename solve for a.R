#' Calculate Population Q/B from Individual Parameters
#'
#' This function simulates a stable age-structured population based on
#' life-history parameters and calculates the emergent, population-level
#' consumption-to-biomass ratio (Q/B).
#'
#' @param a The allometric consumption coefficient (ration of a 1g fish, g/g/day).
#' @param b The allometric consumption exponent (dimensionless).
#' @param W_inf Asymptotic weight (g) from the von Bertalanffy Growth Function (VBGF).
#' @param k VBGF growth coefficient (year^-1).
#' @param t0 Theoretical age at zero length/weight (years).
#' @param M Instantaneous natural mortality rate (year^-1).
#' @param max_age The maximum age to simulate (e.g., 10 or 20 years).
#'
#' @return The calculated annual population Q/B (year^-1).
calculate_Q_B <- function(a, b, W_inf, k, t0, M, max_age = 20) {
  
  # Create a vector of ages from 1 to max_age
  ages <- 1:max_age
  
  # 1. Calculate population structure (Abundance-at-age)
  N_at_age <- 1 * exp(-M * (ages - 1))
  
  # 2. Calculate growth (Weight-at-age)
  W_at_age <- W_inf * (1 - exp(-k * (ages - t0)))^3
  
  # 3. Calculate individual-level specific daily ration (R_spec,i)
  R_spec_at_age <- a * (W_at_age ^ b)
  
  # 4. Calculate total population biomass (B)
  total_B <- sum(N_at_age * W_at_age, na.rm = TRUE)
  
  # 5. Calculate total annual population consumption (Q)
  total_Q <- sum(N_at_age * W_at_age * R_spec_at_age * 365, na.rm = TRUE)
  
  # 6. Calculate final population Q/B
  pop_Q_B <- total_Q / total_B
  
  return(pop_Q_B)
}


#' Solve for Allometric Coefficient 'a'
#'
#' This function iteratively finds the allometric coefficient 'a' that
#' causes the `calculate_Q_B` function to produce a specific `target_Q_B`.
#' It uses R's `uniroot` function to find the root (zero) of the error.
#'
#' @param target_Q_B The desired population Q/B (year^-1) you want to match.
#' @param b The *assumed* allometric consumption exponent (e.g., -0.27).
#' @param W_inf Asymptotic weight (g).
#' @param k VBGF growth coefficient (year^-1).
#' @param t0 Theoretical age at zero length/weight (years).
#' @param M Instantaneous natural mortality rate (year^-1).
#' @param max_age Maximum age for the simulation.
#'
#' @return The solved value for the allometric coefficient 'a'.
solve_for_a <- function(target_Q_B, b, W_inf, k, t0, M, max_age = 20) {
  
  # --- Input Validation ---
  required_params <- c(target_Q_B, b, W_inf, k, t0, M)
  if (any(is.na(required_params)) || any(is.infinite(required_params))) {
    warning(paste("One or more required parameters are NA or Inf. Returning NA."))
    return(NA)
  }
  
  # Define the error function for the root finder.
  error_function <- function(a_val) {
    
    calculated_Q_B <- calculate_Q_B(a = a_val, 
                                    b = b, 
                                    W_inf = W_inf, 
                                    k = k, 
                                    t0 = t0, 
                                    M = M, 
                                    max_age = max_age)
    
    return(calculated_Q_B - target_Q_B)
  }
  
  # A more robust search interval check:
  low_a_check <- try(error_function(0.0001), silent = TRUE)
  high_a_check <- try(error_function(2.0), silent = TRUE)
  
  if (inherits(low_a_check, "try-error") || inherits(high_a_check, "try-error") ||
      is.na(low_a_check) || is.na(high_a_check) || 
      sign(low_a_check) == sign(high_a_check)) {
    
    warning("Root-finding interval may not contain the solution. 
            Target Q/B might be outside the physiologically possible range 
            for the given life-history parameters. Returning NA.")
    return(NA)
  }
  
  solution <- try(uniroot(
    f = error_function,
    interval = c(0.0001, 2.0), # Search interval for 'a'
    tol = 1e-6 # Set a reasonable tolerance for precision
  ), silent = TRUE)
  
  if (inherits(solution, "try-error")) {
    warning("uniroot failed to find a solution. Returning NA.")
    return(NA)
  }
  
  return(solution$root)
}

# --- EXAMPLE USAGE ---
#
# This section demonstrates how to use the `solve_for_a` function
# with a manually defined database of species traits.
# All rfishbase dependencies have been removed.

# 1. Load Required Libraries
library(dplyr)

# 2. Define Full Species Database Manually
species_db <- data.frame(
  species = c(
    # Small Pelagics (High M, High Q/B)
    "Engraulis encrasicolus", "Sardina pilchardus", "Clupea harengus", 
    "Sprattus sprattus", "Brevoortia tyrannus", "Sardinops sagax",
    # Demersal - Medium
    "Merluccius merluccius", "Melanogrammus aeglefinus", "Pollachius virens",
    "Urophycis chuss", "Stenotomus chrysops", "Centropristis striata",
    # Demersal - Large
    "Gadus morhua", "Anoplopoma fimbria", "Sebastes fasciatus",
    # Flatfish
    "Pleuronectes platessa", "Solea solea", "Hippoglossus hippoglossus",
    # Reef Fish
    "Lutjanus campechanus", "Epinephelus morio", "Mycteroperca microlepis",
    # Large Pelagics
    "Scomber scombrus", "Sarda sarda", "Coryphaena hippurus",
    "Thunnus albacares", "Thunnus thynnus", "Katsuwonus pelamis",
    # Elasmobranchs (low M, slow growth)
    "Raja clavata", "Squalus acanthias", "Galeorhinus galeus"
  ),
  target_Q_B = c(
    12.0, 10.0, 9.0, 11.0, 9.5, 10.5, # Small Pelagics
    4.0, 4.2, 4.5, 5.0, 6.0, 5.5,    # Demersal - Medium
    3.5, 2.8, 2.5,                   # Demersal - Large
    4.5, 4.0, 2.0,                   # Flatfish
    3.8, 3.0, 2.8,                   # Reef Fish
    7.0, 6.5, 8.0, 5.0, 4.0, 6.0,    # Large Pelagics
    2.0, 1.8, 1.9                    # Elasmobranchs
  ),
  assumed_b = rep(-0.27, 30),
  L_inf = c(
    20, 25, 38, 16, 35, 30,           # Small Pelagics
    80, 70, 100, 65, 40, 45,          # Demersal - Medium
    130, 100, 50,                    # Demersal - Large
    70, 60, 250,                     # Flatfish
    90, 110, 120,                    # Reef Fish
    60, 80, 150, 200, 300, 100,      # Large Pelagics
    100, 130, 180                    # Elasmobranchs
  ),
  k = c(
    0.8, 0.6, 0.5, 0.9, 0.7, 0.6,    # Small Pelagics
    0.25, 0.3, 0.28, 0.35, 0.4, 0.38, # Demersal - Medium
    0.2, 0.15, 0.1,                  # Demersal - Large
    0.2, 0.25, 0.12,                 # Flatfish
    0.22, 0.18, 0.15,                # Reef Fish
    0.6, 0.5, 0.9, 0.4, 0.18, 0.7,   # Large Pelagics
    0.1, 0.08, 0.07                  # Elasmobranchs
  ),
  t0 = c(
    -0.1, -0.2, -0.3, -0.1, -0.2, -0.2, # Small Pelagics
    -0.5, -0.4, -0.5, -0.3, -0.2, -0.2, # Demersal - Medium
    -0.4, -1.0, -1.5,                # Demersal - Large
    -0.6, -0.5, -0.8,                # Flatfish
    -0.5, -0.6, -0.7,                # Reef Fish
    -0.2, -0.3, -0.1, -0.2, -0.3, -0.2, # Large Pelagics
    -2.0, -2.5, -2.2                 # Elasmobranchs
  ),
  lw_a = c(
    0.005, 0.006, 0.005, 0.004, 0.007, 0.006, # Small Pelagics
    0.008, 0.007, 0.006, 0.009, 0.010, 0.012, # Demersal - Medium
    0.008, 0.005, 0.010,             # Demersal - Large
    0.009, 0.008, 0.004,             # Flatfish
    0.015, 0.012, 0.010,             # Reef Fish
    0.010, 0.011, 0.009, 0.018, 0.015, 0.017, # Large Pelagics
    0.006, 0.005, 0.004              # Elasmobranchs
  ),
  lw_b = c(
    3.0, 3.0, 3.1, 3.0, 3.0, 3.0,    # Small Pelagics
    3.0, 3.05, 3.0, 3.0, 3.0, 2.95,  # Demersal - Medium
    3.05, 3.1, 3.0,                  # Demersal - Large
    3.0, 3.0, 3.1,                   # Flatfish
    2.95, 3.0, 3.05,                 # Reef Fish
    3.0, 3.0, 3.0, 2.9, 2.9, 2.95,   # Large Pelagics
    3.1, 3.1, 3.1                    # Elasmobranchs
  ),
  M = c(
    1.2, 1.0, 0.9, 1.3, 1.0, 1.0,    # Small Pelagics
    0.4, 0.45, 0.4, 0.5, 0.6, 0.55,  # Demersal - Medium
    0.4, 0.25, 0.2,                  # Demersal - Large
    0.35, 0.3, 0.15,                 # Flatfish
    0.3, 0.25, 0.2,                  # Reef Fish
    0.8, 0.7, 0.9, 0.6, 0.3, 0.8,    # Large Pelagics
    0.15, 0.1, 0.1                   # Elasmobranchs
  ),
  stringsAsFactors = FALSE
)


# 3. Calculate W_inf from L_inf
#    W_inf = a * L_inf^b
species_db <- species_db %>%
  mutate(
    W_inf = lw_a * (L_inf ^ lw_b)
  ) %>%
  # Select only the columns needed for the solver
  select(species, target_Q_B, assumed_b, W_inf, k, t0, M)

# 4. Placeholder: Load Your Own Database
#
# species_db_path <- "path/to/your/fish_traits.csv"
# species_db <- read.csv(species_db_path, stringsAsFactors = FALSE)
#
# # Then, run steps 3 (Calculate W_inf) and 5 (Solve for 'a')
#

# 5. Iteratively Solve for 'a' for Each Species
#    We use `mapply` (a multivariate version of `sapply`) to "loop"
#    through each row of the database and pass its parameters to `solve_for_a`.

message("Running reconciliation for each species...")

solved_a_values <- mapply(
  solve_for_a,
  target_Q_B = species_db$target_Q_B,
  b = species_db$assumed_b,
  W_inf = species_db$W_inf,
  k = species_db$k,
  t0 = species_db$t0,
  M = species_db$M
)

# 6. Combine Results
#    Add the newly solved 'a' coefficients as a new column
#    in the original database.
species_db$solved_a <- solved_a_values

# 7. VERIFICATION: Calculate Q/B using the new 'solved_a'
#    This column should be (nearly) identical to 'target_Q_B'
message("Verifying results by recalculating Q/B...")

species_db$calculated_Q_B <- mapply(
  calculate_Q_B,
  a = species_db$solved_a,
  b = species_db$assumed_b,
  W_inf = species_db$W_inf,
  k = species_db$k,
  t0 = species_db$t0,
  M = species_db$M
)

# 8. View Final, Reconciled Database
message("Reconciliation complete. Final parameters:")
# Round the numeric columns for cleaner printing
species_db_rounded <- species_db %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

print(species_db_rounded)

