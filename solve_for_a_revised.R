# solve_for_a.R
#
# Functions to solve for the allometric consumption coefficient (a) that
# reconciles an observed population Q/B with individual-scale bioenergetics,
# given species life-history parameters (W_inf, k, t0, M) and an assumed
# allometric exponent (b).
#
# Accompanies: Woodstock & Bigman — "From population consumption (Q/B) to
# individual daily ration: a workflow for allometric parameterization in
# fisheries models"
#
# GitHub: https://github.com/mattwoodstock/ScaleQBtoDR

library(dplyr)

# =============================================================================
# calculate_Q_B()
# =============================================================================
#' Calculate population Q/B from individual allometric parameters.
#'
#' Simulates a stable, age-structured population using von Bertalanffy growth
#' and exponential mortality, then returns the emergent population-level
#' consumption-to-biomass ratio (Q/B, year^-1).
#'
#' @param a   Allometric consumption coefficient (g food g fish^-1 day^-1
#'            for a 1-g fish).
#' @param b   Allometric consumption exponent (dimensionless; typically -0.2
#'            to -0.4 for teleosts).
#' @param W_inf Asymptotic weight (g) from the von Bertalanffy Growth Function.
#' @param k   VBGF growth coefficient (year^-1).
#' @param t0  Theoretical age at zero length/weight (years).
#' @param M   Instantaneous natural mortality rate (year^-1).
#' @param max_age Maximum age class to simulate (years; default 20).
#'
#' @return Scalar. Annual population Q/B (year^-1).
calculate_Q_B <- function(a, b, W_inf, k, t0, M, max_age = 20) {

  ages <- 1:max_age

  # Abundance-at-age (exponential mortality from recruitment cohort of 1)
  N_at_age <- exp(-M * (ages - 1))

  # Weight-at-age via VBGF
  W_at_age <- W_inf * (1 - exp(-k * (ages - t0)))^3

  # Specific daily ration for each age class (g/g/day)
  R_spec_at_age <- a * (W_at_age ^ b)

  # Total population biomass
  total_B <- sum(N_at_age * W_at_age, na.rm = TRUE)

  # Total annual population consumption
  total_Q <- sum(N_at_age * W_at_age * R_spec_at_age * 365, na.rm = TRUE)

  return(total_Q / total_B)
}


# =============================================================================
# solve_for_a()
# =============================================================================
#' Solve for the allometric coefficient 'a' consistent with a target Q/B.
#'
#' Uses numerical root-finding (uniroot) to identify the value of 'a' that,
#' when propagated through the age-structured population model, reproduces a
#' specified population Q/B. The allometric exponent b must be assumed a priori;
#' a reasonable default for teleosts is -0.27 (Christensen & Pauly, 1992;
#' Palomares & Pauly, 1998).
#'
#' @param target_Q_B Target population Q/B (year^-1) to match.
#' @param b   Assumed allometric exponent (e.g., -0.27).
#' @param W_inf Asymptotic weight (g).
#' @param k   VBGF growth coefficient (year^-1).
#' @param t0  Theoretical age at zero length/weight (years).
#' @param M   Instantaneous natural mortality rate (year^-1).
#' @param max_age Maximum age for the simulation (default 20).
#'
#' @return Scalar. The solved value of 'a', or NA if no solution is found.
solve_for_a <- function(target_Q_B, b, W_inf, k, t0, M, max_age = 20) {

  # Input validation
  required_params <- c(target_Q_B, b, W_inf, k, t0, M)
  if (any(is.na(required_params)) || any(is.infinite(required_params))) {
    warning("One or more required parameters are NA or Inf. Returning NA.")
    return(NA_real_)
  }

  # Error function: difference between calculated and target Q/B
  error_fn <- function(a_val) {
    calculate_Q_B(a = a_val, b = b, W_inf = W_inf,
                  k = k, t0 = t0, M = M, max_age = max_age) - target_Q_B
  }

  # Check that the solution bracket contains a sign change
  low_val  <- try(error_fn(0.0001), silent = TRUE)
  high_val <- try(error_fn(2.0),    silent = TRUE)

  if (inherits(low_val,  "try-error") || inherits(high_val, "try-error") ||
      is.na(low_val)  || is.na(high_val) ||
      sign(low_val) == sign(high_val)) {
    warning(paste0(
      "Root-finding bracket [0.0001, 2.0] does not contain a sign change. ",
      "The target Q/B (", round(target_Q_B, 3), ") may be outside the ",
      "physiologically feasible range for the supplied life-history parameters. ",
      "Returning NA."
    ))
    return(NA_real_)
  }

  solution <- try(
    uniroot(f = error_fn, interval = c(0.0001, 2.0), tol = 1e-6),
    silent = TRUE
  )

  if (inherits(solution, "try-error")) {
    warning("uniroot failed to converge. Returning NA.")
    return(NA_real_)
  }

  return(solution$root)
}


# =============================================================================
# Example usage
# =============================================================================
# The species database below uses manually specified life-history parameters.
# To use your own data, replace species_db with a data.frame containing columns:
#   species, target_Q_B, assumed_b, W_inf (or L_inf + lw_a + lw_b), k, t0, M

species_db <- data.frame(
  species = c(
    # Small pelagics
    "Engraulis encrasicolus", "Sardina pilchardus", "Clupea harengus",
    "Sprattus sprattus", "Brevoortia tyrannus", "Sardinops sagax",
    # Demersal — medium
    "Merluccius merluccius", "Melanogrammus aeglefinus", "Pollachius virens",
    "Urophycis chuss", "Stenotomus chrysops", "Centropristis striata",
    # Demersal — large
    "Gadus morhua", "Anoplopoma fimbria", "Sebastes fasciatus",
    # Flatfish
    "Pleuronectes platessa", "Solea solea", "Hippoglossus hippoglossus",
    # Reef fish
    "Lutjanus campechanus", "Epinephelus morio", "Mycteroperca microlepis",
    # Large pelagics
    "Scomber scombrus", "Sarda sarda", "Coryphaena hippurus",
    "Thunnus albacares", "Thunnus thynnus", "Katsuwonus pelamis",
    # Elasmobranchs
    "Raja clavata", "Squalus acanthias", "Galeorhinus galeus"
  ),
  target_Q_B = c(
    12.0, 10.0,  9.0, 11.0,  9.5, 10.5,
     4.0,  4.2,  4.5,  5.0,  6.0,  5.5,
     3.5,  2.8,  2.5,
     4.5,  4.0,  2.0,
     3.8,  3.0,  2.8,
     7.0,  6.5,  8.0,  5.0,  4.0,  6.0,
     2.0,  1.8,  1.9
  ),
  assumed_b = rep(-0.27, 30),
  L_inf = c(
     20,  25,  38,  16,  35,  30,
     80,  70, 100,  65,  40,  45,
    130, 100,  50,
     70,  60, 250,
     90, 110, 120,
     60,  80, 150, 200, 300, 100,
    100, 130, 180
  ),
  k = c(
    0.80, 0.60, 0.50, 0.90, 0.70, 0.60,
    0.25, 0.30, 0.28, 0.35, 0.40, 0.38,
    0.20, 0.15, 0.10,
    0.20, 0.25, 0.12,
    0.22, 0.18, 0.15,
    0.60, 0.50, 0.90, 0.40, 0.18, 0.70,
    0.10, 0.08, 0.07
  ),
  t0 = c(
    -0.1, -0.2, -0.3, -0.1, -0.2, -0.2,
    -0.5, -0.4, -0.5, -0.3, -0.2, -0.2,
    -0.4, -1.0, -1.5,
    -0.6, -0.5, -0.8,
    -0.5, -0.6, -0.7,
    -0.2, -0.3, -0.1, -0.2, -0.3, -0.2,
    -2.0, -2.5, -2.2
  ),
  lw_a = c(
    0.005, 0.006, 0.005, 0.004, 0.007, 0.006,
    0.008, 0.007, 0.006, 0.009, 0.010, 0.012,
    0.008, 0.005, 0.010,
    0.009, 0.008, 0.004,
    0.015, 0.012, 0.010,
    0.010, 0.011, 0.009, 0.018, 0.015, 0.017,
    0.006, 0.005, 0.004
  ),
  lw_b = c(
    3.00, 3.00, 3.10, 3.00, 3.00, 3.00,
    3.00, 3.05, 3.00, 3.00, 3.00, 2.95,
    3.05, 3.10, 3.00,
    3.00, 3.00, 3.10,
    2.95, 3.00, 3.05,
    3.00, 3.00, 3.00, 2.90, 2.90, 2.95,
    3.10, 3.10, 3.10
  ),
  M = c(
    1.20, 1.00, 0.90, 1.30, 1.00, 1.00,
    0.40, 0.45, 0.40, 0.50, 0.60, 0.55,
    0.40, 0.25, 0.20,
    0.35, 0.30, 0.15,
    0.30, 0.25, 0.20,
    0.80, 0.70, 0.90, 0.60, 0.30, 0.80,
    0.15, 0.10, 0.10
  ),
  stringsAsFactors = FALSE
)

# Derive W_inf from length-weight regression: W_inf = a * L_inf^b
species_db <- species_db %>%
  mutate(W_inf = lw_a * (L_inf ^ lw_b)) %>%
  select(species, target_Q_B, assumed_b, W_inf, k, t0, M)

# Solve for 'a' for each species
message("Solving for allometric coefficient 'a' for each species...")

species_db$solved_a <- mapply(
  solve_for_a,
  target_Q_B = species_db$target_Q_B,
  b          = species_db$assumed_b,
  W_inf      = species_db$W_inf,
  k          = species_db$k,
  t0         = species_db$t0,
  M          = species_db$M
)

# Verify: recalculate Q/B using solved 'a'; should match target_Q_B closely
message("Verifying results...")

species_db$verified_Q_B <- mapply(
  calculate_Q_B,
  a     = species_db$solved_a,
  b     = species_db$assumed_b,
  W_inf = species_db$W_inf,
  k     = species_db$k,
  t0    = species_db$t0,
  M     = species_db$M
)

message("Complete. Final parameter table:")
print(
  species_db %>% mutate(across(where(is.numeric), ~ round(.x, 4))),
  row.names = FALSE
)
