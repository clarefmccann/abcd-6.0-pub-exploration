## Cleaned Pipeline for Longitudinal Measurement Invariance Analysis
# This script implements the streamlined strategy of testing for partial metric
# and partial scalar invariance between adjacent waves for multiple groups.

# -----------------------------------------------------------------------------
# 1. SETUP
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(readr)
library(lavaan)
library(semTools)

set.seed(90025)

# Define file paths (please update if necessary)
pub_data_root <- file.path(
  "/u/projects/silvers/data/ABCD/ABCD-release-6.0/physical-health/puberty"
)

# -----------------------------------------------------------------------------
# 2. DATA LOADING AND PREPARATION
# -----------------------------------------------------------------------------

# Helper function to read, rename, and select columns
load_and_prep <- function(file_path, sex_code, reporter) {
  
  if (reporter == "parent") {
    name_map <- c(pete_p = "fpete", peta_p = "peta", petb_p = "petb", petc_p = "petc", petd_p = "petd", pdss_p = "PDSS")
    if (sex_code == 0) name_map["pete_p"] <- "mpete" # Use mpete for males
  } else {
    name_map <- c(pete_y = "fpete", peta_y = "peta", petb_y = "petb", petc_y = "petc", petd_y = "petd", pdss_y = "PDSS")
    if (sex_code == 0) name_map["pete_y"] <- "mpete"
  }
  
  read_csv(file_path, show_col_types = FALSE) %>%
    rename(any_of(name_map)) %>%
    select(id, wave, age, starts_with("pet"), starts_with("pdss")) %>%
    mutate(sex = sex_code)
}

# Load all four datasets
pub_f_p <- load_and_prep(file.path(pub_data_root, "filtered_parent_tannerstages_f.csv"), 1, "parent")
pub_m_p <- load_and_prep(file.path(pub_data_root, "filtered_parent_tannerstages_m.csv"), 0, "parent")
pub_f_y <- load_and_prep(file.path(pub_data_root, "filtered_youth_tannerstages_f.csv"), 1, "youth")
pub_m_y <- load_and_prep(file.path(pub_data_root, "filtered_youth_tannerstages_m.csv"), 0, "youth")

# Define wave mapping
map_wave <- c("ses-00A" = "bl", "ses-01A" = "fu1", "ses-02A" = "fu2", 
              "ses-03A" = "fu3", "ses-04A" = "fu4", "ses-05A" = "fu5", "ses-06A" = "fu6")

# Combine, join, and pivot to wide format in one pipeline
wide_data <- bind_rows(pub_f_p, pub_m_p) %>%
  full_join(bind_rows(pub_f_y, pub_m_y), by = c("id", "wave", "sex", "age")) %>%
  mutate(wave = recode(wave, !!!map_wave)) %>%
  pivot_wider(
    id_cols = c(id, sex),
    names_from = wave,
    values_from = c(starts_with(c("age", "pet")), starts_with("pdss")),
    names_glue = "{.value}_{wave}"
  ) %>% 
  select(-contains(c("m.x","f.x", "m.y", "f.y")))
  

# Create data subsets for males and females, converting all puberty items to numeric
females_wide <- wide_data %>%
  filter(sex == 1) %>%
  mutate(across(starts_with("peta") | starts_with("petb") | starts_with("petc") | 
                starts_with("petd") | starts_with("pete"), as.numeric))
males_wide <- wide_data %>%
  filter(sex == 0) %>%
  mutate(across(starts_with("peta") | starts_with("petb") | starts_with("petc") | 
                starts_with("petd") | starts_with("pete"), as.numeric))


# -----------------------------------------------------------------------------
# 3. GENERALIZED INVARIANCE TESTING FUNCTION
# -----------------------------------------------------------------------------
# This function automates the entire process for any pair of adjacent waves,
# for any specified reporter (_p or _y).

test_invariance_pair <- function(data, wave1_suffix, wave2_suffix, reporter_suffix) {
  
  cat(paste("\n\n--- TESTING INVARIANCE BETWEEN", wave1_suffix, "AND", wave2_suffix, "---\n"))
  
  # Dynamically create variable lists based on the reporter
  base_vars <- paste0(c("peta", "petb", "petc", "petd", "pete"), reporter_suffix)
  wave1_vars <- paste0(base_vars, "_", wave1_suffix)
  wave2_vars <- paste0(base_vars, "_", wave2_suffix)
  
  # Dynamically create model syntax
  model_config <- paste(
    paste("puberty_", wave1_suffix, " =~ ", paste(wave1_vars, collapse = " + ")),
    paste("puberty_", wave2_suffix, " =~ ", paste(wave2_vars, collapse = " + ")),
    sep = "\n"
  )
  
  print(model_config)
  
  model_metric <- paste(
    paste0("puberty_", wave1_suffix, " =~ 1*", wave1_vars[1], " + l2*", wave1_vars[2], " + l3*", wave1_vars[3], " + l4*", wave1_vars[4], " + l5*", wave1_vars[5]),
    paste0("puberty_", wave2_suffix, " =~ 1*", wave2_vars[1], " + l2*", wave2_vars[2], " + l3*", wave2_vars[3], " + l4*", wave2_vars[4], " + l5*", wave2_vars[5]),
    sep = "\n"
  )
  
  print(model_metric)
  
  fit_config <- cfa(model_config, data = data, estimator = "MLR", missing = "fiml")
  fit_metric <- cfa(model_metric, data = data, estimator = "MLR", missing = "fiml")
  
  cfi_config <- fitmeasures(fit_config, "cfi.scaled")
  cfi_metric <- fitmeasures(fit_metric, "cfi.scaled")
  
  # --- Establish Baseline Metric Model ---
  successful_metric_fit <- NULL
  successful_metric_syntax <- NULL
  metric_level <- "None"
  
  if (cfi_config - cfi_metric <= 0.01) {
    cat("\nFull metric invariance ACHIEVED.\n")
    successful_metric_fit <- fit_metric
    successful_metric_syntax <- model_metric
    metric_level <- "Full Metric"
  } else {
    cat("\nFull metric invariance FAILED. Testing for partial metric invariance...\n")
    mod_indices <- modindices(fit_metric) %>% filter(op == "=~") %>% arrange(desc(mi))
    problem_item_base_name <- gsub("_[^_]+$", "", mod_indices$rhs[1])
    problem_item_index <- which(base_vars == problem_item_base_name)
    problem_item_label <- paste0("l", problem_item_index, "_free")
    cat(paste("Most non-invariant item appears to be:", problem_item_base_name, "\n"))
    
    metric_lines <- strsplit(model_metric, "\n")[[1]]
    metric_lines[2] <- gsub(paste0("l", problem_item_index, "*"), paste0(problem_item_label, "*"), metric_lines[2], fixed = TRUE)
    model_partial_metric <- paste(metric_lines, collapse = "\n")
    print(model_partial_metric)
    
    fit_partial_metric <- cfa(model_partial_metric, data = data, estimator = "MLR", missing = "fiml")
    cfi_partial_metric <- fitmeasures(fit_partial_metric, "cfi.scaled")
    
    if (cfi_config - cfi_partial_metric <= 0.01) {
      cat("Partial metric invariance ACHIEVED.\n")
      successful_metric_fit <- fit_partial_metric
      successful_metric_syntax <- model_partial_metric
      metric_level <- "Partial Metric"
    } else {
      cat("Partial metric invariance FAILED. Stopping analysis for this pair.\n")
      return()
    }
  }
  
  # --- Test for Full Scalar Invariance ---
  cat("\n--- Testing for Scalar Invariance ---\n")
  intercept_constraints <- sapply(base_vars, function(item) {
    item_index <- which(base_vars == item)
    paste0(item, "_", wave1_suffix, " ~ i", item_index, "*1; ", item, "_", wave2_suffix, " ~ i", item_index, "*1")
  })
  
  model_scalar_syntax <- paste(successful_metric_syntax, paste(intercept_constraints, collapse = "\n"), sep = "\n")
  fit_scalar <- cfa(model_scalar_syntax, data = data, estimator = "MLR", missing = "fiml")
  
  cfi_scalar <- fitmeasures(fit_scalar, "cfi.scaled")
  cfi_metric_baseline <- fitmeasures(successful_metric_fit, "cfi.scaled")
  
  # --- Final Conclusion for the Pair ---
  cat(paste("\n--- FINAL CONCLUSION for", wave1_suffix, "vs.", wave2_suffix, "---\n"))
  cat(paste("Highest level of metric invariance achieved:", metric_level, "\n"))
  
  if (cfi_metric_baseline - cfi_scalar <= 0.01) {
    cat("Full scalar invariance ACHIEVED.\n")
    cat("Implication: Latent mean comparisons are valid for this transition.\n")
  } else {
    cat("Scalar invariance FAILED.\n")
    cat("Implication: Latent mean comparisons are NOT valid for this transition.\n")
  }
}

# -----------------------------------------------------------------------------
# 4. MAIN ANALYSIS: LOOP THROUGH GROUPS AND WAVE PAIRS
# -----------------------------------------------------------------------------

# Define the groups to analyze
analysis_groups <- list(
  list(name = "Females - Parent Report", data = females_wide, reporter_suffix = "_p"),
  list(name = "Males - Parent Report",   data = males_wide,   reporter_suffix = "_p"),
  list(name = "Females - Youth Report",  data = females_wide,  reporter_suffix = "_y"),
  list(name = "Males - Youth Report",    data = males_wide,    reporter_suffix = "_y")
)

# Define the wave pairs to test
wave_pairs <- list(
  c("bl", "fu1"),
  c("fu1", "fu2"),
  c("fu2", "fu3"),
  c("fu3", "fu4"),
  c("fu4", "fu5"),
  c("fu5", "fu6")
)

# Main loop
for (group in analysis_groups) {
  cat(paste("\n\n\n========================================================\n"))
  cat(paste("    STARTING ANALYSIS FOR:", group$name, "\n"))
  cat(paste("========================================================\n"))
  
  for (pair in wave_pairs) {
    # Use a try-catch block to prevent one failed model from stopping the whole loop
    tryCatch({
      test_invariance_pair(
        data = group$data, 
        wave1_suffix = pair[1], 
        wave2_suffix = pair[2],
        reporter_suffix = group$reporter_suffix
      )
    }, error = function(e) {
      cat(paste("\nERROR during analysis for", group$name, pair[1], "vs.", pair[2], ":", e$message, "\n"))
    })
  }
}

# -----------------------------------------------------------------------------
# 5. OVERALL CONCLUSION
# -----------------------------------------------------------------------------
# Review the complete output above for the specific results of each group and
# wave-pair comparison. This automated analysis provides a comprehensive test
# of both metric (weak) and scalar (strong) invariance.

# FINAL INTERPRETATION GUIDE:
# For each group and wave transition, the output will tell you:
# 1. If full or partial metric invariance was achieved. This is the minimum
#    requirement for comparing correlations or regression paths over time.
# 2. If scalar invariance was achieved. This is required for comparing the
#    average (latent mean) of the construct over time.

# Based on previous results, you will likely find a lack of scalar invariance
# for most (if not all) groups and transitions. This is a common and important
# finding for developmental scales, highlighting that the scale functions
# differently as children age.
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# 5. OVERALL CONCLUSION
# -----------------------------------------------------------------------------
# Review the output above for the specific results of each wave-pair comparison.
# The analysis now provides a comprehensive test of both metric (weak) and
# scalar (strong) invariance for each adjacent time point.

# KEY FINDINGS (Based on previous runs, to be confirmed by new output):
# 1. Metric Invariance: Partial metric invariance is generally achievable, often
#    requiring the loading for 'pete_p' to be freed. This allows for the comparison
#    of correlations and regression paths across most time points.
# 2. Scalar Invariance: Achieving scalar invariance is more challenging. The
#    function will report whether full, partial, or no scalar invariance was found
#    for each transition.

# IMPLICATIONS:
# The level of invariance achieved determines the types of longitudinal
# comparisons that are psychometrically valid. If scalar invariance (full or
# partial with at least two invariant items) is established for a transition,
# you can meaningfully compare the latent means of the puberty construct. If only
# metric invariance is established, you are restricted to comparing correlations
# and regression coefficients.
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# 2. DEFINE THE SECOND-ORDER GROWTH MODEL
# -----------------------------------------------------------------------------
# This is a complex model with three main parts:
#   Part 1: The first-order measurement models for each reporter at each wave.
#   Part 2: The second-order growth factors (intercept, linear, quadratic)
#           defined by the first-order latent factors.
#   Part 3: The correlations between parent and youth reports.



multirater_model_syntax <- '
  
  pub_bl  =~ 1*peta_p_bl + l2p*petb_p_bl + l3p*petc_p_bl + l4p*petd_p_bl + l5p*pete_p_bl  + l1y*peta_y_bl + l2y*petb_y_bl + l3y*petc_y_bl + l4y*petd_y_bl + l5y*pete_y_bl
  pub_fu1 =~ 1*peta_p_fu1 + l2p*petb_p_fu1 + l3p*petc_p_fu1 + l4p*petd_p_fu1 + l5p_free*pete_p_fu1 + l1y*peta_y_fu1 + l2y*petb_y_fu1 + l3y*petc_y_fu1 + l4y*petd_y_fu1 + l5y_free*pete_y_fu1
  pub_fu2 =~ 1*peta_p_fu2 + l2p*petb_p_fu2 + l3p*petc_p_fu2 + l4p*petd_p_fu2 + l5p_free*pete_p_fu2 + l1y*peta_y_fu2 + l2y*petb_y_fu2 + l3y*petc_y_fu2 + l4y*petd_y_fu2 + l5y_free*pete_y_fu2
  pub_fu3 =~ 1*peta_p_fu3 + l2p*petb_p_fu3 + l3p*petc_p_fu3 + l4p*petd_p_fu3 + l5p_free*pete_p_fu3 + l1y*peta_y_fu3 + l2y*petb_y_fu3 + l3y*petc_y_fu3 + l4y*petd_y_fu3 + l5y_free*pete_y_fu3
  pub_fu4 =~ 1*peta_p_fu4 + l2p*petb_p_fu4 + l3p*petc_p_fu4 + l4p*petd_p_fu4 + l5p_free*pete_p_fu4 + l1y*peta_y_fu4 + l2y*petb_y_fu4 + l3y*petc_y_fu4 + l4y*petd_y_fu4 + l5y_free*pete_y_fu4
  pub_fu5 =~ 1*peta_p_fu5 + l2p*petb_p_fu5 + l3p*petc_p_fu5 + l4p*petd_p_fu5 + l5p_free*pete_p_fu5 + l1y*peta_y_fu5 + l2y*petb_y_fu5 + l3y*petc_y_fu5 + l4y*petd_y_fu5 + l5y_free*pete_y_fu5
 # pub_fu6 =~ 1*peta_p_fu6 + l2p*petb_p_fu6 + l3p*petc_p_fu6 + l4p*petd_p_fu6 + l5p_free*pete_p_fu6 + l1y*peta_y_fu6 + l2y*petb_y_fu6 + l3y*petc_y_fu6 + l4y*petd_y_fu6 + l5y_free*pete_y_fu6
  
  
  # --- Parent Report Factors ---
  p_bl  =~ 1*peta_p_bl + r2p*petb_p_bl + r3p*petc_p_bl + r4p*petd_p_bl + r5p*pete_p_bl
  p_fu1 =~ 1*peta_p_fu1 + r2p*petb_p_fu1 + r3p*petc_p_fu1 + r4p*petd_p_fu1 + r5p_free*pete_p_fu1
  p_fu2 =~ 1*peta_p_fu2 + r2p*petb_p_fu2 + r3p*petc_p_fu2 + r4p*petd_p_fu2 + r5p_free*pete_p_fu2
  p_fu3 =~ 1*peta_p_fu3 + r2p*petb_p_fu3 + r3p*petc_p_fu3 + r4p*petd_p_fu3 + r5p_free*pete_p_fu3
  p_fu4 =~ 1*peta_p_fu4 + r2p*petb_p_fu4 + r3p*petc_p_fu4 + r4p*petd_p_fu4 + r5p_free*pete_p_fu4
  p_fu5 =~ 1*peta_p_fu5 + r2p*petb_p_fu5 + r3p*petc_p_fu5 + r4p*petd_p_fu5 + r5p_free*pete_p_fu5
 # p_fu6 =~ 1*peta_p_fu6 + r2p*petb_p_fu6 + r3p*petc_p_fu6 + r4p*petd_p_fu6 + r5p_free*pete_p_fu6

  # --- Youth Report Factors ---
  y_bl  =~ 1*peta_y_bl + r2y*petb_y_bl + r3y*petc_y_bl + r4y*petd_y_bl + r5y*pete_y_bl
  y_fu1 =~ 1*peta_y_fu1 + r2y*petb_y_fu1 + r3y*petc_y_fu1 + r4y*petd_y_fu1 + r5y_free*pete_y_fu1
  y_fu2 =~ 1*peta_y_fu2 + r2y*petb_y_fu2 + r3y*petc_y_fu2 + r4y*petd_y_fu2 + r5y_free*pete_y_fu2
  y_fu3 =~ 1*peta_y_fu3 + r2y*petb_y_fu3 + r3y*petc_y_fu3 + r4y*petd_y_fu3 + r5y_free*pete_y_fu3
  y_fu4 =~ 1*peta_y_fu4 + r2y*petb_y_fu4 + r3y*petc_y_fu4 + r4y*petd_y_fu4 + r5y_free*pete_y_fu4
  y_fu5 =~ 1*peta_y_fu5 + r2y*petb_y_fu5 + r3y*petc_y_fu5 + r4y*petd_y_fu5 + r5y_free*pete_y_fu5
 # y_fu6 =~ 1*peta_y_fu6 + r2y*petb_y_fu6 + r3y*petc_y_fu6 + r4y*petd_y_fu6 + r5y_free*pete_y_fu6
  
  # variance of latent factors 
  
  pub_bl  ~~ 1*pub_bl 
  pub_fu1 ~~ 1*pub_fu1
  pub_fu2 ~~ 1*pub_fu2
  pub_fu3 ~~ 1*pub_fu3
  pub_fu4 ~~ 1*pub_fu4
  pub_fu5 ~~ 1*pub_fu5
 # pub_fu6 ~~ 1*pub_fu6
  
  p_bl  ~~ 1*p_bl 
  p_fu1 ~~ 1*p_fu1
  p_fu2 ~~ 1*p_fu2
  p_fu3 ~~ 1*p_fu3
  p_fu4 ~~ 1*p_fu4
  p_fu5 ~~ 1*p_fu5
 # p_fu6 ~~ 1*p_fu6
  
  y_bl  ~~ 1*y_bl 
  y_fu1 ~~ 1*y_fu1
  y_fu2 ~~ 1*y_fu2
  y_fu3 ~~ 1*y_fu3
  y_fu4 ~~ 1*y_fu4
  y_fu5 ~~ 1*y_fu5
#  y_fu6 ~~ 1*y_fu6
  
  ## covariances among latent factors 
  pub_bl  ~~ 0*p_bl  + 0*y_bl 
  pub_fu1 ~~ 0*p_fu1 + 0*y_fu1
  pub_fu2 ~~ 0*p_fu2 + 0*y_fu2
  pub_fu3 ~~ 0*p_fu3 + 0*y_fu3
  pub_fu4 ~~ 0*p_fu4 + 0*y_fu4
  pub_fu5 ~~ 0*p_fu5 + 0*y_fu5
  pub_fu6 ~~ 0*p_fu6 + 0*y_fu6


'


fit_multirater <- sem(multirater_model_syntax,
                  data = females_wide, 
                  missing = "fiml",
                  estimator = "MLR")

summary(fit_multirater, fit.measures = TRUE)

lavInspect(fit_multirater, "cov.lv")














second_order_lgm_syntax <- '
  # ===========================================================================
  # PART 1: FIRST-ORDER MEASUREMENT MODELS (METRIC INVARIANCE)
  # ===========================================================================
  # We define the latent factor for each reporter at each wave.
  # Based on your findings, we will implement partial metric invariance for the
  # pete item from FU1 onwards for both reporters.

  # --- Parent Report Factors ---
  pub_p_bl  =~ 1*peta_p_bl + l2p*petb_p_bl + l3p*petc_p_bl + l4p*petd_p_bl + l5p*pete_p_bl
  pub_p_fu1 =~ 1*peta_p_fu1 + l2p*petb_p_fu1 + l3p*petc_p_fu1 + l4p*petd_p_fu1 + l5p_free*pete_p_fu1
  pub_p_fu2 =~ 1*peta_p_fu2 + l2p*petb_p_fu2 + l3p*petc_p_fu2 + l4p*petd_p_fu2 + l5p_free*pete_p_fu2
  pub_p_fu3 =~ 1*peta_p_fu3 + l2p*petb_p_fu3 + l3p*petc_p_fu3 + l4p*petd_p_fu3 + l5p_free*pete_p_fu3
  pub_p_fu4 =~ 1*peta_p_fu4 + l2p*petb_p_fu4 + l3p*petc_p_fu4 + l4p*petd_p_fu4 + l5p_free*pete_p_fu4
  pub_p_fu5 =~ 1*peta_p_fu5 + l2p*petb_p_fu5 + l3p*petc_p_fu5 + l4p*petd_p_fu5 + l5p_free*pete_p_fu5
  pub_p_fu6 =~ 1*peta_p_fu6 + l2p*petb_p_fu6 + l3p*petc_p_fu6 + l4p*petd_p_fu6 + l5p_free*pete_p_fu6

  # --- Youth Report Factors ---
  pub_y_bl  =~ 1*peta_y_bl + l2y*petb_y_bl + l3y*petc_y_bl + l4y*petd_y_bl + l5y*pete_y_bl
  pub_y_fu1 =~ 1*peta_y_fu1 + l2y*petb_y_fu1 + l3y*petc_y_fu1 + l4y*petd_y_fu1 + l5y_free*pete_y_fu1
  pub_y_fu2 =~ 1*peta_y_fu2 + l2y*petb_y_fu2 + l3y*petc_y_fu2 + l4y*petd_y_fu2 + l5y_free*pete_y_fu2
  pub_y_fu3 =~ 1*peta_y_fu3 + l2y*petb_y_fu3 + l3y*petc_y_fu3 + l4y*petd_y_fu3 + l5y_free*pete_y_fu3
  pub_y_fu4 =~ 1*peta_y_fu4 + l2y*petb_y_fu4 + l3y*petc_y_fu4 + l4y*petd_y_fu4 + l5y_free*pete_y_fu4
  pub_y_fu5 =~ 1*peta_y_fu5 + l2y*petb_y_fu5 + l3y*petc_y_fu5 + l4y*petd_y_fu5 + l5y_free*pete_y_fu5
  pub_y_fu6 =~ 1*peta_y_fu6 + l2y*petb_y_fu6 + l3y*petc_y_fu6 + l4y*petd_y_fu6 + l5y_free*pete_y_fu6

  # ===========================================================================
  # PART 2: SECOND-ORDER GROWTH FACTORS (LINEAR ONLY)
  # ===========================================================================
  # We define the intercept (i) and linear slope (s) for each reporter.
  # Time is coded 0, 1, 2, 3, 4, 5, 6.

  # --- Parent Report Growth Factors ---
  i_p =~ 1*pub_p_bl + 1*pub_p_fu1 + 1*pub_p_fu2 + 1*pub_p_fu3 + 1*pub_p_fu4 + 1*pub_p_fu5 + 1*pub_p_fu6
  s_p =~ 0*pub_p_bl + 1*pub_p_fu1 + 2*pub_p_fu2 + 3*pub_p_fu3 + 4*pub_p_fu4 + 5*pub_p_fu5 + 6*pub_p_fu6

  # --- Youth Report Growth Factors ---
  i_y =~ 1*pub_y_bl + 1*pub_y_fu1 + 1*pub_y_fu2 + 1*pub_y_fu3 + 1*pub_y_fu4 + 1*pub_y_fu5 + 1*pub_y_fu6
  s_y =~ 0*pub_y_bl + 1*pub_y_fu1 + 2*pub_y_fu2 + 3*pub_y_fu3 + 4*pub_y_fu4 + 5*pub_y_fu5 + 6*pub_y_fu6

  # ===========================================================================
  # PART 3: MODEL SPECIFICATION & CORRELATIONS
  # ===========================================================================
  # Constrain residuals of first-order factors to 0, as is standard
  pub_p_bl ~~ 0*pub_p_bl; pub_p_fu1 ~~ 0*pub_p_fu1; pub_p_fu2 ~~ 0*pub_p_fu2; pub_p_fu3 ~~ 0*pub_p_fu3; pub_p_fu4 ~~ 0*pub_p_fu4; pub_p_fu5 ~~ 0*pub_p_fu5; pub_p_fu6 ~~ 0*pub_p_fu6
  pub_y_bl ~~ 0*pub_y_bl; pub_y_fu1 ~~ 0*pub_y_fu1; pub_y_fu2 ~~ 0*pub_y_fu2; pub_y_fu3 ~~ 0*pub_y_fu3; pub_y_fu4 ~~ 0*pub_y_fu4; pub_y_fu5 ~~ 0*pub_y_fu5; pub_y_fu6 ~~ 0*pub_y_fu6

  # Correlate the growth factors between reporters
  i_p ~~ i_y
  s_p ~~ s_y
  
  # Correlate the intercepts and slopes within each reporter
  i_p ~~ s_p
  i_y ~~ s_y
  
  # ===========================================================================
  # PART 4: ADDING COVARIATES (NEW)
  # ===========================================================================
  # To control for variables, we regress the growth factors on them.
  # Here, we control for age at baseline (`age_bl`) as an example.

  # --- Parent Report Growth Factors ---
  i_p ~ age_bl
  s_p ~ age_bl
  
  # --- Youth Report Growth Factors ---
  i_y ~ age_bl
  s_y ~ age_bl
'

# -----------------------------------------------------------------------------
# 3. FIT AND SUMMARIZE THE MODEL
# -----------------------------------------------------------------------------
# We will fit the model for the females data as an example.

# Ensure the covariate is numeric and available in the dataset
females_wide <- females_wide %>% mutate(age_bl = as.numeric(age_bl))

# FIT THE MODEL
# Note the addition of `fixed.x = FALSE` to handle missing data on covariates.
fit_lgm_females <- sem(second_order_lgm_syntax, 
                       data = females_wide, 
                       estimator = "MLR", 
                       missing = "fiml",
                       fixed.x = FALSE)

# Print a summary of the results
summary(fit_lgm_females, fit.measures = TRUE, standardized = TRUE)

# -----------------------------------------------------------------------------
# 4. INTERPRETATION OF COVARIATES
# -----------------------------------------------------------------------------
# In the summary output, look for the "Regressions" section.
# - The path from `age_bl` to `i_p` tells you if baseline age predicts the
#   parent-reported starting point of puberty.
# - The path from `age_bl` to `s_p` tells you if baseline age predicts the
#   parent-reported linear rate of change.
# - ... and so on for all six regression paths.

# -----------------------------------------------------------------------------
# 5. DERIVE INDIVIDUAL GROWTH ESTIMATES (NEW)
# -----------------------------------------------------------------------------
# After fitting the model, we can get the estimated factor scores for each
# participant's growth factors using the lavPredict() function.

# Predict the scores for all latent variables in the model
predicted_scores_f <- as.data.frame(lavPredict(fit_lgm_females))

# Select just the growth factor scores we are interested in
individual_growth_estimates_f <- predicted_scores_f %>%
  select(i_p, s_p, i_y, s_y)

# Combine these scores with the participant IDs from the original data
# This creates a new dataframe you can use for plotting or other analyses
final_scores_df_f <- females_wide %>%
  select(id) %>%
  bind_cols(individual_growth_estimates_f)

# View the first few rows of the final scores
head(final_scores_df_f)

# You can now use `final_scores_df` for further analysis. For example:
# - Plot a histogram of the linear slopes: hist(final_scores_df$s_p)
# - Use the scores as predictors in another regression model.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 4. INTERPRETATION
# -----------------------------------------------------------------------------
# In the summary output, you will now have estimates for:
#
# - Latent Variable Means (for i_p, s_p, q_p, i_y, s_y, q_y):
#   - ~1 for 'i_p' and 'i_y': The average starting puberty level.
#   - ~1 for 's_p' and 's_y': The average linear rate of change.
#   - ~1 for 'q_p' and 'q_y': The average quadratic curve (acceleration/deceleration).
#
# - Covariances (for the growth factors):
#   - i_p ~~ i_y: The correlation between parent and youth reports of initial status.
#   - s_p ~~ s_y: The correlation between their reports of the rate of change.
#   - q_p ~~ q_y: The correlation between their reports of the curve.
#
# - Variances (for the growth factors):
#   - i_p ~~ i_p: The amount of individual differences in the starting level.
#
# This model correctly accounts for the lack of scalar invariance while still
# allowing you to model the developmental trajectory and compare reporters.
# -----------------------------------------------------------------------------


piecewise_lgm_syntax <- '
  # ===========================================================================
  # PART 1: FIRST-ORDER MEASUREMENT MODELS (METRIC INVARIANCE)
  # ===========================================================================
  # This part remains the same, defining the latent factor at each wave.

  # --- Parent Report Factors ---
  pub_p_bl  =~ 1*peta_p_bl + l2p*petb_p_bl + l3p*petc_p_bl + l4p*petd_p_bl + l5p*pete_p_bl
  pub_p_fu1 =~ 1*peta_p_fu1 + l2p*petb_p_fu1 + l3p*petc_p_fu1 + l4p*petd_p_fu1 + l5p_free*pete_p_fu1
  pub_p_fu2 =~ 1*peta_p_fu2 + l2p*petb_p_fu2 + l3p*petc_p_fu2 + l4p*petd_p_fu2 + l5p_free*pete_p_fu2
  pub_p_fu3 =~ 1*peta_p_fu3 + l2p*petb_p_fu3 + l3p*petc_p_fu3 + l4p*petd_p_fu3 + l5p_free*pete_p_fu3
  pub_p_fu4 =~ 1*peta_p_fu4 + l2p*petb_p_fu4 + l3p*petc_p_fu4 + l4p*petd_p_fu4 + l5p_free*pete_p_fu4
  pub_p_fu5 =~ 1*peta_p_fu5 + l2p*petb_p_fu5 + l3p*petc_p_fu5 + l4p*petd_p_fu5 + l5p_free*pete_p_fu5
  pub_p_fu6 =~ 1*peta_p_fu6 + l2p*petb_p_fu6 + l3p*petc_p_fu6 + l4p*petd_p_fu6 + l5p_free*pete_p_fu6

  # --- Youth Report Factors ---
  pub_y_bl  =~ 1*peta_y_bl + l2y*petb_y_bl + l3y*petc_y_bl + l4y*petd_y_bl + l5y*pete_y_bl
  pub_y_fu1 =~ 1*peta_y_fu1 + l2y*petb_y_fu1 + l3y*petc_y_fu1 + l4y*petd_y_fu1 + l5y_free*pete_y_fu1
  pub_y_fu2 =~ 1*peta_y_fu2 + l2y*petb_y_fu2 + l3y*petc_y_fu2 + l4y*petd_y_fu2 + l5y_free*pete_y_fu2
  pub_y_fu3 =~ 1*peta_y_fu3 + l2y*petb_y_fu3 + l3y*petc_y_fu3 + l4y*petd_y_fu3 + l5y_free*pete_y_fu3
  pub_y_fu4 =~ 1*peta_y_fu4 + l2y*petb_y_fu4 + l3y*petc_y_fu4 + l4y*petd_y_fu4 + l5y_free*pete_y_fu4
  pub_y_fu5 =~ 1*peta_y_fu5 + l2y*petb_y_fu5 + l3y*petc_y_fu5 + l4y*petd_y_fu5 + l5y_free*pete_y_fu5
  pub_y_fu6 =~ 1*peta_y_fu6 + l2y*petb_y_fu6 + l3y*petc_y_fu6 + l4y*petd_y_fu6 + l5y_free*pete_y_fu6

  # ===========================================================================
  # PART 2: SECOND-ORDER GROWTH FACTORS (PIECEWISE)
  # ===========================================================================
  # We define an intercept (i), a first slope (s1), and a second slope (s2).
  # The knot is at fu2 (time point 2).
  # s1 loadings: 0, 1, 2, 2, 2, 2, 2
  # s2 loadings: 0, 0, 0, 1, 2, 3, 4

  # --- Parent Report Growth Factors ---
  i_p  =~ 1*pub_p_bl + 1*pub_p_fu1 + 1*pub_p_fu2 + 1*pub_p_fu3 + 1*pub_p_fu4 + 1*pub_p_fu5 + 1*pub_p_fu6
  s1_p =~ 0*pub_p_bl + 1*pub_p_fu1 + 2*pub_p_fu2 + 2*pub_p_fu3 + 2*pub_p_fu4 + 2*pub_p_fu5 + 2*pub_p_fu6
  s2_p =~ 0*pub_p_bl + 0*pub_p_fu1 + 0*pub_p_fu2 + 1*pub_p_fu3 + 2*pub_p_fu4 + 3*pub_p_fu5 + 4*pub_p_fu6

  # --- Youth Report Growth Factors ---
  i_y  =~ 1*pub_y_bl + 1*pub_y_fu1 + 1*pub_y_fu2 + 1*pub_y_fu3 + 1*pub_y_fu4 + 1*pub_y_fu5 + 1*pub_y_fu6
  s1_y =~ 0*pub_y_bl + 1*pub_y_fu1 + 2*pub_y_fu2 + 2*pub_y_fu3 + 2*pub_y_fu4 + 2*pub_y_fu5 + 2*pub_y_fu6
  s2_y =~ 0*pub_y_bl + 0*pub_y_fu1 + 0*pub_y_fu2 + 1*pub_y_fu3 + 2*pub_y_fu4 + 3*pub_y_fu5 + 4*pub_y_fu6

  # ===========================================================================
  # PART 3: MODEL SPECIFICATION & CORRELATIONS
  # ===========================================================================
  pub_p_bl ~~ 0*pub_p_bl; pub_p_fu1 ~~ 0*pub_p_fu1; pub_p_fu2 ~~ 0*pub_p_fu2; pub_p_fu3 ~~ 0*pub_p_fu3; pub_p_fu4 ~~ 0*pub_p_fu4; pub_p_fu5 ~~ 0*pub_p_fu5; pub_p_fu6 ~~ 0*pub_p_fu6
  pub_y_bl ~~ 0*pub_y_bl; pub_y_fu1 ~~ 0*pub_y_fu1; pub_y_fu2 ~~ 0*pub_y_fu2; pub_y_fu3 ~~ 0*pub_y_fu3; pub_y_fu4 ~~ 0*pub_y_fu4; pub_y_fu5 ~~ 0*pub_y_fu5; pub_y_fu6 ~~ 0*pub_y_fu6

  # Correlate the growth factors between reporters
  i_p ~~ i_y; s1_p ~~ s1_y; s2_p ~~ s2_y
  
  # Correlate the growth factors within each reporter
  i_p ~~ s1_p; i_p ~~ s2_p; s1_p ~~ s2_p
  i_y ~~ s1_y; i_y ~~ s2_y; s1_y ~~ s2_y
  
  # ===========================================================================
  # PART 4: ADDING COVARIATES
  # ===========================================================================
  # We regress the new growth factors on the covariates.
  i_p ~ age_bl; s1_p ~ age_bl; s2_p ~ age_bl
  i_y ~ age_bl; s1_y ~ age_bl; s2_y ~ age_bl
'

fit_lgm_females <- sem(piecewise_lgm_syntax, 
                       data = females_wide, 
                       estimator = "MLR", 
                       missing = "fiml",
                       fixed.x = FALSE)

# Predict the scores for all latent variables in the model
predicted_scores_f <- as.data.frame(lavPredict(fit_lgm_females))

# Select just the growth factor scores we are interested in
individual_growth_estimates_f <- predicted_scores_f %>%
  select(i_p, s1_p, i_y, s1_y, s2_p, s2_y)

# Combine these scores with the participant IDs from the original data
# This creates a new dataframe you can use for plotting or other analyses
final_scores_df_f <- females_wide %>%
  select(id) %>%
  bind_cols(individual_growth_estimates_f)

# View the first few rows of the final scores
head(final_scores_df_f)

