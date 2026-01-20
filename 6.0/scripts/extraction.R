# load libraries
library(nlme)
library(dplyr)
library(purrr) # loaded but not used in this version, which is fine!

# set seed for reproducibility
set.seed(90025)

# define the root path for data
pub_data_root = "/u/projects/silvers/data/ABCD/ABCD-release-6.0/physical-health/puberty/"
proj_root = "/u/home/c/clarefmc/projects/puberty/"

# --- 1. model and helper functions (your original functions are great) ---

lme_ctrl <- lmeControl(opt = "nlminb", msMaxIter = 200, returnObject = TRUE)
nlme_ctrl <- nlmeControl(maxIter = 1000, pnlsTol = 1e-6, msMaxIter = 1000, 
                            msVerbose = FALSE, returnObject = TRUE)

# mean functions for nonlinear models
logistic_mu <- function(age, lower, upper, alpha, lambda) {
  lower + (upper - lower) / (1 + exp(-alpha * (age - lambda)))
}

gompertz_mu <- function(age, lower, upper, alpha, lambda) {
  lower + (upper - lower) * exp(-exp(-alpha * (age - lambda)))
}

logistic_model <- function(age, lower, upper, alpha, lambda) logistic_mu(age, lower, upper, alpha, lambda)
gompertz_model <- function(age, lower, upper, alpha, lambda) gompertz_mu(age, lower, upper, alpha, lambda)

# function to split data by participant id
split_by_id <- function(dat, p_train = 0.8) {
  ids <- unique(dat$id)
  tr_ids <- sample(ids, size = floor(length(ids) * p_train))
  list(
    train = dplyr::filter(dat, id %in% tr_ids),
    test  = dplyr::filter(dat, !id %in% tr_ids)
  )
}


# --- 2. refactored data loading and filtering function (new) ---
# this function avoids repeating the cleaning steps in the loop.

# --- 2. load + filter (patched) ---
load_and_filter_data <- function(sex, pub_data_root, proj_root) {
  filepath <- paste0(pub_data_root, "parent_tannerstages_", sex, ".csv")
  regressorspath <- paste0(proj_root, "pub_exp_output/regressors/6.0/regression_env_pds_", sex, "_p.csv")
  
  raw_data <- readr::read_csv(filepath, show_col_types = FALSE) |>
    dplyr::rename_with(tolower) |>
    # drop any unnamed/garbage columns like ...1, ...2, etc.
    dplyr::select(-tidyselect::starts_with("..."), dplyr::everything()) |>
    dplyr::mutate(
      id = as.character(id),
      rel_family_id = as.character(rel_family_id),
      age = suppressWarnings(as.numeric(age)),
      pdss = suppressWarnings(as.numeric(pdss))
    ) |>
    # lose non-finite rows early so nlme doesn't choke later
    dplyr::filter(is.finite(age), is.finite(pdss))
  
  regressors <- readr::read_csv(regressorspath, show_col_types = FALSE) |>
    dplyr::rename_with(tolower) |>
    dplyr::transmute(id = as.character(id)) |>
    dplyr::distinct()
  
  cleaned <- dplyr::anti_join(raw_data, regressors, by = "id")
  
  # keep ids with > 3 observations
  sampled <- cleaned |>
    dplyr::add_count(id, name = "n_obs") |>
    dplyr::filter(n_obs > 3) |>
    dplyr::select(-n_obs)
  
  # pick one sibling per family
  id_counts <- sampled |>
    dplyr::count(rel_family_id, id, name = "n_obs")
  
  preferred_ids <- id_counts |>
    dplyr::group_by(rel_family_id) |>
    dplyr::slice_max(order_by = n_obs, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(rel_family_id, id)
  
  filtered_data <- sampled |>
    dplyr::semi_join(preferred_ids, by = c("rel_family_id", "id"))
  
  stopifnot(all(c("id", "pdss", "age") %in% names(filtered_data)))
  filtered_data
}

# --- 3. model fitting and evaluation functions (with minor style changes) ---

# --- 3. fit functions (patched starts) ---
fit_all_nlme <- function(train_df, sex, lower = 1, upper = 5) {
  
  m_lin  <- lme(pdss ~ age,                 random = ~ 1 | id, data = train_df, method = "ML", control = lme_ctrl)
  m_quad <- lme(pdss ~ age + age2,          random = ~ 1 | id, data = train_df, method = "ML", control = lme_ctrl)
  m_cub  <- lme(pdss ~ age + age2 + age3,   random = ~ 1 | id, data = train_df, method = "ML", control = lme_ctrl)
  
  # nlme likes a named numeric vector for start more than that list-wrapper
  start_vals <- if (sex == "f") c(alpha = 0.95, lambda = 11.4) else c(alpha = 0.95, lambda = 12.5)
  
  m_log <- nlme(pdss ~ logistic_model(age, lower, upper, alpha, lambda),
                fixed  = alpha + lambda ~ 1, random = lambda ~ 1,
                groups = ~ id, data = train_df, start = start_vals,
                method = "ML", control = nlme_ctrl,
                na.action = na.omit)
  
  m_gmp <- nlme(pdss ~ gompertz_model(age, lower, upper, alpha, lambda),
                fixed  = alpha + lambda ~ 1, random = lambda ~ 1,
                groups = ~ id, data = train_df, start = start_vals,
                method = "ML", control = nlme_ctrl,
                na.action = na.omit)
  
  list(linear = m_lin, quadratic = m_quad, cubic = m_cub, logistic = m_log, gompertz = m_gmp)
}


# --- 4. main wrapper function ---

compare_models_nlme <- function(sex, pub_data_root, proj_root) {
  df <- load_and_filter_data(sex, pub_data_root, proj_root) %>%
    dplyr::mutate(
      age = as.numeric(age),
      age2 = age^2,
      age3 = age^3
    )
  
  if (!all(c("id","pdss","age") %in% names(df))) {
    stop("Data missing required columns: have ", paste(names(df), collapse=", "))
  }
  
  sp <- split_by_id(df, 0.8)
  train_df <- sp$train
  test_df  <- sp$test
  
  fits <- fit_all_nlme(train_df, sex)
  metrics <- purrr::map_dfr(fits, ~ nlme_metrics(.x, train_df, test_df), .id = "model")
  
  list(fits = fits, metrics = dplyr::arrange(metrics, aic))
}

# --- 5. execution loop (updated to store results) ---

# initialize an empty list to store results for each sex
results_by_sex <- list()

# --- 5. execution loop (patched end block) ---
for (sex in c("f", "m")) {
  model_comparison <- compare_models_nlme(sex, pub_data_root, proj_root)
  
  best_model_name <- model_comparison$metrics$model[1]
  best_fit_cv <- model_comparison$fits[[best_model_name]]
  
  cat("\n--- Results for sex:", toupper(sex), "---\n")
  cat("Best model determined by AIC:", best_model_name, "\n\n")
  print(model_comparison$metrics)
  
  df_all <- load_and_filter_data(sex, pub_data_root, proj_root) |>
    dplyr::mutate(age2 = age^2, age3 = age^3)
  
  # refit on full data; keep same structure but allow fresh groups
  final_fit <- tryCatch(
    update(best_fit_cv, data = df_all, control = if (inherits(best_fit_cv, "lme")) lme_ctrl else nlme_ctrl),
    error = function(e) {
      # if update is grumpy, refit from scratch with the chosen formula
      if (best_model_name == "linear") {
        lme(pdss ~ age, random = ~1|id, data = df_all, method = "ML", control = lme_ctrl)
      } else if (best_model_name == "quadratic") {
        lme(pdss ~ age + age2, random = ~1|id, data = df_all, method = "ML", control = lme_ctrl)
      } else if (best_model_name == "cubic") {
        lme(pdss ~ age + age2 + age3, random = ~1|id, data = df_all, method = "ML", control = lme_ctrl)
      } else if (best_model_name == "logistic") {
        start_vals <- if (sex == "f") c(alpha = 0.95, lambda = 11.4) else c(alpha = 0.95, lambda = 12.5)
        nlme(pdss ~ logistic_model(age, 1, 5, alpha, lambda),
             fixed = alpha + lambda ~ 1, random = lambda ~ 1,
             groups = ~ id, data = df_all, start = start_vals,
             method = "ML", control = nlme_ctrl, na.action = na.omit)
      } else { # gompertz
        start_vals <- if (sex == "f") c(alpha = 0.95, lambda = 11.4) else c(alpha = 0.95, lambda = 12.5)
        nlme(pdss ~ gompertz_model(age, 1, 5, alpha, lambda),
             fixed = alpha + lambda ~ 1, random = lambda ~ 1,
             groups = ~ id, data = df_all, start = start_vals,
             method = "ML", control = nlme_ctrl, na.action = na.omit)
      }
    }
  )
  
  re_tbl <- tibble::as_tibble(ranef(final_fit), rownames = "id")
  if ("lambda" %in% names(re_tbl)) {
    individual_params <- dplyr::rename(re_tbl, timing = lambda)
  } else {
    individual_params <- re_tbl
  }
  
  results_by_sex[[sex]] <- list(
    metrics = model_comparison$metrics,
    best_model_name = best_model_name,
    final_fit_on_all_data = final_fit,
    individual_params = individual_params   # <- fixed name
  )
}


# you can now access the results like so:
# names(results_by_sex)
# head(results_by_sex$f$individual_params)
# results_by_sex$m$metrics