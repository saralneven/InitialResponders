# =============================================================================
# Script:  05_nonlinear_analysis_compare.R
# Purpose: Compare linear, log, and natural spline (df=2,3) transformations
#          for each continuous predictor using BIC, both in isolation and in
#          full multivariable context. A robust GLMM fitter with allFit rescue
#          is used throughout.
# Author:  Sara Neven
# Date:    2025
# =============================================================================

# Set working directory to the location of this script
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# ---- libraries ----
library(lme4)
library(readr)
library(dplyr)
library(tidyr)
library(splines)

# ---- paths ----
info_dir    <- "../data"
results_dir <- "../outputs"
data_path   <- file.path(info_dir, "filtered_observations.csv")

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ---- load & factors ----
data <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(
    Species = factor(Species),
    Dep_ID  = factor(Dep_ID)
  )
data$Species <- relevel(data$Species, ref = "BrownChromis")

# ---- safe logs + scaling ----
safe_log <- function(x) {
  x <- as.numeric(x)
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  off <- if (is.finite(min_pos)) min_pos/2 else 1
  log(x + off)
}

data <- data %>%
  mutate(
    # convert to "raw" units per your pipeline (divide by 1000 for mm->m)
    Distance_to_Coral_raw    = Distance_to_Coral / 1000,
    Distance_to_Stimulus_raw = Distance_to_Stimulus / 1000,
    Neighbor_Proximity_raw        = Neighbor_Proximity,
    Loom_Speed_raw           = Loom_Speed,
    Orientation_Angle_raw    = Orientation_Angle,
    Viewing_Angle_raw        = Viewing_Angle,
    Trial_Event_raw          = Trial_Event,

    # safe logs, then numeric scale
    log_Loom_Speed           = as.numeric(scale(safe_log(Loom_Speed_raw))),
    log_Orientation_Angle    = as.numeric(scale(safe_log(Orientation_Angle_raw))),
    log_Viewing_Angle        = as.numeric(scale(safe_log(Viewing_Angle_raw))),
    log_Distance_to_Coral    = as.numeric(scale(safe_log(Distance_to_Coral_raw))),
    log_Distance_to_Stimulus = as.numeric(scale(safe_log(Distance_to_Stimulus_raw))),
    log_Trial_Event          = as.numeric(scale(safe_log(Trial_Event_raw))),
    log_Neighbor_Proximity        = as.numeric(scale(safe_log(Neighbor_Proximity_raw))),

    # numeric scaled linear versions
    Loom_Speed           = as.numeric(scale(Loom_Speed_raw)),
    Orientation_Angle    = as.numeric(scale(Orientation_Angle_raw)),
    Viewing_Angle        = as.numeric(scale(Viewing_Angle_raw)),
    Distance_to_Coral    = as.numeric(scale(Distance_to_Coral_raw)),
    Distance_to_Stimulus = as.numeric(scale(Distance_to_Stimulus_raw)),
    Trial_Event          = as.numeric(scale(Trial_Event_raw)),
    Neighbor_Proximity        = as.numeric(scale(Neighbor_Proximity_raw))
  ) %>%
  drop_na(Response_Binary, Species, Dep_ID)


# ---- robust fitter ----
ctrl <- glmerControl(
  optimizer = "bobyqa",
  optCtrl   = list(maxfun = 200000)   # higher budget to help convergence
)

# Suppress warnings during fitting so allFit rescue logic can cleanly re-evaluate
trySuppressWarnings <- function(expr) {
  withCallingHandlers(try(expr, silent = TRUE), warning = function(w) invokeRestart("muffleWarning"))
}

fit_glmer_robust <- function(formula, dat) {
  mod <- trySuppressWarnings(
    glmer(formula, data = dat, family = binomial(link = "logit"), control = ctrl)
  )
  if (inherits(mod, "try-error")) {
    return(list(model = NULL, bic = NA_real_, note = "fit_error"))
  }
  note <- character(0)
  if (isSingular(mod, tol = 1e-4)) note <- c(note, "singular")
  oi <- try(summary(mod)$optinfo$conv$lme4$messages, silent = TRUE)
  if (!inherits(oi, "try-error") && !is.null(oi)) note <- c(note, paste(oi, collapse = "; "))

  # optional: try allFit rescue if warnings detected
  if (length(note)) {
    af <- try(lme4::allFit(mod), silent = TRUE)
    if (!inherits(af, "try-error")) {
      oks <- sapply(af, function(x) !inherits(x, "try-error") && isTRUE(x@optinfo$conv$opt == 0))
      if (any(oks)) {
        cand <- af[oks]
        devs <- sapply(cand, deviance)
        mod  <- cand[[which.min(devs)]]
        note <- c(note, "allFit_rescue")
      }
    }
  }

  list(
    model = mod,
    bic   = tryCatch(BIC(mod), error = function(e) NA_real_),
    note  = if (length(note)) paste(unique(note), collapse = " | ") else NA_character_
  )
}

# ---- variables to test ----
nonlinear_variables <- c(
  "Loom_Speed",
  "Orientation_Angle",
  "Viewing_Angle",
  "Distance_to_Coral",
  "Distance_to_Stimulus",
  "Trial_Event",
  "Neighbor_Proximity"
)

# ---- Baseline linear model (all linear predictors) — saved for reference ----
base_linear_fit <- fit_glmer_robust(
  Response_Binary ~ Loom_Speed + Orientation_Angle + Viewing_Angle +
    Distance_to_Coral + Distance_to_Stimulus + Trial_Event + Neighbor_Proximity +
    Species + (1 | Dep_ID),
  data
)
base_linear_model <- base_linear_fit$model
if (!is.null(base_linear_model)) {
  saveRDS(base_linear_model, file.path(results_dir, "X_baseline_linear_model.rds"))
}

# =========================================================
# 1) Single-variable transformation comparison (isolation)
# =========================================================
compare_transformations <- function(var) {
  others <- setdiff(nonlinear_variables, var)
  base_terms <- paste(others, collapse = " + ")
  common <- paste0(base_terms, " + Species + (1 | Dep_ID)")

  forms <- list(
    Linear  = var,
    Log     = paste0("log_", var),
    Spline2 = paste0("ns(", var, ", 2)"),
    Spline3 = paste0("ns(", var, ", 3)")
  )

  do.call(rbind, lapply(names(forms), function(name) {
    form_txt <- paste("Response_Binary ~", forms[[name]], "+", common)
    res <- fit_glmer_robust(as.formula(form_txt), data)
    data.frame(
      Variable = var, Form = name, BIC = res$bic, Note = res$note,
      stringsAsFactors = FALSE
    )
  }))
}

all_results <- bind_rows(lapply(nonlinear_variables, compare_transformations)) %>%
  filter(!is.na(BIC)) %>%
  arrange(Variable, BIC)

# extract best form per variable
best_iso <- all_results %>%
  group_by(Variable) %>%
  slice_min(BIC, with_ties = FALSE) %>%
  ungroup()

# turn into a named list compatible with the contextual reassessment
best_forms <- setNames(
  as.list(best_iso$Form),
  best_iso$Variable
)

# Map the simple names to actual RHS pieces
map_form_to_rhs <- function(var, form_name) {
  switch(form_name,
         "Linear"  = var,
         "Log"     = paste0("log_", var),
         "Spline2" = paste0("ns(", var, ", 2)"),
         "Spline3" = paste0("ns(", var, ", 3)"))
}

best_forms_rhs <- lapply(names(best_forms), function(v) map_form_to_rhs(v, best_forms[[v]]))
names(best_forms_rhs) <- names(best_forms)

# =========================================================
# 2) Multivariable reassessment (contextual)
# =========================================================
reassess_forms_in_context <- function(var, data, best_forms_rhs) {
  alt_forms <- list(
    Linear  = var,
    Log     = paste0("log_", var),
    Spline2 = paste0("ns(", var, ", 2)"),
    Spline3 = paste0("ns(", var, ", 3)")
  )
  do.call(rbind, lapply(names(alt_forms), function(form_name) {
    updated <- best_forms_rhs
    updated[[var]] <- alt_forms[[form_name]]
    formula_rhs <- paste(unlist(updated), collapse = " + ")
    full_formula <- as.formula(paste("Response_Binary ~", formula_rhs, "+ Species + (1 | Dep_ID)"))
    res <- fit_glmer_robust(full_formula, data)
    data.frame(
      Variable = var, Form = form_name, BIC = res$bic, Note = res$note,
      stringsAsFactors = FALSE
    )
  }))
}

final_context_results <- bind_rows(
  lapply(nonlinear_variables, reassess_forms_in_context, data = data, best_forms_rhs = best_forms_rhs)
) %>%
  filter(!is.na(BIC)) %>%
  arrange(Variable, BIC)

write_csv(final_context_results, file.path(results_dir, "S4_BestNonlinearForms_Context.csv"))

# ---- also save the final "best in context" set ----
best_ctx <- final_context_results %>%
  group_by(Variable) %>%
  slice_min(BIC, with_ties = FALSE) %>%
  ungroup()

write_csv(best_ctx, file.path(results_dir, "X_BestForms_Selected_in_Context.csv"))
print(best_ctx)
