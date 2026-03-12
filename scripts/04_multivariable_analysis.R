# =============================================================================
# Script:  04_multivariable_analysis.R
# Purpose: Fit multivariable GLMMs with alternative social context variables
#          (Neighbor_Proximity, N, NAS), extract fixed and random effects, and
#          compute a correlation matrix among predictors.
# Author:  Sara Neven
# Date:    2025
# =============================================================================

# Set working directory to the location of this script
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# ---- Setup ----
library(lme4)
library(dplyr)
library(readr)

data_path <- "../data/filtered_observations.csv"
out_dir   <- "../outputs"

# ---- Load data ----
data <- read_csv(data_path, show_col_types = FALSE)

# Ensure factors + reference level
data <- data %>%
  mutate(
    Species = factor(Species),
    Dep_ID  = factor(Dep_ID)
  )
data$Species <- relevel(data$Species, ref = "BrownChromis")

# ---- Unit conversion (mm → m) ----
data <- data %>%
  mutate(
    Distance_to_Coral    = Distance_to_Coral / 1000,
    Distance_to_Stimulus = Distance_to_Stimulus / 1000
  )

# ---- Scale relevant variables ----
data <- data %>%
  mutate(
    Loom_Speed           = as.numeric(scale(Loom_Speed)),
    Orientation_Angle    = as.numeric(scale(Orientation_Angle)),
    Viewing_Angle        = as.numeric(scale(Viewing_Angle)),
    Distance_to_Coral    = as.numeric(scale(Distance_to_Coral)),
    Distance_to_Stimulus = as.numeric(scale(Distance_to_Stimulus)),
    Trial_Event          = as.numeric(scale(Trial_Event)),
    Neighbor_Proximity        = as.numeric(scale(Neighbor_Proximity)),
    N                    = as.numeric(scale(N)),
    NAS                  = as.numeric(scale(NAS))
  )

# ---- Helper: fit model for a given social variable ----
fit_model <- function(social_var, data) {
  formula_str <- paste0(
    "Response_Binary ~ ",
    "Species + Loom_Speed + Orientation_Angle + Viewing_Angle + ",
    "Distance_to_Coral + Distance_to_Stimulus + Trial_Event + ",
    social_var, " + (1 | Dep_ID)"
  )
  form <- as.formula(formula_str)

  control_options <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
  model <- glmer(form, data = data, family = binomial(link = "logit"), control = control_options)

  # Extract and save fixed effects
  fe <- as.data.frame(coef(summary(model)))
  fe <- fe %>%
    tibble::rownames_to_column("Variable") %>%
    rename(Estimate = Estimate, Std_Error = `Std. Error`, z_value = `z value`, p_value = `Pr(>|z|)`)
  write_csv(fe, file.path(out_dir, paste0("S3_MultivariableNormalized_", social_var, ".csv")))

  # Extract and save random effects
  re_tab <- as.data.frame(VarCorr(model))
  re_df <- data.frame(
    Group    = re_tab$grp,
    Term     = re_tab$var1,
    Variance = re_tab$vcov,
    Std_Dev  = sqrt(re_tab$vcov)
  )
  write_csv(re_df, file.path(out_dir, paste0("X_RandomEffectSummary_", social_var, ".csv")))

  # Return summary metrics for comparison
  tibble::tibble(
    Social_Var = social_var,
    AIC  = AIC(model),
    BIC  = BIC(model)
  )
}

# ---- Run models for Neighbor_Proximity, N, NAS ----
social_vars <- c("Neighbor_Proximity", "N", "NAS")
model_summaries <- lapply(social_vars, fit_model, data = data) %>% bind_rows()
print(model_summaries)

# ---- Correlation matrix ----
data_correlation <- data %>%
  select(
    Loom_Speed, Orientation_Angle, Viewing_Angle,
    Distance_to_Coral, Distance_to_Stimulus,
    Neighbor_Proximity, Trial_Event, N, NAS
  )
correlation_matrix <- cor(data_correlation, use = "pairwise.complete.obs")
write.csv(correlation_matrix, file.path(out_dir, "S8_CorrelationMatrix.csv"), row.names = TRUE)
