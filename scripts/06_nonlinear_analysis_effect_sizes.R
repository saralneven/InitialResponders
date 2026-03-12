# =============================================================================
# Script:  06_nonlinear_analysis_effect_sizes.R
# Purpose: Fit the final multivariable GLMM using the best-performing
#          transformation for each predictor (from script 5), extract fixed
#          effects, VIF, random effects, and a correlation matrix.
# Author:  Sara Neven
# Date:    2025
# =============================================================================

# Set working directory to the location of this script
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(lme4)
library(dplyr)
library(car)
library(readr)
library(splines)

# Load the dataset
data <- read_csv("../data/filtered_observations.csv", show_col_types = FALSE)

# Ensure 'Species' and 'Dep_ID' are factors
data <- data %>%
  mutate(
    Species = as.factor(Species),
    Dep_ID = as.factor(Dep_ID)
  )

data$Species <- relevel(data$Species, ref = "BrownChromis")

# Convert relevant variables from mm to meters
data <- data %>%
  mutate(
    Distance_to_Coral = Distance_to_Coral / 1000,
    Distance_to_Stimulus = Distance_to_Stimulus / 1000
  )

# LOG-transform Distance_to_Coral and Distance_to_Stimulus before scaling
data <- data %>%
  mutate(
    log_Distance_to_Coral = log1p(Distance_to_Coral),
    log_Distance_to_Stimulus = log1p(Distance_to_Stimulus)
  )

# Scale all variables (including log-transformed)
data <- data %>%
  mutate(
    Loom_Speed = as.numeric(scale(Loom_Speed)),
    Orientation_Angle = as.numeric(scale(Orientation_Angle)),
    Viewing_Angle = as.numeric(scale(Viewing_Angle)),
    log_Distance_to_Coral = as.numeric(scale(log_Distance_to_Coral)),
    log_Distance_to_Stimulus = as.numeric(scale(log_Distance_to_Stimulus)),
    Trial_Event = as.numeric(scale(Trial_Event)),
    Neighbor_Proximity = as.numeric(scale(Neighbor_Proximity))
  )

# Define final formula using best forms
final_formula <- Response_Binary ~ 
  Species +
  ns(Loom_Speed, 3) +
  Orientation_Angle +
  Viewing_Angle +
  log_Distance_to_Coral +
  log_Distance_to_Stimulus +
  ns(Trial_Event, 3) +
  Neighbor_Proximity +
  (1 | Dep_ID)

# Fit the model
control_options <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
final_model <- glmer(
  final_formula, 
  data = data, 
  family = binomial(link = "logit"),
  control = control_options
)

# Extract fixed effects
fixed_effects <- summary(final_model)$coefficients %>%
  as.data.frame() %>%
  mutate(Variable = rownames(.)) %>%
  select(Variable, Estimate, `Std. Error`, `z value`, `Pr(>|z|)`) %>%
  rename(
    Std_Error = `Std. Error`,
    z_value = `z value`,
    p_value = `Pr(>|z|)`
  )

write_csv(fixed_effects, "../outputs/S5_Multivariable_WithNonLinear.csv")

# Calculate VIF (from an lm without random effect)
lm_for_vif <- lm(
  update(final_formula, . ~ . - (1 | Dep_ID)), 
  data = data
)
vif_values <- vif(lm_for_vif)
write.csv(vif_values, "../outputs/X_MultivariableVIF_WithNonLinear.csv")

# Random effects summary
random_effects <- as.data.frame(VarCorr(final_model))
random_effects_df <- data.frame(
  Group = random_effects$grp,
  Term = random_effects$var1,
  Variance = random_effects$vcov,
  Std_Dev = sqrt(random_effects$vcov)
)
write_csv(random_effects_df, "../outputs/X_RandomEffectSummary_WithNonLinear.csv")

# Correlation matrix of scaled predictors before spline/log-transform
correlation_data <- data %>%
  select(Loom_Speed, Orientation_Angle, Viewing_Angle,
         Distance_to_Coral, Distance_to_Stimulus,
         Trial_Event, Neighbor_Proximity)

correlation_matrix <- cor(correlation_data)
write.csv(correlation_matrix, "../outputs/X_CorrelationMatrix_WithNonLinear.csv")
