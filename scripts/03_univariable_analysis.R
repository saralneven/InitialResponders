# =============================================================================
# Script:  03_univariable_analysis.R
# Purpose: Fit univariable GLMMs for each predictor and extract fixed-effect
#          estimates and marginal R2 (relative to a random-intercept null model).
#          Results are saved to the outputs folder.
# Author:  Sara Neven
# Date:    2025
# =============================================================================

# Set working directory to the location of this script
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(lme4)
library(dplyr)
library(readr)
library(MuMIn)

# Load the dataset
data <- read_csv("../data/filtered_observations.csv", show_col_types = FALSE)

# Ensure 'Species' and 'Dep_ID' are factors
data <- data %>%
  mutate(
    Species = as.factor(Species),
    Dep_ID = as.factor(Dep_ID)
  )

data$Species <- relevel(data$Species, ref = "BrownChromis")

# Convert relevant variables from mm to m
data <- data %>%
  mutate(
    Distance_to_Coral = Distance_to_Coral / 1000,
    Distance_to_Stimulus = Distance_to_Stimulus / 1000
  )

# Scale relevant variables
data <- data %>%
  mutate(
    Loom_Speed = as.numeric(scale(Loom_Speed)),
    Orientation_Angle = as.numeric(scale(Orientation_Angle)),
    Viewing_Angle = as.numeric(scale(Viewing_Angle)),
    Distance_to_Coral = as.numeric(scale(Distance_to_Coral)),
    Distance_to_Stimulus = as.numeric(scale(Distance_to_Stimulus)),
    Trial_Event = as.numeric(scale(Trial_Event)),
    Neighbor_Proximity = as.numeric(scale(Neighbor_Proximity)),
    N = as.numeric(scale(N)),
    NAS = as.numeric(scale(NAS))
  )

# Define predictors
predictors <- c("Orientation_Angle", "Viewing_Angle", "Distance_to_Coral", 
                "Distance_to_Stimulus", "Loom_Speed", "Neighbor_Proximity", 
                "Species", "Trial_Event", "N", "NAS")

# Define the null model with random effect only
null_model <- glmer(Response_Binary ~ 1 + (1 | Dep_ID),
                    data = data,
                    family = binomial(link = "logit"),
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

# Fit univariable models and calculate R2 values
univariable_results <- lapply(predictors, function(var) {
  formula <- as.formula(paste("Response_Binary ~", var, "+ (1 | Dep_ID)"))
  tryCatch({
    # Fit the univariable model
    model <- glmer(formula, data = data, family = binomial(link = "logit"),
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
    model  # Return the model object
  }, error = function(e) {
    message(paste("Error with variable:", var))
    NULL
  })
})

# Extract key statistics including R2 values
univariable_summary <- lapply(seq_along(predictors), function(i) {
  model <- univariable_results[[i]]
  if (!is.null(model)) {
    coef_summary <- coef(summary(model))
    r2_values <- r.squaredGLMM(model, null = null_model)
    data.frame(
      Variable = predictors[i],
      Estimate = coef_summary[2, "Estimate"],
      Std_Error = coef_summary[2, "Std. Error"],
      z_value = coef_summary[2, "z value"],
      p_value = coef_summary[2, "Pr(>|z|)"],
      R2_marginal = r2_values[1]
    )
  } else {
    data.frame(
      Variable = predictors[i],
      Estimate = NA,
      Std_Error = NA,
      z_value = NA,
      p_value = NA,
      R2_marginal = NA
    )
  }
})

# Combine results into a single data frame
univariable_summary_df <- do.call(rbind, univariable_summary)

# Save to CSV
write_csv(univariable_summary_df, "../outputs/S2_UnivariableNormalized.csv")
