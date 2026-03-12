# =============================================================================
# Script:  07_test_species.R
# Purpose: Test whether species differences in escape response can be explained
#          by other predictors by sequentially adding each to a species-only
#          GLMM. Also visualizes the species effect before/after controlling for
#          distance to coral, and tests species differences in coral distance
#          with a Wilcoxon test.
# Author:  Sara Neven
# Date:    2025
# =============================================================================

# Set working directory to the location of this script
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(lme4)
library(readr)
library(dplyr)
library(broom.mixed)
library(splines)

# Load the dataset
data <- read_csv("../data/filtered_observations.csv", show_col_types = FALSE)

# Factors and reference level
data <- data %>%
  mutate(
    Species = as.factor(Species),
    Dep_ID  = as.factor(Dep_ID)
  )
data$Species <- relevel(data$Species, ref = "BrownChromis")

# Convert relevant variables from mm to meters
data <- data %>%
  mutate(
    Distance_to_Coral    = Distance_to_Coral / 1000,
    Distance_to_Stimulus = Distance_to_Stimulus / 1000
  )

# --- Safe logs BEFORE any scaling ---
# Use log1p to handle zeros safely (log1p(x) = log(1 + x))
data <- data %>%
  mutate(
    log_Distance_to_Coral    = log1p(Distance_to_Coral),
    log_Distance_to_Stimulus = log1p(Distance_to_Stimulus)
  )

# Scale continuous predictors used as linear terms
data <- data %>%
  mutate(
    Loom_Speed        = as.numeric(scale(Loom_Speed)),
    Orientation_Angle = as.numeric(scale(Orientation_Angle)),
    Viewing_Angle     = as.numeric(scale(Viewing_Angle)),
    log_Distance_to_Coral    = as.numeric(scale(log_Distance_to_Coral)),
    log_Distance_to_Stimulus = as.numeric(scale(log_Distance_to_Stimulus)),
    Trial_Event       = as.numeric(scale(Trial_Event)),
    Neighbor_Proximity     = as.numeric(scale(Neighbor_Proximity))
  )

### Test if species differences can be explained by other variables ###

# Predictors to test (note: splines are included directly as function calls in the formula)
extra_predictors <- c(
  "log_Distance_to_Coral",
  "log_Distance_to_Stimulus",
  "Viewing_Angle",
  "ns(Loom_Speed, 3)",
  "ns(Trial_Event, 3)",
  "Neighbor_Proximity",
  "Orientation_Angle"
)

# Base model
base_model <- glmer(Response_Binary ~ Species + (1 | Dep_ID),
                    data = data,
                    family = binomial(link = "logit"),
                    control = glmerControl(optimizer = "bobyqa",
                                           optCtrl = list(maxfun = 100000)))

# Function to fit models and extract stats
compare_models <- function(predictor) {
  fml <- as.formula(paste("Response_Binary ~ Species +", predictor, "+ (1 | Dep_ID)"))
  model <- glmer(fml, data = data, family = binomial(link = "logit"),
                 control = glmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 100000)))
  coef_summary <- tidy(model, effects = "fixed") %>%
    filter(term == "SpeciesBicolorDamselfish") %>%
    select(estimate, std.error, p.value) %>%
    rename(species_estimate = estimate, species_se = std.error, species_p = p.value)

  tibble(
    predictor = predictor,
    AIC = AIC(model),
    species_estimate = coef_summary$species_estimate,
    species_se = coef_summary$species_se,
    species_p = coef_summary$species_p
  )
}

# Run the comparison
comparison_results <- bind_rows(lapply(extra_predictors, compare_models))

# Add base model for reference
base_stats <- tidy(base_model, effects = "fixed") %>%
  filter(term == "SpeciesBicolorDamselfish") %>%
  transmute(
    predictor = "None",
    AIC = AIC(base_model),
    species_estimate = estimate,
    species_se = std.error,
    species_p = p.value
  )

comparison_results <- bind_rows(base_stats, comparison_results)

# Save results to CSV
write_csv(comparison_results, "../outputs/S6_Species_Comparison_Results.csv")

print(comparison_results)

### Test relationship Species & Distance to Coral ###
# Wilcoxon test
wilcox_test <- wilcox.test(Distance_to_Coral ~ Species, data = data)
print(wilcox_test)

wilcox_df <- data.frame(
  W_statistic = wilcox_test$statistic,
  p_value     = wilcox_test$p.value,
  method      = wilcox_test$method,
  alternative = wilcox_test$alternative
)
write_csv(wilcox_df, "../outputs/X_Wilcoxon_DCo_Species.csv")
