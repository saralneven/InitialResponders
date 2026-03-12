# =============================================================================
# Script:  02_random_variable.R
# Purpose: Test whether trial session ID (Dep_ID) should be included as a
#          random effect by comparing a mixed-effects GLMM to a fixed-effects
#          GLM using a likelihood ratio test and BIC.
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

# Scale continuous predictors
data <- data %>%
  mutate(
    Loom_Speed = as.numeric(scale(Loom_Speed)),
    Orientation_Angle = as.numeric(scale(Orientation_Angle)),
    Viewing_Angle = as.numeric(scale(Viewing_Angle)),
    Distance_to_Coral = as.numeric(scale(Distance_to_Coral)),
    Distance_to_Stimulus = as.numeric(scale(Distance_to_Stimulus)),
    Trial_Event = as.numeric(scale(Trial_Event)),
    Neighbor_Proximity = as.numeric(scale(Neighbor_Proximity))
  )

### FIT MIXED AND FIXED MODELS ###

# Mixed model: all main effects with Dep_ID as random effect
base_formula <- Response_Binary ~ 
  Species + Loom_Speed + Orientation_Angle + Viewing_Angle +
  Distance_to_Coral + Distance_to_Stimulus + Trial_Event + 
  Neighbor_Proximity + (1 | Dep_ID)

# Fit the base model
base_model <- glmer(
  formula = base_formula,
  data = data,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

# Fixed-effects-only model (no random effect) for comparison
fixed_model <- glm(
  Response_Binary ~ Species + Loom_Speed + Orientation_Angle + 
                    Viewing_Angle + Distance_to_Coral + Distance_to_Stimulus + 
                    Trial_Event + Neighbor_Proximity,
  data = data,
  family = binomial
)

logLik_fixed <- logLik(fixed_model)
logLik_base <- logLik(base_model)

# Compute the chi-squared statistic
chisq_value <- -2 * (as.numeric(logLik_fixed) - as.numeric(logLik_base))
df <- attr(logLik_base, "df") - attr(logLik_fixed, "df")

# Compute the p-value
p_value <- pchisq(chisq_value, df, lower.tail = FALSE)

cat("Chi-squared statistic:", chisq_value, "\n")
cat("Degrees of freedom:", df, "\n")
cat("P-value:", p_value, "\n")

# Save LRT results to CSV
lrt_results <- data.frame(
  Chi_Sq    = chisq_value,
  df        = df,
  p_value   = p_value,
  BIC_mixed = BIC(base_model),
  BIC_fixed = BIC(fixed_model),
  Delta_BIC = BIC(fixed_model) - BIC(base_model)
)
write_csv(lrt_results, "../outputs/S1_RandomEffect_Test.csv")