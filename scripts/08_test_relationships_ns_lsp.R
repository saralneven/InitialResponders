# =============================================================================
# Script:  08_test_relationships_ns_lsp.R
# Purpose: Test all pairwise two-way interactions among ns(Loom_Speed,3),
#          Neighbor_Proximity, log(Distance_to_Coral), and Species using likelihood
#          ratio tests (LRT). Each interaction is compared to a reduced model
#          with main effects only. Results are saved per interaction and in a
#          master summary with BH-adjusted p-values.
# Author:  Sara Neven
# Date:    2025
# =============================================================================

# Set working directory to the location of this script
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Libraries
library(readr)
library(dplyr)
library(tidyr)
library(lme4)
library(splines)

# ---- INPUT / OUTPUT ----
in_path    <- "../data/filtered_observations.csv"
master_out <- "../outputs/S7_All_Interaction_LRT_Summary.csv"

# ---- LOAD RAW DATA ----
data <- read_csv(in_path, show_col_types = FALSE) %>%
  mutate(
    Species = as.factor(Species),
    Dep_ID  = as.factor(Dep_ID)
  )
data$Species <- stats::relevel(data$Species, ref = "BrownChromis")

# ---- UNITS + SAFE LOGS ----
# Convert mm to meters; apply log1p transforms before scaling
data <- data %>%
  mutate(
    Distance_to_Coral        = Distance_to_Coral / 1000,
    Distance_to_Stimulus     = Distance_to_Stimulus / 1000,
    log_Distance_to_Coral    = log1p(Distance_to_Coral),
    log_Distance_to_Stimulus = log1p(Distance_to_Stimulus)
  )

# ---- SCALE LINEAR CONTINUOUS TERMS ----

data <- data %>%
  mutate(
    Loom_Speed         = as.numeric(scale(Loom_Speed)),
    Neighbor_Proximity      = as.numeric(scale(Neighbor_Proximity)),
    log_Distance_to_Coral = as.numeric(scale(log_Distance_to_Coral)),
    log_Distance_to_Stimulus = as.numeric(scale(log_Distance_to_Stimulus)),
    Orientation_Angle = as.numeric(scale(Orientation_Angle)),
    Viewing_Angle = as.numeric(scale(Viewing_Angle)),
    Trial_Event = as.numeric(scale(Trial_Event))
  )

# ---- PRECOMPUTE SPLINE BASIS ON SCALED LOOM SPEED ----
# ns() returns a matrix with 3 columns (df = 3). We then scale each column.
ns_mat <- ns(data$Loom_Speed, df = 3)
colnames(ns_mat) <- c("nsLS1", "nsLS2", "nsLS3")

# Bind and scale each spline column to improve conditioning
data <- bind_cols(data, as.data.frame(ns_mat)) %>%
  mutate(
    nsLS1 = as.numeric(scale(nsLS1)),
    nsLS2 = as.numeric(scale(nsLS2)),
    nsLS3 = as.numeric(scale(nsLS3))
  )

# ---- CONSISTENT COMPLETE-CASE DATA ----
vars_needed <- c(
  "Response_Binary", "Dep_ID", "Species",
  "Neighbor_Proximity", "log_Distance_to_Coral",
  "nsLS1", "nsLS2", "nsLS3", "log_Distance_to_Stimulus",
  "Orientation_Angle", "Viewing_Angle", "Trial_Event"
)

analysis_data <- data %>%
  select(all_of(vars_needed)) %>%
  drop_na()

cat("Rows in analysis_data:", nrow(analysis_data), "\n")

# ---- MODEL SPEC HELPERS ----
# Grouped representation of the spline term (sum of its basis columns)
SPL_GROUP <- "(nsLS1 + nsLS2 + nsLS3)"

# Main effects (as strings to paste into formulas)
main_effects <- c(
  SPL_GROUP,                 # ns(Loom_Speed, 3) via precomputed, scaled basis
  "Neighbor_Proximity",
  "log_Distance_to_Coral",
  "Species",
  "log_Distance_to_Stimulus",
  "Orientation_Angle",
  "Viewing_Angle",
  "ns(Trial_Event, 3)"
)

# Base (reduced) model: no interactions
base_formula_txt <- paste0("Response_Binary ~ ",
                           paste(main_effects, collapse = " + "),
                           " + (1 | Dep_ID)")
base_formula <- as.formula(base_formula_txt)

# Control
ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))

# ---- FIT REDUCED MODEL ----
reduced_model <- glmer(base_formula, data = analysis_data, family = binomial, control = ctrl)

conv_note <- function(m) {
  msg <- try(m@optinfo$conv$lme4$messages, silent = TRUE)
  if (!inherits(msg, "try-error") && !is.null(msg)) paste(msg, collapse = "; ") else NA_character_
}
reduced_note <- conv_note(reduced_model)

# ---- BUILD ALL 2-WAY INTERACTIONS ----
# We'll generate the 6 pairwise interactions among the 4 terms,
# being careful to use the grouped spline term in parentheses.
terms_for_pairs <- c(SPL_GROUP, "Neighbor_Proximity", "log_Distance_to_Coral", "Species")

pairs <- combn(terms_for_pairs, 2, simplify = FALSE)

# Pretty label for console/outputs
pretty_label <- function(term) {
  if (term == SPL_GROUP) return("ns(Loom_Speed,3)")
  term
}

# ---- LOOP OVER INTERACTIONS ----
lrt_rows <- list()

for (pair in pairs) {
  var1 <- pair[[1]]
  var2 <- pair[[2]]

  # Nice readable name
  name <- paste0(pretty_label(var1), "_x_", pretty_label(var2))

  cat("Testing interaction:", name, "\n")

  # Build the interaction using grouped spline parentheses
  interaction_term <- paste0("(", var1, "):(", var2, ")")
  full_formula <- as.formula(paste(base_formula_txt, "+", interaction_term))

  # Fit full model
  full_model <- try(
    glmer(full_formula, data = analysis_data, family = binomial, control = ctrl),
    silent = TRUE
  )

  if (inherits(full_model, "try-error")) {
    cat("  -> Skipping due to fit error\n")
    lrt_rows[[name]] <- data.frame(
      Interaction = name,
      Chisq = NA_real_, Df = NA_real_, p_value = NA_real_,
      AIC_full = NA_real_, AIC_base = AIC(reduced_model), dAIC = NA_real_,
      Note_full = "fit_error", Note_base = reduced_note,
      stringsAsFactors = FALSE
    )
    next
  }

  # LRT (same rows guaranteed via analysis_data)
  lrt <- anova(reduced_model, full_model, test = "Chisq")
  lrt_df <- as.data.frame(lrt)

  Chisq   <- suppressWarnings(lrt_df$Chisq[2])
  Df      <- suppressWarnings(lrt_df$Df[2])
  p_value <- suppressWarnings(lrt_df$`Pr(>Chisq)`[2])

  note_full <- conv_note(full_model)

  lrt_rows[[name]] <- data.frame(
    Interaction = name,
    Chisq = Chisq, Df = Df, p_value = p_value,
    AIC_full = AIC(full_model), AIC_base = AIC(reduced_model),
    dAIC = AIC(full_model) - AIC(reduced_model),
    Note_full = note_full, Note_base = reduced_note,
    stringsAsFactors = FALSE
  )
}

# ---- MASTER SUMMARY + BH ----
lrt_summary <- bind_rows(lrt_rows) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

write_csv(lrt_summary, master_out)
print(lrt_summary)
