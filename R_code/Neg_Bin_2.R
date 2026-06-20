

library(dplyr)
library(tidyverse)
library(countrycode)

#  LOAD & PREPARE DATA
df <- read.csv('dataset_thesis_new_hdi_v8.csv')

# FIXED CRITICAL FLAW: Retain metadata ('iso_code_mapped') before filtering NAs
df_pca_clean <- df %>%
  mutate(
    log_GDP = log(GDP_merged + 1),
    log_pop = log(population + 1)
  ) %>%
  # Keep the metadata columns along with your PCA variables
  select(iso_code_mapped, year, log_GDP, log_pop, perc_athletic_prime, Team_Size, gdi, hdi) %>%
  drop_na() %>%
  # Automatically generate continents using the ISO country codes
  mutate(continent = countrycode(iso_code_mapped, origin = "iso3c", destination = "continent")) %>%
  # Fill any unresolved regions (like custom codes) as 'Other'
  mutate(continent = ifelse(is.na(continent), "Other", continent))

#  RUN PCA ON NUMERIC MATRIX ONLY
final_pca_vars <- c('log_GDP', 'log_pop', 'perc_athletic_prime', 'Team_Size', 'gdi', 'hdi')
X <- df_pca_clean[, final_pca_vars]

pca_res <- prcomp(X, center = TRUE, scale. = TRUE)

#  EXTRACT COORDINATES & DEFINE HIGHLIGHT LOGIC
df_plot <- df_pca_clean %>%
  mutate(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    # Create text labels ONLY for USA and China
    highlight_label = ifelse(iso_code_mapped %in% c("USA", "CHN"), iso_code_mapped, ""),
    # Make highlighted points larger than the others
    point_size = ifelse(iso_code_mapped %in% c("USA", "CHN"), 4, 1.5)
  )

# Extract variable loading vectors (arrows)
loadings <- as.data.frame(pca_res$rotation) %>%
  mutate(Variable = rownames(.))

# Expansion factor to stretch arrows to fit the scale of the country points
arrow_scale <- 3

#  PLOT HIGH-QUALITY ggplot2 BIPLOT
ggplot() +
  # Plot all country points colored by continent with custom sizes
  geom_point(data = df_plot, aes(x = PC1, y = PC2, color = continent, size = point_size), alpha = 0.5) +
  # Overlay bold labels specifically for USA and CHN positions
  geom_text(data = df_plot %>% filter(highlight_label != ""), 
            aes(x = PC1, y = PC2, label = highlight_label), 
            fontface = "bold", vjust = -1.2, size = 4.5, color = "black") +
  # Draw loading arrows for structural variables
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
               arrow = arrow(length = unit(0.2, "cm")), color = "darkred", lwd = 0.8) +
  # Label loading vectors
  geom_text(data = loadings, aes(x = PC1 * arrow_scale * 1.1, y = PC2 * arrow_scale * 1.1, label = Variable),
            color = "darkred", fontface = "bold", size = 3.5) +
  # Final formatting
  scale_size_identity() + 
  theme_minimal(base_size = 12) +
  labs(
    title = "Principal Component Analysis: Global Olympic Profiles",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 1), "% Variance)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 1), "% Variance)"),
    color = "Continent"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

######################_____________Negative Binomial ________________################### 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, MASS, broom, car, performance)

df_model <- df %>%
  arrange(iso_code_mapped, year) %>%
  group_by(iso_code_mapped) %>%
  mutate(
    #  Lagging (t-1)
    lag_GDP      = lag(GDP_merged),
    lag_pop      = lag(population),
    lag_gdi      = lag(gdi),
    lag_hdi     = lag(hdi),
    lag_team     = lag(Team_Size),
    lag_prime    = lag(perc_athletic_prime),
    # Log Transformations
    log_lag_GDP  = log(lag_GDP + 1),
    log_lag_pop  = log(lag_pop + 1)
  ) %>%
  ungroup()

#  CRITICAL FIX: RUS/BLR 2024 

df_model <- df_model %>%
  mutate(lag_team = case_when(
    iso_code_mapped == "RUS" & year == 2024 ~ 335, # ROC 2021 Size
    iso_code_mapped == "BLR" & year == 2024 ~ 101, # BLR 2021 Size
    TRUE ~ lag_team
  ))

#ADJUSTEMENT OF THE DATA 
if (!require("glmmTMB")) install.packages("glmmTMB")
library(glmmTMB)
library(car)
library(broom.mixed)


df_model <- df_model %>%
  mutate(
    # Create GDP per Capita to kill VIF
    lag_gdp_pc = lag_GDP / lag_pop,
    log_lag_gdp_pc = log(lag_gdp_pc + 1),
    log_lag_pop = log(lag_pop + 1)
  )

# RUN ZERO-INFLATED NEGATIVE BINOMIAL (ZINB)
# We use log_lag_pop in the 'ziformula' because small population 
# is the primary driver of 'structural zeros'.
zinb_fit <- glmmTMB(
  Total ~ log_lag_gdp_pc + log_lag_pop + lag_gdi + lag_hdi + lag_team + Is_Host, 
  ziformula = ~ log_lag_pop, 
  family = nbinom2, 
  data = df_model
)

#  RESULTS & DIAGNOSTICS
summary(zinb_fit)

# Check VIF again - it should be MUCH lower now (< 5 for most)

performance::check_collinearity(zinb_fit)

# INTERPRETATION (IRR)

results_table <- tidy(zinb_fit, effects = "fixed", exponentiate = TRUE, conf.int = TRUE)
print(results_table)

#OTHER MODEL

zinb2_fit = glmmTMB(Total ~ log_lag_pop + lag_gdi + lag_hdi + lag_team + Is_Host,
                    ziformula = ~ log_lag_pop, 
                    family = nbinom2, 
                    data = df_model
)
summary(zinb2_fit)

performance::check_collinearity(zinb2_fit)
# #  RUN NEGATIVE BINOMIAL MODEL
# 
# nb_fit <- glm.nb(Total ~ log_lag_GDP + log_lag_pop + lag_gdi + 
#                    lag_hdi + lag_team + Is_Host + lag_prime, 
#                  data = df_model)
# 
# #  OUTPUT RESULTS
# summary(nb_fit)
# 
# # Calculate Incidence Rate Ratios (IRR) and 95% CIs
# # IRR = exp(Estimate). It represents the multiplicative change in medal count.
# results_table <- tidy(nb_fit, exponentiate = TRUE, conf.int = TRUE)
# print(results_table)
# 
# 
# # Check for Multicollinearity (VIF > 5-10 is problematic)
# vif(nb_fit)
# 
# # Check for Overdispersion
# #check_overdispersion(nb_fit)
# 
# #plotting
# plot(nb_fit$fitted.values, nb_fit$residuals)

#OTHER CHECK QUALITY MODEL
library(DHARMa)
simulateResiduals(zinb2_fit, plot = T)
plotResiduals(zinb2_fit)
#WE SEE THAT THERE ARE SOME OUTLIERS IN THE MODEL AND THE CURVES OF THE PLOT ARE NOT SUPER STRAIGHT. THAT SUGGESTS NON-LINEARITY. LET'S BUILD THE GAM TO SEE IF IT MAKES ANY DIFFERENCE

#PSEUDO R^2 FOR ZINB
library(performance)

r2_values <- r2(zinb2_fit)
print(r2_values)

#It seems veery high, too good to be true

cor(df_model$Total, df_model$lag_team, use = "complete.obs")

# This aligns predictions and actuals perfectly
preds <- predict(zinb2_fit, type = "response")
actuals <- model.frame(zinb2_fit)$Total  

plot(preds, actuals, 
     xlab = "Predicted Medals", 
     ylab = "Actual Medals", 
     main = "ZINB: Predicted vs. Actual",
     pch = 16, col = rgb(0, 0, 1, 0.3)) 
abline(a = 0, b = 1, col = "red", lwd = 2) 

#Look at the values here and see where are they in the pca

#The R^2 of 0.99 was too "fake", but the correlation of 0.87 and the non perfect plot are telling me thaet the model is not perfect. (good news: no data leakage).
#Should I only keep the R^2 of the GAM?


library(mgcv)
gam_model <- gam(Total ~ s(log_lag_pop) + s(lag_gdi) + s(lag_hdi) + 
                   s(lag_team) + Is_Host, 
                 family = nb(), 
                 data = df_model, 
                 method = "REML")

summary(gam_model)

plot(gam_model, pages = 1, shade = TRUE, residuals = TRUE)

#Use Picewise linear model for interpretation of team size starting from the GAM plot
#PICEWISE LINEAR MODEL
#  knot value based on your GAM visual inspection (when the curve of team size reaches a plateau)
knot_value <- 150

# Piecewise Variables

df_model <- df_model %>%
  mutate(
    # The base slope 
    team_base = lag_team, 
    
    # The "penalty" or changed slope (only triggers for athletes AFTER the knot)
    # pmax(0, ...) ensures this is exactly 0 for teams smaller than the knot
    team_extra = pmax(0, lag_team - knot_value) 
  )

#  Piecewise ZINB
piecewise_fit <- glmmTMB(
  Total ~ log_lag_pop + lag_gdi + lag_hdi + team_base + team_extra + Is_Host, 
  ziformula = ~ log_lag_pop, 
  family = nbinom2, 
  data = df_model
)

summary(piecewise_fit)

############__________DOING XGBOOST AS WELL_______________########################
# 1. LOAD REQUIRED PACKAGES
if (!require("xgboost")) install.packages("xgboost")
if (!require("pdp")) install.packages("pdp")
library(tidyverse)
library(xgboost)
library(pdp)
library(glmmTMB)

#  PREPARE THE METADATA AND CLEAN OBSERVATIONS
# We must use rows that have complete cases for our covariates
model_vars <- c("log_lag_pop", "lag_gdi", "lag_hdi", "lag_team", "Is_Host")

df_xgb_clean <- df_model %>%
  select(Total, all_of(model_vars)) %>%
  drop_na()

# Convert features into a numeric matrix required by XGBoost
X_matrix <- as.matrix(df_xgb_clean[, model_vars])
y_vector <- df_xgb_clean$Total

# TRAIN THE EXPLORATORY XGBOOST MODEL
# Hyperparameters optimized for regularizing smaller panel datasets
xgb_fit <- xgboost(
  data = X_matrix,
  label = y_vector,
  nrounds = 100,
  objective = "count:poisson", # Suitable baseline objective for raw medal counts
  max_depth = 4,
  eta = 0.1,
  subsample = 0.8,
  verbose = 0
)

#  GENERATE AND PLOT PARTIAL DEPENDENCE (PDP) FOR TEAM SIZE
# This isolates the marginal effect of lag_team on predicted medals
pdp_team <- pdp::partial(
  xgb_fit,
  pred.var = "lag_team",
  train = X_matrix,
  type = "regression"
)

# Render the plot to visually confirm the flattening threshold (the knot)
plotPartial(pdp_team, smooth = TRUE, lwd = 2, col = "darkblue",
            xlab = "Lagged Team Size", ylab = "Partial Dependence (Predicted Medals)")

#  DEFINE PIECEWISE VARIABLES BASED ON THE XGBOOST THRESHOLD
# Replace 100 with the exact point where the PDP plot curve flattens out
xgb_knot <- 150 

df_final_model <- df_model %>%
  mutate(
    team_base  = lag_team,
    team_extra = pmax(0, lag_team - xgb_knot)
  )

#  ESTIMATE THE FINAL PARAMETRIC PIECEWISE ZINB MODEL
piecewise_zinb <- glmmTMB(
  Total ~ log_lag_pop + lag_gdi + lag_hdi + team_base + team_extra + Is_Host, 
  ziformula = ~ log_lag_pop, 
  family = nbinom2, 
  data = df_final_model
)

summary(piecewise_zinb)

############___________robust standard errors_________________-#################

#  Load the clubSandwich package for robust clustered standard errors
library(clubSandwich)

#  Generate the robust variance-covariance matrix (returns the 7 count model parameters)
robust_vcov <- vcovCR(piecewise_zinb_clean, cluster = df_complete$iso_code_mapped, type = "CR2")

# Extract ONLY the conditional model coefficients to align with the vcov matrix
estimates_cond <- fixef(piecewise_zinb_clean)$cond

# Extract the parameter names from the robust matrix to ensure a perfect match
target_names <- colnames(robust_vcov)

# ilter and align the estimates to match the robust matrix dimensions exactly
matched_estimates <- estimates_cond[target_names]

# Extract robust standard errors from the diagonal
robust_se <- sqrt(diag(robust_vcov))

#  Compute Z-statistics and P-values using matched dimensions
z_stat <- matched_estimates / robust_se
p_vals <- 2 * (1 - pnorm(abs(z_stat)))

#  Bind into the final, bulletproof diagnostic table
robust_results_table <- data.frame(
  Variable  = target_names,
  Estimate  = round(matched_estimates, 5),
  Robust_SE = round(robust_se, 5),
  Z_value   = round(z_stat, 3),
  P_value   = format.pval(p_vals, digits = 4, eps = 0.0001),
  row.names = NULL
)

# Print the final robust table
print(robust_results_table)
