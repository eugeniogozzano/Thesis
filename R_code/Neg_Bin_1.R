
######################_____________DOING PCA ________________################### 

# Load necessary libraries
library(tidyverse)
library(ggcorrplot) 
library(dplyr)
#Load data
df <- read.csv('dataset_thesis_new_hdi_v8.csv')


pca_vars <- c('GDP_merged', 'population', 'perc_athletic_prime', 'Team_Size', 'gdi', 'hdi')

df_pca <- df %>%
  dplyr::select(all_of(pca_vars)) %>%
  drop_na() %>%
  mutate(
    log_GDP = log(GDP_merged + 1),
    log_pop = log(population + 1)
  )


final_pca_vars <- c('log_GDP', 'log_pop', 'perc_athletic_prime', 'Team_Size', 'gdi', 'hdi')
X <- df_pca[, final_pca_vars]

pca_res <- prcomp(X, center = TRUE, scale. = TRUE)



explained_variance <- pca_res$sdev^2 / sum(pca_res$sdev^2)
cumulative_variance <- cumsum(explained_variance)
scree_data <- data.frame(
  PC = 1:length(explained_variance),
  Individual = explained_variance,
  Cumulative = cumulative_variance
)

loadings <- as.data.frame(pca_res$rotation)

# --- Plotting ---


p1 <- ggplot(scree_data, aes(x = PC)) +
  geom_bar(aes(y = Individual), stat = "identity", fill = "skyblue", alpha = 0.7) +
  geom_line(aes(y = Cumulative, group = 1), color = "orange", size = 1) +
  geom_point(aes(y = Cumulative), color = "orange") +
  labs(title = "Scree Plot", x = "Principal Components", y = "Variance Explained") +
  theme_minimal()


loadings_subset <- loadings[, 1:3]
p2 <- ggcorrplot(loadings_subset, lab = TRUE) +
  labs(title = "Variable Loadings (PC1, PC2, PC3)") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", limit = c(-1,1))


print(p1)
print(p2)


cat("Explained Variance Ratio:\n")
print(explained_variance)

cat("\nLoadings Table:\n")
print(loadings_subset)

cat("\nNumber of observations used for PCA:", nrow(df_pca), "\n")


######################_____________Negative Binomial ________________################### 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, MASS, broom, car, performance)

df_model <- df %>%
  arrange(iso_code_mapped, year) %>%
  group_by(iso_code_mapped) %>%
  mutate(
    # Strict Lagging (t-1)
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

#  CRITICAL FIX: RUS/BLR 2024 Data Recovery

df_model <- df_model %>%
  mutate(lag_team = case_when(
    iso_code_mapped == "RUS" & year == 2024 ~ 335, # ROC 2021 Size
    iso_code_mapped == "BLR" & year == 2024 ~ 101, # BLR 2021 Size
    TRUE ~ lag_team
  ))

#  RUN NEGATIVE BINOMIAL MODEL

nb_fit <- glm.nb(Total ~ log_lag_GDP + log_lag_pop + lag_gdi + 
                   lag_hdi + lag_team + Is_Host + lag_prime, 
                 data = df_model)

#  OUTPUT RESULTS
summary(nb_fit)

# Calculate Incidence Rate Ratios (IRR) and 95% CIs
# IRR = exp(Estimate). It represents the multiplicative change in medal count.
results_table <- tidy(nb_fit, exponentiate = TRUE, conf.int = TRUE)
print(results_table)


# Check for Multicollinearity (VIF > 5-10 is problematic)
vif(nb_fit)

# Check for Overdispersion
#check_overdispersion(nb_fit)
