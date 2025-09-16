#################
# MODEL FITTING #
#################
# Author: Giovanni Colavizza

# Set your own working directory here
setwd("~/Dropbox/db_projects/Odoma_projects/2024_PLOS/analysis")

options(scipen=999) # prevents excessive use of scientific notation

require(reshape2)
require(ggplot2)
require(GGally)
require(dplyr)
require(readr)

##############################################################
### Data loading and preparation: SKIP TO START HERE BELOW ###
##############################################################

# load dataset
DATASET <- read_delim(
  file = "../data/FOSM_OA_Crossref_dump_v2.csv.gz", 
  delim = ";", 
  col_types = cols(
    .default = col_guess(),  # Guess other columns
  )
)
# load dataset: OSI comparator dataset
# PLEASE NOTE that some functions below need to be adjusted to use this dataset, which is much smaller than the other one
#DATASET <- read_delim(
#  file = "../data/FOSM_OA_Crossref_dump_v2_OSI_overlap.csv.gz", 
#  delim = ";", 
#  col_types = cols(
#    .default = col_guess(),  # Guess other columns
#  )
#)

# convert language to en and fr only, keep a spare category for the rest (lack of data points)
DATASET <- DATASET %>%
  mutate(
    lang_reduce = case_when(
      lang %in% c('en', 'fr') ~ lang,
      TRUE ~ 'other'
    )
  )

# make transformations as necessary
DATASET$preprint_match <- factor(DATASET$preprint_match)
DATASET$software_created <- factor(DATASET$software_created)
DATASET$software_shared <- factor(DATASET$software_shared)
DATASET$data_created <- factor(DATASET$data_created)
DATASET$data_shared <- factor(DATASET$data_shared)
DATASET$is_oa <- factor(DATASET$is_oa)
DATASET$bso_classification <- factor(DATASET$bso_classification, ordered = FALSE)

#Extra interesting fields
DATASET$genre <- factor(DATASET$genre, ordered = FALSE)
DATASET$oa_host_type <- factor(DATASET$oa_host_type, ordered = FALSE)
# convert all NA to closed from Unpaywall status data
DATASET$unpaywall_oa_status[is.na(DATASET$unpaywall_oa_status)] <- "closed"
DATASET$unpaywall_oa_status <- factor(DATASET$unpaywall_oa_status, ordered = FALSE)
DATASET$lang_reduce <- factor(DATASET$lang_reduce, ordered = FALSE)

DATASET$software_used <- factor(DATASET$software_used)
DATASET$data_used <- factor(DATASET$data_used)

# filter for time if necessary
df_filtered <- DATASET[DATASET$year<2025,]

# log-transform citation counts (add 1 to bound between zero and infinity)
df_filtered$n_cit_tot_log <- df_filtered$cited_by_count + 1
df_filtered$n_cit_tot_log <- df_filtered$n_cit_tot_log %>% replace(is.na(.), 1)
df_filtered$n_cit_tot_log <- sapply(df_filtered$n_cit_tot_log,log)
df_filtered$n_authors_log <- df_filtered$n_authors + 1
df_filtered$n_authors_log <- df_filtered$n_authors_log %>% replace(is.na(.), 1)
df_filtered$n_authors_log <- sapply(df_filtered$n_authors_log,log)
df_filtered$n_references_log <- df_filtered$n_references + 1
df_filtered$n_references_log <- df_filtered$n_references_log %>% replace(is.na(.), 1)
df_filtered$n_references_log <- sapply(df_filtered$n_references_log,log)

# select the dataset which will be used in regressions
DATASET <- df_filtered

# save
write.csv(DATASET, "dataset/DATASET.csv", row.names=FALSE)

##################
### START HERE ###
##################
# read
DATASET <- read.csv("dataset/DATASET.csv")

# TABLE controls and dependent
summary(DATASET[c("cited_by_count", "n_authors", "n_references", "year", "month")])

summary(DATASET[c('preprint_match', 'software_used', 'software_created', 'software_shared', 'data_used', 'data_created', 'data_shared', 'is_oa')])

###############
# MODELS: OLS #
###############
# https://stats.idre.ucla.edu/r/dae/robust-regression/

require(MASS)
require(stargazer)

# BASE MODEL:
summary(m_ols <- lm(n_cit_tot_log ~ n_authors_log + n_references_log + year + month + C(preprint_match) 
                    + C(software_created) + C(software_shared) + C(data_created) + C(data_shared)
                  , data = DATASET))

# ROBUST Base model: 
summary(m_rols <- rlm(n_cit_tot_log ~ n_authors_log + n_references_log + year + month + C(preprint_match) 
                      + C(software_created) + C(software_shared) + C(data_created) + C(data_shared)
                      , data = DATASET))

###
# Output in LaTeX
###
stargazer(m_ols, m_rols, title="Results", align=TRUE, mean.sd = FALSE)

# check for multicollinearities
library(car)
vif(m_ols)
alias(m_ols)

# Check only journal articles
summary(m_ols <- lm(n_cit_tot_log ~ n_authors_log + n_references_log + year + month + C(preprint_match) 
                    + C(software_created) + C(software_shared) + C(data_created) + C(data_shared)
                    , data = DATASET %>% filter(genre == "journal-article")))

# FULL MODELS #
# Set the reference categories
DATASET$bso_classification <- relevel(DATASET$bso_classification, ref = "Medical research")
DATASET$genre <- relevel(DATASET$genre, ref = "journal-article")
DATASET$unpaywall_oa_status <- relevel(DATASET$unpaywall_oa_status, ref = "closed")
DATASET$oa_host_type <- relevel(DATASET$oa_host_type, ref = "closed")
DATASET$lang_reduce <- relevel(DATASET$lang_reduce, ref = "en")

# BSO classes
summary(m_ols <- lm(n_cit_tot_log ~ n_authors_log + n_references_log + year + month + C(preprint_match) 
                    + C(software_created) + C(software_shared) + C(data_created) + C(data_shared) + C(is_oa)
                    + C(bso_classification), data = DATASET))

# Add all fields
summary(m_ols <- lm(n_cit_tot_log ~ n_authors_log + n_references_log + year + month + C(genre) 
                    + C(preprint_match) 
                    + C(software_used) + C(software_created) + C(software_shared) + C(data_used) + C(data_created) + C(data_shared) 
                    + C(oa_host_type) + C(unpaywall_oa_status) #+ C(is_oa) # is_oa is collinear with oa_host and oa_status
                    + C(bso_classification) , data = DATASET))

# FULL MODEL:
summary(m_ols <- lm(n_cit_tot_log ~ n_authors_log + n_references_log + year + month + C(genre) 
                    + C(preprint_match) 
                    + C(software_used) + C(software_created) + C(software_shared) + C(data_used) + C(data_created) + C(data_shared) 
                    + C(lang_reduce) + C(is_oa) #+ C(unpaywall_oa_status) # swap is_oa and unpaywall_os_status as needed
                    + C(bso_classification) , data = DATASET))
# ROBUST Full model:
summary(m_rols <- rlm(n_cit_tot_log ~ n_authors_log + n_references_log + year + month + C(genre) 
                      + C(preprint_match) 
                      + C(software_used) + C(software_created) + C(software_shared) + C(data_used) + C(data_created) + C(data_shared) 
                      + C(lang_reduce) + C(is_oa) #+ C(unpaywall_oa_status) # swap is_oa and unpaywall_os_status as needed
                      + C(bso_classification) , data = DATASET))

###
# Output in LaTeX
###
stargazer(m_ols, m_rols, title="Results", align=TRUE, mean.sd = FALSE)

# Check only journal articles
summary(m_ols <- lm(n_cit_tot_log ~ n_authors_log + n_references_log + year + month
                    + C(preprint_match) 
                    + C(software_used) + C(software_created) + C(software_shared) + C(data_used) + C(data_created) + C(data_shared) 
                    + C(oa_host_type) + C(unpaywall_oa_status) #+ C(is_oa) 
                    + C(bso_classification) , data = DATASET %>% filter(genre == "journal-article")))

summary(m_ols <- lm(n_cit_tot_log ~ n_authors_log + n_references_log + year + month
                    + C(preprint_match) 
                    + C(software_used) + C(software_created) + C(software_shared) + C(data_used) + C(data_created) + C(data_shared) 
                    + C(is_oa) + C(lang_reduce)
                    + C(bso_classification) , data = DATASET %>% filter(genre == "journal-article")))

###
# Robust OLS
###

# Base model: 
summary(m_rols <- rlm(n_cit_tot_log ~ n_authors_log + n_references_log + year + month + C(preprint_match) 
                      + C(software_created) + C(software_shared) + C(data_created) + C(data_shared)
                      , data = DATASET))

# Full model:
summary(m_rols <- rlm(n_cit_tot_log ~ n_authors_log + n_references_log + year + month + C(genre) 
                      + C(preprint_match) 
                      + C(software_used) + C(software_created) + C(software_shared) + C(data_used) + C(data_created) + C(data_shared) 
                      + C(is_oa)
                      + C(bso_classification) , data = DATASET))

###
# Output in LaTeX
###
stargazer(m_ols, m_rols, title="Results", align=TRUE, mean.sd = FALSE)

###
# ORDINAL REGRESSION ON QUARTILES
###

# Create a quartile variable for citation counts
DATASET <- DATASET %>%
  mutate(
    cit_quartile = ntile(cited_by_count, 4),
    cit_quartile = factor(cit_quartile, 
                          ordered = TRUE, 
                          levels = 1:4, 
                          labels = c("Q1", "Q2", "Q3", "Q4"))
  )

# Check the distribution
table(DATASET$cit_quartile, useNA = "ifany")

# Run the ordinal logistic regression
# Note: The response variable must be a factor, ordered from lowest to highest
DATASET$cit_quartile <- factor(DATASET$cit_quartile, ordered = TRUE, levels = c("Q1", "Q2", "Q3", "Q4"))

# Run the model
m_ord <- polr(
  cit_quartile ~ n_authors_log + n_references_log + year + month + 
    genre + preprint_match + software_used + software_created + software_shared +
    data_used + data_created + data_shared + is_oa + lang_reduce + bso_classification,
  data = DATASET,
  Hess = TRUE
)

# Get the coefficient summary from polr
ctable <- coef(summary(m_ord))

# Calculate p-values for t statistics
pvals <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

# Format the p-values nicely (e.g., scientific notation, 3 significant digits)
formatted_pvals <- format.pval(pvals, digits = 3, eps = .001)

# Combine into a data frame
summary_df <- as.data.frame(ctable)
summary_df$p.value <- formatted_pvals

# View nicely
print(summary_df)

# Even nicer: use knitr::kable for markdown/html output
if (requireNamespace("knitr", quietly = TRUE)) {
  knitr::kable(summary_df, digits = 3, caption = "Ordinal Regression Coefficient Table")
}

###
# SANITY CHECKS
###

# Calculate residuals
residuals <- resid(m_ols)
# Open a new plotting window
#dev.new()
# Generate the Q-Q plot of the residuals
qqnorm(residuals)
qqline(residuals, col = "red")

# Calculate residuals
residuals <- resid(m_rols)
# Open a new plotting window
#dev.new()
# Generate the Q-Q plot of the residuals
qqnorm(residuals)
qqline(residuals, col = "red")

# check residuals
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(m_ols, las = 1)

# Check leverage of data points
# Calculate Cook's Distance
cooks.distance <- cooks.distance(m_ols)
# Plot Cook's Distance
plot(cooks.distance, type="h", main="Cook's Distance")
abline(h=4/length(cooks.distance), col="red")

# Sort Cook's Distance in descending order and get the top 10
top_influential <- order(cooks.distance, decreasing = TRUE)[1:10]
# Print the indices of the top 10 influential observations
print(top_influential)
# If you want to see the Cook's Distance values as well
top_cooks_values <- sort(cooks.distance, decreasing = TRUE)[1:10]
print(top_cooks_values)
# To retrieve the actual observations
DATASET[top_influential, c("cited_by_count", "n_authors", "n_references", "year")]

### PLOT table with percentage changes for each BSO category ###

# Use ONLY BSO class X and genre journal article
# List of BSO categories:
unique(DATASET$bso_classification)
"
[1] Computer and  information sciences         Biology (fond.)                           
[3] Chemistry                                  Medical research                          
[5] Physical sciences, Astronomy               Humanities                                
[7] Social sciences                            Earth, Ecology, Energy and applied biology
[9] Mathematics                                Engineering                               
[11] unknown 
"

summary(m_ols <- lm(n_cit_tot_log ~ n_authors_log + n_references_log + year + month
                    + C(preprint_match) 
                    + C(software_used) + C(software_created) + C(software_shared) + C(data_used) + C(data_created) + C(data_shared) 
                    + C(is_oa)
                    , data = DATASET %>% filter(genre == "journal-article" & bso_classification == "Engineering")))


# Load required libraries
library(dplyr)
library(broom)

# Initialize an empty list to store results
results_list <- list()

# Loop over unique bso_classification values
# Exclude the 'unknown' category from the unique BSO classifications
unique_bso_classes <- unique(DATASET$bso_classification) %>% 
  setdiff("unknown")  # Remove 'unknown'

for (bso_class in unique_bso_classes) {
  # Filter the dataset for the current bso_class
  subset_data <- DATASET %>% 
    filter(genre == "journal-article", bso_classification == bso_class)
  
  # Fit the regression model
  m_ols <- lm(n_cit_tot_log ~ n_authors_log + n_references_log + year + month
              + C(preprint_match)
              + C(software_used) + C(software_created) + C(software_shared)
              + C(data_used) + C(data_created) + C(data_shared)
              + C(is_oa), 
              data = subset_data)
  
  # Get the coefficients and p-values
  model_summary <- tidy(m_ols) %>% as_tibble() # Convert to tibble
  
  # Transform coefficients to percentage change
  model_summary <- model_summary %>%
    mutate(
      percentage_change = (exp(estimate) - 1) * 100
    ) %>%
    mutate(bso_classification = bso_class)  # Add classification info
  
  # Append to the results list
  results_list[[bso_class]] <- model_summary
}

# Combine results into a single data frame
final_results <- bind_rows(results_list)

# Display the final results
print(final_results)

# Optionally save the results to a file
# write.csv(final_results, "regression_results.csv", row.names = FALSE)

# Load required library
library(gt)
library(dplyr)
library(tidyr)

# Filter final_results to include only the specified variables
filtered_results <- final_results %>%
  filter(term %in% c(
    "C(preprint_match)TRUE", "C(software_used)TRUE", "C(software_created)TRUE", 
    "C(software_shared)TRUE", "C(data_used)TRUE", "C(data_created)TRUE", 
    "C(data_shared)TRUE", "C(is_oa)TRUE"
  )) %>%
  select(bso_classification, term, percentage_change, p.value)

# Add bold formatting for significant values (p-value ≤ 0.05)
filtered_results <- filtered_results %>%
  mutate(
    percentage_change = ifelse(
      p.value <= 0.05,
      paste0("**", round(percentage_change, 1), "**"),  # Add bold markdown
      as.character(round(percentage_change, 1))  # Keep non-significant values plain
    )
  )

# Pivot the data to make bso_classification columns
table_data <- filtered_results %>%
  select(bso_classification, term, percentage_change) %>%
  pivot_wider(
    names_from = bso_classification, 
    values_from = percentage_change
  )

# Create a gt table
gt_table <- table_data %>%
  gt(rowname_col = "term") %>%
  # Enable Markdown rendering for bold text
  fmt_markdown(columns = everything()) %>%
  # Add title
  tab_header(
    title = "Percentage Changes by BSO Classification",
    subtitle = "Significant values (p ≤ 0.05) are bold"
  ) %>%
  # Format columns
  cols_label(term = "Variable")

# Print the gt table
print(gt_table)
