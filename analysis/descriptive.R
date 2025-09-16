########################
# DESCRIPTIVE ANALYSIS #
########################
# Author: Giovanni Colavizza

# Set your own working directory here
setwd("~/Dropbox/db_projects/Odoma_projects/2024_PLOS/analysis")

options(scipen=999) # prevents excessive use of scientific notation

require(reshape2)
require(ggplot2)
require(GGally)
require(dplyr)
require(readr)

# read
DATASET <- read.csv("dataset/DATASET.csv")

summary(DATASET)

# correlations
corr <- round(cor(DATASET[, c("cited_by_count", "year", "month")], method = "pearson", use="complete.obs"), 2)
upper <- corr
upper[upper.tri(corr, diag = TRUE)] <- ""
upper <- as.data.frame(upper)
upper
#ggpairs(DATASET[, c("cited_by_count", "year", "month")])

corr <- round(cor(DATASET[, c("cited_by_count", "year", "month", 'preprint_match', 'software_created', 'software_shared', 'data_created', 'data_shared', 'is_oa')], method = "pearson", use="complete.obs"), 2)
upper <- corr
upper[upper.tri(corr, diag = TRUE)] <- ""
upper <- as.data.frame(upper)
upper
#ggpairs(DATASET[, c('cited_by_count', 'preprint_match', 'software_used', 'software_created', 'software_shared', 'data_used', 'data_created', 'data_shared')])

# check for lognormal distribution (and compare vs Pareto)
qqnorm(DATASET$n_cit_tot_log)
qex <- function(x) qexp((rank(x)-.375)/(length(x)+.25))
plot(qex(DATASET$n_cit_tot),DATASET$n_cit_tot_log)

# check value counts

mat <- stack(table(DATASET$bso_classification)) # USE, filtering low value counts
mat <- mat[order(mat$values, decreasing = TRUE), ]
tail(mat,20)

mat <- stack(table(DATASET$genre)) # USE, filtering low value counts
mat <- mat[order(mat$values, decreasing = TRUE), ]
tail(mat,20)

mat <- stack(table(DATASET$oa_host_type)) # USE, filtering low value counts
mat <- mat[order(mat$values, decreasing = TRUE), ]
tail(mat,20)

mat <- stack(table(DATASET$unpaywall_oa_status)) # USE, filtering low value counts
mat <- mat[order(mat$values, decreasing = TRUE), ]
tail(mat,20)

mat <- stack(table(DATASET$lang)) # USE, filtering low value counts
mat <- mat[order(mat$values, decreasing = TRUE), ]
head(mat,10)

# DESCRIPTIVE PLOTS

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# 0: Publications over time
ggplot(DATASET, aes(x = factor(year))) +
  geom_bar() +
  labs(x = "Year", y = "Number of Publications", title = "Number of Publications per Year") +
  theme_minimal()

# 1: % of OSI over time

# Calculate the percentage of 1s for each variable by year
DATASET_aggregated <- DATASET %>%
  group_by(year) %>%
  summarise(across(c(preprint_match,software_used, software_created, software_shared, data_used, data_created, data_shared), ~mean(.x, na.rm = TRUE) * 100)) # Calculate the mean and convert to percentage

# Reshape the data from wide to long format for plotting
DATASET_long <- reshape2::melt(DATASET_aggregated, id.vars = "year", variable.name = "variable", value.name = "percentage")

# Plotting the data
ggplot(DATASET_long, aes(x = year, y = percentage, linetype = variable)) +
  geom_line(aes(color = variable)) + # Drawing the lines
  scale_color_manual(values = rep("black", 7)) + # Set the colors to black
  theme_minimal(base_size = 16) + # Minimal theme
  labs(x = "Year", y = "Percentage of publications", title = "Adoption of OSI over time") +
  theme(legend.title = element_blank()) + # Remove legend title
  scale_linetype_manual(values = c("solid", "dotted", "twodash", "longdash", "dotdash", "dashed", "solid")) # Custom line types

# With colors
ggplot(DATASET_long, aes(x = year, y = percentage, linetype = variable, color = variable)) +
  geom_line(linewidth = 1.2) + # Drawing the lines
  scale_color_manual(values = c(
    "preprint_match" = "blue",
    "software_used" = "green",
    "software_created" = "green",
    "software_shared" = "green",
    "data_used" = "purple",
    "data_created" = "purple",
    "data_shared" = "purple"
  )) + # Custom colors by group
  scale_linetype_manual(values = c(
    "preprint_match" = "solid",
    "software_used" = "dotted",
    "software_created" = "twodash",
    "software_shared" = "solid",
    "data_used" = "dotted",
    "data_created" = "twodash",
    "data_shared" = "solid"
  )) + # Custom line types
  scale_x_continuous(breaks = c(2020, 2021, 2022)) + # Specify x-axis breaks
  theme_minimal(base_size = 16) + # Minimal theme
  labs(x = "Year", y = "Percentage of publications", title = "Adoption of OSI over time") +
  theme(legend.title = element_blank()) # Remove legend title

# 2: OSI by BSO class

division <- "Medical research"
# Step 1: Filter by BSO class
division_1_data <- DATASET %>% filter(bso_classification == division)

# Step 2: Calculate percentages
percentages <- division_1_data %>%
  summarise(across(c(preprint_match,software_used, software_created, software_shared, data_used, data_created, data_shared), ~mean(.x, na.rm = TRUE) * 100)) %>%
  pivot_longer(cols = c(preprint_match,software_used, software_created, software_shared, data_used, data_created, data_shared), names_to = "variable", values_to = "percentage")

# Step 3: Plot
ggplot(percentages, aes(x = variable, y = percentage, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal(base_size = 12) +
  labs(x = "Variable", y = "Percentage of 1s", title = paste("Percentage of 1s in ",division)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ALL categories

# Step 1: Calculate percentages for each BSO category
percentages_all <- DATASET %>%
  group_by(bso_classification) %>%
  summarise(across(c(preprint_match, software_used, software_created, software_shared, data_used, data_created, data_shared), ~mean(.x, na.rm = TRUE) * 100)) %>%
  pivot_longer(cols = c(preprint_match, software_used, software_created, software_shared, data_used, data_created, data_shared), 
               names_to = "variable", values_to = "percentage")

# Step 2: Plot with facets for each BSO category
ggplot(percentages_all, aes(x = variable, y = percentage, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal(base_size = 12) +
  labs(x = "Variable", y = "Percentage of 1s", title = "Percentage of 1s by BSO Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ bso_classification) # Creates a separate panel for each BSO category

# 3: citation counts x OA status

library(ggplot2)

ggplot(DATASET, aes(x = cited_by_count, fill = as.factor(is_oa))) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_x_continuous(trans = 'log1p') +   # log scale if highly skewed, optional
  scale_fill_manual(values = c("gray70", "dodgerblue"), name = "Open Access") +
  labs(x = "Citation count", y = "Density", 
       title = "Citation Distributions by Open Access Status") +
  theme_minimal()

ggplot(
  data = tidyr::drop_na(DATASET, cited_by_count, is_oa),
  aes(x = as.factor(is_oa), y = cited_by_count, fill = as.factor(is_oa))
) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  scale_y_continuous(trans = 'log1p') +
  scale_fill_manual(values = c("gray70", "dodgerblue"), name = "Open Access") +
  labs(
    x = "Open Access",
    y = "Citation count",
    title = "Citation Distributions by OA Status"
  ) +
  theme_minimal()

