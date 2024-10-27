library(dplyr)
library(tidyr)
library(MatchIt)

file_path <- "../cpet_clinical_data_v2.tsv"

##################################
##### Add difference columns #####
##################################

# Specify columns that are the same between days
non_duplicate_cols <- c("matched_pair", "sex", "phenotype", "race", "age", 
                        "height_in", "weight_lb", "bmi", "bas_score", 
                        "q_education", "q_reclined", "q_sleeprefreshing", 
                        "q_hoursinbed", "test_site")

# Read in the data
df <- read.delim(file_path, sep = "\t")

# Exclude participants that didn't meet criteria for max, only for max time point
ids_to_exclude <- c('PI-012', 'PI-043', 'PI-147', 'PI-008', 'PI-018', 'PI-029', 'PI-057', 'PI-082', 'PI-091', 'PI-128')
df <- df[!(df$Time_Point == "max" & df$ParticipantID %in% ids_to_exclude), ]

# Convert specified columns to numeric, coercing errors to NA, suppressing warnings
df <- df %>%
  mutate(across(all_of(setdiff(colnames(df), c('ParticipantID', 'Study_Visit', 'Time_Point', 'sex', 'phenotype'))), 
                ~ suppressWarnings(as.numeric(as.character(.)))))

# Convert these columns to factors
df <- df %>%
  mutate(across(all_of(c('Study_Visit', 'Time_Point', 'sex', 'phenotype')), 
                ~ as.factor(.)))


# Function to create a dataframe with differences between D2 and D1 for specified columns
calculate_day_diffs <- function(df, diff_cols) {
  
  # Step 1: Calculate differences
  df_diffs <- df %>%
    group_by(ParticipantID, Time_Point) %>%
    arrange(ParticipantID, Time_Point, Study_Visit) %>%
    summarise(
      across(!all_of(c(diff_cols, "Study_Visit")), first),
      across(all_of(diff_cols), 
             list(
               diff_abs = ~if_else(any(is.na(.)), NA_real_, last(.) - first(.)),
               diff_pct = ~if_else(any(is.na(.)) | first(.) == 0, NA_real_, (last(.) - first(.))/first(.) * 100)
             )
      ),
      .groups = "drop"
    ) %>%
    pivot_wider(
      id_cols = !matches("_diff_(abs|pct)$"),
      names_from = Time_Point,
      values_from = matches("_diff_(abs|pct)$"),
      names_glue = "{Time_Point}_{.value}"
    )
  
  # Step 2: Pivot wider to get single day metrics (D1, D2) for specified columns
  df_single <- df %>%
    pivot_wider(
      id_cols = ParticipantID,
      names_from = c(Time_Point, Study_Visit),
      values_from = all_of(diff_cols),
      names_glue = "{Study_Visit}_{Time_Point}_{.value}"
    )
  
  # Step 3: Merge single day metrics with the diffs
  diff_df <- left_join(df_single, df_diffs, by = "ParticipantID")
  
  return(diff_df)
}

# Columns where we don't need to calculate difference between days
non_diff_cols <- union(non_duplicate_cols, c("ParticipantID", "Study_Visit", "Time_Point"))

# Columns where we are going to calculate difference between days
diff_cols <- setdiff(colnames(df), non_diff_cols)

# Run the difference calculation on the specified columns and create the full cohort dataframe
diff_df <- calculate_day_diffs(df, diff_cols)


################################
### Create matching pairs ######
################################


# Find matches. Set caliper values to max SD that rows in matched pairs can differ. 
m.out3 <- matchit(phenotype ~ bmi + age + D1_AT_wkld,
                  data = diff_df,
                  method = "nearest", distance = "mahalanobis",
                  exact = ~ sex,
                  caliper = c(bmi = 1, age = 1, D1_AT_wkld = 0.5))

# Display summary of matching
summary(m.out3)

# Plots to visualize matching
plot(summary(m.out3))

plot(m.out3, type = "density", interactive = FALSE,
     which.xs = ~bmi + age + D1_AT_wkld)

# Create dataframe with matched pairs
m.data <- match.data(m.out3)

# Remove weights column
m.data <- subset(m.data, select = -weights)

# Rename matched pair column
colnames(m.data)[colnames(m.data) == 'subclass'] <- 'new_matched_pair'

# Final matched dataset
full_matched <- m.data

# Export full dataset CSV
write.csv(full_matched, 'matched_dataset.csv', row.names = FALSE)


#################################
##### Statistical tests #########
#################################


only_diffs_df <- full_matched %>% 
  select(contains('diff'), 'ParticipantID', 'phenotype')

check_distinct_values <- function(df, column_name) {
  # Get unique non-NA values
  distinct_values <- unique(na.omit(df[[column_name]]))
  # Check if we have at least 2
  length(distinct_values) >= 2
}

results <- data.frame(feature=character(),
                      p_value=double(),
                      stringsAsFactors = FALSE)

for (col in setdiff(colnames(only_diffs_df), c('ParticipantID', 'phenotype'))) {
  
  # Skip if not enough distinct values
  if (!check_distinct_values(only_diffs_df, col)) {
    message(sprintf("Skipping %s - fewer than 2 distinct values", col))
    next
  }
  
  # Remove rows that have a null for the current column
  nulls_removed <- only_diffs_df[complete.cases(only_diffs_df[col]), ]
  
  # Create the formula dynamically using as.formula()
  formula <- as.formula(paste(col, "~ phenotype"))
  
  test_result <- wilcox.test(formula, data = nulls_removed)
  result_df <- data.frame(feature=col,
                       p_value=test_result$p.value,
                       stringsAsFactors = FALSE)  
  
  results <- rbind(results, result_df)
}

# Add BH-corrected p-values as a new column
results$p_adjusted <- p.adjust(results$p_value, method = "BH")

# Sort by adjusted p-value
results <- results[order(results$p_value, results$p_adjusted), ]


###############################
### Output results as image ###
###############################

results_table <- tableGrob(results,
                           rows = NULL
)
h = grid::convertHeight(sum(results_table$heights), "in", TRUE)
w = grid::convertWidth(sum(results_table$widths), "in", TRUE)
ggplot2::ggsave("results.png", results_table, width=w, height=h)
