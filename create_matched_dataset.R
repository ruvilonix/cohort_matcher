library(dplyr)
library(tidyr)

# Specify non-duplicate columns
non_duplicate_cols <- c("matched_pair", "sex", "phenotype", "race", "age", 
                        "height_in", "weight_lb", "bmi", "bas_score", 
                        "q_education", "q_reclined", "q_sleeprefreshing", 
                        "q_hoursinbed", "test_site")

# Read in the data
df <- read.delim("../cpet_clinical_data_v2.tsv", sep = "\t")

ids_to_exclude <- c('PI-012', 'PI-043', 'PI-147', 'PI-008', 'PI-018', 'PI-029', 'PI-057', 'PI-082', 'PI-091', 'PI-128')
df <- df[!(df$Time_Point == "max" & df$ParticipantID %in% ids_to_exclude), ]

# First convert specified columns to numeric, coercing errors to NA, suppressing warnings
df <- df %>%
  mutate(across(all_of(setdiff(colnames(df), c('ParticipantID', 'Study_Visit', 'Time_Point', 'sex', 'phenotype'))), 
                ~ suppressWarnings(as.numeric(as.character(.)))))

# Convert these columns to factors
df <- df %>%
  mutate(across(all_of(c('Study_Visit', 'Time_Point', 'sex', 'phenotype')), 
                ~ as.factor(.)))


# Function to create a dataframe with differences between D2 and D1 for specified columns
calculate_day_diffs <- function(df, diff_cols) {
  # Convert specified columns to numeric, coercing errors to NA, suppressing warnings
  df <- df %>%
    mutate(across(all_of(diff_cols), ~suppressWarnings(as.numeric(as.character(.)))))
  
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
  df_combined <- left_join(df_single, df_diffs, by = "ParticipantID")
  
  return(diff_df)
}



# Columns where we don't need to calculate difference between days
non_diff_cols <- union(non_duplicate_cols, c("ParticipantID", "Study_Visit", "Time_Point"))

# Columns where we are going to calculate difference between days
diff_cols <- setdiff(colnames(df), non_diff_cols)

# Run the difference calculation on the specified columns
diff_df <- calculate_day_diffs(df, diff_cols)




