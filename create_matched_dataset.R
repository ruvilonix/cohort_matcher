library(dplyr)
library(tidyr)

MATCHING_COLUMNS <- c("sex", "bmi", "age", "wkld_D1_AT")

# Define function to reshape participant data
reshape_participant_data <- function(df, non_duplicate_cols = NULL) {
  # Define the columns to be duplicated
  duplicate_cols <- setdiff(names(df), c("ParticipantID", "Time_Point", "Study_Visit", non_duplicate_cols))
  
  # Pivot the data to a wide format
  df_reshaped <- df %>%
    pivot_wider(
      id_cols = c(ParticipantID, all_of(non_duplicate_cols)),
      names_from = c(Study_Visit, Time_Point),
      values_from = all_of(duplicate_cols),
      names_sep = "_"
    )
  
  return(df_reshaped)
}

# Read in the data
df <- read.delim("../cpet_clinical_data_v2.tsv", sep = "\t")

# Specify non-duplicate columns
non_duplicate_cols <- c("matched_pair", "sex", "phenotype", "race", "age", 
                        "height_in", "weight_lb", "bmi", "bas_score", 
                        "q_education", "q_reclined", "q_sleeprefreshing", 
                        "q_hoursinbed")

# Reshape the data
df_merged <- reshape_participant_data(df, non_duplicate_cols)

# Select necessary columns for propensity score analysis (assuming PROPENSITY_COLUMNS is defined)
# Replace 'PROPENSITY_COLUMNS' with actual column names
prop_data <- df_merged %>%
  select(ParticipantID, phenotype, all_of(MATCHING_COLUMNS))

# Display final data
prop_data
