COLUMNS_TO_MATCH <- c("sex", "bmi", "age", "D1_AT_wkld")

# Create dataframe only with columns for matching
prop_data <- diff_df %>%
  select(ParticipantID, phenotype, all_of(COLUMNS_TO_MATCH))

# Display final data
prop_data

head(prop_data)

m.out3 <- matchit(phenotype ~ bmi + age + D1_AT_wkld,
                  data = prop_data,
                  method = "nearest", distance = "mahalanobis",
                  exact = ~ sex,
                  caliper = c(bmi = .5, age = .5, D1_AT_wkld = .5))

summary(m.out3)

plot(summary(m.out3))

plot(m.out3, type = "density", interactive = FALSE,
     which.xs = ~bmi + age + D1_AT_wkld)

# Create dataframe with matched pairs
m.data <- match.data(m.out3)

# Remove weights column
m.data <- subset(m.data, select = -weights)

# Rename matched pair column
colnames(m.data)[colnames(m.data) == 'subclass'] <- 'new_matched_pair'

# Filter `diff_df` to include only rows where `ParticipantID` is in `m.data`
full_matched <- diff_df[diff_df$ParticipantID %in% m.data$ParticipantID, ]

# Add new_matched_pair column to final df
full_matched <- merge(diff_df, m.data[c('ParticipantID', 'new_matched_pair')])
