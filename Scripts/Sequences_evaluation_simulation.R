library(dplyr)


# Read the CSV file
data <- read.csv("C:/Users/Folashadé/Desktop/Masterarbeit/Datasets/simulation_data/simulation_test_1/dada2/DADA2_table.tsv",
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)



# Select the "sequence" column and the sample columns
selected_columns <- c("ASV_ID", "sequence", "sample_1", "sample_2", "sample_3", "sample_4", "sample_5", 
                      "sample_6", "sample_7", "sample_8", "sample_9")

# Create a new data frame with the selected columns
dada2_data <- data %>% select(all_of(selected_columns))


sample_columns <- setdiff(colnames(dada2_data), c("ASV_ID", "sequence"))

# Convert the selected sample columns to 0 or 1 based on whether the value is greater than 0
dada2_data[, sample_columns] <- lapply(dada2_data[, sample_columns], function(x) ifelse(x > 0, 1, 0))



# Define a function to merge, replace NA values with 0, and verify the length of unique values
merge_and_verify <- function(dada2_data, expected_sequences_file) {
  exp_sequences <- read.csv(expected_sequences_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  merged_data <- merge(dada2_data, exp_sequences, by.x = "sequence", by.y = "exp_sequences", all = TRUE)
  merged_data[is.na(merged_data)] <- 0
  return(merged_data)
}

# File paths for expected sequences files
file <- c(
  "C:/Users/Folashadé/Desktop/Masterarbeit/Datasets/simulation_data/abundances_mock_12_simulation_data.txt"
)

# Merge each file with dada2_data and verify
merged_data_list <- lapply(file, function(file) merge_and_verify(dada2_data, file))

# Access each merged data separately
merged_data <- merged_data_list[[1]]


# Extract .x and .y columns from each merged data
extract_x_y_columns <- function(merged_data) {
  x_columns <- grep("\\.x$", colnames(merged_data), value = TRUE)
  y_columns <- grep("\\.y$", colnames(merged_data), value = TRUE)
  extracted_data <- merged_data[, c("sequence", "ASV_ID", x_columns, y_columns)]
  return(extracted_data)
}
# Extract .x and .y columns from each merged data
extracted_mock_12 <- extract_x_y_columns(merged_data)

# Combine the extracted data into one table based on the "sequence" column
combined_data <- extracted_mock_12

# Replace NA values with 0
combined_data[is.na(combined_data)] <- 0


# Define your sample columns
sample_columns <- c(
  "sample_1", "sample_2", "sample_3", "sample_4", "sample_5", 
  "sample_6", "sample_7", "sample_8", "sample_9"
)

# .x for the output and .y for the expected. I will have to rename everything with the good names

# Iterate through the sample columns and create the "result" column
for (col in sample_columns) {
  combined_data <- combined_data %>%
    mutate(!!paste0(col, ".result") := case_when(
      .data[[paste0(col, '.x')]] == 1 & .data[[paste0(col, '.y')]] > 0 ~ "TP",
      .data[[paste0(col, '.x')]] == 0 & .data[[paste0(col, '.y')]] > 0 ~ "FN",
      .data[[paste0(col, '.x')]] == 1 & .data[[paste0(col, '.y')]] == 0 ~ "FP",
      .data[[paste0(col, '.x')]] == 0 & .data[[paste0(col, '.y')]] == 0 ~ "TN"
    ))
}

# Create an empty data frame to store the results
results_table <- data.frame(sample = character(0), recall = numeric(0), precision = numeric(0), tp = numeric(0), fn = numeric(0), fp = numeric(0))

# Iterate through the sample columns
for (col in sample_columns) {
  # Calculate TP, FN, FP, TN based on the result column
  tp <- sum(combined_data[[paste0(col, ".result")]] == "TP")
  fn <- sum(combined_data[[paste0(col, ".result")]] == "FN")
  fp <- sum(combined_data[[paste0(col, ".result")]] == "FP")
  
  recall <- tp / (tp + fn)
  precision <- tp / (tp + fp)
  
  f1_score <- 2 * (recall * precision) / (recall + precision)
  
  # Create a new row for the results table
  new_row <- data.frame(sample = col, recall = recall, precision = precision, f1_score = f1_score, tp = tp, fn = fn, fp = fp)
  
  # Bind the new row to the results table
  results_table <- bind_rows(results_table, new_row)
  
}
write.table(results_table, sep = "\t", file = "C:/Users/Folashadé/Desktop/Masterarbeit/Datasets/simulation_data/performance.tsv", row.names = FALSE)
