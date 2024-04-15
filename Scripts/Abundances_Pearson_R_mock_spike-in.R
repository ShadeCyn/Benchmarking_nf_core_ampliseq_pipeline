library(dplyr)
library(ggpubr)
library(ggplot2)

# Read the CSV file
data <- read.csv("C:/Users/Folashadé/Desktop/Masterarbeit/Datasets/test_1_v2/dada2/DADA2_table.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)



# Select the "sequence" column and the sample columns
selected_columns <- c("ASV_ID", "sequence", "SRR3832217", "SRR3832216", "SRR3832221", "SRR3985723", "SRR3985733", "mock_16", "SRR3985725", "SRR3985724", "SRR3985728", "SRR3985729", "SRR3985740", "SRR3985731", "SRR3985732", "SRR3985748")

# Create a new data frame with the selected columns
dada2_data <- data %>% select(all_of(selected_columns))


sample_columns <- setdiff(colnames(dada2_data), c("ASV_ID", "sequence"))

# Define a function to merge, replace NA values with 0, and verify the length of unique values
merge_and_verify <- function(dada2_data, expected_sequences_file) {
  exp_sequences <- read.csv(expected_sequences_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  merged_data <- merge(dada2_data, exp_sequences, by.x = "sequence", by.y = "exp_sequence", all = TRUE)
  merged_data[is.na(merged_data)] <- 0
  return(merged_data)
}

# File paths for expected sequences files
files <- c(
  "C:/Users/Folashadé/Desktop/Masterarbeit/Datasets/test_1_output/expected_sequences_end_files/expected_sequences_all_mix3_filtered_abundances.txt",
  "C:/Users/Folashadé/Desktop/Masterarbeit/Datasets/test_1_output/expected_sequences_end_files/mock_16_abundances_percentages.txt",
  "C:/Users/Folashadé/Desktop/Masterarbeit/Datasets/test_1_output/expected_sequences_end_files/expected_sequences_Q_filtered_abundances.txt"
)

# Merge each file with dada2_data and verify
merged_data_list <- lapply(files, function(file) merge_and_verify(dada2_data, file))

# Access each merged data separately
merged_data_all_mix3 <- merged_data_list[[1]]
merged_data_mock_16 <- merged_data_list[[2]]
merged_data_Q <- merged_data_list[[3]]


# Extract .x and .y columns from each merged data
extract_x_y_columns <- function(merged_data) {
  x_columns <- grep("\\.x$", colnames(merged_data), value = TRUE)
  y_columns <- grep("\\.y$", colnames(merged_data), value = TRUE)
  extracted_data <- merged_data[, c("sequence", "ASV_ID", x_columns, y_columns)]
  return(extracted_data)
}

# Extract .x and .y columns from each merged data
extracted_all_mix3 <- extract_x_y_columns(merged_data_all_mix3)
extracted_mock_16 <- extract_x_y_columns(merged_data_mock_16)
extracted_Q <- extract_x_y_columns(merged_data_Q)

# Combine the extracted data into one table based on the "sequence" column
combined_data <- Reduce(function(x, y) merge(x, y, by = c("sequence", "ASV_ID"), all = TRUE), list(extracted_all_mix3, extracted_mock_16, extracted_Q))

# Replace NA values with 0
combined_data[is.na(combined_data)] <- 0


# Define your sample columns
sample_columns <- c(
  "SRR3832217", "SRR3832216", "SRR3832221", "SRR3985723", "SRR3985733", "mock_16",
  "SRR3985725", "SRR3985724", "SRR3985728", "SRR3985729", "SRR3985740", "SRR3985731",
  "SRR3985732", "SRR3985748"
)


# Iterate through the sample columns and create the "result" column and compute correlation
correlation_columns <- c()
for (col in sample_columns) {
  x_column <- paste0(col, ".x")
  y_column <- paste0(col, ".y")
  result_column <- paste0(col, ".result")
  correlation_column <- paste0(col, ".correlation")
  
  combined_data <- combined_data %>%
    mutate(!!result_column := case_when(
      .data[[x_column]] > 0 & .data[[y_column]] > 0 ~ "TP",
      .data[[x_column]] == 0 & .data[[y_column]] > 0 ~ "FN",
      .data[[x_column]] > 0 & .data[[y_column]] == 0 ~ "FP",
      .data[[x_column]] == 0 & .data[[y_column]] == 0 ~ "TN"
    ),
    !!correlation_column := cor(.data[[x_column]], .data[[y_column]], method = "pearson"))
  
  correlation_columns <- c(correlation_columns, correlation_column)
}

# Create a data frame to store the counts
result_table <- data.frame(Sample = character(), TP = numeric(), FN = numeric(), FP = numeric(), stringsAsFactors = FALSE)

# Iterate through the sample columns and compute counts
for (col in sample_columns) {
  result_column <- paste0(col, ".result")
  
  result_table <- result_table %>%
    add_row(
      Sample = col,
      TP = sum(combined_data[[result_column]] == "TP"),
      FN = sum(combined_data[[result_column]] == "FN"),
      FP = sum(combined_data[[result_column]] == "FP")
    )
}

# Define a function to convert non-zero values to percentages
convert_to_percentage <- function(x) {
  return(ifelse(x != 0, (x / sum(x)) * 100, 0))
}

# Apply the function to the relevant columns
columns_to_convert <- c("SRR3832217.x", "SRR3832221.x", "SRR3985723.x", "SRR3985733.x", "SRR3985725.x",
                        "SRR3985728.x", "SRR3985729.x", "SRR3985740.x", "SRR3985731.x", "SRR3985732.x",
                        "SRR3985748.x", "SRR3832217.y", "SRR3832221.y", "SRR3985723.y", "SRR3985733.y",
                        "SRR3985725.y", "SRR3985728.y", "SRR3985729.y", "SRR3985740.y", "SRR3985731.y",
                        "SRR3985732.y", "SRR3985748.y", "mock_16.x", "mock_16.y", "SRR3832216.x",
                        "SRR3985724.x", "SRR3832216.y", "SRR3985724.y")

# Initialize an empty data frame to store TP correlations
tp_correlation_results <- data.frame()

# Initialize an empty list to store TP data for each sample
tp_data_list <- list()

# Iterate through the sample columns and compute correlations for TP
for (col in sample_columns) {
  tp_data <- combined_data %>% filter(get(paste0(col, ".result")) == "TP")
  
  # Check if there are TP records
  if (nrow(tp_data) > 0) {
    
    # Convert TP values to percentages
    tp_data[, paste0(col, ".x")] <- convert_to_percentage(tp_data[[paste0(col, ".x")]])
    tp_data[, paste0(col, ".y")] <- convert_to_percentage(tp_data[[paste0(col, ".y")]])
    
    # Calculate correlation for TP
    tp_cor <- cor(tp_data[[paste0(col, ".x")]], tp_data[[paste0(col, ".y")]], method = "pearson")
    
    tp_correlation_results <- rbind(tp_correlation_results, data.frame(
      Sample = col,
      Correlation = tp_cor
    ))
    
    # Store TP data for each sample
    tp_data_list[[col]] <- tp_data
  }
}

output_dir <- "C:/Users/Folashadé/Desktop/Masterarbeit/Datasets/test_1_output/scatterplots"
dir.create(output_dir, showWarnings = FALSE)  # Create the directory if it doesn't exist

# Iterate through each sample and create scatterplots
for (sample in names(tp_data_list)) {
  tp_data <- tp_data_list[[sample]]
  
  # Create scatterplot
  scatter_plot <-ggscatter(tp_data, x = paste0(sample, ".x"), y = paste0(sample, ".y"),
                            add = "reg.line", conf.int = TRUE,
                            cor.coef = TRUE, cor.method = "pearson",
                            xlab = "computed_abundance", ylab = "exp_abundance",
                            ggtheme = theme_pubr()) +
                           labs(title = paste("Scatterplot for", sample))
   
  
  # Save the scatterplot
  plot_filename <- file.path(output_dir, paste0(sample, "_scatterplot.png"))
  ggsave(plot_filename, plot = scatter_plot, width = 6, height = 4)
}


# Create a directory to save plots
output_dir <- "C:/Users/Folashadé/Desktop/Masterarbeit/Datasets/test_1_output/jitterplots"
dir.create(output_dir, showWarnings = FALSE)


# Iterate through each sample
for (col in sample_columns) {
  # Filter data for TP, FN, and FP
  tp_data <- combined_data %>% filter(get(paste0(col, ".result")) == "TP")
  fn_data <- combined_data %>% filter(get(paste0(col, ".result")) == "FN")
  fp_data <- combined_data %>% filter(get(paste0(col, ".result")) == "FP")
  
  # Create a data frame for TP
  tp_df <- data.frame(Type = "TP", x = tp_data[[paste0(col, ".x")]], y = tp_data[[paste0(col, ".y")]])
  
  # Create a data frame for FN
  fn_df <- if (nrow(fn_data) > 0) {
    data.frame(Type = "FN", x = ifelse(fn_data[[paste0(col, ".x")]] == 0, 0, fn_data[[paste0(col, ".x")]]), y = fn_data[[paste0(col, ".y")]])
  } else {
    fn_df <- data.frame(Type = "FN", x = numeric(0), y = numeric(0))
  }
  
  # Create a data frame for FP or set to NULL if no FP records
  fp_df <- if (nrow(fp_data) > 0) {
    data.frame(Type = "FP", x = fp_data[[paste0(col, ".x")]], y = fp_data[[paste0(col, ".y")]])
  } else {
    fp_df <- NULL
  }
  
  # Combine data frames
  result_df <- dplyr::bind_rows(tp_df, fn_df, fp_df)
  
  # Add a new column for sample
  result_df$Sample <- col
  
  # Create a jitter plot for each sample
  if (nrow(result_df) > 0) {  # Check if the data frame is not empty
    jitter_plot <- ggplot(result_df, aes(x = x, y = y, color = Type)) +
      geom_jitter(shape= 18, size=3.3, alpha = 1.0) +
      labs(title = paste("Jitter Plot for", col),
           x = "Calculated Abundance", y = "Expected Abundance") +
      scale_color_manual(values = c("black", "red", "green")) +
      scale_x_log10(limits = c(0.001, max(result_df$x, na.rm = TRUE))) +
      scale_y_log10(limits = c(0.001, max(result_df$y, na.rm = TRUE))) +
      theme_minimal() +
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
    
    # Save the jitter plot
    plot_filename <- file.path(output_dir, paste0("jitter_plot_abundance_", col, ".jpeg"))
    ggsave(plot_filename, plot = jitter_plot, width = 6, height = 4)
  }
}




