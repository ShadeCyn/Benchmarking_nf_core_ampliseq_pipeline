library(dplyr)
library(ggpubr)
library(ggplot2)

# Read the CSV file
data <- read.csv("C:/Users/Folashadé/Desktop/Masterarbeit/Datasets/simulation_data/simulation_test_1/dada2/DADA2_table.tsv", 
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)



# Select the "sequence" column and the sample columns
selected_columns <- c("ASV_ID", "sequence", "sample_1", "sample_2", "sample_3", "sample_4", "sample_5", 
                      "sample_6", "sample_7", "sample_8", "sample_9")

# Create a new data frame with the selected columns
dada2_data <- data %>% select(all_of(selected_columns))


sample_columns <- setdiff(colnames(dada2_data), c("ASV_ID", "sequence"))


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
columns_to_convert <- c("sample_1.x", "sample_2.x", "sample_3.x", "sample_4.x", "sample_5.x", 
                        "sample_6.x", "sample_7.x", "sample_8.x", "sample_9.x",
                        "sample_1.y", "sample_2.y", "sample_3.y", "sample_4.y", "sample_5.y", 
                        "sample_6.y", "sample_7.y", "sample_8.y", "sample_9.y")

combined_data[, columns_to_convert] <- apply(combined_data[, columns_to_convert], 2, convert_to_percentage)


# Initialize an empty data frame to store TP correlations
tp_correlation_results <- data.frame()

# Initialize an empty list to store TP data for each sample
tp_data_list <- list()

# Iterate through the sample columns and compute correlations for TP
for (col in sample_columns) {
  tp_data <- combined_data %>% filter(get(paste0(col, ".result")) == "TP")
  
  # Check if there are TP records
  if (nrow(tp_data) > 0) {
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


output_dir <- "C:/Users/Folashadé/Desktop/Masterarbeit/Datasets/simulation_data/scatterplots"
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
output_dir <- "C:/Users/Folashadé/Desktop/Masterarbeit/Datasets/simulation_data/jitterplots"
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
      scale_color_manual(values = c("red", "green", "black")) +
      scale_x_log10(limits = c(0.001, max(result_df$x, na.rm = TRUE))) +
      scale_y_log10(limits = c(0.001, max(result_df$y, na.rm = TRUE))) +
      theme_minimal() +
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
    
    # Save the jitter plot
    plot_filename <- file.path(output_dir, paste0("jitter_plot_abundance_", col, ".jpeg"))
    ggsave(plot_filename, plot = jitter_plot, width = 6, height = 4)
  }
}






