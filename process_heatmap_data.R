library(tidyverse)

# Step 1: List all CSV files
csv_files <- list.files("./data", pattern = "^tt_.*\\.csv$")
csv_files <- paste0("./data/", csv_files)

# Step 2: Function to parse filename and extract parameters
parse_filename <- function(filename) {
  parts <- str_split(basename(str_remove(filename, ".csv")), "_")[[1]]
  output <- parts[2]
  R0 <- as.numeric(sub("R0_", "", parts[4]))
  gammaa <- as.numeric(sub("gammaa_", "", parts[6]))
  gammap <- as.numeric(sub("gammap_", "", parts[8]))
  return(tibble(output = output, R0 = R0, gammaa = gammaa, gammap = gammap, filename = filename))
}

# Step 3: Parse all filenames
file_info <- map_dfr(csv_files, parse_filename)

# Step 4: Function to read and reshape a single CSV
read_and_reshape <- function(filename, output, R0, gammaa, gammap) {
  # Read CSV
  data <- read.csv(filename)[-1] # Remove row names column
  
  # Assume 8 columns for delays 0 to 7
  colnames(data) <- seq(0, 7, length.out = ncol(data))
  # Assume 21 rows for qq from 0 to 1 by 0.05
  qq_values <- seq(0, 1, length.out = nrow(data))
  
  # Add qq column
  data$qq <- qq_values
  
  # Reshape to long format
  data_long <- data %>%
    pivot_longer(cols = -qq, names_to = "delay", values_to = "inc") %>%
    mutate(
      R0 = R0,
      inc_type = output,
      gammaa = gammaa,
      gammap = gammap,
      delay = as.numeric(delay),
      AWT = delay, # Numeric delay for AWT
      AWT2 = paste(delay, "days") # String for AWT2
    )
  
  return(data_long)
}

# Step 5: Process all CSVs
all_data <- file_info %>%
  filter(output %in% c("cinc", "cincres", "clinc")) %>% # Only relevant outputs
  pmap_dfr(~ read_and_reshape(..5, ..1, ..2, ..3, ..4))

# Step 6: Pivot wider to combine outputs
combined_data <- all_data %>%
  pivot_wider(
    id_cols = c(R0, gammaa, gammap, qq, delay, AWT, AWT2),
    names_from = inc_type,
    values_from = inc
  ) %>%
  select(R0, qq, gammaa, gammap, cincres, cinc, clinc, AWT, AWT2)

# Step 7: Write to HeatMap_data.csv
write.csv(combined_data, "HeatMap_data.csv", row.names = FALSE)

# Optional: Print summary
cat("Generated HeatMap_data.csv with", nrow(combined_data), "rows\n")