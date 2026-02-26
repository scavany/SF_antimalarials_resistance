library(tidyverse)

# Load data (replace with your data frame loading method if needed)
dataMQ <- read_csv("HeatMap_data.csv") # Uncomment if loading from CSV
# Assume dataMQ has columns: R0, qq, cincres, cinc, clinc, AWT, AWT2, gammaa, gammap

# Define unique R0 values and outputs
r0_values <- unique(dataMQ$R0)
outputs <- c("cincres", "cinc", "clinc")

# Function to create a 3x3 contour plot for a given R0 and output
create_contour_plot <- function(r0_val, output, data) {
  # Filter data for specific R0
  plot_data <- filter(data, R0 == r0_val)
  
  # Set plot title and contour variable
  title <- paste("Filled Contours for", output, "at R0 =", r0_val)
  fill_label <- switch(output,
                       "cincres" = "Cumulative Resistance Cases",
                       "cinc" = "Incidence Cases",
                       "clinc" = "Infected Clinical (per 1000)")
  
  # Create plot
  p <- ggplot(plot_data, aes(x = qq, y = AWT, z = .data[[output]])) +
    geom_contour_filled(bins = 10) + # Adjust bins for desired contour levels
    facet_grid(gammaa ~ gammap, labeller = label_both) +
    scale_fill_brewer(name = fill_label, palette = "YlOrRd") + # Color scheme
    labs(
      title = title,
      x = "API (% Active Pharmaceutical Ingredient)",
      y = "Delay to Treatment (days)"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 8, face = "bold"),
      axis.text = element_text(size = 7, face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 8, face = "bold")
    )
  
  # Save plot
  filename <- paste0("contour_", output, "_R0_", r0_val, ".png")
  ggsave(filename, p, width = 10, height = 8, dpi = 300)
  
  return(p)
}

# Generate all plots
for (r0 in r0_values) {
  for (out in outputs) {
    create_contour_plot(r0, out, dataMQ)
  }
}

# Print confirmation
cat("Generated", length(r0_values) * length(outputs), "contour plots as PNG files.\n")
