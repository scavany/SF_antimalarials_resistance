library(tidyverse)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(grid)
library(cowplot)   # <<-- add this

# ... (your data loading / preprocessing unchanged) ...

create_heatmap <- function(r0_val, output, data, include_title = TRUE, include_legend = TRUE,
                           show_strip_x = TRUE, show_strip_y = TRUE) {
  plot_data <- filter(data, R0 == r0_val)

  fill_label <- switch(output,
                       "cincres" = "Cumulative resistant cases",
                       "cinc"    = "Cumulative cases",
                       "clinc"   = "Cumulative clinical cases")
  title <- paste("R0 =", r0_val)

  # Build the ggplot heatmap as before but adjust the guide and legend theme
  p <- ggplot(plot_data, aes(x = qq, y = AWT, fill = .data[[output]])) +
    geom_tile() +
    facet_grid(gammaa ~ gammap) +
    scale_fill_distiller(
      palette = "Blues",
      direction = -1,
      name = fill_label,
      # Use a guide_colorbar so we can control the width/height (prevents overlap)
      guide = guide_colorbar(title.position = "top",
                             title.hjust = 0.5,
                             barwidth = unit(10, "cm"),   # increase width to stop overlap
                             barheight = unit(0.6, "cm"))
    ) +
    labs(
      title = if (include_title) title else NULL,
      subtitle = if (include_title) "Partner-drug resistance" else NULL,
      x = "API (% Active Pharmaceutical Ingredient)",
      y = "Delay to Treatment (days)"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      strip.text.x = if (show_strip_x) element_text(face = "bold", size = 10) else element_blank(),
      strip.text.y = if (show_strip_y) element_text(face = "bold", size = 10) else element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, face = "bold", size = 10),
      axis.text.y = element_text(face = "bold", size = 10),
      axis.title = element_text(face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),
      # legend at bottom, horizontal, more spacing and bigger keys
      legend.position = if (include_legend) "bottom" else "none",
      legend.direction = "horizontal",
      legend.text = element_text(face = "bold", size = 10),
      legend.title = element_text(face = "bold", size = 11),
      legend.key.width = unit(2.0, "cm"),
      legend.key.height = unit(0.6, "cm"),
      legend.spacing.x = unit(0.5, "cm"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      # enlarge bottom margin a bit to keep legend from getting clipped
      plot.margin = unit(c(6, 40, 6, 6), "pt")  # top, right, bottom, left (increase right margin for our right-side label)
    )

  return(p)
}

# Save loop: use cowplot to add the right-side vertical label and then save.
for (r0 in r0_values) {
  for (out in outputs) {
    base_plot <- create_heatmap(r0, out, dataMQ, include_title = TRUE, include_legend = TRUE)

    # Compose final figure: draw the base plot and add a right-side label outside plot area.
    final_plot <- ggdraw() +
      draw_plot(base_plot, x = 0, y = 0, width = 1, height = 1) +
      # place the vertical label near the right edge (x close to 1), centered vertically
      draw_label("Artemisinin resistance",
                 x = 0.99, y = 0.5,
                 angle = -90,
                 vjust = 0.5, hjust = 0.5,
                 fontface = "bold",
                 size = 14)

    filename <- paste0("figs/heatmap_", out, "_R0_", r0, ".png")
    # optionally increase width to give more room for the wide legend
    ggsave(filename, final_plot, width = 12, height = 8, dpi = 300)
  }
}
