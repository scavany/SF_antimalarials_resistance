# selective_window_png.R

png("selective_window.png",
    width = 3200,
    height = 2400,
    res = 300,
    bg = "white")

# Plot limits
xlim <- c(0, 10)
ylim <- c(0, 1.05)

# Selective-window positions
x_left  <- 3
x_right <- 7

# Exponential decay curve
A <- 0.95
k <- 0.3
offset <- 0.03

xs <- seq(xlim[1], xlim[2], length.out = 1000)
ys <- A * exp(-k * xs) + offset

# Exact intersection heights
y_at_left  <- A * exp(-k * x_left)  + offset
y_at_right <- A * exp(-k * x_right) + offset

par(mar = c(6, 6, 4, 2) + 0.1)

plot(NA, xlim = xlim, ylim = ylim,
     xaxs = "i", yaxs = "i",
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "",
     bty = "n")

# Axes (black arrows)
arrows(0.5, 0.03, 0.5, 1.02, length = 0.12, lwd = 3)
arrows(0.5, 0.03, 9.8, 0.03, length = 0.12, lwd = 3)

# Solid curve
lines(xs, ys, lwd = 4)

# Dashed vertical lines (stop at curve)
segments(x_left, 0.03, x_left, 1, lty = 2, lwd = 2)
segments(x_right, 0.03, x_right, 1, lty = 2, lwd = 2)

# Dashed horizontal guidelines (meet curve exactly)
segments(0.5, y_at_left, x_left, y_at_left, lty = 2, lwd = 2)
segments(x_right, y_at_right, 0.5, y_at_right, lty = 2, lwd = 2)

# Intersection points
points(c(x_left, x_right),
       c(y_at_left, y_at_right),
       pch = 19, cex = 1.3)

# Selective window arrow
select_y <- min(y_at_left, y_at_right) - 0.08
arrows(x_left + 0.05, select_y,
       x_right - 0.05, select_y,
       code = 3, angle = 20,
       length = 0.12, lwd = 3)

text((x_left + x_right)/2,
     select_y + 0.06,
     "Selective window",
     cex = 1.3)

# MIC labels at x-axis where dashed verticals meet
text(-0.15,
     y_at_left,
     expression(MIC[res]),
     xpd = NA,
     cex = 1.8)

text(-0.15,
     y_at_right,
     expression(MIC[sens]),
     xpd = NA,
     cex = 1.8)

# Region labels
text((0.5 + x_left)/2, 0.92,
     "Prevents both", cex = 1.2)

text((x_left + x_right)/2, 0.92,
     "Prevents sensitive only", cex = 1.2)

text((x_right + 9.8)/2, 0.92,
     "Prevents neither", cex = 1.2)

# Axis labels
mtext(expression(italic("Concentration")),
      side = 2, line = 3.5, cex = 1.3)

mtext(expression(italic("Time since treatment")),
      side = 1, line = 4, cex = 1.3)

# X-axis math labels
text(x = c(x_left, x_right),
     y = par("usr")[3] - 0.07,
     labels = c(expression(frac(1, nu[D[res]])),
                expression(frac(1, nu[D[sens]]))),
     xpd = NA,
     cex = 1.3)

dev.off()