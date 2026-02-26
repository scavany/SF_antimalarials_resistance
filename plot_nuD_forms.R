nuDmin <- 365/27
nuDmax <- 365/4
nupmax <- 365 / 2
Ds <- 12
dose <- seq(0, 1, 0.01) * Ds
nup <- (dose/Ds)*nupmax            # clearance rate on treatment day 4+ (partner-drug) - SHOULDN'T THIS TAKE A SIMILAR FORM to nuD?

par(mfrow=c(2,1))
plot(dose, nup, type = 'l', lwd = 2)#, ylim = range(nuD1, nuD2, nuD3))
#lines(dose, nuD2, lty = 2, lwd = 2)
#lines(dose, nuD3, col = "red", lwd = 2)
# legend("top", legend = 1:3, lty = c(1,2,1), lwd = 2)
plot(dose, nup, type = 'l', lwd = 2, log = "y")#, ylim = range(nuD1, nuD2, nuD3), log = "y")
# lines(dose, nuD2, lty = 2, lwd = 2)
# lines(dose, nuD3, col = "red", lwd = 2)
# legend("top", legend = 1:3, lty = c(1,2,1), lwd = 2, col = c("black", "black", "red"))
nup[length(nup)]
