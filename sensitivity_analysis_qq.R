# packages
library(lhs)
library(pbapply)
library(deSolve)
library(ppcor)
library(parallel)

rm(list = ls())

# ============================================================
# 1) Parameter ranges — qq is EXCLUDED from LHS
# ============================================================
params <- list(
  wait_treat = c(0, 7),
  R0         = c(2.0, 3.5),
  gammaa     = c(0.25, 0.75),
  gammap     = c(0.25, 0.75),
  nuDmin     = c(365 / 30, 365 / 10),
  nupmax     = c(365 / 30, 365 / 1)
)

k <- length(params)
n <- 500

# qq sweep grid (separate from LHS)
qq_grid <- seq(0.01, 1, length.out = 30)

# 2) Create LHS samples (without qq)
lhs_mat <- randomLHS(n, k)
param_mat <- sapply(1:k, function(i) {
  lo <- params[[i]][1]; hi <- params[[i]][2]
  lhs_mat[, i] * (hi - lo) + lo
})
colnames(param_mat) <- names(params)

# ============================================================
# MedQual ODE function (unchanged except cpmax bugfix)
# ============================================================
MedQual <- function(t, X, parameters) {
  with(as.list(c(X, parameters)), {
    S   <- X[Sindex]
    T1  <- X[T1index]
    T2  <- X[T2index]
    T3  <- X[T3index]
    Tp  <- X[Tpindex]
    Tr  <- X[Trindex]
    IC1 <- X[IC1index]
    IA1 <- X[IA1index]
    IU1 <- X[IU1index]
    P   <- X[Pindex]
    R   <- X[Rindex]
    CumInc         <- X[CumIncindex]
    Fail           <- X[Failindex]
    positiveDay3up <- X[positiveDay3upindex]
    positiveDay1up <- X[positiveDay1upindex]
    CumIC1         <- X[CumIC1index]
    CliFail        <- X[CliFailindex]

    N    <- S + T1 + T2 + T3 + Tp + Tr + IC1 + IA1 + IU1 + P + R
    seas <- 1 + amp * cos(2 * pi * (t - phi))
    nu   <- 1 / ((1 / nuC) + (1 / nuA) + (1 / nuU))
    I    <- T1 + T2 + T3 + Tp + Tr + IC1 + IA1 + IU1
    beta <- R0 * (muo + nu) * seas
    lam  <- beta * (IC1 + kA * (IA1 + T1 + T2 + T3 + Tp) + kU * (IU1 + Tr)) / N
    nuD  <- nuDmin / (dose / parameters$Ds + 0.01 / 365)
    nup  <- (dose / parameters$Ds) * nupmax
    theta <- 1 / (1 + (3 / 365) * nuC)
    treat <- (t >= t_treat) * 365 / (wait_treat)
    tf   <- 7 / 365
    rf   <- 52
    fw   <- ((t > tf) * rf) * f

    r <- c(0.96, 0.92, 0.76, 0.5)
    africanD1 <- sum(exp(-5.25), 3.39, 1.7, 0)
    africanD2 <- sum(exp(-7.35), 2.28, 1.43, 0)
    africanD3 <- sum(exp(-8.63), 1.85, 4.19, 0)
    africanD4 <- sum(exp(-9.00), 1.00, 5.00, 0)
    k.africa  <- c(africanD1, africanD2, africanD3, africanD4)

    p.africa <- matrix(0, nrow = length(dose), ncol = 4)
    c.africa <- matrix(0, nrow = length(dose), ncol = 4)

    for (i in 1:4) {
      p.africa[, i] <- (k.africa[i] * r[i]^dose) / (1 + (k.africa[i] * r[i]^dose))
    }

    c.africa[, 1] <- 1 - p.africa[, 1]
    c.africa[, 2] <- 1 - p.africa[, 2] / p.africa[, 1]
    c.africa[, 3] <- 1 - p.africa[, 3] / p.africa[, 2]
    c.africa[, 4] <- 1 - p.africa[, 4] / p.africa[, 3]

    c1max[1] <- c.africa[, 1]
    c2max[1] <- c.africa[, 2]
    c3max[1] <- c.africa[, 3]
    # BUGFIX: was p.africa[,4]/p.africa[,4] (always 0), should be day4/day3
    cpmax[1] <- c.africa[, 4]

    prec <- (1 - dose / Ds) * precmax + dose / Ds * precmin

    c1 <- dose / parameters$Ds * c1max
    c2 <- dose / parameters$Ds * c2max
    c3 <- dose / parameters$Ds * c3max
    cp <- dose / parameters$Ds * cpmax
    c1[2] <- (1 - gammaa) * c1[1]
    c2[2] <- (1 - gammaa) * c2[1]
    c3[2] <- (1 - gammaa) * c3[1]
    cp[2] <- (1 - gammap) * cp[1]

    # ---- INFLUENCES ----
    # Part A: cross-immunity
    ps[1] <- ps[1] * (1 - R[2] / N[2]) + pr * (R[2] / N[2])
    ps[2] <- ps[2] * (1 - R[1] / N[1]) + pr * (R[1] / N[1])

    # Part B+D
    lam[2] <- lam[2] * (1 - (1 / N[2]) * ((nuD[1] / (365 + nuD[1])) * T1[1] +
      (2 * nuD[1] / (365 + 2 * nuD[1])) * T2[1] +
      (3 * nuD[1] / (365 + 3 * nuD[1])) * T3[1] +
      IC1[1] + IA1[1] + IU1[1] + Tr[1] + (nuD[1] / nuD[2] * P[1])))
    kRes2Sens[1] <- 0
    kRes2Sens[2] <- ((S[1] + R[1]) / N[1]) * lam[1]

    # Part C+D
    lam[1] <- lam[1] * (1 - (1 / N[1]) * ((365 / (365 + nuD[1])) * T1[2] +
      (365 / (365 + 2 * nuD[1])) * T2[2] +
      (365 / (365 + 3 * nuD[1])) * T3[2] + Tp[2] + P[2]))
    kSens2Res[1] <- ((S[2] + R[2]) / N[2]) * lam[2]
    kSens2Res[2] <- 0

    # Part E: novel mutations
    jSens2Res1[1] <- -365 * (365 / (365 + nuD[1])) * m1 * T1[1] / N[2] * (S[2] + R[2])
    jSens2Res1[2] <-  365 * (365 / (365 + nuD[1])) * m1 * T1[1] / N[2] * (S[2] + R[2])
    jSens2Res2[1] <- -365 * (365 / (365 + 2 * nuD[1])) * m2 * T2[1] / N[2] * (S[2] + R[2])
    jSens2Res2[2] <-  365 * (365 / (365 + 2 * nuD[1])) * m2 * T2[1] / N[2] * (S[2] + R[2])
    jSens2Res3[1] <- -365 * (365 / (365 + 3 * nuD[1])) * m3 * T3[1] / N[2] * (S[2] + R[2])
    jSens2Res3[2] <-  365 * (365 / (365 + 3 * nuD[1])) * m3 * T3[1] / N[2] * (S[2] + R[2])
    jSens2Resp[1] <- -365 * mp * Tp[1] / N[2] * (S[2] + R[2])
    jSens2Resp[2] <-  365 * mp * Tp[1] / N[2] * (S[2] + R[2])
    jSens2ResR[1] <- 0
    jSens2ResR[2] <- -365 * (365 / (365 + nuD[1])) * m1 * T1[1] / N[2] -
      365 * (365 / (365 + 2 * nuD[1])) * m2 * T2[1] / N[2] -
      365 * (365 / (365 + 3 * nuD[1])) * m3 * T3[1] / N[2] -
      365 * mp * Tp[1] / N[2]
    jSens2ResS[1] <- 0
    jSens2ResS[2] <- -365 * (365 / (365 + nuD[1])) * m1 * T1[1] / N[2] -
      365 * (365 / (365 + 2 * nuD[1])) * m2 * T2[1] / N[2] -
      365 * (365 / (365 + 3 * nuD[1])) * m3 * T3[1] / N[2] -
      365 * mp * Tp[1] / N[2]
    jSens2ResP[1] <- 365 * (365 / (365 + nuD[1])) * m1 * T1[1] / N[2] * (S[2] + R[2]) +
      365 * (365 / (365 + 2 * nuD[1])) * m2 * T2[1] / N[2] * (S[2] + R[2]) +
      365 * (365 / (365 + 3 * nuD[1])) * m3 * T3[1] / N[2] * (S[2] + R[2]) +
      365 * mp * Tp[1] / N[2] * (S[2] + R[2])
    jSens2ResP[2] <- 0

    # ---- ODEs ----
    dS   <- mui * N + omega * R - lam * S - muo * S + jSens2ResS * S
    dT1  <- treat * IC1 - nu1 * T1 - (nuD[1] / (365 + nuD[1])) * kRes2Sens * T1 -
      (365 / (365 + nuD[1])) * kSens2Res * T1 - muo * T1 + jSens2Res1
    dT2  <- (1 - c1) * nu1 * T1 - nu2 * T2 - (2 * nuD[1] / (365 + 2 * nuD[1])) * kRes2Sens * T2 -
      (365 / (365 + 2 * nuD[1])) * kSens2Res * T2 - muo * T2 + jSens2Res2
    dT3  <- (1 - c2) * nu2 * T2 - nu3 * T3 - (3 * nuD[1] / (365 + 3 * nuD[1])) * kRes2Sens * T3 -
      (365 / (365 + 3 * nuD[1])) * kSens2Res * T3 - muo * T3 + jSens2Res3
    dTp  <- (1 - theta) * (1 - c3) * nu3 * T3 - (nuD + nup) * Tp - kSens2Res * Tp - muo * Tp + jSens2Resp
    dTr  <- prec * c1 * nu1 * T1 + prec * c2 * nu2 * T2 + prec * c3 * nu3 * T3 +
      prec * nuD * Tp - rho * Tr - kRes2Sens * Tr - muo * Tr
    dIC1 <- ps * lam * S + pr * lam * R + theta * (1 - c3) * nu3 * T3 + rho * Tr -
      treat * IC1 - nuC * IC1 - kRes2Sens * IC1 - muo * IC1
    dIA1 <- lam * (1 - ps) * S + lam * (1 - pr) * R + (1 - prec) * nuD * Tp +
      nuC * IC1 - nuA * IA1 - kRes2Sens * IA1 - muo * IA1
    dIU1 <- nuA * IA1 - nuU * IU1 - kRes2Sens * IU1 - muo * IU1
    dP   <- (1 - prec) * c1 * nu1 * T1 + (1 - prec) * c2 * nu2 * T2 +
      (1 - prec) * c3 * nu3 * T3 + nup * Tp - nuD * P - muo * P + jSens2ResP
    dR   <- nuU * IU1 + nuD * P - omega * R - lam * R +
      kRes2Sens * (IC1 + IA1 + IU1 + Tr) +
      ((nuD[1] / (365 + nuD[1])) * kRes2Sens * T1 +
       (2 * nuD[1] / (365 + 2 * nuD[1])) * kRes2Sens * T2 +
       (3 * nuD[1] / (365 + 3 * nuD[1])) * kRes2Sens * T3) +
      kSens2Res * Tp +
      (365 / (365 + nuD[1])) * kSens2Res * T1 +
      (365 / (365 + 2 * nuD[1])) * kSens2Res * T2 +
      (365 / (365 + 3 * nuD[1])) * kSens2Res * T3 -
      muo * R + jSens2ResR * R

    dCumInc         <- treat * IC1
    dFail           <- fw * Tp + sensC * fw * IC1 + sensA * fw * IA1 + sensU * fw * IU1
    dpositiveDay3up <- T3 + Tp
    dpositiveDay1up <- T1 + T2 + T3 + Tp
    dCumIC1         <- ps * lam * S + pr * lam * R + theta * (1 - c3) * nu3 * T3 + rho * Tr
    dCliFail        <- theta * (1 - c3) * nu3 * T3 + rho * Tr

    list(c(dS, dT1, dT2, dT3, dTp, dTr, dIC1, dIA1, dIU1,
           dP, dR, dCumInc, dFail, dpositiveDay3up, dpositiveDay1up, dCumIC1, dCliFail))
  })
}

# ============================================================
# run_model: unchanged from original except bugfix applied above
# ============================================================
run_model <- function(pars) {
  A <- 2
  maxt      <- 35
  starttime <- -10
  dt        <- 1 / 12
  times     <- seq(starttime, maxt, by = dt)

  Sindex             <- 1:A
  T1index            <- (A + 1):(2 * A)
  T2index            <- (2 * A + 1):(3 * A)
  T3index            <- (3 * A + 1):(4 * A)
  Tpindex            <- (4 * A + 1):(5 * A)
  Trindex            <- (5 * A + 1):(6 * A)
  IC1index           <- (6 * A + 1):(7 * A)
  IA1index           <- (7 * A + 1):(8 * A)
  IU1index           <- (8 * A + 1):(9 * A)
  Pindex             <- (9 * A + 1):(10 * A)
  Rindex             <- (10 * A + 1):(11 * A)
  CumIncindex        <- (11 * A + 1):(12 * A)
  Failindex          <- (12 * A + 1):(13 * A)
  positiveDay3upindex <- (13 * A + 1):(14 * A)
  positiveDay1upindex <- (14 * A + 1):(15 * A)
  CumIC1index        <- (15 * A + 1):(16 * A)
  CliFailindex       <- (16 * A + 1):(17 * A)

  parameters <- list(
    N = 1000, mui = 1 / 50, muo = 1 / 50,
    R0 = rep(pars$R0, A),
    ps = c(0.9, 0.9), pr = 0.1,
    wait_treat = pars$wait_treat,
    omega = 1 / 2, nuC = 365 / 10, nuA = 365 / 60, nuU = 365 / 120,
    t_treat = 5,
    c1max = c(0.99, 0.00), c2max = c(0.99, 0.00),
    c3max = c(0.99, 0.00), cpmax = c(0.99, 0.60),
    nupmax = pars$nupmax,
    rho = 365 / 20,
    nu1 = 365 / 1, nu2 = 365 / 1, nu3 = 365 / 1,
    precmax = 0.5, precmin = 0.05,
    thetamax = 0.77,
    nuDmin = c(pars$nuDmin, pars$nuDmin * 27 / 10),
    amp = 0.7, phi = 0,
    sensC = 0.95, sensA = 0.50, sensU = 0.00,
    kA = 0.7, kU = 0.3,
    m1 = 1 / 10^6, m2 = 1 / 10^6, m3 = 1 / 10^6, mp = 1 / 10^6,
    pn = 1, Ds = 12, f = 0, dose = 0,
    kRes2Sens = rep(0, A), kSens2Res = rep(0, A),
    jSens2Res1 = rep(0, A), jSens2Res2 = rep(0, A),
    jSens2Res3 = rep(0, A), jSens2Resp = rep(0, A),
    jSens2ResP = rep(0, A), jSens2ResR = rep(0, A), jSens2ResS = rep(0, A),
    lam = rep(0, A),
    gammaa = pars$gammaa, gammap = pars$gammap,
    Sindex = Sindex, T1index = T1index, T2index = T2index, T3index = T3index,
    Tpindex = Tpindex, Trindex = Trindex, IC1index = IC1index,
    IA1index = IA1index, IU1index = IU1index, Pindex = Pindex,
    Rindex = Rindex, CumIncindex = CumIncindex, Failindex = Failindex,
    positiveDay3upindex = positiveDay3upindex,
    positiveDay1upindex = positiveDay1upindex,
    CumIC1index = CumIC1index, CliFailindex = CliFailindex
  )

  initS   <- rep(0.9 * parameters$N, A)
  initIA1 <- rep(0.05 * parameters$N, A)
  initR   <- rep(parameters$N, A) - initS - initIA1
  zeros   <- rep(0, A)
  X <- c(initS, zeros, zeros, zeros, zeros, zeros, zeros, initIA1, zeros,
         zeros, initR, zeros, zeros, zeros, zeros, zeros, zeros)

  parameters$dose <- pars$qq * parameters$Ds

  out <- ode(y = X, times = times, func = MedQual, parms = parameters,
             method = "vode", maxsteps = 5000)

  # Post-process
  maxt_idx <- nrow(out)
  CumInc  <- out[maxt_idx, 1 + CumIncindex]
  CumIC1  <- out[maxt_idx, 1 + CumIC1index]
  CliFail <- out[maxt_idx, 1 + CliFailindex]

  cinc     <- sum(CumInc) - prod(CumInc) / ((maxt - parameters$t_treat) * 1000)
  cincres  <- as.numeric(CumInc[2])
  clinc    <- sum(CumIC1) - prod(CumIC1) / ((maxt - parameters$t_treat) * 1000)
  clincres <- as.numeric(CumIC1[2])
  cliFail  <- sum(CliFail) - prod(CliFail) / ((maxt - parameters$t_treat) * 1000)

  return(c(cinc = cinc, cincres = cincres, clinc = clinc,
           clincres = clincres, cliFail = cliFail))
}

# ============================================================
# 3) For each LHS sample, sweep qq and find max resistance
# ============================================================
sweep_qq <- function(pars_vec) {
  pars <- as.list(pars_vec)

  best_clincres <- -Inf
  best_cincres  <- -Inf
  best_qq_clinc <- NA
  best_qq_cinc  <- NA
  best_all      <- NULL   # full output at max clincres

  for (qq in qq_grid) {
    pars$qq <- qq
    res <- tryCatch(run_model(pars), error = function(e) rep(NA, 5))
    if (any(is.na(res))) next

    if (res["clincres"] > best_clincres) {
      best_clincres  <- res["clincres"]
      best_qq_clinc  <- qq
      best_all       <- res
    }
    if (res["cincres"] > best_cincres) {
      best_cincres <- res["cincres"]
      best_qq_cinc <- qq
    }
  }

  if (is.null(best_all)) return(rep(NA, 7))

  return(c(best_all,
           qq_max_clincres = best_qq_clinc,
           qq_max_cincres  = best_qq_cinc))
}

cat("Running LHS samples with qq sweep (n =", n, ", qq grid =", length(qq_grid), "points)...\n")

# ============================================================
# 4) PRCC on remaining parameters vs worst-case resistance
# ============================================================
cl <- makeCluster(detectCores() - 1)

# Export everything the workers need
clusterExport(cl, c("sweep_qq", "run_model", "MedQual",
                    "qq_grid"))  # add any globals used inside

clusterEvalQ(cl, library(deSolve))  # load packages on workers

results <- pbapply(param_mat, 1, sweep_qq, cl = cl)
stopCluster(cl)
results_df <- as.data.frame(t(results))

combined <- cbind(as.data.frame(param_mat), results_df)
combined[] <- lapply(combined, as.numeric)

# Drop rows where sweep failed
combined <- combined[complete.cases(combined), ]
cat("Successful samples:", nrow(combined), "/", n, "\n")

param_names   <- colnames(param_mat)
# Outcomes of interest: max resistance burden + the qq that produced it
outcome_names <- c("clincres", "cincres", "qq_max_clincres", "qq_max_cincres",
                   "cinc", "clinc", "cliFail")
k <- length(param_names)

prcc_est <- matrix(NA, nrow = k, ncol = length(outcome_names),
                   dimnames = list(param_names, outcome_names))
prcc_p   <- matrix(NA, nrow = k, ncol = length(outcome_names),
                   dimnames = list(param_names, outcome_names))

for (oj in seq_along(outcome_names)) {
  oname <- outcome_names[oj]
  if (!(oname %in% colnames(combined))) next
  y <- combined[[oname]]
  for (i in seq_len(k)) {
    x <- combined[[param_names[i]]]
    other_idx <- setdiff(seq_len(k), i)
    Z <- as.matrix(combined[, param_names[other_idx], drop = FALSE])
    test <- tryCatch(
      pcor.test(x, y, Z, method = "spearman"),
      error = function(e) NULL
    )
    if (!is.null(test)) {
      prcc_est[i, oj] <- test$estimate
      prcc_p[i, oj]   <- test$p.value
    }
  }
}

prcc_est_df <- as.data.frame(prcc_est)
prcc_p_df   <- as.data.frame(prcc_p)

cat("\n=== PRCC for max clinical resistance incidence (clincres) ===\n")
print(data.frame(PRCC = prcc_est_df$clincres, p = prcc_p_df$clincres))

cat("\n=== PRCC for max cumulative resistance incidence (cincres) ===\n")
print(data.frame(PRCC = prcc_est_df$cincres, p = prcc_p_df$cincres))

cat("\n=== PRCC for qq at max clincres (which dose fraction maximises resistance?) ===\n")
print(data.frame(PRCC = prcc_est_df$qq_max_clincres, p = prcc_p_df$qq_max_clincres))

save(prcc_est_df, prcc_p_df, combined, file = "./prcc_qq_sweep_results.RData")