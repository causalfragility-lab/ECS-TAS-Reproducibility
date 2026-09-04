# Evidence-Calibration-Stability (ECS)
# Purpose:
#   1. Reproduce the numerical results reported in the manuscript and supplement.
#   2. Print the principal results directly to the R console.
#   3. Display the four manuscript figures in the R graphics device.
#
# Required external file:
#   ECS_diabetes_data.csv

rm(list = ls())
graphics.off()
options(width = 160)
set.seed(20260826)

print("ECS analysis started")
flush.console()

ALPHA       <- 0.05
POWER_REQ   <- 0.80
DELTA_D     <- 0.30
DELTA_M     <- 0.30
B0_DEFAULT  <- 0.10
R0_DEFAULT  <- 0.10
S_DEFAULT   <- 1.00
TOL         <- 1e-10

SAVE_FIGURES <- FALSE

# Optional Monte Carlo stress tests are computationally slower.
RUN_STRESS <- FALSE
B_STRESS   <- 30000


# -----------------------------------------------------------------------------
# File location
# -----------------------------------------------------------------------------

get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)

  if (length(file_arg) > 0L) {
    p <- sub("^--file=", "", file_arg[1L])
    return(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
  }

  frames <- sys.frames()
  if (length(frames)) {
    for (i in rev(seq_along(frames))) {
      ofile <- frames[[i]]$ofile
      if (!is.null(ofile) && nzchar(ofile)) {
        return(dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
      }
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

SCRIPT_DIR <- get_script_dir()

resolve_diabetes_file <- function() {
  env_path <- Sys.getenv("ECS_DIABETES_DATA", unset = "")

  candidates <- c(
    if (nzchar(env_path)) env_path else character(0),
    file.path(SCRIPT_DIR, "ECS_diabetes_data.csv"),
    file.path(getwd(), "ECS_diabetes_data.csv"),
    path.expand("~/Downloads/ECS_diabetes_data.csv")
  )

  candidates <- unique(candidates)
  hit <- candidates[file.exists(candidates)]

  if (length(hit)) {
    return(normalizePath(hit[1L], winslash = "/", mustWork = TRUE))
  }

  if (interactive()) {
    message("Select ECS_diabetes_data.csv")
    return(normalizePath(file.choose(), winslash = "/", mustWork = TRUE))
  }

  stop(
    "ECS_diabetes_data.csv was not found. Place it beside this R script or in the working directory.",
    call. = FALSE
  )
}


# -----------------------------------------------------------------------------
# functions
# -----------------------------------------------------------------------------

exact_t_power <- function(n, d, alpha = ALPHA) {
  df <- n - 1
  crit <- qt(1 - alpha, df = df)
  1 - pt(crit, df = df, ncp = d * sqrt(n))
}

calibration_row <- function(n,
                            d_star = DELTA_D,
                            test_alpha = ALPHA,
                            alpha_target = ALPHA,
                            power_req = POWER_REQ) {
  pow <- exact_t_power(n, d_star, test_alpha)
  size_actual <- test_alpha
  size_margin <- alpha_target - size_actual
  power_margin <- pow - power_req
  pass <- (size_actual <= alpha_target + TOL) &&
          (pow >= power_req - TOL)

  data.frame(
    n = n,
    size_actual = size_actual,
    size_margin = size_margin,
    power_at_Delta_D = pow,
    power_margin = power_margin,
    C = min(size_margin, power_margin),
    calibration_pass = pass
  )
}

ecs_one_sample <- function(x,
                           alpha = ALPHA,
                           b0 = B0_DEFAULT,
                           r0 = R0_DEFAULT) {
  n <- length(x)
  xb <- mean(x)
  sx <- sd(x)
  se <- sx / sqrt(n)
  crit <- qt(1 - alpha, df = n - 1)
  t_obs <- xb / se
  E <- t_obs - crit
  m0 <- xb - crit * se
  a <- c(b0, crit * se * r0)
  S <- if (m0 > 0) m0 / sqrt(sum(a^2)) else 0

  list(
    n = n,
    mean = xb,
    sd = sx,
    se = se,
    critical = crit,
    t = t_obs,
    p_one_sided = 1 - pt(t_obs, n - 1),
    E = E,
    mean_scale_margin = m0,
    stability_S = S,
    perturbation_loading = a
  )
}

ellipsoidal_stability <- function(m0, a, W = diag(length(a))) {
  if (m0 <= 0) return(0)
  Winv_a <- solve(W, a)
  m0 / sqrt(drop(crossprod(a, Winv_a)))
}

closest_reversal <- function(m0, a, W = diag(length(a))) {
  Winv_a <- solve(W, a)
  drop(m0 / drop(crossprod(a, Winv_a))) * Winv_a
}

exact_stable_probability <- function(n,
                                     mu,
                                     b0 = B0_DEFAULT,
                                     r0 = R0_DEFAULT,
                                     s_star = S_DEFAULT,
                                     alpha = ALPHA) {
  df <- n - 1
  crit <- qt(1 - alpha, df)
  sig_xbar <- 1 / sqrt(n)

  integrand <- function(v) {
    s <- sqrt(v / df)
    se <- s / sqrt(n)
    threshold <- crit * se +
      s_star * sqrt(b0^2 + (crit * se * r0)^2)

    p <- pnorm(
      (threshold - mu) / sig_xbar,
      lower.tail = FALSE
    )

    p * dchisq(v, df)
  }

  lo <- qchisq(1e-10, df)
  hi <- qchisq(1 - 1e-10, df)

  integrate(
    integrand,
    lo,
    hi,
    rel.tol = 2e-9,
    subdivisions = 1000L
  )$value
}

exact_ES_correlation <- function(n,
                                 mu,
                                 b0 = B0_DEFAULT,
                                 r0 = R0_DEFAULT,
                                 alpha = ALPHA) {
  df <- n - 1
  crit <- qt(1 - alpha, df)
  sig <- 1 / sqrt(n)

  pieces <- function(v) {
    s <- sqrt(v / df)
    if (s <= 0) return(rep(0, 6))

    q <- crit * s / sqrt(n)
    z <- (q - mu) / sig
    P <- pnorm(z, lower.tail = FALSE)
    ph <- dnorm(z)
    M1 <- mu * P + sig * ph
    M2 <- (mu^2 + sig^2) * P + sig * (mu + q) * ph

    AE <- sqrt(n) / s
    BE <- -crit

    den <- sqrt(b0^2 + (crit * s / sqrt(n) * r0)^2)
    AS <- 1 / den
    BS <- -q / den

    EI <- AE * M1 + BE * P
    SI <- AS * M1 + BS * P
    E2I <- AE^2 * M2 + 2 * AE * BE * M1 + BE^2 * P
    S2I <- AS^2 * M2 + 2 * AS * BS * M1 + BS^2 * P
    ESI <- AE * AS * M2 +
      (AE * BS + AS * BE) * M1 +
      BE * BS * P

    c(P, EI, SI, E2I, S2I, ESI) * dchisq(v, df)
  }

  vals <- sapply(1:6, function(j) {
    integrate(
      function(v) {
        vapply(v, function(x) pieces(x)[j], numeric(1))
      },
      0,
      Inf,
      rel.tol = 5e-8,
      subdivisions = 1000L
    )$value
  })

  P <- vals[1]
  mE <- vals[2] / P
  mS <- vals[3] / P
  vE <- vals[4] / P - mE^2
  vS <- vals[5] / P - mS^2
  covES <- vals[6] / P - mE * mS

  covES / sqrt(vE * vS)
}


# -----------------------------------------------------------------------------
# 1. Calibration
# -----------------------------------------------------------------------------

calibration <- do.call(
  rbind,
  lapply(c(40, 71, 100), calibration_row)
)

print("CALIBRATION")
print(calibration, row.names = FALSE, digits = 7)
flush.console()


# -----------------------------------------------------------------------------
# 2. Table 1 and E-S dependence
# -----------------------------------------------------------------------------

n_grid <- c(40, 71, 100)
mu_grid <- c(0.00, 0.15, 0.30, 0.50)

rows <- list()
k <- 1L

for (n in n_grid) {
  pow <- exact_t_power(n, DELTA_D)
  calpass <- pow >= POWER_REQ - TOL
  crit <- qt(1 - ALPHA, n - 1)

  for (mu in mu_grid) {
    reject <- 1 - pt(
      crit,
      n - 1,
      ncp = mu * sqrt(n)
    )

    stable <- exact_stable_probability(n, mu)

    rows[[k]] <- data.frame(
      n = n,
      true_mean = mu,
      power_at_030 = pow,
      conventional_rejection = reject,
      stable_among_rejections = stable / reject,
      ECS_criterion_rate = if (calpass) stable else 0,
      corr_E_S_given_rejection = exact_ES_correlation(n, mu)
    )

    k <- k + 1L
  }
}

Table1 <- do.call(rbind, rows)

print("TABLE 1")
print(Table1, row.names = FALSE, digits = 7)
write.csv(
  Table1,
  "Table1_exact_normal_R.csv",
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 3. Threshold sensitivity
# -----------------------------------------------------------------------------

b0_grid <- c(.05, .10, .15, .20)
s_grid <- c(.5, 1, 1.5, 2)

threshold_sensitivity <- expand.grid(
  b0 = b0_grid,
  s_star = s_grid
)

threshold_sensitivity$ecs_criterion_probability <- mapply(
  function(b, ss) {
    exact_stable_probability(
      71,
      .30,
      b0 = b,
      r0 = .10,
      s_star = ss
    )
  },
  threshold_sensitivity$b0,
  threshold_sensitivity$s_star
)

print("THRESHOLD SENSITIVITY")
print(threshold_sensitivity, row.names = FALSE, digits = 7)
flush.console()

write.csv(
  threshold_sensitivity,
  "threshold_sensitivity_R.csv",
  row.names = FALSE
)

r0_grid <- c(0, .05, .10, .20, .30)

r0_sensitivity <- expand.grid(
  b0 = b0_grid,
  r0 = r0_grid
)

r0_sensitivity$ecs_criterion_probability <- mapply(
  function(b, rr) {
    exact_stable_probability(
      71,
      .30,
      b0 = b,
      r0 = rr,
      s_star = 1
    )
  },
  r0_sensitivity$b0,
  r0_sensitivity$r0
)

write.csv(
  r0_sensitivity,
  "r0_sensitivity_R.csv",
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 4. Claim-relative calculation
# -----------------------------------------------------------------------------

claim_n <- c(40, 100, 500, 2000)
claim_mu <- .20

claim_rows <- lapply(claim_n, function(n) {
  df <- n - 1
  crit <- qt(1 - ALPHA, df)
  power_n <- exact_t_power(n, DELTA_D)
  calibrated <- power_n >= POWER_REQ - TOL

  reject0 <- 1 - pt(
    crit,
    df,
    ncp = claim_mu * sqrt(n)
  )

  criterion0 <- if (calibrated) {
    exact_stable_probability(n, claim_mu)
  } else {
    0
  }

  rejectM <- 1 - pt(
    crit,
    df,
    ncp = (claim_mu - DELTA_M) * sqrt(n)
  )

  data.frame(
    n = n,
    reject_H0_mu_le_0 = reject0,
    ECS_criterion_positivity = criterion0,
    reject_H0_mu_le_030 = rejectM
  )
})

claim_results <- do.call(rbind, claim_rows)

print("CLAIM-RELATIVE RESULTS")
print(claim_results, row.names = FALSE, digits = 10)

write.csv(
  claim_results,
  "claim_relative_exact_R.csv",
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 5. Diabetes specification analysis
# -----------------------------------------------------------------------------

DIABETES_DATA_FILE <- resolve_diabetes_file()
diabetes <- read.csv(DIABETES_DATA_FILE)

base_vars <- c("age", "sex", "bmi")
optional <- c("bp", "s1", "s2", "s3", "s4", "s5", "s6")

fit_spec <- function(add_vars) {
  rhs <- c(base_vars, add_vars)
  f <- as.formula(
    paste("target ~", paste(rhs, collapse = " + "))
  )

  fit <- lm(f, data = diabetes)
  sm <- summary(fit)$coefficients

  b <- unname(sm["age", "Estimate"])
  se <- unname(sm["age", "Std. Error"])
  tt <- unname(sm["age", "t value"])
  p1 <- pt(
    tt,
    df = df.residual(fit),
    lower.tail = FALSE
  )

  data.frame(
    k = length(add_vars),
    controls = if (length(add_vars)) {
      paste(add_vars, collapse = "+")
    } else {
      "(none)"
    },
    coef = b,
    se = se,
    t = tt,
    p_one_sided = p1,
    support = (tt > 0 && p1 < ALPHA),
    stringsAsFactors = FALSE
  )
}

specs <- list()
ii <- 1L

for (mask in 0:(2^length(optional) - 1)) {
  keep <- as.logical(
    intToBits(mask)[seq_along(optional)]
  )

  specs[[ii]] <- fit_spec(optional[keep])
  ii <- ii + 1L
}

spec_results <- do.call(rbind, specs)
spec_results <- spec_results[
  order(spec_results$k, spec_results$controls),
]

benchmark <- spec_results[
  spec_results$k == 0,
][1, ]

min_loss <- min(
  spec_results$k[!spec_results$support]
)

min_sign_reversal <- min(
  spec_results$k[spec_results$coef < 0]
)

print("DIABETES: 128 SPECIFICATIONS")
print(spec_results, row.names = FALSE, digits = 8)
flush.console()

print("DIABETES SUMMARY")
print(
  data.frame(
    benchmark_age_coef = benchmark$coef,
    benchmark_se = benchmark$se,
    benchmark_t = benchmark$t,
    benchmark_p_one_sided = benchmark$p_one_sided,
    supporting_specs = sum(spec_results$support),
    total_specs = nrow(spec_results),
    nearest_loss = min_loss,
    nearest_sign_reversal = min_sign_reversal
  ),
  row.names = FALSE,
  digits = 8
)

write.csv(
  spec_results,
  "diabetes_specifications_R.csv",
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 6. Student sleep data
# -----------------------------------------------------------------------------

data(sleep)

x1 <- sleep$extra[sleep$group == 1]
x2 <- sleep$extra[sleep$group == 2]
d_sleep <- x2 - x1

sleep025 <- ecs_one_sample(
  d_sleep,
  b0 = .25,
  r0 = .10
)

sleep050 <- ecs_one_sample(
  d_sleep,
  b0 = .50,
  r0 = .10
)

sleep_m0 <- sleep025$mean_scale_margin
sleep_q <- sleep025$critical * sleep025$se * .10

b0_cross <- sqrt(
  sleep_m0^2 - sleep_q^2
)

additive_only_reversal <- sleep_m0
reversal_after_10pct_se <- sleep_m0 - sleep_q

sleep_results <- data.frame(
  t = sleep025$t,
  p_one_sided = sleep025$p_one_sided,
  E = sleep025$E,
  S_b0_025 = sleep025$stability_S,
  S_b0_050 = sleep050$stability_S,
  b0_crossing = b0_cross,
  additive_only_reversal = additive_only_reversal,
  reversal_after_10pct_se = reversal_after_10pct_se
)

print("STUDENT SLEEP RESULTS")
print(sleep_results, row.names = FALSE, digits = 8)
flush.console()

sleep_profile <- data.frame(
  b0 = seq(.05, 1.10, length.out = 251)
)

sleep_profile$S <- sleep_m0 /
  sqrt(sleep_profile$b0^2 + sleep_q^2)

write.csv(
  sleep_profile,
  "sleep_stability_profile_R.csv",
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 7. Smooth nonlinear lower bound
# -----------------------------------------------------------------------------

nonlinear_stability_lower_bound <- function(m0,
                                            g,
                                            W = diag(length(g)),
                                            K = 0,
                                            R = Inf) {
  if (m0 <= 0 || R == 0) return(0)

  G <- sqrt(
    drop(
      crossprod(
        g,
        solve(W, g)
      )
    )
  )

  if (G <= 0) {
    bound <- if (K > 0) sqrt(2 * m0 / K) else Inf
  } else if (K == 0) {
    bound <- m0 / G
  } else {
    bound <- (
      sqrt(G^2 + 2 * K * m0) - G
    ) / K
  }

  min(R, bound)
}


# -----------------------------------------------------------------------------
# 8. out-of-class stress tests
# -----------------------------------------------------------------------------

r_standardized_lognormal <- function(n, log_sd = .75) {
  raw <- rlnorm(n, 0, log_sd)
  m <- exp(log_sd^2 / 2)
  sd0 <- sqrt(
    (exp(log_sd^2) - 1) *
      exp(log_sd^2)
  )

  (raw - m) / sd0
}

r_error <- function(n, dgp) {
  switch(
    dgp,
    Normal = rnorm(n),
    t3 = rt(n, 3) / sqrt(3),
    Lognormal = r_standardized_lognormal(n),
    stop("Unknown data-generating process.")
  )
}

stress_cell <- function(n,
                        mu,
                        dgp,
                        B = B_STRESS) {
  pow <- exact_t_power(n, DELTA_D)
  cal <- pow >= POWER_REQ - TOL

  rej <- logical(B)
  stable <- logical(B)

  for (i in seq_len(B)) {
    z <- ecs_one_sample(
      mu + r_error(n, dgp)
    )

    rej[i] <- z$E > 0
    stable[i] <- z$stability_S > S_DEFAULT
  }

  data.frame(
    dgp = dgp,
    n = n,
    mu = mu,
    rejection_rate = mean(rej),
    stable_given_rejection = if (any(rej)) {
      mean(stable[rej])
    } else {
      NA_real_
    },
    ecs_criterion_rate = mean(
      rej & stable & cal
    )
  )
}

if (RUN_STRESS) {
  set.seed(20260826)

  stress <- do.call(
    rbind,
    lapply(
      c("Normal", "t3", "Lognormal"),
      function(dgp) {
        do.call(
          rbind,
          lapply(
            n_grid,
            function(n) {
              do.call(
                rbind,
                lapply(
                  mu_grid,
                  function(mu) {
                    stress_cell(
                      n,
                      mu,
                      dgp
                    )
                  }
                )
              )
            }
          )
        )
      }
    )
  )

  print("OPTIONAL OUT-OF-CLASS STRESS TESTS")
  print(stress, row.names = FALSE, digits = 7)

  write.csv(
    stress,
    "supplement_out_of_class_stress_R.csv",
    row.names = FALSE
  )
}


# -----------------------------------------------------------------------------
# Figure helpers
# -----------------------------------------------------------------------------


save_png <- function(filename, plot_fun, width = 1800, height = 1400) {
  dir.create("figures", showWarnings = FALSE)
  png(
    file.path("figures", filename),
    width = width,
    height = height,
    res = 200
  )
  plot_fun()
  dev.off()
}

save_pdf <- function(filename, plot_fun, width = 8, height = 6.5) {
  dir.create("figures", showWarnings = FALSE)
  pdf(
    file.path("figures", filename),
    width = width,
    height = height
  )
  plot_fun()
  dev.off()
}


# -----------------------------------------------------------------------------
# FIGURE 1: exact operating characteristics and threshold sensitivity
# -----------------------------------------------------------------------------

plot_figure1 <- function() {
  oldpar <- par(
    mfrow = c(2, 2),
    mar = c(4.2, 4.4, 2.8, 1.2),
    oma = c(0, 0, 1, 0)
  )
  on.exit(par(oldpar))

  nvals <- sort(unique(Table1$n))
  muvals <- sort(unique(Table1$true_mean))

  mat_reject <- sapply(
    nvals,
    function(nn) {
      Table1$conventional_rejection[
        Table1$n == nn
      ]
    }
  )

  mat_stable <- sapply(
    nvals,
    function(nn) {
      Table1$stable_among_rejections[
        Table1$n == nn
      ]
    }
  )

  mat_ecs <- sapply(
    nvals,
    function(nn) {
      Table1$ECS_criterion_rate[
        Table1$n == nn
      ]
    }
  )

  matplot(
    muvals,
    mat_reject,
    type = "b",
    lty = 1,
    pch = 1:length(nvals),
    ylim = c(0, 1),
    xlab = expression(mu),
    ylab = "Probability",
    main = "A. Conventional rejection"
  )
  legend(
    "bottomright",
    legend = paste("n =", nvals),
    lty = 1,
    pch = 1:length(nvals),
    bty = "n",
    cex = .8
  )

  matplot(
    muvals,
    mat_stable,
    type = "b",
    lty = 1,
    pch = 1:length(nvals),
    ylim = c(0, 1),
    xlab = expression(mu),
    ylab = "Conditional proportion",
    main = "B. Stable among rejections"
  )

  matplot(
    muvals,
    mat_ecs,
    type = "b",
    lty = 1,
    pch = 1:length(nvals),
    ylim = c(0, 1),
    xlab = expression(mu),
    ylab = "Probability",
    main = "C. ECS criterion rate"
  )

  x <- sort(unique(threshold_sensitivity$b0))
  svals <- sort(unique(threshold_sensitivity$s_star))

  ymat <- sapply(
    svals,
    function(ss) {
      threshold_sensitivity$ecs_criterion_probability[
        threshold_sensitivity$s_star == ss
      ]
    }
  )

  matplot(
    x,
    ymat,
    type = "b",
    lty = 1,
    pch = 1:length(svals),
    ylim = c(0, 1),
    xlab = expression(b[0]),
    ylab = "Probability",
    main = expression(paste(
      "D. Threshold sensitivity, n=71, ",
      mu,
      "=0.30"
    ))
  )
  legend(
    "topright",
    legend = paste("s* =", svals),
    lty = 1,
    pch = 1:length(svals),
    bty = "n",
    cex = .8
  )
}

plot_figure1()
if (SAVE_FIGURES) {
  save_png("Figure1_ExactRates_ThresholdSensitivity.png", plot_figure1)
  save_pdf("Figure1_ExactRates_ThresholdSensitivity.pdf", plot_figure1)
}


# -----------------------------------------------------------------------------
# FIGURE 2: claim-relative exact probabilities
# -----------------------------------------------------------------------------

plot_figure2 <- function() {
  plot(
    claim_results$n,
    claim_results$reject_H0_mu_le_0,
    type = "b",
    log = "x",
    ylim = c(0, 1),
    xlab = "Sample size n (log scale)",
    ylab = "Probability",
    pch = 1,
    lty = 1,
    main = "Claim-relative exact probabilities"
  )

  lines(
    claim_results$n,
    claim_results$ECS_criterion_positivity,
    type = "b",
    pch = 2,
    lty = 2
  )

  lines(
    claim_results$n,
    claim_results$reject_H0_mu_le_030,
    type = "b",
    pch = 3,
    lty = 3
  )

  legend(
    "right",
    legend = c(
      expression(H[0]: mu <= 0),
      "ECS criterion: positivity",
      expression(H[0]: mu <= 0.30)
    ),
    lty = 1:3,
    pch = 1:3,
    bty = "n",
    cex = .85
  )
}

plot_figure2()


if (SAVE_FIGURES) {
  save_png("Figure2_ClaimRelative_Exact.png", plot_figure2)
  save_pdf("Figure2_ClaimRelative_Exact.pdf", plot_figure2)
}


# -----------------------------------------------------------------------------
# FIGURE 3: diabetes specification stability
# -----------------------------------------------------------------------------

plot_figure3 <- function() {
  xj <- jitter(
    spec_results$k,
    amount = .10
  )

  plot(
    xj,
    spec_results$coef,
    pch = ifelse(spec_results$support, 19, 1),
    xlab = "Number of added controls",
    ylab = "Estimated age coefficient",
    main = "Diabetes specification stability",
    xaxt = "n"
  )

  axis(
    1,
    at = 0:7
  )

  abline(
    h = 0,
    lty = 2
  )

  points(
    0,
    benchmark$coef,
    pch = 19,
    cex = 1.2
  )

  legend(
    "topright",
    legend = c(
      "Supports positive-age claim",
      "Does not support claim"
    ),
    pch = c(19, 1),
    bty = "n",
    cex = .8
  )
}


plot_figure3()

if (SAVE_FIGURES) {
  save_png("Figure3_Diabetes_SpecificationStability.png", plot_figure3)
  save_pdf("Figure3_Diabetes_SpecificationStability.pdf", plot_figure3)
}


# -----------------------------------------------------------------------------
# FIGURE 4: Student sleep stability-scale profile
# -----------------------------------------------------------------------------

plot_figure4 <- function() {
  plot(
    sleep_profile$b0,
    sleep_profile$S,
    type = "l",
    lwd = 2,
    xlab = expression(b[0]~"(hours)"),
    ylab = "Stability S",
    main = "Student sleep stability-scale profile"
  )

  abline(
    h = 1,
    lty = 2
  )

  abline(
    v = b0_cross,
    lty = 3
  )

  points(
    c(.25, .50),
    c(
      sleep025$stability_S,
      sleep050$stability_S
    ),
    pch = 19
  )

  text(
    b0_cross,
    1,
    labels = sprintf(
      "  crossing = %.3f h",
      b0_cross
    ),
    pos = 4,
    cex = .85
  )
}

plot_figure4()

if (SAVE_FIGURES) {
  save_png("Figure4_Sleep_StabilityProfile.png", plot_figure4)
  save_pdf("Figure4_Sleep_StabilityProfile.pdf", plot_figure4)
}


