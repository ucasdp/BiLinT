# Post-run diagnostics for BiLinT sampler
# Generates trace plots and quantitative metrics to demonstrate posterior exploration.
# This file must not alter sampling behavior; it only reads saved objects (Trace, MCMC_par, beta_log_like).

BiLinT_extract_trace_df <- function(Trace) {
  if (length(Trace) == 0) {
    return(data.frame())
  }

  get_scalar <- function(x, idx = 1L) {
    if (is.null(x)) return(NA_real_)
    if (length(x) < idx) return(NA_real_)
    as.numeric(x[[idx]])
  }

  n <- length(Trace)

  df <- data.frame(
    iter = seq_len(n),
    likelihood = vapply(Trace, function(s) as.numeric(s$likelihood), numeric(1)),
    lik_crispr = vapply(Trace, function(s) as.numeric(s$lik_crispr), numeric(1)),
    lik_traits = vapply(Trace, function(s) as.numeric(s$lik_traits), numeric(1)),
    l = vapply(Trace, function(s) as.numeric(s$l), numeric(1)),
    c = vapply(Trace, function(s) as.numeric(s$c), numeric(1)),
    alpha = vapply(Trace, function(s) as.numeric(s$alpha), numeric(1)),
    stringsAsFactors = FALSE
  )

  # A few representative coordinates for vector parameters (plots stay readable)
  df$mu1 <- vapply(Trace, function(s) get_scalar(s$mu, 1L), numeric(1))
  df$theta1 <- vapply(Trace, function(s) get_scalar(s$theta, 1L), numeric(1))
  df$sigma1 <- vapply(Trace, function(s) get_scalar(s$sigma, 1L), numeric(1))
  df$S1 <- vapply(Trace, function(s) get_scalar(s$S, 1L), numeric(1))

  df
}

BiLinT_compute_metrics <- function(df) {
  if (nrow(df) == 0) {
    return(list(summary = data.frame(), ess = data.frame(), acf = list()))
  }

  # Use coda for ESS. For ACF we use stats::acf.
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for diagnostics (effectiveSize).")
  }

  metrics_vars <- c(
    "likelihood", "lik_crispr", "lik_traits",
    "l", "c", "alpha", "mu1", "theta1", "sigma1", "S1"
  )
  metrics_vars <- metrics_vars[metrics_vars %in% colnames(df)]

  ess_vals <- vapply(metrics_vars, function(v) {
    x <- df[[v]]
    x <- x[is.finite(x)]
    if (length(x) < 3) return(NA_real_)
    as.numeric(coda::effectiveSize(x))
  }, numeric(1))

  ess_df <- data.frame(
    variable = metrics_vars,
    ess = ess_vals,
    ess_per_sample = ess_vals / nrow(df),
    stringsAsFactors = FALSE
  )

  summ_df <- data.frame(
    variable = metrics_vars,
    mean = vapply(metrics_vars, function(v) mean(df[[v]], na.rm = TRUE), numeric(1)),
    sd = vapply(metrics_vars, function(v) stats::sd(df[[v]], na.rm = TRUE), numeric(1)),
    q025 = vapply(metrics_vars, function(v) stats::quantile(df[[v]], 0.025, na.rm = TRUE, names = FALSE), numeric(1)),
    q50 = vapply(metrics_vars, function(v) stats::quantile(df[[v]], 0.50, na.rm = TRUE, names = FALSE), numeric(1)),
    q975 = vapply(metrics_vars, function(v) stats::quantile(df[[v]], 0.975, na.rm = TRUE, names = FALSE), numeric(1)),
    stringsAsFactors = FALSE
  )

  # ACF values (returned for potential further use)
  acf_list <- lapply(metrics_vars, function(v) {
    x <- df[[v]]
    x <- x[is.finite(x)]
    if (length(x) < 3) return(NULL)
    stats::acf(x, plot = FALSE)
  })
  names(acf_list) <- metrics_vars

  list(summary = summ_df, ess = ess_df, acf = acf_list)
}

BiLinT_write_diagnostics <- function(
  foldername,
  dataname,
  Trace,
  MCMC_par = NULL,
  beta_log_like = NULL,
  file_prefix = "diagnostics"
) {
  dir.create(foldername, showWarnings = FALSE, recursive = TRUE)

  df <- BiLinT_extract_trace_df(Trace)
  metrics <- BiLinT_compute_metrics(df)

  # Write CSV metrics
  base <- file.path(foldername, paste0(dataname, "_", file_prefix))
  utils::write.csv(metrics$summary, paste0(base, "_summary.csv"), row.names = FALSE)
  utils::write.csv(metrics$ess, paste0(base, "_ess.csv"), row.names = FALSE)

  # Swap/tempering metrics if available
  if (!is.null(beta_log_like) && is.matrix(beta_log_like)) {
    # beta_log_like rows = chains, cols = saved samples
    bl <- beta_log_like
    wbic_like_mean <- mean(bl[nrow(bl), , drop = TRUE], na.rm = TRUE)
    utils::write.csv(
      data.frame(wbic_chain_index = nrow(bl), wbic_like_mean = wbic_like_mean),
      paste0(base, "_tempering.csv"),
      row.names = FALSE
    )
  }

  if (!is.null(MCMC_par)) {
    # dump key run settings (does not include large objects)
    run_df <- data.frame(
      burnin = if (!is.null(MCMC_par$burnin)) MCMC_par$burnin else NA,
      burnin_used = if (!is.null(MCMC_par$burnin_used)) MCMC_par$burnin_used else NA,
      Nsamp = if (!is.null(MCMC_par$Nsamp)) MCMC_par$Nsamp else NA,
      Nsamp_used = if (!is.null(MCMC_par$Nsamp_used)) MCMC_par$Nsamp_used else NA,
      Nchain = if (!is.null(MCMC_par$Nchain)) MCMC_par$Nchain else NA,
      swap_interval = if (!is.null(MCMC_par$swap_interval)) MCMC_par$swap_interval else NA,
      stringsAsFactors = FALSE
    )
    utils::write.csv(run_df, paste0(base, "_run.csv"), row.names = FALSE)
  }

  # PDF trace plots
  grDevices::pdf(paste0(base, "_trace.pdf"), width = 8.5, height = 11)
  on.exit(grDevices::dev.off(), add = TRUE)

  if (nrow(df) == 0) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No samples in Trace; diagnostics skipped.")
    return(invisible(list(df = df, metrics = metrics)))
  }

  plot_series <- function(y, main) {
    graphics::plot(df$iter, y, type = "l", xlab = "Saved sample", ylab = main, main = main)
  }

  plot_series(df$likelihood, "log-likelihood (total)")
  plot_series(df$lik_crispr, "log-likelihood (crispr)")
  plot_series(df$lik_traits, "log-likelihood (traits)")
  plot_series(df$l, "l")
  plot_series(df$c, "c")
  plot_series(df$alpha, "alpha")
  plot_series(df$mu1, "mu[1]")
  plot_series(df$theta1, "theta[1]")
  plot_series(df$sigma1, "sigma[1]")
  plot_series(df$S1, "S[1]")

  # ACF plots
  acf_vars <- c("likelihood", "l", "c", "alpha", "mu1", "theta1")
  acf_vars <- acf_vars[acf_vars %in% colnames(df)]
  for (v in acf_vars) {
    x <- df[[v]]
    x <- x[is.finite(x)]
    if (length(x) >= 3) {
      stats::acf(x, main = paste0("ACF: ", v))
    }
  }

  invisible(list(df = df, metrics = metrics))
}
