################################################################################
#
# this file is used to specify model parameters
#
########################################################################################


########################################################
############## MODEL PARAMETERS ############################
########################################################

set.seed(myseed)

MCMC_par <- Params <- list()

#### number of cells
Params$M <-  M

#### number of targets
Params$N <- N

#### number of dimension of traits
Params$K <- dim(obs_traits)[2]

#### number of states
Params$St <- length(state)

#### hyper parameter of scaring rate S
Params$S_exp <- 2

##### hyper parameter of loss rate l
Params$l_meanlog <- -5
Params$l_sdlog <- 1

##### hyper parameter of clock rate c
Params$c_meanlog <- -5
Params$c_sdlog <- 1


##### hyper parameter of proposal sd
Params$pro_sd <- 0.01

# Optional overrides (useful for SLURM runs)
if (nzchar(Sys.getenv("BILINT_PRO_SD"))) {
	Params$pro_sd <- as.numeric(Sys.getenv("BILINT_PRO_SD"))
}

##### hyper parameter of tree
# Params$tree_br_mean <- -3.5  ## the mean of LogNormal distribution for birth rate
# Params$tree_br_sd <- 0.5     ## the sd of LogNormal distribution for birth rate
#
# Params$tree_dr_mean <- -7  ## the mean of LogNormal distribution for death rate
# Params$tree_dr_sd <- 0.5     ## the sd of LogNormal distribution for death rate
#
# Params$tree_sp_a <- 4     ## the parameter a of Beta distribution for sampling proportion
# Params$tree_sp_b <- 8     ## the parameter b of Beta distribution for sampling proportion

Params$tree_brith_rate <- 0.18 # tree brith rate
Params$tree_death_rate <- 0.03 # tree death rate
Params$tree_sampling_rate <- 0.5 # tree sampling rate

########################################################
############  specify MCMC parameters ##################
########################################################

MCMC_par$burnin <- 2000  # burnin sample size
MCMC_par$Nsamp <- 2000  # number of samples for inference
# MCMC_par$Ntune <- 300  # number of samples used for adaptive parameter tuning

# Allow common settings to be overridden via environment variables (useful for parameter tuning when submitting jobs with SLURM sbatch)
if (nzchar(Sys.getenv("BILINT_BURNIN"))) {
	MCMC_par$burnin <- as.integer(Sys.getenv("BILINT_BURNIN"))
}
if (nzchar(Sys.getenv("BILINT_NSAMP"))) {
	MCMC_par$Nsamp <- as.integer(Sys.getenv("BILINT_NSAMP"))
}

########################################################
############  optional auto-stop settings ##############
########################################################
# Automatic stopping: determine convergence based on whether the log-likelihood
# of chain 1 is approximately stable within a sliding window.
# - Controls when burn-in ends (transition to sampling phase)
# - Determines whether sampling can terminate early
# Enabled by default. To use a fixed number of iterations instead,
# set the environment variable: BILINT_AUTOSTOP=0
MCMC_par$auto_stop <- TRUE
MCMC_par$auto_stop_window <- 20        # window size (iterations) used for stability check
MCMC_par$auto_stop_check_every <- 20   # check interval (iterations)
MCMC_par$auto_stop_stable <- 1         # require this many consecutive stable checks (very permissive)
MCMC_par$auto_stop_tol_sd <- 1.0       # mean change threshold as multiple of window sd (very permissive)

# Set minimum iterations/samples for auto-stop based on problem size (number of cells).
# Rule:
# - Params$M < 100      -> 300
# - 100 <= Params$M <= 1000 -> 500
# - Params$M > 1000     -> 1000
.auto_stop_min_by_M <- function(M) {
  M <- as.numeric(M)
  if (!is.finite(M)) return(500L)
  if (M < 100) return(300L)
  if (M <= 1000) return(500L)
  1000L
}
.auto_min <- .auto_stop_min_by_M(Params$M)
MCMC_par$auto_stop_burnin_min <- .auto_min  # minimum burnin iterations before allowing burnin to stop
MCMC_par$auto_stop_samp_min <- .auto_min    # minimum saved samples before allowing early stop

MCMC_par$swap_interval <- 3 # make Matroplis Hastings move in every how many samples
MCMC_par$Nchain <- 4 # number of paralel chains
if (nzchar(Sys.getenv("BILINT_NCHAIN"))) {
  MCMC_par$Nchain <- as.integer(Sys.getenv("BILINT_NCHAIN"))
}
MCMC_par$delta_T <- 0.001


if (nzchar(Sys.getenv("BILINT_DELTA_T"))) {
	MCMC_par$delta_T <- as.numeric(Sys.getenv("BILINT_DELTA_T"))
}
if (nzchar(Sys.getenv("BILINT_SWAP_INTERVAL"))) {
	MCMC_par$swap_interval <- as.integer(Sys.getenv("BILINT_SWAP_INTERVAL"))
}

Temperature <- seq(1, by=MCMC_par$delta_T, length.out = MCMC_par$Nchain)  # temperatures


########################################################
############  optional inference for t2 ################
########################################################
# If the user sets `t2 <- NULL` (e.g., in BiLinT_app_original.R) before sourcing
# this file, we will infer t2; otherwise t2 is treated as fixed.

Params$sample_t2 <- exists("t2", inherits = TRUE) && is.null(t2)

if (isTRUE(Params$sample_t2)) {
  if (!exists("t1", inherits = TRUE) || !is.numeric(t1) || length(t1) != 1) {
    stop("To infer t2, t1 must exist and be a numeric scalar.")
  }
  if (!exists("t_tol", inherits = TRUE) || !is.numeric(t_tol) || length(t_tol) != 1) {
    stop("To infer t2, t_tol must exist and be a numeric scalar.")
  }

  Params$t2_lower <- as.numeric(t1)
  Params$t2_upper <- as.numeric(t_tol)
  if (!is.finite(Params$t2_lower) || !is.finite(Params$t2_upper) || Params$t2_upper <= Params$t2_lower) {
    stop("Invalid t2 prior bounds: require finite t1 < t_tol.")
  }

  # Prior: t2 ~ Uniform(t1, t_tol)
  # Proposal SD (in absolute time units). Keep it reasonably large by default.
  Params$t2_pro_sd <- 1
  if (nzchar(Sys.getenv("BILINT_T2_PRO_SD"))) {
    Params$t2_pro_sd <- as.numeric(Sys.getenv("BILINT_T2_PRO_SD"))
  }

  # Choose an initial numeric value for t2.
  # Prefer simulator-provided value if available; otherwise use midpoint of the prior.
  t2_init <- NULL
  if (exists("muti_data", inherits = TRUE) && !is.null(muti_data$seq_para_data$t2)) {
    t2_init <- as.numeric(muti_data$seq_para_data$t2)
  }
  if (is.null(t2_init) || !is.finite(t2_init)) {
    t2_init <- 0.5 * (Params$t2_lower + Params$t2_upper)
  }
  # Clamp to the support (avoid boundary issues)
  eps <- 1e-3
  t2 <- min(max(t2_init, Params$t2_lower + eps), Params$t2_upper - eps)
} else {
  # Fixed t2 mode (backward compatible). Keep Params fields defined for downstream checks.
  Params$t2_lower <- NULL
  Params$t2_upper <- NULL
  Params$t2_pro_sd <- NULL
}


