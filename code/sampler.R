################################################################
# this file contains sampling procedure procedure of model
##################################################################

print_red <- function(text) {
  cat(sprintf("\033[31m%s\033[0m\n", text))
}

######################################################
######## initial value for parameters  ###############
######################################################

New_par <- vector('list',MCMC_par$Nchain)

for(i in 1:MCMC_par$Nchain){
  tree_init <- generate_tree(seq_data_mt, obs_traits_original)
  New_par[[i]] <- init_gen_fun(Params, tree_init)
}

tree_info <- tree_transform2(New_par[[1]]$tree)
model.OU.BM <- BM_OU_obj(New_par[[1]], tree_info, Params)

for(i in 1:MCMC_par$Nchain){
  l <- New_par[[i]]$l
  Sc <- New_par[[i]]$S * New_par[[i]]$c
  r <- l + sum(Sc)

  tree <- New_par[[i]]$tree
  tree_Cas9 <- tree_transform1(tree, t1, New_par[[i]]$t2)
  tree_PostOrder <- c(reorder_tree_edges(tree_Cas9, root_to_tips = FALSE)$edge[,2], length(tree_Cas9$tip.label)+1)
  childlist <- ChildList(tree_Cas9, tree_PostOrder)
  dis <- get_all_distances_to_root(tree_Cas9, as_edge_count=FALSE)
  dis[which(dis==0)] <- 10^(-5)

  lik_crispr <- tree_likehood(tree_Cas9, tree_PostOrder, dis, childlist, seq_data_mt, state, t1, New_par[[i]]$t2, r, l, Sc)
  while(is.infinite(lik_crispr) && lik_crispr < 0){
    tree_init <- generate_tree(seq_data_mt, obs_traits_original)
    New_par[[i]] <- init_gen_fun(Params, tree_init)
    l <- New_par[[i]]$l
    Sc <- New_par[[i]]$S * New_par[[i]]$c
    r <- l + sum(Sc)

    tree <- New_par[[i]]$tree
    tree_Cas9 <- tree_transform1(tree, t1, New_par[[i]]$t2)
    tree_PostOrder <- c(reorder_tree_edges(tree_Cas9, root_to_tips = FALSE)$edge[,2], length(tree_Cas9$tip.label)+1)
    childlist <- ChildList(tree_Cas9, tree_PostOrder)
    dis <- get_all_distances_to_root(tree_Cas9, as_edge_count=FALSE)
    dis[which(dis==0)] <- 10^(-5)
    lik_crispr <- tree_likehood(tree_Cas9, tree_PostOrder, dis, childlist, seq_data_mt, state, t1, New_par[[i]]$t2, r, l, Sc)
  }

  New_par[[i]]$lik_crispr <- lik_crispr
  tree_info <- tree_transform2(New_par[[i]]$tree)
  model.OU.BM2 <- BM_OU_obj_change(New_par[[i]], tree_info, model.OU.BM, Params)
  New_par[[i]]$lik_traits <- lik_obs_traits(New_par[[i]], tree_info, model.OU.BM2, obs_traits)
  New_par[[i]]$likelihood <- New_par[[i]]$lik_crispr + New_par[[i]]$lik_traits
  #New_par[[i]]$lik_crispr <- tree_likehood_R(tree_Cas9, seq_data_mt, state, t1, t2,  r, l, Sc)
}

#################################################################
################### MCMC sampling  #####################
################################################################

start_t <- Sys.time()

burnin_max <- MCMC_par$burnin
nsamp_max <- MCMC_par$Nsamp
Nrep_max <- burnin_max + nsamp_max # maximum total iterations

auto_stop <- isTRUE(MCMC_par$auto_stop)
if (nzchar(Sys.getenv("BILINT_AUTOSTOP")) && Sys.getenv("BILINT_AUTOSTOP") == "0") {
  auto_stop <- FALSE
}

auto_verbose <- nzchar(Sys.getenv("BILINT_AUTOSTOP_VERBOSE")) && Sys.getenv("BILINT_AUTOSTOP_VERBOSE") == "1"

auto_window <- if (!is.null(MCMC_par$auto_stop_window)) MCMC_par$auto_stop_window else 50
auto_check_every <- if (!is.null(MCMC_par$auto_stop_check_every)) MCMC_par$auto_stop_check_every else 25
auto_stable_required <- if (!is.null(MCMC_par$auto_stop_stable)) MCMC_par$auto_stop_stable else 3
auto_tol_sd <- if (!is.null(MCMC_par$auto_stop_tol_sd)) MCMC_par$auto_stop_tol_sd else 0.1
burnin_min <- if (!is.null(MCMC_par$auto_stop_burnin_min)) MCMC_par$auto_stop_burnin_min else max(50, floor(0.2 * burnin_max))
nsamp_min <- if (!is.null(MCMC_par$auto_stop_samp_min)) MCMC_par$auto_stop_samp_min else max(100, floor(0.4 * nsamp_max))

if (nzchar(Sys.getenv("BILINT_AUTOSTOP_WINDOW"))) {
  auto_window <- as.integer(Sys.getenv("BILINT_AUTOSTOP_WINDOW"))
}
if (nzchar(Sys.getenv("BILINT_AUTOSTOP_CHECK_EVERY"))) {
  auto_check_every <- as.integer(Sys.getenv("BILINT_AUTOSTOP_CHECK_EVERY"))
}
if (nzchar(Sys.getenv("BILINT_AUTOSTOP_STABLE"))) {
  auto_stable_required <- as.integer(Sys.getenv("BILINT_AUTOSTOP_STABLE"))
}
if (nzchar(Sys.getenv("BILINT_AUTOSTOP_TOL_SD"))) {
  auto_tol_sd <- as.numeric(Sys.getenv("BILINT_AUTOSTOP_TOL_SD"))
}
if (nzchar(Sys.getenv("BILINT_AUTOSTOP_BURNIN_MIN"))) {
  burnin_min <- as.integer(Sys.getenv("BILINT_AUTOSTOP_BURNIN_MIN"))
}
if (nzchar(Sys.getenv("BILINT_AUTOSTOP_SAMP_MIN"))) {
  nsamp_min <- as.integer(Sys.getenv("BILINT_AUTOSTOP_SAMP_MIN"))
}

auto_window <- max(5L, as.integer(auto_window))
auto_check_every <- max(1L, as.integer(auto_check_every))
auto_stable_required <- max(1L, as.integer(auto_stable_required))
auto_tol_sd <- as.numeric(auto_tol_sd)
burnin_min <- max(1L, as.integer(burnin_min))
nsamp_min <- max(1L, as.integer(nsamp_min))

burnin_min <- min(burnin_min, burnin_max)
nsamp_min <- min(nsamp_min, nsamp_max)

# Reproducibility switch: enabled by default (ensures identical results across
# repeated runs on the same machine with the same number of cores/chains).
# For more dynamic load balancing, set the environment variable: BILINT_REPRO=0
reproducible <- !(nzchar(Sys.getenv("BILINT_REPRO")) && Sys.getenv("BILINT_REPRO") == "0")

tune_samp <- vector('list',MCMC_par$Nchain)
Trace <- vector('list',nsamp_max)
beta_log_like <- matrix(0, MCMC_par$Nchain, nsamp_max)  # keep beta chain likelihood for WBIC
lik_chain1_saved <- rep(NA_real_, nsamp_max)
lik_chain1_iter <- rep(NA_real_, Nrep_max)
# Trace <- vector('list',Nrep)
# beta_log_like <- matrix(0,MCMC_par$Nchain, Nrep)

# stanmodel <- stan_model(file = 'OU_BM_sampling2.stan')

# Windows does not support fork; therefore, parallel::mclapply(mc.cores > 1)
# will fail on Windows.
# To preserve parallel execution:
# - Windows: use a PSOCK cluster with parLapply
# - macOS/Linux: use the default fork-based mclapply (can be switched to a
#   persistent cluster via an environment variable to improve performance)
.get_allocated_cores <- function() {
  if (nzchar(Sys.getenv("BILINT_CORES"))) {
    n <- suppressWarnings(as.integer(Sys.getenv("BILINT_CORES")))
    if (!is.na(n) && n > 0) return(n)
  }
  n <- NA_integer_
  if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK"))) {
    n <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
  }
  if (!is.na(n) && n > 0) return(n)
  parallel::detectCores()
}
.get_parallel_backend <- function() {
  # Values:
  # - "mclapply" (default on non-Windows)
  # - "psock"    (persistent PSOCK cluster)
  # - "fork"     (persistent FORK cluster; Unix only)
  # - "auto"     (Windows -> psock; others -> mclapply)
  val <- Sys.getenv("BILINT_PARALLEL_BACKEND")
  if (!nzchar(val)) val <- Sys.getenv("BILINT_BACKEND")
  if (!nzchar(val)) val <- "auto"
  val <- tolower(trimws(val))
  if (val %in% c("auto", "mclapply", "psock", "fork")) return(val)
  # Backward-compatible booleans
  if (val %in% c("1", "true", "cluster")) return("psock")
  "auto"
}
.use_cluster_backend <- function() {
  if (.Platform$OS.type == "windows") return(TRUE)
  backend <- .get_parallel_backend()
  if (backend == "auto") return(FALSE)
  backend != "mclapply"
}
.cluster_type <- function() {
  backend <- .get_parallel_backend()
  if (backend == "fork" && .Platform$OS.type != "windows") return("FORK")
  "PSOCK"
}
.get_project_root <- function() {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) {
    return(normalizePath(file.path(dirname(ofile), ".."), winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
project_root <- .get_project_root()
code_dir <- file.path(project_root, "code")

cl <- NULL
if (.use_cluster_backend()) {
  nWorkers <- min(MCMC_par$Nchain, .get_allocated_cores())
  cl <- parallel::makeCluster(nWorkers, type = .cluster_type())
  doParallel::registerDoParallel(cl)

  parallel::clusterExport(cl, c("code_dir"), envir = environment())

  # Load required packages and function definitions on each worker
  # (prevents errors such as 'tree_transform2' or 'samp_mu_all' not found)
  parallel::clusterEvalQ(cl, {
    suppressMessages(library(ape))
    suppressMessages(library(phangorn))
    suppressMessages(library(PCMBase))
    suppressMessages(library(mclust))
    suppressMessages(library(PCMBaseCpp))
    suppressMessages(library(msm))
    suppressMessages(library(castor))
    suppressMessages(library(TreeTools))
    suppressMessages(library(TreePar))

    source(file.path(code_dir, "Crispr_Cas9_lik.R"))
    if (file.exists(file.path(code_dir, "Crispr_Cas9_cpplik.R"))) {
      source(file.path(code_dir, "Crispr_Cas9_cpplik.R"))
    }
    source(file.path(code_dir, "model.OU.BM_pro.R"))
    source(file.path(code_dir, "par_samp_para.R"))
    source(file.path(code_dir, "main_fun.R"))

    BiLinT_worker_fun <- function(x, h) {
      MCMC_onesamp_II(h, Params, seq_data_mt, x$state, obs_traits, Temperature[x$i])
    }
    NULL
  })

  # Export static objects once (these are required by the parallel chains)
  parallel::clusterExport(
    cl,
    c(
      "Params", "seq_data_mt", "obs_traits", "Temperature",
      "t1", "t2", "t_tol", "state",
      "model.OU.BM",
      "obs_traits_original"
    ),
    envir = environment()
  )

  if (exists("myseed")) {
    parallel::clusterSetRNGStream(cl, as.integer(myseed))
  }
}
sawp_time <- 0
burnin_used <- burnin_max
stable_count <- 0
ct <- 0

stop_reason <- NULL
stop_iter <- NA_integer_
last_iter <- 0L

auto_last_burnin_iter <- NA_integer_
auto_last_burnin_diff <- NA_real_
auto_last_burnin_thresh <- NA_real_
auto_last_burnin_sd <- NA_real_
auto_last_burnin_stable_count <- 0L
auto_max_stable_burnin <- 0L

auto_last_samp_iter <- NA_integer_
auto_last_samp_saved <- NA_integer_
auto_last_samp_diff <- NA_real_
auto_last_samp_thresh <- NA_real_
auto_last_samp_sd <- NA_real_
auto_last_samp_stable_count <- 0L
auto_max_stable_samp <- 0L

for(h in 1:Nrep_max){
  # print_red(paste('the current iteration count:', h, sep=' '))

  last_iter <- h

  ######## perform MCMC on each chahn ##########
  # for (z in 1:MCMC_par$Nchain){
  #   temp <- MCMC_onesamp_II(h, Params, seq_data_mt, New_par[[z]], obs_traits, temper=Temperature[z])
  #   New_par[[z]] <- temp
  # }
  # export_vars <- c("Params","seq_data_mt", "New_par", "obs_traits", "Temperature", 't1', 't2')
  # New_par <- foreach(i=1:4,.packages=c('ape', 'phangorn', 'PCMBase', 'ggtree', 'mclust', 'Rcpp', 'PCMBaseCpp', 'msm',
  #                                      'castor', 'Rcpp2doParallel', 'TreeTools', 'TreePar'),.noexport = 'tree_likehood') %do% {
  #                                        if (!exists("tree_likehood")) {
  #                                          source("Crispr_Cas9_cpplik.R")
  #                                        }
  #                                        MCMC_onesamp_II(h, Params, seq_data_mt, New_par[[i]], obs_traits, Temperature[i])
  #                                      }
  if (!is.null(cl)) {
    task_input <- Map(
      function(i, state) list(i = i, state = state),
      seq_len(MCMC_par$Nchain),
      New_par
    )
    wrapper_fun <- function(x, h) {
      tryCatch(
        BiLinT_worker_fun(x, h),
        error = function(e) {
          structure(
            list(
              i = x$i,
              h = h,
              message = conditionMessage(e),
              call = tryCatch(paste(deparse(conditionCall(e)), collapse = ""), error = function(e2) NA_character_),
              sys_calls = tryCatch(vapply(sys.calls(), function(z) paste(deparse(z), collapse = ""), character(1)), error = function(e2) character(0))
            ),
            class = "BiLinT_worker_error"
          )
        }
      )
    }
    if (reproducible) {
      New_par <- parallel::parLapply(cl, task_input, wrapper_fun, h = h)
    } else {
      New_par <- parallel::parLapplyLB(cl, task_input, wrapper_fun, h = h)
    }
  } else {
    mc_cores <- min(MCMC_par$Nchain, .get_allocated_cores())
    New_par <- parallel::mclapply(
      1:MCMC_par$Nchain,
      function(i) MCMC_onesamp_II(h, Params, seq_data_mt, New_par[[i]], obs_traits, Temperature[i]),
      mc.cores = mc_cores,
      mc.preschedule = reproducible
    )
  }

  ############### switch between chains ###########################

  if(MCMC_par$Nchain>1){
    if(h %% MCMC_par$swap_interval == 0)
    {
      a <- sample(c(1:(MCMC_par$Nchain-1)),1) # randomly select a chain
      b <- a+1
      #log  of numerater
      logp_num <- New_par[[a]]$likelihood/Temperature[b] + New_par[[b]]$likelihood/Temperature[a]
      logp_num <-  logp_num + log_prior_all(New_par[[a]],Params,Temperature[b]) +
        log_prior_all(New_par[[b]],Params,Temperature[a])

      # logp_num <-  log_prior_all(New_par[[a]],Params,1) +
      #   log_prior_all(New_par[[b]],Params,1)

      #log of denominator
      logp_denom <-  New_par[[b]]$likelihood/Temperature[b] + New_par[[a]]$likelihood/Temperature[a]
      logp_denom <-  logp_denom + log_prior_all(New_par[[a]],Params,Temperature[a]) +
        log_prior_all(New_par[[b]],Params,Temperature[b])

      # logp_denom <-  log_prior_all(New_par[[a]],Params,1) +
      #   log_prior_all(New_par[[b]],Params,1)
      acc_prob <- min(1,exp(logp_num-logp_denom))  # probability of accepting swap

      if(runif(1)<acc_prob)
      {
        sawp_time <- sawp_time + 1
        temp <- New_par[[a]]
        New_par[[a]] <- New_par[[b]]
        New_par[[b]] <- temp
        # message('samples switched')
      }
    }
  }

  # Record the log-likelihood trajectory of chain 1 (used for automatic burn-in and early stopping decisions)
  lik_chain1_iter[h] <- New_par[[1]]$likelihood

  in_burnin <- h <= burnin_used

  # Automatic burn-in detection: terminate burn-in early if the log-likelihood
  # stabilizes within a sliding window
  if (auto_stop && in_burnin) {
    if (h >= burnin_min && (h %% auto_check_every == 0) && h >= 2 * auto_window) {
      w1 <- lik_chain1_iter[(h - auto_window + 1):h]
      w0 <- lik_chain1_iter[(h - 2 * auto_window + 1):(h - auto_window)]
      s <- sd(c(w0, w1), na.rm = TRUE)
      s <- if (is.finite(s) && s > 0) s else 1
      mean_diff <- abs(mean(w1, na.rm = TRUE) - mean(w0, na.rm = TRUE))
      thresh <- auto_tol_sd * s
      stable_now <- mean_diff <= thresh
      stable_count <- if (stable_now) stable_count + 1 else 0

      auto_last_burnin_iter <- h
      auto_last_burnin_diff <- mean_diff
      auto_last_burnin_thresh <- thresh
      auto_last_burnin_sd <- s
      auto_last_burnin_stable_count <- stable_count
      auto_max_stable_burnin <- max(auto_max_stable_burnin, stable_count)

      if (auto_verbose) {
        cat(
          sprintf(
            "[AutoBurninCheck] iter=%d diff=%g thresh=%g (sd=%g, tol_sd=%g) stable=%s stable_count=%d/%d\n",
            h, mean_diff, thresh, s, auto_tol_sd, if (stable_now) "TRUE" else "FALSE", stable_count, auto_stable_required
          )
        )
      }
      if (stable_count >= auto_stable_required) {
        old_burnin_used <- burnin_used
        burnin_used <- h
        stable_count <- 0
        if (!isTRUE(old_burnin_used == burnin_used)) {
          cat(sprintf("[AutoBurnin] burnin ended early at iter=%d (burnin_used=%d)\n", h, burnin_used))
        }
      }
    }
  }

  ######## after burnin stage: save new samples into Trace ###############
  #### only keep chain 1
  if (h == burnin_used + 1) {
    # Reset the stability counter upon entering the sampling phase to prevent
    # burn-in stability counts from carrying over to early stopping
    stable_count <- 0
  }

  if (h > burnin_used) {
    ct <- ct + 1
    if (ct <= nsamp_max) {
      Trace[[ct]] <- New_par[[1]]
      lik_chain1_saved[ct] <- New_par[[1]]$likelihood
      beta_log_like[, ct] <- vapply(New_par, function(p) p$likelihood, numeric(1))
    }

    if (auto_stop && ct >= nsamp_min && (ct %% auto_check_every == 0) && ct >= 2 * auto_window) {
      w1 <- lik_chain1_saved[(ct - auto_window + 1):ct]
      w0 <- lik_chain1_saved[(ct - 2 * auto_window + 1):(ct - auto_window)]
      s <- sd(c(w0, w1), na.rm = TRUE)
      s <- if (is.finite(s) && s > 0) s else 1
      mean_diff <- abs(mean(w1, na.rm = TRUE) - mean(w0, na.rm = TRUE))
      thresh <- auto_tol_sd * s
      stable_now <- mean_diff <= thresh
      stable_count <- if (stable_now) stable_count + 1 else 0

      auto_last_samp_iter <- h
      auto_last_samp_saved <- ct
      auto_last_samp_diff <- mean_diff
      auto_last_samp_thresh <- thresh
      auto_last_samp_sd <- s
      auto_last_samp_stable_count <- stable_count
      auto_max_stable_samp <- max(auto_max_stable_samp, stable_count)

      if (auto_verbose) {
        cat(
          sprintf(
            "[AutoStopCheck] iter=%d saved=%d diff=%g thresh=%g (sd=%g, tol_sd=%g) stable=%s stable_count=%d/%d\n",
            h, ct, mean_diff, thresh, s, auto_tol_sd, if (stable_now) "TRUE" else "FALSE", stable_count, auto_stable_required
          )
        )
      }
      if (stable_count >= auto_stable_required) {
        stop_reason <- "auto_stop_sampling_stable"
        stop_iter <- h
        cat(sprintf("[AutoStop] sampling stopped early at iter=%d (saved=%d/%d)\n", h, ct, nsamp_max))
        break
      }
    }

    if (ct >= nsamp_max) {
      stop_reason <- "reached_nsamp_max"
      stop_iter <- h
      break
    }
  }


  ########## display progress ###############

  if (h %% 100 == 0) {
    ave_speed <- difftime(Sys.time(), start_t, units = "mins")/h
    target_total <- burnin_used + nsamp_max
    target_total <- min(target_total, Nrep_max)
    phase <- if (h <= burnin_used) "burnin" else "sampling"
    rem_iters <- max(0, target_total - h)
    cat(
      sprintf(
        "[Progress] iter=%d/%d (target=%d) phase=%s burnin_used=%d saved=%d/%d; time remaining: %s\n",
        h, Nrep_max, target_total, phase, burnin_used, ct, nsamp_max, as.character(rem_iters * ave_speed)
      )
    )
  }
}

if (is.null(stop_reason)) {
  stop_reason <- "reached_max_iterations"
  stop_iter <- last_iter
}

cat(sprintf("[Done] reason=%s; iter=%d; burnin_used=%d; saved=%d/%d\n", stop_reason, stop_iter, burnin_used, ct, nsamp_max))

if (auto_stop) {
  burnin_iter_txt <- if (is.na(auto_last_burnin_iter)) "NA" else as.character(auto_last_burnin_iter)
  samp_iter_txt <- if (is.na(auto_last_samp_iter)) "NA" else as.character(auto_last_samp_iter)
  samp_saved_txt <- if (is.na(auto_last_samp_saved)) "NA" else as.character(auto_last_samp_saved)
  cat(
    sprintf(
      "[AutoStopSummary] burnin_last(iter=%s diff=%g thresh=%g max_stable=%d/%d); samp_last(iter=%s saved=%s diff=%g thresh=%g max_stable=%d/%d)\n",
      burnin_iter_txt, auto_last_burnin_diff, auto_last_burnin_thresh, auto_max_stable_burnin, auto_stable_required,
      samp_iter_txt, samp_saved_txt, auto_last_samp_diff, auto_last_samp_thresh, auto_max_stable_samp, auto_stable_required
    )
  )
}

if (!is.null(cl)) {
  parallel::stopCluster(cl)
}

cat('completed. Time consumed:',difftime(Sys.time(),start_t,units = "mins"),'mins \n \n')

# Truncate to the actual number of collected samples
if (ct == 0) {
  Trace <- list()
  beta_log_like <- beta_log_like[, 0, drop = FALSE]
} else if (ct < nsamp_max) {
  Trace <- Trace[seq_len(ct)]
  beta_log_like <- beta_log_like[, seq_len(ct), drop = FALSE]
}

# Record the actual burn-in and Nsamp used (the original fields still store the maximum values)
MCMC_par$burnin_used <- burnin_used
MCMC_par$Nsamp_used <- ct




##################################################################
########### save data ############################################
##################################################################

# cur_file=paste(foldername,'/seed',myseed,'_N', N, '_M', M, '_K', K, '_', rep, '.Rdata',sep='')
cur_file=paste(foldername,'/', dataname, '.Rdata',sep='')
save(Params, MCMC_par, obs_traits_original, Trace, myseed, beta_log_like, file = cur_file)
