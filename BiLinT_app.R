
rm(list=ls())


# ---- Reproducibility guardrails ----
# Even with the same `set.seed()`, multi-threaded BLAS/OpenMP can introduce tiny
# floating-point differences that may flip MH accept/reject and lead to different traces.
# To make repeated `source('BiLinT_app.R')` runs stable, default to 1 thread unless the
# user already set these env vars.
.setenv_if_unset <- function(name, value) {
  if (!nzchar(Sys.getenv(name))) do.call(Sys.setenv, setNames(list(value), name))
}
.setenv_if_unset("OMP_NUM_THREADS", "1")
.setenv_if_unset("OPENBLAS_NUM_THREADS", "1")
.setenv_if_unset("MKL_NUM_THREADS", "1")
.setenv_if_unset("VECLIB_MAXIMUM_THREADS", "1")  # macOS Accelerate
.setenv_if_unset("NUMEXPR_NUM_THREADS", "1")

################################################################
############# load packages ##################################
################################################################

require(tidyr)
require(ggplot2)
require(dplyr)
require(coda)
require(Rcpp)
require(RcppArmadillo)
require(LaplacesDemon)
require(PCMBase)
require(ape)
require(ggtree)
require(castor)
require(phangorn)
require(TreeSearch)
require(TreePar)
require(TreeSim)
require(msm) ##rtmsm
require(Rfast)
require(stringr)
require(treeio)
require(ggstance)
require(TreeTools)
require(phytools)
require(PCMBaseCpp)
require(rdist)
require(foreach)
require(doParallel)
require(Rcpp2doParallel)


########################################################
######## load functions ################################
########################################################

source('code/Crispr_Cas9_cpplik.R')
source('code/Crispr_Cas9_lik.R')
source('code/main_fun.R')
source('code/par_samp_para.R')
source('code/assist_fun.R')
source('code/model.OU.BM_pro.R')



#################################################################
########## MODEL INPUT ############################################
#################################################################

myseed <-  1               # set random seed
RNGkind("L'Ecuyer-CMRG")
set.seed(myseed)
foldername <-  "temp_out_toy"          # set output foldername
dataname <- 'toy_model'

# Optional overrides for batch/Slurm runs (defaults above remain if unset)
if (nzchar(Sys.getenv("BILINT_FOLDER"))) foldername <- Sys.getenv("BILINT_FOLDER")
if (nzchar(Sys.getenv("BILINT_DATANAME"))) dataname <- Sys.getenv("BILINT_DATANAME")

dir.create(foldername, showWarnings = FALSE)  # folder where outputs are saved

muti_data <- readRDS('toy_model_data.rds')
seq_data_mt <- muti_data$seq_data_mt ### barcode data, where row represents cells, columns represent targets
M <- dim(muti_data$seq_data_mt)[1]
N <- dim(muti_data$seq_data_mt)[2]
state <- sort(setdiff(unique(as.vector(seq_data_mt)), c(0, -1))) ### 0 denotes the unedited state, whereas -1 denotes the silenced state.

obs_traits_original <- muti_data$obs_traits ### gene expression data, where row represents cells, columns represent genes
obs_traits <- prcomp(obs_traits_original, rank. = 5)$x
K <- dim(obs_traits)[2]

t1 <- 0  # scaring start time
t2 <- muti_data$seq_para_data$t2 # scaring end time (set to NULL to infer t2)
# t2 <- NULL
t_tol <- muti_data$seq_para_data$t_tol # total experiment time


##############################################
######## load parameter file ################
############################################
source('code/specify_pars.R')

# Ensure deterministic parallel behavior unless user overrides.
if (!nzchar(Sys.getenv("BILINT_REPRO"))) Sys.setenv(BILINT_REPRO = "1")
# Fix worker count so `detectCores()` changes (or scheduler env) don't change the RNG stream.
if (!nzchar(Sys.getenv("BILINT_CORES"))) Sys.setenv(BILINT_CORES = as.character(MCMC_par$Nchain))


##############################################
######## sampling ##########################
############################################
source("code/sampler.R")

if (!(nzchar(Sys.getenv("BILINT_SKIP_POST")) && Sys.getenv("BILINT_SKIP_POST") == "1")) {
  ##############################################
  ######## diagnostics (trace/ESS/ACF) ########
  ##############################################
  source("code/diagnostics.R")
  BiLinT_write_diagnostics(
    foldername = foldername,
    dataname = dataname,
    Trace = Trace,
    MCMC_par = MCMC_par,
    beta_log_like = beta_log_like
  )

  ########################################################
  ########## get tree estimates ###############
  ########################################################
  tree_est <- tryCatch(
    tree_estimation(Trace),
    error = function(e) {
      warning(sprintf("tree_estimation skipped: %s", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(tree_est)) {
    ggtree(tree_est, size = 1,  layout = "circular", branch.length = "none") +
      geom_tiplab(size = 5) +
      geom_point2(aes(subset = (node %in% c(65:127))), size = 2) +
      theme_tree() + theme(legend.position = "none")
  }
}





