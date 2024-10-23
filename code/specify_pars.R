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

MCMC_par$burnin <- 1000  # burnin sample size
MCMC_par$Nsamp <- 1000  # number of samples for inference
# MCMC_par$Ntune <- 300  # number of samples used for adaptive parameter tuning

MCMC_par$swap_interval <- 30 # make Matroplis Hastings move in every how many samples
MCMC_par$Nchain <- 4 # number of paralel chains
MCMC_par$delta_T <- 0.35

Temperature <- seq(1, by=MCMC_par$delta_T, length.out = MCMC_par$Nchain)  # temperatures


