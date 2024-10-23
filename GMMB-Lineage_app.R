
rm(list=ls())
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

source('Crispr_Cas9_cpplik.R')
source('Crispr_Cas9_lik.R')
source('main_fun.R')
source('par_samp_para.R')
source('assist_fun.R')
source('model.OU.BM_pro.R')



#################################################################
########## MODEL INPUT ############################################
#################################################################

myseed <-  1               # set random seed
foldername <-  "temp_out"          # set output foldername
dir.create(foldername)  # folder where outputs are saved

muti_data <- readRDS('toy_model_data.rds')
seq_data_mt <- muti_data$seq_data_mt ### barcode data, where row represents cells, columns represent targets
M <- dim(muti_data$seq_data_mt)[1]
N <- dim(muti_data$seq_data_mt)[2]

obs_traits_original <- muti_data$obs_traits ### normalized gene expression data, where row represents cells, columns represent genes
obs_traits <- prcomp(obs_traits_original, rank. = 5)$x
K <- dim(obs_traits)[2]
state <- sort(setdiff(unique(as.vector(seq_data_mt)), c(0, -1)))

t1 <- 0
t2 <- muti_data$seq_para_data$t2
t_tol <- muti_data$seq_para_data$t_tol

##############################################
######## load parameter file ################
############################################
myseed <- 1
source('specify_pars.R')


##############################################
######## sampling ##########################
############################################
source("sampler.R")



########################################################
########## get tree estimates ###############
########################################################
tree_est <- tree_estimate(Trace)








