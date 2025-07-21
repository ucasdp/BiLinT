
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
foldername <-  "temp_out"          # set output foldername
dataname <- 'toy_model'
dir.create(foldername)  # folder where outputs are saved

muti_data <- readRDS('toy_model_data.rds')
seq_data_mt <- muti_data$seq_data_mt ### barcode data, where row represents cells, columns represent targets
M <- dim(muti_data$seq_data_mt)[1]
N <- dim(muti_data$seq_data_mt)[2]

obs_traits_original <- muti_data$obs_traits ### normalized gene expression data, where row represents cells, columns represent genes
obs_traits <- prcomp(obs_traits_original, rank. = 5)$x
K <- dim(obs_traits)[2]
state <- sort(setdiff(unique(as.vector(seq_data_mt)), c(0, -1)))

t1 <- 0  # scaring start time
t2 <- muti_data$seq_para_data$t2 # scaring end time
t_tol <- muti_data$seq_para_data$t_tol # total experiment time

##############################################
######## load parameter file ################
############################################
source('specify_pars.R')


##############################################
######## sampling ##########################
############################################
source("sampler.R")

########################################################
########## get tree estimates ###############
########################################################
tree_est <- tree_estimation(Trace)

ggtree(tree_est, size = 1,  layout = "circular",branch.length = "none") +
  geom_tiplab(size=5) +
  geom_point2(aes(subset=(node %in% c(65:127))),size=2)+
  theme_tree()+ theme(legend.position = "none")





