# GMMB-Lineage
GMMB-Lineage: A Guassian model-based Bayesian lineage tree inference with multi-modal traits

## Software dependencies

GMMB-Lineage is written with `R` and `C++`. Before implementing our software, please install the following packages in `R`:

 `tidyr`, `dplyr`, `coda`, `Rcpp`, `RcppArmadillo`, `ggplot2`, `PCMBase`, `ape`, `ggtree`, `castor`, `phangorn`, `TreeSearch`, `TreePar`, `TreeSim`, `msm`, `Rfast`,  `stringr`, `treeio`, `TreePar`, `ggstance`, `TreeTools`, `phytools`, `PCMBaseCpp`, `rdist`, `foreach`, `doParallel`, `Rcpp2doParallel`





## Usage

To use GMMB-Lineage, please set `R` working directory to `GMMB-Lineage-master` after downloading this repository. Please make sure you have installed all the dependencies correctly, and then open source code `GMMB-Lineage_app.R` to execute the commands line by line as following.

* In *Model Input* section, first set random seed `myseed`, and then specify the folder `foldername` to save output files (a new folder will be created if it does not exist). For example: 
  ```
  #################################################################
  ###################### MODEL INPUT ##############################
  #################################################################

  myseed <-  1               # set random seed
  foldername <-  "temp_out"          # set output foldername
  dir.create(foldername)  # folder where outputs are saved
  ```

  Then input the barcode matrix, `seq_data_mt`, the gene expression matrix `obs_traits_original` and the scaring window for barocde data, including `t1`, `t2` and `t_tol`. For the given example data:
  
  ```
  muti_data <- readRDS('toy_model_data.rds')
  seq_data_mt <- muti_data$seq_data_mt # barcode data, where row represents cells, columns represent targets
  M <- dim(muti_data$seq_data_mt)[1]
  N <- dim(muti_data$seq_data_mt)[2]

  obs_traits_original <- muti_data$obs_traits # normalized gene expression data, where row represents cells, columns represent genes
  obs_traits <- prcomp(obs_traits_original, rank. = 5)$x
  K <- dim(obs_traits)[2]
  state <- sort(setdiff(unique(as.vector(seq_data_mt)), c(0, -1)))

  t1 <- 0 # scaring start time
  t2 <- muti_data$seq_para_data$t2 # scaring end time
  t_tol <- muti_data$seq_para_data$t_tol # total experiment time
  ```
   

* Next, assign Bayesian sampling parameters in `specify_pars.R`. For most of the parameters, BiTSC2 works just fine with default values. Some of the parameters you can change are:
  ```
  MCMC_par$burnin <- 500  # burnin sample size
  MCMC_par$Nsamp <- 500   # number of samples for inference
  ```

* After that, execute `sampler.R` to perform MCMC sampling, then the samples used for inference are stored in the `.Rdata` files;
* finally, get the tree estimation by
  ```
  tree_est <- tree_estimate(Trace)
  ```
