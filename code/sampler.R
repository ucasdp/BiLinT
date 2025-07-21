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
  New_par[[i]] <- init_gen_fun(Params)
}

tree_info <- tree_transform2(New_par[[1]]$tree)
model.OU.BM <- BM_OU_obj(New_par[[1]], tree_info, Params)

for(i in 1:MCMC_par$Nchain){
  l <- New_par[[i]]$l
  Sc <- New_par[[i]]$S * New_par[[i]]$c
  r <- l + sum(Sc)

  tree <- New_par[[i]]$tree
  tree_Cas9 <- tree_transform1(tree, t1, t2)
  tree_PostOrder <- c(reorder_tree_edges(tree_Cas9, root_to_tips = FALSE)$edge[,2], length(tree_Cas9$tip.label)+1)
  childlist <- ChildList(tree_Cas9, tree_PostOrder)
  dis <- get_all_distances_to_root(tree_Cas9, as_edge_count=FALSE)
  dis[which(dis==0)] <- 10^(-5)

  lik_crispr <- tree_likehood(tree_Cas9, tree_PostOrder, dis, childlist, seq_data_mt, state, t1, t2, r, l, Sc)
  while(is.infinite(lik_crispr) && lik_crispr < 0){
    New_par[[i]] <- init_gen_fun(Params)
    l <- New_par[[i]]$l
    Sc <- New_par[[i]]$S * New_par[[i]]$c
    r <- l + sum(Sc)

    tree <- New_par[[i]]$tree
    tree_Cas9 <- tree_transform1(tree, t1, t2)
    tree_PostOrder <- c(reorder_tree_edges(tree_Cas9, root_to_tips = FALSE)$edge[,2], length(tree_Cas9$tip.label)+1)
    childlist <- ChildList(tree_Cas9, tree_PostOrder)
    dis <- get_all_distances_to_root(tree_Cas9, as_edge_count=FALSE)
    dis[which(dis==0)] <- 10^(-5)
    lik_crispr <- tree_likehood(tree_Cas9, tree_PostOrder, dis, childlist, seq_data_mt, state, t1, t2, r, l, Sc)
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

Nrep <- MCMC_par$burnin+ MCMC_par$Nsamp #+ MCMC_par$Ntune # number of total samples

tune_samp <- vector('list',MCMC_par$Nchain)
Trace <- vector('list',MCMC_par$Nsamp)
beta_log_like <- matrix(0,MCMC_par$Nchain ,MCMC_par$Nsamp)  # keep beta chain likelihood for WBIC
# Trace <- vector('list',Nrep)
# beta_log_like <- matrix(0,MCMC_par$Nchain, Nrep)

# stanmodel <- stan_model(file = 'OU_BM_sampling2.stan')
cl <- makeCluster(4) # Use 4 cores, adjust as needed
registerDoParallel(cl)
sawp_time <- 0
for(h in 1:Nrep){
  print_red(paste('the current iteration count:', h, sep=' '))

  ######## perform MCMC on each chahn ##########
  # for (z in 1:MCMC_par$Nchain){
  #   temp <- MCMC_onesamp_II(h, Params, seq_data_mt, New_par[[z]], obs_traits, temper=Temperature[z])
  #   New_par[[z]] <- temp
  # }
  # export_vars <- c("Params","seq_data_mt", "New_par", "obs_traits", "Temperature", 't1', 't2')
  New_par <- foreach(i=1:4,.packages=c('ape', 'phangorn', 'PCMBase', 'ggtree', 'mclust', 'Rcpp', 'PCMBaseCpp', 'msm',
                                       'castor', 'Rcpp2doParallel', 'TreeTools', 'TreePar'),.noexport = 'tree_likehood') %do% {
                                         if (!exists("tree_likehood")) {
                                           source("Crispr_Cas9_cpplik.R")
                                         }
                                         MCMC_onesamp_II(h, Params, seq_data_mt, New_par[[i]], obs_traits, Temperature[i])
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
        message('samples switched')
      }
    }
  }

  ######## after burnin stage: save new samples into Trace ###############
  #### only keep chain 1
  ct <- h - MCMC_par$burnin #- MCMC_par$Ntune  # count number of saved samples
  if(ct > 0) {
    Trace[[ct]] <- New_par[[1]]
    for(z in 1:length(New_par)){
      beta_log_like[z,ct] <-  New_par[[z]]$likelihood # calculate Ln for all chains
    }
  }


  ########## display progress ###############

  if( h %% 100 == 0){
    ave_speed <- difftime(Sys.time(),start_t,units = "mins")/h
    cat("accomplished",100*h/(MCMC_par$Nsamp+MCMC_par$burnin), #+MCMC_par$Ntune
        "%; time remaining: ",(MCMC_par$Nsamp+MCMC_par$burnin-h)*ave_speed,"mins \n") #+MCMC_par$Ntune
  }
}

stopCluster(cl)

cat('completed. Time consumed:',difftime(Sys.time(),start_t,units = "mins"),'mins \n \n')




##################################################################
########### save data ############################################
##################################################################

# cur_file=paste(foldername,'/seed',myseed,'_N', N, '_M', M, '_K', K, '_', rep, '.Rdata',sep='')
cur_file=paste(foldername,'/', dataname, '.Rdata',sep='')
save(Params, MCMC_par, obs_traits_original, Trace, myseed, beta_log_like, file = cur_file)
