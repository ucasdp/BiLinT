##################################################
##### parameter random initiation ###############
##################################################
init_gen_fun <- function(Params){

  out <- list()

  M <- Params$M
  N <- Params$N
  K <- Params$K

  out$alpha <- rexp(1, 2)
  out$mu <- as.vector(rmvnorm(1, rep(0, K), diag(rep(1, K))))
  out$theta <- as.vector(rmvnorm(1, rep(0, K), diag(rep(1, K))))

  out$sigma <- LaplacesDemon::rhalft(3, scale = 0.25, nu=1)

  # initial scaring rate
  out$S <- rexp(Params$St, rate = 2)
  out$S[length(out$S)] <- 1

  # initial loss rate
  out$l <- rlnorm(1, meanlog = Params$l_meanlog, sdlog = Params$l_sdlog)

  # initial clock rate
  out$c <- rlnorm(1, meanlog = Params$c_meanlog, sdlog = Params$c_sdlog)

######### using the cluster tree to initial tree

  seq_data_cluster <- unique(seq_data_mt)
  rownames(seq_data_cluster) <- c(1:nrow(seq_data_cluster))
  dis_cluster <- pegas::dist.hamming(seq_data_cluster)

  # Create a hierarchical clustering object using UPGMA
  hc <- hclust(as.dist(dis_cluster), method = "average")  # "average" is UPGMA
  # Plot the tree
  plot(as.dendrogram(hc))

  cluster_tree <- as.phylo(hc)
  cluster_tree$edge.length <- cluster_tree$edge.length / max(diag(vcv(cluster_tree))) * t2

  # RNA_data <- muti_data$obs_traits
  RNA_data <- obs_traits_original

  joint_tree <- cluster_tree
  ggtree(cluster_tree) + geom_text2(aes(label=label),hjust=-.3,color="red")

  for(i in 1:nrow(seq_data_cluster)){
    cells_icluster <- NULL

    for(j in 1:nrow(seq_data_mt)){
      if(sum(abs(seq_data_mt[j, ]- seq_data_cluster[i,]))==0){
        cells_icluster <- c(cells_icluster, j)
      }
    }

    if(length(cells_icluster)==1){
      RNA_tree <- ape::rtree(1)
      RNA_tree$edge.length[1] <- (t_tol - t2)
      RNA_tree$tip.label <- paste('cell', cells_icluster, sep='')
      joint_tree <- bind.tree(joint_tree, RNA_tree, where = which(joint_tree$tip.label==i))
    }else{
      RNA_ibranch <- obs_traits_original[cells_icluster, ]
      rownames(RNA_ibranch) <- cells_icluster
      dis_RNA <- rdist::pdist(RNA_ibranch, 'correlation')
      dis_RNA <- dis_RNA * dim(RNA_data)[2]
      rownames(dis_RNA) <- paste('cell', cells_icluster, sep='')
      colnames(dis_RNA) <- paste('cell', cells_icluster, sep='')
      RNA_tree <- hclust(as.dist(dis_RNA), method = "average")
      RNA_tree <- as.phylo(RNA_tree)
      RNA_tree <- compute.brlen(RNA_tree, method = "equal")
      RNA_tree$edge.length <- RNA_tree$edge.length / max(diag(vcv(RNA_tree))) * (t_tol - t2)
      joint_tree <- bind.tree(joint_tree, RNA_tree, where = which(joint_tree$tip.label==i))
      # ggtree(joint_tree) + geom_text2(aes(label=node),hjust=-.3,color="red")
    }
  }

  leaves_label <- joint_tree$tip.label
  real_cell_bumber <- as.numeric(gsub("\\D", "", joint_tree$tip.label))
  real_tree_edge <- joint_tree$edge
  Nleaves <- length(joint_tree$tip.label)
  for(i in 1:nrow(real_tree_edge)){
    if(real_tree_edge[i,2] <= Nleaves){
      real_tree_edge[i, 2] <- real_cell_bumber[real_tree_edge[i, 2]]
    }
  }
  joint_tree$edge <- real_tree_edge
  joint_tree$tip.label <- paste('cell', c(1:Nleaves), sep='')
  out$tree <- joint_tree

  out
}


RNA_dis <- function(obs_traits){
  d1 <- dist(obs_traits , method="euclidean")
  return(d1)
}



#######################################################
########## acquire one MCMC sample ####################
#######################################################

MCMC_onesamp_II <- function(h, Params, obs_seq, init, obs_traits, temper=1)
{ # init: last round MCMC sample
  # temper: temperature
  # adap: adaptive tuning parameters

  # init <- New_par[[z]]
  out <- init

  tree_info <- tree_transform2(init$tree)

  ########## update mu #######
  out <- samp_mu_all(out, obs_traits, model.OU.BM, tree_info, Params, temper)
  # print('complete mu sampling!')

  ########## update alpha #######
  out <- samp_alpha(out, obs_traits, model.OU.BM, tree_info, Params, temper)
  # print('complete alpha sampling!')

  ########## update theta #######
  out <- samp_theta_all(out, obs_traits, model.OU.BM, tree_info, Params, temper)
  # print('complete theta sampling!')

  ########## update sigma #######
  out <- samp_sigma_all(out, obs_traits, model.OU.BM, tree_info, Params, temper)
  # print('complete sigma sampling!')

  ########## update S #######

  out <- samp_S_all(out, obs_seq, state, Params, temper)
  # print('complete S sampling!')

  ########## update c #######
  out <- samp_c(out, obs_seq, state, Params, temper)
  # print('complete c sampling!')

  ########## update l #######
  out <- samp_l(out, obs_seq, state, Params, temper)
  # print('complete l sampling!')

  ########## update tree #######
  # if(h%%20==0){
  #
  # }
  out <- samp_tree(out, obs_seq, state, Params, obs_traits, temper)
  # print('complete tree sampling!')

  out
}
