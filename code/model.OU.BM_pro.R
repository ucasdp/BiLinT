
BM_OU_obj <- function(temp, tree_info, Params){
  ## tree: the original tree
  tree <- temp$tree
  # tree_info <- tree_transform2(tree)
  tree.ab <- tree_info$tree
  regime <- tree_info$regime
  N <- ape::Ntip(tree)
  K <- Params$K
  modeltypes <- c("BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
                  "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x")
  mapping_vec <- c(1, 2, 1)
  names(mapping_vec) <- regime

  model.OU.BM <- MixedGaussian( k = K,
                                modelTypes = modeltypes,
                                className = 'MixedGaussian_BMOU',
                                mapping = mapping_vec, Sigmae_x = structure(
                                  0,
                                  class = c("MatrixParameter", "_Omitted", description = "upper triangular factor of the non-phylogenetic variance-covariance")))
  model.OU.BM <- PCMApplyTransformation(model.OU.BM)
  model.OU.BM
}


BM_OU_obj_change <- function(temp, tree_info, model.OU.BM, Params){
  tree <- temp$tree
  N <- ape::Ntip(tree)
  K <- Params$K
  # tree_info <- tree_transform2(tree)
  tree.ab <- tree_info$tree
  model.OU.BM[["X0"]][] <- temp$mu
  model.OU.BM[[2]]$Sigma_x[,,1] <- UpperTriFactor(diag(rep(temp$sigma[1], K)))
  # selection strength matrix for the OU process
  alpha <- temp$alpha
  model.OU.BM[[3]]$H[,,1] <- diag(rep(alpha, K))
  # long term optimum vector for the OU process:
  model.OU.BM[[3]]$Theta[] <- temp$theta
  # random drift matrix for the OU process
  # V <- matrix(model_para[grep("V\\[", names(model_para))], K, K)
  V <- diag(rep(temp$sigma[2], K))
  # R <- 2 * alpha * V
  model.OU.BM[[3]]$Sigma_x[,,1][] <- UpperTriFactor(V)

  model.OU.BM[[4]]$Sigma_x[,,1] <- UpperTriFactor(diag(rep(temp$sigma[3], K)))

  names(model.OU.BM)[2:length(model.OU.BM)] <- tree_info$regime
  model.OU.BM
}


lik_obs_traits <- function(temp, tree_info, model.OU.BM, obs_traits){
  #tree <- temp$tree
  # tree_info <- tree_transform2(tree)
  tree.ab <- tree_info$tree
  PCMBase_lh <- PCMLik(t(obs_traits), tree.ab, model.OU.BM, metaI = PCMInfoCpp)
  return(PCMBase_lh)
}

