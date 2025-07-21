
#####################################################
########## sampling of mu #########################
#######################################################

samp_mu <- function(samp, i, obs_traits, model.OU.BM, tree_info, Params, temper)
{
  mu0 <- samp$mu[i]
  MU <- samp$mu
  Mu <- rnorm(1, mean = mu0, sd = Params$pro_sd)
  # while (Mu < 0) {
  #   Mu <- rnorm(1, mean = mu0, sd = Params$pro_sd)
  # }
  MU[i] <- Mu

  # lpriors <- dmvnorm(samp$mu, rep(0, Params$K), diag(rep(1, Params$K)), log = TRUE)/temper
  # lpriors2 <- dmvnorm(MU, rep(0, Params$K), diag(rep(1, Params$K)), log = TRUE)/temper
  lpriors <- sum(mclust::dmvnorm(samp$mu, 0, 1, log = TRUE)/temper)
  lpriors2 <- sum(mclust::dmvnorm(MU, 0, 1, log = TRUE)/temper)

  p0 <- samp$lik_traits/temper

  temp <- samp
  temp$mu <- MU

  model.OU.BM2 <- BM_OU_obj_change(temp, tree_info, model.OU.BM, Params)

  p1 <- lik_obs_traits(temp, tree_info, model.OU.BM2, obs_traits)
  temp$lik_traits <- p1
  temp$likelihood <- temp$lik_crispr + temp$lik_traits
  p1 <- p1/temper

  if(runif(n=1L, min=0, max=1) <= exp(p1 + lpriors2 - p0 - lpriors)){
    # print('mu changed')
    return(temp)
  }else{
    return(samp)
  }
}

samp_mu_all <- function(samp, obs_traits, model.OU.BM, tree_info, Params, temper){
  for(i in 1:length(samp$mu)){
    samp <- samp_mu(samp, i, obs_traits, model.OU.BM, tree_info, Params, temper)
  }
  return(samp)
}


#########################################################
#################### sampling of alpha ######################
#########################################################

samp_alpha <- function(samp, obs_traits, model.OU.BM, tree_info, Params, temper)
{
  alpha0 <- samp$alpha
  alpha <- rtnorm(1, mean = alpha0, sd = Params$pro_sd, lower = 0, upper = Inf)

  lpriors <- dexp(alpha0, rate = 2, log = TRUE)/temper
  lpriors2 <- dexp(alpha, rate = 2, log = TRUE)/temper

  temp <- samp
  temp$alpha <- alpha
  model.OU.BM2 <- BM_OU_obj_change(temp, tree_info, model.OU.BM, Params)

  p0 <- samp$lik_traits/temper
  p1 <- lik_obs_traits(temp, tree_info, model.OU.BM2, obs_traits)
  temp$lik_traits <- p1
  temp$likelihood <- temp$lik_crispr + temp$lik_traits
  p1 <- p1/temper

  if(runif(n=1L, min=0, max=1) <= exp(p1 + lpriors2 - p0 - lpriors)){
    # print('alpha changed')
    return(temp)
  }else{
    return(samp)
  }
}


#####################################################
########## sampling of theta #########################
#######################################################

samp_theta <- function(samp, i, obs_traits, model.OU.BM, tree_info, Params, temper)
{
  theta0 <- samp$theta[i]
  THETA <- samp$theta
  Theta <- rnorm(1, mean = theta0, sd = Params$pro_sd)
  THETA[i] <- Theta

  # lpriors <- dmvnorm(samp$theta, rep(0, Params$K), diag(rep(1, Params$K)), log = TRUE)/temper
  # lpriors2 <- dmvnorm(THETA, rep(0, Params$K), diag(rep(1, Params$K)), log = TRUE)/temper

  lpriors <- sum(mclust::dmvnorm(samp$theta, 0, 1, log = TRUE)/temper)
  lpriors2 <- sum(mclust::dmvnorm(THETA, 0, 1, log = TRUE)/temper)

  p0 <- samp$lik_traits/temper

  temp <- samp
  temp$theta <- THETA

  model.OU.BM2 <- BM_OU_obj_change(temp, tree_info, model.OU.BM, Params)

  p1 <- lik_obs_traits(temp, tree_info, model.OU.BM2, obs_traits)
  temp$lik_traits <- p1
  temp$likelihood <- temp$lik_crispr + temp$lik_traits
  p1 <- p1/temper

  if(runif(n=1L, min=0, max=1) <= exp(p1 + lpriors2 - p0 - lpriors)){
    # print('theta changed')
    return(temp)
  }else{
    return(samp)
  }
}

samp_theta_all <- function(samp, obs_traits, model.OU.BM, tree_info, Params, temper){
  for(i in 1:length(samp$theta)){
    samp <- samp_theta(samp, i, obs_traits, model.OU.BM, tree_info, Params, temper)
  }
  return(samp)
}



#####################################################
########## sampling of sigma #########################
#######################################################

samp_sigma <- function(samp, i, obs_traits, model.OU.BM, tree_info, Params, temper)
{
  sigma0 <- samp$sigma[i]
  # print(samp$sigma)
  SIGMA <- samp$sigma
  Sigma <- rtnorm(1, mean = sigma0, sd = Params$pro_sd, lower = 0, upper = Inf)
  SIGMA[i] <- Sigma
  # print(SIGMA)

  lpriors <-  LaplacesDemon::dhalft(sigma0, scale= 0.25, nu=1, log = TRUE)/temper
  lpriors2 <- LaplacesDemon::dhalft(Sigma, scale= 0.25, nu=1, log = TRUE)/temper

  p0 <- samp$lik_traits/temper
  # print(p0)

  temp <- samp
  temp$sigma <- SIGMA

  model.OU.BM2 <- BM_OU_obj_change(temp, tree_info, model.OU.BM, Params)
  p1 <- lik_obs_traits(temp, tree_info, model.OU.BM2, obs_traits)
  if(is.na(p1)){
    return(samp)
  }else{
    temp$lik_traits <- p1
    temp$likelihood <- temp$lik_crispr + temp$lik_traits
    p1 <- p1/temper
    # print(p1)
    if(runif(n=1L, min=0, max=1) <= exp(p1 + lpriors2 - p0 - lpriors)){
      # print('sigma changed')
      return(temp)
    }else{
      return(samp)
    }
  }
}

samp_sigma_all <- function(samp, obs_traits, model.OU.BM, tree_info, Params, temper){
  for(i in 1:length(samp$sigma)){
    samp <- samp_sigma(samp, i, obs_traits, model.OU.BM, tree_info, Params, temper)
  }
  return(samp)
}


#####################################################
########## sampling of S #########################
#######################################################

samp_S <- function(samp, i, obs_seq, state, Params, temper)
{
  S0 <- samp$S[i]
  SS <- samp$S
  S <- rtnorm(1, mean = S0, sd = Params$pro_sd, lower = 0, upper = Inf)
  SS[i] <- S

  Sc <- SS * samp$c
  r <- samp$l + sum(Sc)

  lpriors <- dexp(S0, rate = Params$S_exp, log = TRUE)/temper
  lpriors2 <- dexp(S, rate = Params$S_exp, log = TRUE)/temper

  p0 <- samp$lik_crispr/temper

  temp <- samp
  temp$S <- SS

  tree <- samp$tree
  tree_Cas9 <- tree_transform1(tree, t1, t2)
  tree_PostOrder <- c(reorder_tree_edges(tree_Cas9, root_to_tips = FALSE)$edge[,2], length(tree_Cas9$tip.label)+1)
  childlist <- ChildList(tree_Cas9, tree_PostOrder)
  dis <- get_all_distances_to_root(tree_Cas9, as_edge_count=FALSE)
  dis[which(dis==0)] <- 10^(-5)

  p1 <- tree_likehood(tree_Cas9, tree_PostOrder, dis, childlist, obs_seq, state, t1, t2, r, samp$l, Sc)
  # p1 <- tree_likehood_R(tree, obs_seq, state, t1, t2, r, samp$l, Sc)
  temp$lik_crispr <- p1
  temp$likelihood <- temp$lik_crispr + temp$lik_traits
  p1 <- p1/temper

  if(runif(n=1L, min=0, max=1) <= exp(p1 + lpriors2 - p0 - lpriors)){
    return(temp)
  }else{
    return(samp)
  }
}

samp_S_all <- function(samp, obs_seq, state, Params, temper){
  S0 <- samp$S
  SS <- S0
  for(i in 1:(length(S0)-1)){
    samp <- samp_S(samp, i, obs_seq, state, Params, temper)
  }
  return(samp)
}



#########################################################
#################### sampling of c ######################
#########################################################

samp_c <- function(samp, obs_seq, state, Params, temper)
{
  c0 <- samp$c
  c <- rtnorm(1, mean = c0, sd = Params$pro_sd, lower = 0, upper = Inf)

  Sc <- samp$S * c
  r <- samp$l + sum(Sc)

  lpriors <- sum(dlnorm(c0, meanlog = Params$c_meanlog, sdlog = Params$c_sdlog, log = TRUE))/temper
  lpriors2 <- sum(dlnorm(c, meanlog = Params$c_meanlog, sdlog = Params$c_sdlog, log = TRUE))/temper


  p0 <- samp$lik_crispr/temper
  temp <- samp
  temp$c <- c

  tree <- samp$tree
  tree_Cas9 <- tree_transform1(tree, t1, t2)
  tree_PostOrder <- c(reorder_tree_edges(tree_Cas9, root_to_tips = FALSE)$edge[,2], length(tree_Cas9$tip.label)+1)
  childlist <- ChildList(tree_Cas9, tree_PostOrder)
  dis <- get_all_distances_to_root(tree_Cas9, as_edge_count=FALSE)
  dis[which(dis==0)] <- 10^(-5)

  p1 <- tree_likehood(tree_Cas9, tree_PostOrder, dis, childlist, obs_seq, state, t1, t2, r, samp$l, Sc)
  temp$lik_crispr <- p1
  temp$likelihood <- temp$lik_crispr + temp$lik_traits
  p1 <- p1/temper


  if(runif(n=1L, min=0, max=1) <= exp(p1 + lpriors2 - p0 - lpriors)){
    return(temp)
  }else{
    return(samp)
  }
}




#########################################################
#################### sampling of l ######################
#########################################################

samp_l <- function(samp, obs_seq, state, Params, temper)
{
  l0 <- samp$l
  l <- rtnorm(1, mean = l0, sd = Params$pro_sd, lower = 0, upper = Inf)

  Sc <- samp$S * samp$c
  r <- l + sum(Sc)

  lpriors <- sum(dlnorm(l0, meanlog = Params$l_meanlog, sdlog = Params$l_sdlog, log = TRUE))/temper
  lpriors2 <- sum(dlnorm(l, meanlog = Params$l_meanlog, sdlog = Params$l_sdlog, log = TRUE))/temper

  p0 <- samp$lik_crispr/temper

  temp <- samp
  temp$l <- l

  tree <- samp$tree
  tree_Cas9 <- tree_transform1(tree, t1, t2)
  tree_PostOrder <- c(reorder_tree_edges(tree_Cas9, root_to_tips = FALSE)$edge[,2], length(tree_Cas9$tip.label)+1)
  childlist <- ChildList(tree_Cas9, tree_PostOrder)
  dis <- get_all_distances_to_root(tree_Cas9, as_edge_count=FALSE)
  dis[which(dis==0)] <- 10^(-5)

  p1 <- tree_likehood(tree_Cas9, tree_PostOrder, dis, childlist, obs_seq, state, t1, t2, r, l, Sc)
  temp$lik_crispr <- p1
  temp$likelihood <- temp$lik_crispr + temp$lik_traits
  p1 <- p1/temper

  if(runif(n=1L, min=0, max=1) <= exp(p1 + lpriors2 - p0 - lpriors)){
    return(temp)
  }else{
    return(samp)
  }
}


#####################################################
########## sampling of tree #########################
#######################################################

samp_tree <- function(samp, obs_seq, state, Params, obs_traits, temper){
  tree0 <- samp$tree
  tree <- tree_propose(tree0, t_tol, t2, prob_internal=NULL)
  # tree <- try(u <- tree_propose(tree0, t_tol, t2, prob_internal=NULL), silent=T)
  p0 <- (samp$lik_traits + samp$lik_crispr)/temper

  if(!is.null(tree)){
    # print('tree is not NULL')
    # print(tree$edge.length)
    if(length(which(tree$edge.length==0))>0){
      tree$edge.length[which(tree$edge.length==0)] <- 10^(-5)
    }

    tree_info <- tree_transform2(tree)

    lpriors <- -LikConstantn(lambda = Params$tree_brith_rate, mu = Params$tree_death_rate, sampling=Params$tree_sampling_rate,
                             x = getx(tree0),root = Params$M +1)/temper
    lpriors2 <- -LikConstantn(lambda = Params$tree_brith_rate, mu = Params$tree_death_rate, sampling=Params$tree_sampling_rate,
                              x = getx(tree),root = Params$M + 1)/temper

    temp <- samp
    temp$tree <- tree

    model.OU.BM2 <- BM_OU_obj_change(temp, tree_info, model.OU.BM, Params)
    lik_traits2 <- lik_obs_traits(temp, tree_info, model.OU.BM2, obs_traits)

    Sc <- samp$S * samp$c
    r <-  samp$l + sum(Sc)
    tree_Cas9 <- tree_transform1(tree, t1, t2)
    tree_PostOrder <- c(reorder_tree_edges(tree_Cas9, root_to_tips = FALSE)$edge[,2], length(tree_Cas9$tip.label)+1)
    childlist <- ChildList(tree_Cas9, tree_PostOrder)
    dis <- get_all_distances_to_root(tree_Cas9, as_edge_count=FALSE)
    dis[which(dis==0)] <- 10^(-5)

    lik_crispr2 <- tree_likehood(tree_Cas9, tree_PostOrder, dis, childlist, obs_seq, state, t1, t2, r, temp$l, Sc)
    p1 <- (lik_traits2 + lik_crispr2)/temper

    if(!is.na(lik_traits2)){
      temp$lik_traits <- lik_traits2
      temp$lik_crispr <- lik_crispr2
      temp$likelihood <- temp$lik_crispr + temp$lik_traits
    }

    # print(paste('p1: ', p1, 'lpriors2: ', lpriors2, 'p0: ', p0, 'lpriors: ', lpriors, sep=''))
    if(!is.na(p1) & runif(n=1L, min=0, max=1) <= exp(p1 + lpriors2 - p0 - lpriors)){
      return(temp)
    }else{
      return(samp)
    }
  }else{
    return(samp)
  }
}



tree_propose <- function(tree0, t_tol, t2, prob_internal=NULL){
  # tree0 <- TreeTools::Preorder(tree0)
  # tree0 <- samp$tree
  # if(!is.null(tree)){
  #   tree0 <- tree
  # }
  # tree0$node.label <- c(1:tree0$Nnode)+Ntip(tree0)
  dis <- get_all_distances_to_root(tree0, as_edge_count=FALSE)
  splitNode <- NULL
  for(i in 1:nrow(tree0$edge)){
    node1 <- tree0$edge[i,1]
    node2 <- tree0$edge[i,2]
    if(t2 > dis[node1] && t2 <= dis[node2]){
      splitNode <- c(splitNode, node2) # add point on the tree at time t2 one by one
    }
  }
  CrisprNode <- NULL
  for(i in 1:nrow(tree0$edge)){
    node1 <- tree0$edge[i,1]
    node2 <- tree0$edge[i,2]
    if(dis[node2] < t2){
      CrisprNode <- c(CrisprNode, node2)
    }
  }

  if (is.null(prob_internal)){
    prob_internal <- dunif(1:(tree0$Nnode + length(tree0$tip.label)),1,tree0$Nnode + length(tree0$tip.label))
  }

  ### the pruned subtree
  subtree_swap <- sample(tree0$Nnode + length(tree0$tip.label), 1, prob = prob_internal)
  tree0_root <- ape::Ntip(tree0) + 1
  tree0_leaves <- c(1:ape::Ntip(tree0))
  while(subtree_swap %in% c(splitNode, tree0_root, tree0_leaves)){
    subtree_swap <- sample(tree0$Nnode + length(tree0$tip.label), 1, prob = prob_internal)
  }

  splitnode_Ans <- intersect(Ancestors(tree0, subtree_swap), splitNode) ### the subtree root, i.e., one of the split node

  edge <- tree0$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  edge.length <- tree0$edge.length
  names(edge.length) <- child
  CrisprNodes <- child %in% CrisprNode
  edgeToBreak <- which(child == subtree_swap)
  nEdge <- length(parent)

  brokenEdge <- seq_along(parent) == edgeToBreak
  brokenEdge.parentNode <- parent[edgeToBreak] ## split node parent
  brokenEdge.childNode <- child[edgeToBreak] ## split node
  # edgesCutAdrift <- DescendantEdges(edgeToBreak, parent, child,  nEdge) ## subtree edges
  brokenEdge_Des <- c(getDescendants(tree0, brokenEdge.childNode), brokenEdge.childNode)
  edgesCutAdrift <- child %in% brokenEdge_Des
  edgesOnAdriftSegment <- edgesCutAdrift | brokenEdge
  brokenEdgeParent <- child == brokenEdge.parentNode  ## split node parent's parent
  brokenEdgeSister <- parent == brokenEdge.parentNode & !brokenEdge ## sibling node
  brokenEdgeDaughters <- parent == brokenEdge.childNode ## split node children
  nearBrokenEdge <- brokenEdge | brokenEdgeSister | brokenEdgeParent | brokenEdgeDaughters
  if (breakingRootEdge <- !any(brokenEdgeParent)) {
    brokenRootDaughters <- parent == child[brokenEdgeSister]
    nearBrokenEdge <- nearBrokenEdge | brokenRootDaughters
  }
  if(length(splitnode_Ans)>0){
    #edgesBreakSplitNodeDes <- DescendantEdges(which(child==splitnode_Ans), parent, child,  nEdge)
    edgesBreakSplitNode_Des <- c(getDescendants(tree0, splitnode_Ans), splitnode_Ans)
    edgesBreakSplitNodeDes <- child %in% edgesBreakSplitNode_Des
    mergeEdge <- which((CrisprNodes | edgesBreakSplitNodeDes) & !nearBrokenEdge & !edgesOnAdriftSegment)
  }else{
    mergeEdge <- which(CrisprNodes & !nearBrokenEdge & !edgesOnAdriftSegment)
  }
  if(length(mergeEdge) > 0){
    nCandidates <- length(mergeEdge)
    if (nCandidates > 1) {
      mergeEdge <- SampleOne(mergeEdge, len = nCandidates)
    }
    if(breakingRootEdge) {
      parent[brokenRootDaughters] <- brokenEdge.parentNode
      spareNode <- child[brokenEdgeSister]
      child[brokenEdgeSister] <- child[mergeEdge]
      parent[brokenEdge | brokenEdgeSister] <- spareNode
      child[mergeEdge] <- spareNode

      tree <- tree0
      tree$edge <- cbind(parent, child)
      # tree$node.label <- child
      # tree <- TreeTools::Preorder(tree) #nearBrokenEdge
      newEdgeLength <- tree$edge.length
      names(newEdgeLength) <- child

      oldEdgeLength <- tree0$edge.length
      names(oldEdgeLength) <- tree0$edge[,2]

      ### adjust the original sister subtree
      brokenEdge_siblingTree0 <- Siblings(tree0, brokenEdge.childNode, include.self = FALSE)
      brokenEdgeLength_siblingTree0 <- oldEdgeLength[as.character(brokenEdge_siblingTree0)]
      spareNode_child <- Descendants(tree0, spareNode, "child")
      newEdgeLength[as.character(spareNode_child)] <- newEdgeLength[as.character(spareNode_child)] + brokenEdgeLength_siblingTree0

      ## adjust the new sister subtree
      brokenEdge.childNode_Sibling <- Siblings(tree, brokenEdge.childNode, include.self = FALSE)
      newEdgeLength[as.character(brokenEdge.childNode_Sibling)] <- edge.length[as.character(brokenEdge.childNode_Sibling)]/2
      newEdgeLength[as.character(brokenEdge_siblingTree0)] <- newEdgeLength[as.character(brokenEdge.childNode_Sibling)]
      tree$edge.length <- newEdgeLength

      ## adjust the pruned subtree
      newSubtreeNodes <- c(brokenEdge_siblingTree0, brokenEdge.childNode, Descendants(tree, brokenEdge.childNode, 'all'))
      newSubtreeNodes_noroot <- setdiff(newSubtreeNodes, brokenEdge_siblingTree0)

      NewDis <- get_all_distances_to_root(tree, as_edge_count=FALSE)
      h_newSubtreeNodes <- NewDis[newSubtreeNodes]
      hh <- NewDis[brokenEdge_siblingTree0]
      newEdgeLength[as.character(newSubtreeNodes_noroot)] <- newEdgeLength[as.character(newSubtreeNodes_noroot)]/(max(h_newSubtreeNodes)-min(h_newSubtreeNodes))*(t_tol-hh)
      tree$edge.length <- newEdgeLength

    }else{
      parent[brokenEdgeSister] <- parent[brokenEdgeParent] ### original sister connected to ancestor on upper layer
      parent[brokenEdgeParent] <- parent[mergeEdge] ### connect the split node parent to the sampled node
      parent[mergeEdge] <- brokenEdge.parentNode ### make the sampled node's daughter as the split node's daughter

      tree <- tree0
      tree$edge <- cbind(parent, child)
      # tree$node.label <- child
      # tree <- TreeTools::Preorder(tree) #nearBrokenEdge
      newEdgeLength <- tree$edge.length
      names(newEdgeLength) <- child

      brokenEdge.childNode_Sibling <- Siblings(tree, brokenEdge.childNode, include.self = FALSE)
      brokenEdgeLength <- newEdgeLength[as.character(brokenEdge.parentNode)]
      newEdgeLength[as.character(brokenEdge.parentNode)] <- edge.length[as.character(brokenEdge.childNode_Sibling)]/2
      newEdgeLength[as.character(brokenEdge.childNode_Sibling)] <- edge.length[as.character(brokenEdge.childNode_Sibling)]/2
      tree$edge.length <- newEdgeLength

      ### adjust the edge length of brokenEdgeSister
      brokenEdge_siblingTree0 <- Siblings(tree0, brokenEdge.childNode, include.self = FALSE)
      newEdgeLength[as.character(brokenEdge_siblingTree0)] <- brokenEdgeLength + edge.length[brokenEdgeSister]

      newSubtreeRoot <- parent[mergeEdge]
      newSubtreeNodes <- c(parent[brokenEdgeParent], newSubtreeRoot, Descendants(tree, newSubtreeRoot, 'all'))
      newSubtreeNodes_noroot <- c(newSubtreeRoot, Descendants(tree, newSubtreeRoot, 'all'))

      ## adjust the new whole subtree length
      NewDis <- get_all_distances_to_root(tree, as_edge_count=FALSE)
      h_newSubtreeNodes <- NewDis[newSubtreeNodes]
      h <- dis[parent[brokenEdgeParent]]
      newEdgeLength[as.character(newSubtreeNodes_noroot)] <- newEdgeLength[as.character(newSubtreeNodes_noroot)]/(max(h_newSubtreeNodes)-min(h_newSubtreeNodes))*(t_tol-h)
      tree$edge.length <- newEdgeLength

      ## adjust the subtree length
      NewDis2 <- get_all_distances_to_root(tree, as_edge_count=FALSE)
      brokenEdge.childNode_Siblings <- c(brokenEdge.childNode_Sibling, Descendants(tree, brokenEdge.childNode_Sibling, 'all'))
      h <- NewDis2[newSubtreeRoot]
      h1 <- NewDis2[brokenEdge.childNode_Siblings]
      brokenEdge.childNode_Des <- unique(c(brokenEdge.childNode, Descendants(tree, brokenEdge.childNode, 'all')))
      h2 <- NewDis2[brokenEdge.childNode_Des]
      if(max(h1)<(t_tol-10^(-6))){
        newEdgeLength[as.character(brokenEdge.childNode_Siblings)] <- newEdgeLength[as.character(brokenEdge.childNode_Siblings)]/(max(h1)-h)*(t_tol-h)
        tree$edge.length <- newEdgeLength
      }
      if(max(h2)<(t_tol-10^(-6))){
        newEdgeLength[as.character(brokenEdge.childNode_Des)] <- newEdgeLength[as.character(brokenEdge.childNode_Des)]/(max(h2)-h)*(t_tol-h)
        tree$edge.length <- newEdgeLength
      }
    }
    names(tree$edge.length) <- ""
    tree <- TreeTools::Postorder(tree)
    tree <- TreeTools::Preorder(tree)
    #ggtree(tree) + geom_text2(aes(label=node),hjust=-.3,color="red")
  }else{
    tree <- NULL
    # print('tree not changed')
  }
  return(tree)
}




#####################################################
########## all prior #########################
#######################################################
log_prior_all <- function(samp, Params, temper){
  S_lpriors <- sum(dexp(samp$S[length(samp$S)-1], rate = Params$S_exp, log = TRUE))
  l_lpriors <- sum(dlnorm(samp$l, meanlog = Params$l_meanlog, sdlog = Params$l_sdlog, log = TRUE))
  c_lpriors <- sum(dlnorm(samp$c, meanlog = Params$c_meanlog, sdlog = Params$c_sdlog, log = TRUE))
  tree_lpriors <- -LikConstantn(lambda = Params$tree_brith_rate, mu = Params$tree_death_rate, sampling=Params$tree_sampling_rate,
                                x = getx(samp$tree),root = Params$M + 1)
  alpha_lpriors <- dexp(samp$alpha, rate = 2, log = TRUE)
  # mu_lpriors <- dmvnorm(samp$mu, rep(0, Params$K), diag(rep(1, Params$K)), log = TRUE)
  # theta_lpriors <- dmvnorm(samp$theta, rep(0, Params$K), diag(rep(1, Params$K)), log = TRUE)
  mu_lpriors <- sum(mclust::dmvnorm(samp$mu, 0, 1, log = TRUE))
  theta_lpriors <- sum(mclust::dmvnorm(samp$theta, 0, 1, log = TRUE))
  sigma_lpriors <- sum(LaplacesDemon::dhalft(samp$sigma, scale= 0.25, nu=1, log = TRUE))
  all_lpriors <- S_lpriors + l_lpriors + c_lpriors + tree_lpriors + alpha_lpriors +
    mu_lpriors + theta_lpriors + sigma_lpriors
  return(all_lpriors/temper)
}

