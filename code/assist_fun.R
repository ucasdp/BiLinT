##################################################################
############## get samples for given parameters ##################
##################################################################

get_samp <- function(MCMCout,var_name,subset=NULL)
{
  #subset: a subset of samples

  if(is.null(subset)){subset <- c(1:length(MCMCout))}
  MCMCout <- MCMCout[subset]

  Ns <- length(MCMCout)  # number of samples
  temp <- MCMCout[[1]][[var_name]]
  if(is.vector(temp)){temp <- matrix(temp,nrow = 1)}
  if(is.matrix(temp)){
    if(ncol(temp)==1){
      temp <- t(temp)
    }
  }

  d1 <- nrow(temp)
  d2 <- ncol(temp)

  if(is.null(d1))# no row dimension,then a scalar
  {
    out <- c(1:Ns)
    for(z in 1:Ns)
    {
      temp <- MCMCout[[z]][[var_name]]
      out[z] <- temp
    }
  } else if(d1==1)  # if parameter is vector, get samples as a matrix
  {
    out <- matrix(0,Ns,d2)
    for(z in 1:Ns)
    {
      temp <- MCMCout[[z]][[var_name]]
      out[z,] <- temp
    }
  } else {  # else get samples as an array
    out <- array(0,dim=c(d1,d2,Ns))
    for(z in 1:Ns)
    {
      temp <- MCMCout[[z]][[var_name]]
      out[,,z] <- temp
    }
  }
  out
}


get_samp_tree <- function(MCMCout,subset=NULL){
  if(is.null(subset)){subset <- c(1:length(MCMCout))}
  MCMCout <- MCMCout[subset]

  Ns <- length(MCMCout)  # number of samples
  out <- list()
  for(i in 1:Ns){
    out <- c(out, list(MCMCout[[Ns]]$tree))
  }
  return(out)
}

tree_estimate <- function(Trace){
  temp_traits <- get_samp(Trace,'lik_traits')
  temp_crispr <- get_samp(Trace,'lik_crispr')
  temp <- temp_traits + temp_crispr
  tree <- Trace[[which.max(temp)]]$tree
  tree$edge <- tree$edge[order(tree$edge[,1], tree$edge[,2], decreasing = TRUE), ]
  tree$edge.length <- tree$edge.length[order(tree$edge[,1], tree$edge[,2], decreasing = TRUE)]
  tree
}





