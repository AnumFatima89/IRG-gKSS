library(igraph)
library(doParallel)
library(foreach)
library(Matrix)
library(reticulate)

##do not use if anaconda is working ## use_python("C:/ProgramData/anaconda3/python.exe")
#conda_create("r-reticulate-env")
# Use your conda environment where grakel is installed
#use_condaenv("r-reticulate-env", required = TRUE)

# Alternative (newer version)
Sys.setenv(RETICULATE_PYTHON = "C:/Users/fatima/grakel-env/Scripts/python.exe")

# Import GraKeL and other Python modules
#conda_install("r-reticulate-env", "grakel")
grakel <- import("grakel")

################################################################################
################# Test statistic and related GOF tests ########################
################################################################################

# test_G: An igraph object containing the observed network
# C:  A vector of group labels for all the vertices
# p: The edge probability matrix for the null model
# g.kernel: Kernel to use to compute the test statistic
# alpha: level of significance for the test 
# M: Number of iterations in Monte Carlo setting to compute the null set for test statistic

generate.one.gKSS_grakel=function(G, C, p, g.kernel = wl_kernel){
  V(G)$name <- seq(0,gorder(G)-1,1)
  
  X = as_adjacency_matrix(G, type = "both")
  S = X[upper.tri(X)]
  n = length(S)
  S.t =  p[upper.tri(p)]
  S.t.vec = abs(S - S.t)
  S.mat = S.t.vec %*% t(S.t.vec)
  
  P = compute.transition.list_grakel(X, C)
  K = compute.kernel_grakel(P, g.kernel)
  
  J.kernel = S.mat * K
  
  stats.value = mean(J.kernel)
  
  #Return:
  #stats.value: gKSS^2
  #J.kernel: J kernel matrix for wild bootstrap 
  list(stats.value=stats.value,J.kernel=J.kernel, K=K, S=S.mat)
}


compute.transition.list_grakel = function(X,C){
  P=list()
  counter = 1
  for (i in 2:nrow(X)) {
    for (j in 1:(i - 1)) {
      Y = X
      Y[i,j] = abs(1-Y[i,j])
      Y[j,i] = Y[i,j]
      g = graph_from_adjacency_matrix(Y, mode = "undirected")
      V(g)$label = as.character(C)
      V(g)$name = as.numeric(V(g)$name)
      P[[counter]] = convert_to_grakel_format(g)
      counter = counter + 1
    }
  }
  g = graph_from_adjacency_matrix(X, mode = "undirected")
  V(g)$label = as.character(C)
  V(g)$name = as.numeric(V(g)$name)
  P[[counter]] = convert_to_grakel_format(g)
  return(P)
}


method <- grakel$WeisfeilerLehman
base <- grakel$VertexHistogram

wl_kernel <- method(n_iter = 3L)
wl_kernel$base_kernel <- base()

compute.kernel_grakel = function(P, g.kernel = wl_kernel){
  n = length(P) - 1 
  # Fit & transform to compute similarity matrix
  kernel.matrix <- g.kernel$fit_transform(P)
  #kernel.matrix = g.kernel(P, level)
  rm(P)
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K - outer(K.vec, K.vec, function(x,y)x+y)
  return(K)
}

convert_to_grakel_format <- function(g) {
  # Get the edge list as a list of tuples (convert to integer pairs)
  edges <- as_edgelist(g)
  edge_tuples <- lapply(1:nrow(edges), function(i) tuple(as.integer(edges[i, 1]), as.integer(edges[i, 2])))
  
  
  # Get the node labels as a dictionary (ensure labels are set as node names)
  labels <- setNames(V(g)$label, as.character(V(g)$name))
  
  # Ensure labels are passed as a dictionary: node id -> label
  label_dict <- dict(as.list(labels))
  
  # Return as tuple (edges, labels)
  return(tuple(list(edge_tuples, label_dict)))
}
