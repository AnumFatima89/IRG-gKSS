library(igraph)
library(netrankr)
library(network)
library(graphkernels)
library("randnet")
library(arrangements)
library(HelpersMG)
library(igraphdata)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(doParallel)
library(foreach)
library(qgraph)
library(sjmisc)
library(boot)
{
#install.packages("remotes")
#remotes::install_github("LTLA/bluster")
#library(bluster)
}
{
library(devtools)
library(remotes)
library(GoodFitSBM)
#remotes::install_github("Roy-SR-007/GoodFitSBM")
}
library(ggrepel)
library(digest)
library(Matrix)
library(tictoc)

####Generating an ERMM with vertex groups####

MOAR_LETTERS <- function(n=2) {
  n <- as.integer(n[1L])
  if(!is.finite(n) || !is.integer(n))
    stop("'n' must be a length-1 integer")
  if(n == 1) return(LETTERS)
  
  res <- vector("list", n)
  res[[1]] <- LETTERS
  for(i in 2:n)
    res[[i]] <- c(sapply(res[[i-1L]], function(y) paste0(y, LETTERS)))
  
  unlist(res)
}
ml <- MOAR_LETTERS(4)

################################################################################
########################### Simulating IRGs ####################################
################################################################################

#input
# p : edge probability matrix

sample_IRG = function(p, vlab = MOAR_LETTERS(3)){
  n = nrow(p)
  
  # Creating the n x n adjacency matrix  
  adj <- matrix(0, n, n)
  
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      adj[i, j] = rbinom(1, 1, p[i,j]) # We include the edge with probability p
    }
  }
  
  adjsymm = symmetricize(adj, method = "ld")
  
  # graph from the adjacency matrix
  G = igraph::graph_from_adjacency_matrix(adjsymm, mode = "undirected", weighted = NULL)
  
  # assigning vertex attributes
  #V(G)$name <- vlab[1:n]
  
  return(G)
}

############# Simulating ER with a planted clique #################

network_with_clique_nogroup = function(G, K){
  G_pc = G
  smpl = sample(1:gorder(G_pc),K)
  new_edges = combinations(smpl, k = 2, layout = "column")
  ap = as_ids(E(G_pc)[new_edges[1,] %--% new_edges[2,]])
  to_del = c(1:gsize(G_pc))
  to_del = to_del[!to_del %in% ap]
  del_eid = sample(to_del,ncol(new_edges)-length(ap))
  G_pc = delete_edges(G_pc, del_eid)
  G_pc <- G_pc + edges(as.vector(new_edges))
  G_pc = simplify(G_pc, remove.multiple = TRUE, remove.loops = TRUE)
  return(G_pc)
}

with_clique = function(p_0, K){
  mismatches = 0
  net_rejected = 0
  
  if(K < 2) {
    G_K = sample_IRG(p_0)
    mismatches = 0
  } else{
    counter <- 0
    cond <- FALSE
    while (!cond && counter < 30) {
      counter = counter +1
      G = sample_IRG(p_0)
      cond = gsize(G) >= choose(K,2)
    }
    net_rejected = ifelse(counter >= 30, 30, counter - 1)
    G_K = network_with_clique_nogroup(G, K)
    if(gsize(G_K) != gsize(G)) mismatches = mismatches + 1
  }
  list("Graph" = G_K, "tries" = net_rejected, "mismatches" = mismatches)
}

##########Simulating ERMM with cliques#################

#simulate a network with a clique of size K number of edges equal to the ERMM with given edge probabilities
# K must be at least 2 or greater

network_with_clique = function(G, C, K, max.rep = 20){
  G_pc = G
  E(G_pc)$label = edge.labels(G, C)
  counter = 0
  repeat{
    counter = counter + 1
    smpl = sample(1:gorder(G_pc),K)
    new_edges = combinations(smpl, k = 2, layout = "column")
    
    L= length(unique(C))
    edge.types = g_labels(L)
    new.edge.labels = c()
    for (i in 1:ncol(new_edges)) {
      new.edge.labels[i] = edge.types[C[new_edges[1,i]],C[new_edges[2,i]]]
    }
    new_type = as.data.frame(table(new.edge.labels, dnn = list("type")), responseName = "type_n")
    #new_type = sort(unique(new.edge.labels))
    #new_type_n = as.numeric(table(new.edge.labels))
    
    ap = as_ids(E(G_pc)[new_edges[1,] %--% new_edges[2,]])
    
    if(length(ap) == 0){ 
      to_del_type = new_type
      to_del_type$to_del = to_del_type$type_n
    }
    else{
      ap_type = as.data.frame(table(E(G_pc)[ap]$label, dnn = list("type")), responseName = "type_n")
      #ap_type = sort(unique(E(G_pc)[ap]$label))
      #ap_type_n = as.numeric(table(E(G_pc)[ap]$label))
      to_del_type = merge(new_type, ap_type,by=c("type"), all=TRUE)
      to_del_type[is.na(to_del_type)] <- 0
      to_del_type$to_del = to_del_type$type_n.x - to_del_type$type_n.y
    }
    
    to_del = c(1:gsize(G_pc))
    to_del = to_del[!to_del %in% ap]
    
    ###### Checking if we have enough edges to delete ##########
    cond = c()
    for (j in 1:length(to_del_type$type)) {
      if(length(to_del[E(G_pc)[to_del]$label == to_del_type$type[j]]) < to_del_type$to_del[j]) cond[j] = 1
      else cond[j] = 0
    }
    if(sum(cond) == 0) break
    if(counter > max.rep){ #print("Error: Network too sparse. Not enought edges to delete")
      break}
  }
  
  if (counter > max.rep) {
    print("Error: Network too sparse. Not enought edges to delete.")
    list("network" = make_empty_graph(n = 0, directed = FALSE),"rep" = counter)
  }
  else{
    ####### Sampling the edges to delete from the reduced edge list #######
    
    del_eid = c()
    for (j in 1:length(to_del_type$type)) {
      set = to_del[E(G_pc)[to_del]$label == to_del_type$type[j]]
      if(to_del_type$to_del[j] == 0) del_eid = del_eid
      else{
        if(length(set) == 1 && to_del_type$to_del[j] == 1){del_eid = c(del_eid, set)}
        else del_eid = c(del_eid,sample(set,to_del_type$to_del[j],,replace = FALSE))
      }
    }
    #del_eid 
    G_pc = delete_edges(G_pc, c(ap,del_eid))
    G_pc <- G_pc + edges(as.vector(new_edges))
    
    list("network" = G_pc, "rep" = counter)
  }
}

#simulate network with a planted cliques of size K 
sample_IRG_planted_clique = function(p,K){
  if(K==0 || K==1){G = sample_IRG(p)}
  if(K>1){
    G = sample_IRG(p)
    smpl = sample(V(G)$name,K)
    new_edges = as.vector(combinations(smpl, k = 2, layout = "column"))
    G <- G + edges(new_edges)
    G = simplify(G, remove.multiple = TRUE, remove.loops = TRUE)
  }
  return(G)
}

##################### Simulating ERMMs with planted hubs #######################

planted_hub <- function(G, d_m, k, C) {
  d <- degree(G)
  s_d <- sd(d)
  
  u_max <- which(d == d_m)[1]
  if (is.na(u_max)) return(G)
  
  #range_max <- max(1, ceiling(k * s_d))
  #d_star <- sample(seq(d_m + 1, d_m + range_max), 1)
  d_star <- max(1, d_m + ceiling(k * s_d))
  n_new_nei <- d_star - degree(G, u_max)
  
  target_v <- setdiff(V(G), neighbors(G, u_max))
  target_v <- setdiff(target_v, u_max)
  
  if (length(target_v) == 0 || n_new_nei == 0) return(G)
  
  # Filter target_v by type
  filtered_targets <- c()
  for (v in target_v) {
    nei_v <- neighbors(G, v)
    if (any(C[nei_v] == C[u_max])) {
      filtered_targets <- c(filtered_targets, v)
    }
  }
  
  if (length(filtered_targets) == 0) return(G)
  
  if (length(filtered_targets) <= n_new_nei) {
    pot_nei <- filtered_targets
  } else {
    pot_nei <- sample(filtered_targets, n_new_nei)
  }
  
  #new_edges <- as.vector(rbind(rep(u_max, length(pot_nei)), new_nei))
  
  # Loop to build to_del while avoiding duplicates
  
  to_del <- c()
  new_nei <- c()
  for (v in pot_nei) {
    to_del_temp = to_del
    nei_v <- neighbors(G, v)
    same_type <- nei_v[C[nei_v] == C[u_max]]
    
    if (length(same_type) == 0) next
    
    # Try deleting a unique edge without creating duplicates
    found <- FALSE
    to_sample <- same_type
    while (length(to_sample) > 0 && !found) {
      del_node <- if (length(to_sample) == 1) to_sample else sample(to_sample, 1)
      del_edge <- sort(c(v, del_node))
      to_del_temp <- c(to_del, del_edge)
      
      # Split into edge list and check for duplicates
      edge_list <- split(to_del_temp, ceiling(seq_along(to_del_temp) / 2))
      if (!any(duplicated(edge_list))) {
        to_del <- to_del_temp
        new_nei <- c(new_nei, v)
        found <- TRUE
      } else {
        to_sample <- setdiff(to_sample, del_node)
      }
    }
  }
  
  if (length(to_del) %% 2 != 0) to_del <- to_del[-length(to_del)]
  
  if (length(to_del) == 0) return(G)
  if (length(to_del) != length(new_nei) * 2) {
    message("Mismatch: to_del has ", length(to_del), 
            ", expected ", length(new_nei) * 2)
    return(G)
  }
  
  if (length(to_del) > 0) {
    eid <- get.edge.ids(G, to_del)
    G_hub <- delete_edges(G, eid[eid > 0])
  }
  
  new_edges <- as.vector(rbind(rep(u_max, length(new_nei)), new_nei))
  G_hub <- G_hub + edges(new_edges)
  
  return(G_hub)
}

planted_hubs <- function(G, Rep = 5, k = 3, C) {
  G_h <- G
  d_m = max(degree(G))
  i=1
  while (!is.na(d_m) && i <= Rep) {
    G_h <- planted_hub(G_h, d_m, k, C)
    i = i+1
    d_unique <- sort(unique(degree(G_h)), decreasing = TRUE)
    d_m <- d_unique[i]
  } 
  
  list("Graph" = G_h, "Iter" = i-1)
}

################################################################################
#################### Configuration model with two groups #######################
################################################################################

config.m_net = function(deg1,deg2){
  degs = c(deg1,deg2)
  if (sum(degs) %% 2 != 0){print("Sum of degrees is not even.")}
  else{
    G <- sample_degseq(degs, method = "vl")
    ml <- MOAR_LETTERS(3)
    V(G)$name = ml[1:gorder(G)]
    V(G)$group = c(rep(1,length(deg1)),rep(2,length(deg2)))
    
    G
  }
}

################################################################################
######################### Random Geometric Graphs ##############################
################################################################################


sample.RGG = function(n,r){
  ##Euclidean distance  
  get.dist = function(x)
  {
    sqrt(x[1]^2 + x[2]^2)
  }
  ## Generate matrix with total number of possible edges N(N-1)/2 ##
  get.pairs = combinations(n,2) 
  ##random points in 2D with uniform distribution
  rnd.points = matrix(runif(2 * n), ncol = 2)
  perms = get.pairs
  ###Caculate the difference between the points to calculate euclidien distance between each pair
  Edges = apply(perms, 1, FUN = function(x){
    vec.diff = rnd.points[x[1], ] - rnd.points[x[2], ]
    get.dist(vec.diff)
  })
  
  res = cbind(Edges, perms)
  colnames(res) = c('E', 'V1', 'V2')
  
  matnew<-res[res[, 'E'] <= r,] 
  if(is.matrix(matnew) ==  "FALSE") {
    matnew = matrix(matnew, nrow = 1)
    colnames(matnew) = c('E', 'V1', 'V2')
  }
  rgedge<-cbind(matnew[,'V1'], matnew[,'V2'])
  
  if(length(matnew[,'E']) == 0) G = make_empty_graph(n, directed = FALSE)
  else G = graph_from_edgelist(rgedge, directed=FALSE)
  G
}  

sample.RGG.with.two.groups = function(bs,r){
  n = bs[1]+bs[2]
  ### Euclidian distance ###
  get.dist = function(x)
  {
    sqrt(x[1]^2 + x[2]^2)
  }
  ## Generate matrix with total number of possible edges N(N-1)/2 ##
  v = seq(1,n,1)
  v1 = v[1:bs[1]]
  v2 = v[(bs[1] + 1):(bs[1] + bs[2])]
  get.pairs = combinations(v,2) 
  ##random points in 2D with uniform distribution
  rnd.points = matrix(runif(2 * n), ncol = 2)
  perms = get.pairs
  ###Caculate the difference between the points to calculate euclidien distance between each pair
  Edges = apply(perms, 1, FUN = function(x){
    vec.diff = rnd.points[x[1], ] - rnd.points[x[2], ]
    get.dist(vec.diff)
  })
  
  res = cbind(Edges, perms)
  colnames(res) = c('E', 'V1', 'V2')
  
  edgees = c()
  for (i in 1:nrow(res)) {
    if(res[i,2] %in% v1 && res[i,3] %in% v1 && res[i,1] < r[1]) edgees[i] = TRUE
    else
      if(res[i,2] %in% v1 && res[i,3] %in% v2 && res[i,1] < r[3]) edgees[i] = TRUE
      else
        if(res[i,2] %in% v2 && res[i,3] %in% v2 && res[i,1] < r[2]) edgees[i] = TRUE
        else edgees[i] = FALSE
  }
  
  matnew<-res[edgees,] 
  
  if(is.matrix(matnew) ==  "FALSE") {
    matnew = matrix(matnew, nrow = 1)
    colnames(matnew) = c('E', 'V1', 'V2')
  }
  rgedge<-cbind(matnew[,'V1'], matnew[,'V2'])
  
  if(length(matnew[,'E']) == 0) G = make_empty_graph(n, directed = FALSE)
  else G = graph_from_edgelist(rgedge, directed=FALSE)
  V(G)$name = ml[1:gorder(G)]
  V(G)$group = c(rep(1,bs[1]), rep(2,bs[2]))
  G
}

################################################################################
################# Test statistic and related GOF tests ########################
################################################################################
# Using vvRKHS
# Different kernel choices
# CalculateWLKernel(P,level = 3)
# CalculateGeometricRandomWalkKernel(P,level = 3)
# compute.sp.kernel(P,level = 1)
# CalculateVertexEdgeHistGaussKernel(P,level = 0.1)
# CalculateGraphletKernel(P,level = 3)
# CalculateConnectedGraphletKernel(P,level = 3)

# test_G: An igraph object containing the observed network
# C:  A vector of group labels for all the vertices
# p: The edge probability matrix for the null model
# g.kernel: Kernel to use to compute the test statistic
# alpha: level of significance for the test 
# M: Number of iterations in Monte Carlo setting to compute the null set for test statistic

generate.one.gKSS=function(G, C, p, g.kernel=CalculateWLKernel,level=3){
  V(G)$name = C
  E(G)$label = edge.labels(G, C)
  # as_adjacency_matrix(...) function uses $name attribute to assign row and column names to the matrix 
  # graph_from_adjacency_matrix(...) function use row and column names to assign vertex attributes
  # GetGraphInfo(...) uses first vertex attribute in the list to create vlabels amd first edge attribute to create edge labels
  # graphkernel package use vlabels as node labels and edge labels from 'GetGraphInfo(...)' to compute kernel values
  
  X = as_adjacency_matrix(G, type = "both")
  S = X[upper.tri(X)]
  n = length(S)
  S.t =  p[upper.tri(p)]
  S.t.vec = abs(S - S.t)
  S.mat = S.t.vec %*% t(S.t.vec)
  
  P = compute.transition.list(X)
  K = compute.kernel(P, g.kernel, level)
  
  J.kernel = S.mat * K
  
  W=rep(1/n,n)
  J.kernel.out=J.kernel
  stats.value = as.numeric(t(W)%*%J.kernel%*%W) 

  #Return:
  #stats.value: gKSS^2
  #J.kernel: J kernel matrix for wild bootstrap 
  list(stats.value=stats.value,J.kernel=J.kernel.out, K=K, S=S.mat)
}

GOF_IRG = function(test_G, C , p_0, M = 200, g.kernel= CalculateWLKernel, level = 3, alpha = 0.05){
  
  ######## Test Statistics from observed network ###########
  
  test_stat = generate.one.gKSS(test_G, C, p_0, g.kernel,level)$stats.value
  print(test_stat)
  
  ############## Simulations from null model ERMM-gKSS ##############
  
  cores = detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  statistic = c()
  for (i in 1:M){
    g = sample_IRG(p_0)
    statistic[i] = generate.one.gKSS(g, C, p_0, g.kernel,level)$stats.value
    print(i)
  }
  stopCluster(cl)
  
  critical_value_lower = quantile(statistic,alpha/2)
  critical_value_upper = quantile(statistic,(1-alpha/2))
  #median(statistic)
  
  ######## Decision and P-value #############
  
  pval <- two_tailed_pval(test_stat, statistic)
  
  ################### Output ERMM-gKSS#########################
  list("obs.test.stat" = test_stat, "LQ_e" = critical_value_lower, "UQ_e" = critical_value_upper, "P_value" = pval, "null_set" = statistic)
}

# Re-Sample algorithms 
generate.one.gKSS.sampled=function(G, C, p, sample.index, g.kernel=CalculateWLKernel,level=3){
  V(G)$name = C
  E(G)$label = edge.labels(G, C)
  # as_adjacency_matrix(...) function uses $name attribute to assign row and column names to the matrix 
  # graph_from_adjacency_matrix(...) function use row and column names to assign vertex attributes
  # GetGraphInfo(...) uses first vertex attribute in the list to create vlabels 
  # graphkernel package use vlabels as node labels to compute kernel values
  
  X = as_adjacency_matrix(G, type = "both")
  S = matrix(X,byrow=TRUE)[sample.index]
  S.t =  p[sample.index]
  S.t.vec = abs(S - S.t)
  S.mat = S.t.vec %*% t(S.t.vec)
  
  n = length(sample.index)
  P = compute.sampled.list(X, sample.index)
  K = compute.kernel(P, g.kernel, level)
  
  J.kernel = S.mat * K
  
  W=rep(1/n,n)
  J.kernel.out=J.kernel
  stats.value = as.numeric(t(W)%*%J.kernel%*%W) 
  
  #Return:
  #stats.value: gKSS^2
  #J.kernel: J kernel matrix for wild bootstrap 
  list(stats.value=stats.value,J.kernel=J.kernel.out, K=K, S=S.mat)
}

GOF_IRG_resamp = function(test_G, C , p_0, s_size, M = 200, g.kernel= CalculateWLKernel, level = 3, alpha = 0.05){
  
  ######## Test Statistics from observed network ###########
  
  sample.index = sample_index(n, s.size = s_size, replace = TRUE)
  test_stat = generate.one.gKSS.sampled(test_G, C, p_0, sample.index, g.kernel,level)$stats.value
  
  print(test_stat)
  
  ############## Simulations from null model ERMM-gKSS ##############
  
  statistic = c()
  for (i in 1:M){
    g = sample_IRG(p_0)
    sample.index = sample_index(n, s.size = s_size, replace = TRUE)
    statistic[i] = generate.one.gKSS.sampled(g, C, p_0, sample.index, g.kernel,level)$stats.value
    print(i)
  }
  
  critical_value_lower = quantile(statistic,alpha/2)
  critical_value_upper = quantile(statistic,(1-alpha/2))
  #median(statistic)
  
  ######## Decision and P-value #############
  
  #decision = isTRUE(critical_value_lower >= test_stat ||test_stat >= critical_value_upper)
  #False means do not reject null
  
  pval <- two_tailed_pval(test_stat, statistic)

  ################### Output ERMM-gKSS#########################
  list("obs.test.stat" = test_stat, "LQ_e" = critical_value_lower, "UQ_e" = critical_value_upper, "P_value" = pval, "null_set" = statistic)
}

generate.one.gKSS.cons.Ker=function(X, p, diagonal=1, v.scale=TRUE){ 
  S = X[upper.tri(X)]
  n.sim=length(S)
  S.t =  p[upper.tri(p)]
  S.t.vec = abs(S - S.t)
  S.mat = S.t.vec %*% t(S.t.vec)
  #Ky.mat = (S*2-1)%*%t(S*2-1) 
  
  J.kernel = S.mat
  
  n = length(X[upper.tri(X)])
  
  W=rep(1/n,n)
  J.kernel.out=J.kernel
  v.kernel = var(S)
  stats.value = t(W)%*%J.kernel%*%W 
  if(v.scale)KSD = stats.value*sqrt(v.kernel)
  #Return:
  #stats.value: gKSS^2
  #J.kernel: J kernel matrix for wild bootstrap 
  list(stats.value=stats.value,J.kernel=J.kernel.out, S=S.mat, vs.KSD = KSD)
}

############################ Conditional gKSS #####################################  

ERMM_gKSS = function(G, vgroup, g.kernel = CalculateWLKernel, level = 3){
  C= vgroup
  Y = G[,]
  X = G[,]
  X[lower.tri(X)] <- NA
  diag(X) <- NA
  el = edge.labels(G, C, L)
  edges = which(X[,]==1)
  non_edges = which(X[,]==0)
  
  P=list()
  const = c()
  counter = 1
  for (i in 1:length(edges)) {
    type = el[edges[i]]
    to_add = non_edges[which(el[non_edges] == type)]
    if(length(to_add) == 0) print("Warning: One type in the mixtiure graph is complete")
    for (j in 1:length(to_add)) {
      Y = G[,]
      Y[edges[i]] = abs(1 - Y[edges[i]])
      Y[to_add[j]] = abs(1 - Y[to_add[j]])
      Y = symmetricize(Y ,method = c("ud"), adjacencyList = FALSE)
      G_new = graph_from_adjacency_matrix(Y, mode = "undirected")
      V(G_new)$names = V(G)$name
      V(G_new)$group = V(G)$group
      P[[counter]] = G_new
      const[counter] = 1/(length(to_add))
      counter = counter + 1
    }
    Y = G[,]
    G_last = graph_from_adjacency_matrix(Y, mode = "undirected")
    V(G_last)$names = V(G)$name
    V(G_last)$group = V(G)$group
    P[[counter]] = G_last
  }
  #P
  #const
  #outer(const, const, function(x,y)x*y)
  
  
  kernel.matrix = g.kernel(P, level)
  n = ncol(kernel.matrix)-1
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K - outer(K.vec, K.vec, function(x,y)x+y)
  sum(K*outer(const, const, function(x,y)x*y))/(gsize(G)*gsize(G))
}

GOF_ERMM = function(G, vmem, M, g.kernel = CalculateWLKernel, level = 3){
  
  n_cores = detectCores()
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  stat.list = c()
  V(G)$group = vmem
  stat.list[1] = ERMM_gKSS(G, V(G)$group)
  for (i in 2:M) {
    Y = G[,]                        # get the adjacency matrix from the graph
    X = G[,]                        # get the adjacency matrix from the graph
    X[lower.tri(X)] <- NA           # make the lower triangular of the adjacency matrix empty
    diag(X) <- NA                   # make the diagonal of the adjacency matrix empty
    el = edge.labels(G, V(G)$group, L)       # label the edges according to their type
    edges.ind = which(X[,]==1)      # get the indices in X that corresponds to an edge
    non_edges.ind = which(X[,]==0)  # get the indices in X that corresponds to a non-edge
    
    s = sample(edges.ind,1)
    g_s = el[s]
    to_add = non_edges.ind[which(el[non_edges.ind] == g_s)]
    r = sample(to_add,1)
    g_r = el[r]
    if(g_s != g_r)print("g_s != g_r")
    
    Y[s] = abs(1 - Y[s])
    Y[r] = abs(1 - Y[r])
    Y = symmetricize(Y ,method = c("ud"), adjacencyList = FALSE)
    G_new = graph_from_adjacency_matrix(Y, mode = "undirected")
    V(G_new)$group = V(G)$group             # Re-assign node types
    stat.list[i] = ERMM_gKSS(G_new, V(G_new)$group, g.kernel, level)
    #plot(G_new, vertex.color = V(G_new)$group)
    G <- G_new
    #print(i)
  }
  stopCluster(cl)
  
  r = rank(stat.list, ties.method = "random")
  pvalue <- min(r[1]/(length(r)), (length(r) + 1 - r[1])/(length(r)))
  
  list("statistic" = stat.list[1], "null_set" = stat.list, "p.value" = pvalue)
}

################################################################################
########################### Kernel Computation #################################
################################################################################

##################### WL-IMQ kernel ############################################

WL_RBF <- function(graphs, h_max = 2, sigma = 0.1,  label_list) {
  n <- length(graphs)  # Number of graphs
  K_mat <- Matrix(0, n, n, sparse = TRUE)  # Initialize kernel matrix
  
  # Extract number of nodes and edges per graph
  num_v <- sapply(graphs, vcount)
  num_e <- sapply(graphs, ecount)
  
  # Compute the maximum degree across all graphs
  degree_max <- sapply(graphs, function(g) max(degree(g)))
  
  # Initialize label storage
  #label_list <- vector("list", n)
  
  #for (i in 1:n) {
  #  label_list[[i]] <- rep(1, num_v[i])  # Initial labels are all "1"
  #}
  
  # Compute kernel values based on initial vertex labels
  for (i in 1:n) {
    for (j in i:n) {
      combined_labels <- unique(c(label_list[[i]],label_list[[j]]))
      tab1 = sapply(combined_labels, function(x) sum(label_list[[i]] == x))
      tab2 = sapply(combined_labels, function(x) sum(label_list[[j]] == x))
      K_mat[i, j] <- exp(-(norm(tab1-tab2, type="2")^2)/(2*sigma*sigma))
      K_mat[j, i] <- K_mat[i, j]
    }
  }
  
  # Weisfeiler-Lehman Iterations
  for (h in 1:h_max) {
    new_labels <- vector("list", n)
    
    for (i in 1:n) {
      graph <- graphs[[i]]
      old_labels <- label_list[[i]]
      new_labels[[i]] <- rep("", num_v[i])
      
      for (v in V(graph)) {
        neighbor_labels <- sort(old_labels[neighbors(graph, v)])
        combined_label <- paste0(old_labels[v], ":", paste(neighbor_labels, collapse = ","))
        new_labels[[i]][v] <- digest(combined_label, algo = "md5", serialize = FALSE)
      }
    }
    
    label_list <- new_labels  # Update labels
    
    # Update kernel values
    for (i in 1:n) {
      for (j in i:n) {
        combined_labels <- unique(c(label_list[[i]],label_list[[j]]))
        tab1 = sapply(combined_labels, function(x) sum(label_list[[i]] == x))
        tab2 = sapply(combined_labels, function(x) sum(label_list[[j]] == x))
        K_mat[i, j] <- K_mat[i, j] + exp(-(norm(tab1-tab2, type="2")^2)/(2*sigma*sigma))
        K_mat[j, i] <- K_mat[i, j]
      }
    }
  }
  
  return(K_mat)  # Return the computed WL kernel matrix
}

##################### WL-IMQ kernel ############################################

WL_IMQ_kernel <- function(graphs, h_max = 2, label_list, mag = 1000) {
  n <- length(graphs)  # Number of graphs
  K_mat <- Matrix(0, n, n, sparse = TRUE)  # Initialize kernel matrix
  
  # Extract number of nodes and edges per graph
  num_v <- sapply(graphs, vcount)
  num_e <- sapply(graphs, ecount)
  
  # Compute the maximum degree across all graphs
  degree_max <- sapply(graphs, function(g) max(degree(g)))
  
  # Compute kernel values based on initial vertex labels
  for (i in 1:n) {
    for (j in i:n) {
      combined_labels <- unique(c(label_list[[i]], label_list[[j]]))
      tab1 = sapply(combined_labels, function(x) sum(label_list[[i]] == x))
      tab2 = sapply(combined_labels, function(x) sum(label_list[[j]] == x))
      med = median(((tab1/sum(tab1)) - (tab2/sum(tab2)))^2)
      offset = mean(c(tab1, tab2)) / mag
      if (med == 0) {
        K_mat[i, j] <- 1 / sqrt((norm((tab1 / sum(tab1)) - (tab2 / sum(tab2)), type = "2")^2) + offset)
      } else {
        K_mat[i, j] <- 1 / sqrt((norm((tab1 / sum(tab1)) - (tab2 / sum(tab2)), type = "2")^2) / med + offset)
      }
      K_mat[j, i] <- K_mat[i, j]
    }
  }
  
  # Weisfeiler-Lehman Iterations
  for (h in 1:h_max) {
    new_labels <- vector("list", n)
    
    # For each graph, generate new labels by combining node labels with neighbors
    for (i in 1:n) {
      graph <- graphs[[i]]
      old_labels <- label_list[[i]]
      new_labels[[i]] <- rep(0, num_v[i])  # Initialize new labels as zeros (will be replaced by integers)
      
      # Create a map to store the new labels
      label_map <- list()
      new_label_id <- 1  # Start with ID 1 for new labels
      
      # Update the node labels based on the neighbors
      for (v in V(graph)) {
        # Get the neighbors' labels
        neighbor_labels <- sort(old_labels[neighbors(graph, v)])
        
        # Combine the current node's label with its neighbors' labels
        combined_label <- paste0(old_labels[v], ":", paste(neighbor_labels, collapse = ","))
        
        # Check if this combined label already has a unique ID in the map
        if (!combined_label %in% names(label_map)) {
          label_map[[combined_label]] <- new_label_id
          new_label_id <- new_label_id + 1  # Increment ID for next label
        }
        
        # Assign the new label ID to the node
        new_labels[[i]][v] <- label_map[[combined_label]]
      }
    }
    
    label_list <- new_labels  # Update labels
    
    # Update kernel values based on new labels
    for (i in 1:n) {
      for (j in i:n) {
        combined_labels <- unique(c(label_list[[i]], label_list[[j]]))
        tab1 = sapply(combined_labels, function(x) sum(label_list[[i]] == x))
        tab2 = sapply(combined_labels, function(x) sum(label_list[[j]] == x))
        med = median(((tab1/sum(tab1)) - (tab2/sum(tab2)))^2)
        offset = mean(c(tab1, tab2)) / mag
        if (med == 0) {
          K_mat[i, j] <- K_mat[i, j] + 1 / sqrt((norm((tab1 / sum(tab1)) - (tab2 / sum(tab2)), type = "2")^2) + offset)
        } else {
          K_mat[i, j] <- K_mat[i, j] + 1 / sqrt((norm((tab1 / sum(tab1)) - (tab2 / sum(tab2)), type = "2")^2) / med + offset)
        }
        K_mat[j, i] <- K_mat[i, j]
      }
    }
  }
  
  return(K_mat)  # Return the computed WL kernel matrix
}

####### Computes Kernel inner product #######
compute.kernel=function(P, g.kernel = CalculateWLKernel, level=3){
  n = length(P) - 1 
#  if (identical(g.kernel, CalculateShortestPathKernel)) {
#    kernel.matrix = g.kernel(P)
#  } else {kernel.matrix = g.kernel(P, level)}
  kernel.matrix = g.kernel(P, level)
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K - outer(K.vec, K.vec, '+')
  return(K)
}

#removed from above function
#else{
#  if(identical(g.kernel,WL_IMQ_kernel)) {
#    label_list = vector(mode = "list", length = length(P))
#    label_list = lapply(label_list, function(x) C)
#    kernel.matrix = g.kernel(P, level, label_list)
#  }
#  else{
#    if(identical(g.kernel,WL_RBF)) {
#      label_list = vector(mode = "list", length = length(P))
#      label_list = lapply(label_list, function(x) C)
#      kernel.matrix = g.kernel(P, level, sigma = 0.1, label_list)
#    }

compute.transition.list = function(X){
  n <- nrow(X)
  P <- vector("list", n * (n - 1) / 2 + 1)  # Preallocate list
  
  counter <- 1
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      Y <- X
      Y[i, j] <- abs(1 - X[i, j])
      Y[j, i] <- Y[i, j]
      P[[counter]] <- graph_from_adjacency_matrix(Y, mode = "undirected")
      counter <- counter + 1
    }
  }
  
  # Original graph as the last element
  P[[counter]] <- graph_from_adjacency_matrix(X, mode = "undirected")
  
  return(P)
}

compute.sampled.list = function(X, sample.index){
  P=list()
  l = length(sample.index)
  for (w in 1:l){
    x = X
    x[sample.index[w]] = abs(1 - X[sample.index[w]])
    x = symmetricize(x ,method = c("ud"), adjacencyList = FALSE)
    G = graph_from_adjacency_matrix(x, mode = "undirected")
    P[[w]] = G
  }
  P[[l+1]] = graph_from_adjacency_matrix(X, mode = "undirected")
  P
}

compute.normalized = function(K){
  V = diag(K)
  D.mat = diag(1.0/sqrt(V))
  #  D.mat[sapply(D.mat, is.infinite)] <- 0 ### Need to confirm
  D.mat%*%K%*%D.mat
}

################################################################################
########################## Other helping codes #################################
################################################################################

two_tailed_pval <- function(test_stat, statistic) {
  tstat_set = c(test_stat, statistic)
  p_lower <- mean(tstat_set <= test_stat)
  p_upper <- mean(tstat_set >= test_stat)
  pval <- 2 * min(p_lower, p_upper)
  return(min(pval, 1))  # Cap at 1
}

#stat_set = c(test_stat, statistic)
#r = rank(stat_set, ties.method = "random")
#pval = 2*min(r[1]/(length(r)), (length(r) + 1 - r[1])/(length(r)))

# The following are valid only if the null distribution is roughly symmetric and unimodal.
#pval <- mean(statistic >= as.numeric(test_stat)) + mean(statistic <= (2 * mean(statistic) - as.numeric(test_stat)))
#pval <- mean(abs(statistic - mean(statistic)) >=  abs(as.numeric(test_stat) - mean(statistic)))

perturb = function(X, a){
  ids = seq(1,length(X),1)
  S = sample(ids, length(X)/2, replace = FALSE)
  NS = setdiff(ids, S)
  
  X[S] <- sapply(X[S], function(x) max(x - a, 0))
  X[NS] <- sapply(X[NS], function(x) min(x + a, 1))
  
  return(X)
}

ran_perturb = function(X, q, sd = 0.5){
  ids = seq(1,length(X),1)
  S = sample(ids, length(X)/2, replace = FALSE)
  NS = setdiff(ids, S)
  
  X[S] <- sapply(X[S], function(x) min(max(x - rnorm(1, mean = 0, sd), 0),1))
  
  return(X)
}


## Code to find row id and column id of a specific index in a matrix

# ind: index for which we want to find row and column id
# nr: number of rows of the matrix we are searching
rc_id = function(ind, nr){
  if(ind%%nr == 0){col_id = ind%/%nr
  row_id = nr}
  else 
  {col_id = ind%/%nr+1
  row_id = ind - ind%/%nr*6}
  list("row_id" = row_id, "col_id" = col_id)
}

######## MLE for ERMM parameters with number of observed edges of each type ##############
# v_group is the vector containing vertex group information, vertex should be labeled as (1,2,...,L)

edge_counter_ermm = function(G, C)
{
  L = length(unique(C)) # no of block
  bs = as.numeric(table(C))
  X = as_adjacency_matrix(G, type = "both")
  n_e = matrix(0,nrow = L,ncol = L)
  for (i in 1:L) {
    for (j in i:L) {
      n_e[i,j] = sum(X[C == i, C == j])
    }
  }
  diag(n_e) = diag(n_e)/2
  n_e = symmetricize(n_e, method = "ud") 
  n_e
}

MLE.est = function(G, C)
{
  L = length(unique(C)) # no of block
  bs = as.numeric(table(C))
  X = as_adjacency_matrix(G, type = "both")
  N_ij = matrix(0,nrow = L,ncol = L)
  for (i in 1:L) {
    for (j in 1:L) {
      if(i==j) N_ij[i,i] = (bs[i]*(bs[i]-1)/2)
      else {N_ij[i,j] = bs[i]*bs[j]
      }
    }
  }
  N_ij
  
  
  n_e = matrix(0,nrow = L,ncol = L)
  for (i in 1:L) {
    for (j in i:L) {
      n_e[i,j] = sum(X[C == i, C == j])
    }
  }
  diag(n_e) = diag(n_e)/2
  n_e = symmetricize(n_e, method = "ud") 
  
  p_hat = matrix(0,nrow = L, ncol = L)
  for (i in 1:L) {
    for (j in i:L) {
      if(i==j) p_hat[i,j] = n_e[i,i]/(bs[i]*(bs[i] - 1)/2)
      else p_hat[i,j] = n_e[i,j]/(bs[i]*bs[j]) 
    }
  }
  p_hat = symmetricize(p_hat, method = "ud")

  list("estimates" = p_hat, "E" = n_e, "N_ij" = N_ij )
}

##### Sample upper triangular indices #####
sample_index = function(n, s.size, replace = TRUE){
  Non = matrix(seq(1,n*n,1), nrow = n, ncol = n, byrow = FALSE)
  N_set = Non[upper.tri(Non)]
  sample.index = sample(N_set, size = s.size, replace)
  sample.index
}

############# Potential Edge labels generator #########
# L is the number of groups

g_labels = function(L){
  g_labels = matrix(0,nrow = L,ncol = L)
  counter = L + 1
  for (i in 1:L) {
    for (j in i:L) {
      if(i==j) g_labels[i,i] = i
      else {g_labels[i,j] = counter
      g_labels[j,i] = g_labels[i,j]
      counter = counter + 1}
    }
  }
  g_labels
}

edge.labels = function(G, C){
  L= length(unique(C))
  V(G)$name = C
  e_mat = as_edgelist(G, names = TRUE)
  elabels = c()
  
  if(length(e_mat)== 0) return(elabels)
  
  edge.types = g_labels(L)
  for (i in 1:nrow(e_mat)) {
    elabels[i] = edge.types[e_mat[i,1],e_mat[i,2]]
  }
  return(elabels)
}

###################################################################

community_aware_WL_kernel <- function(graphs, h = 3, C) {
  graph_kernels <- list()
  
  for (g in graphs) {
    V(g)$label <- as.character(C)
    
    # Step 3: Perform Weisfeiler-Lehman label propagation (community-aware)
    for (i in 1:h) {
      new_labels <- sapply(V(g), function(v) {
        # Collect neighbor labels, restricting to the same community
        neighbors <- neighbors(g, v)
        same_community_neighbors <- neighbors[C[neighbors] == C[v]]
        sorted_labels <- sort(unique(c(V(g)$label[v], V(g)$label[same_community_neighbors])))
        return(paste(sorted_labels, collapse = "_"))
      })
      V(g)$label <- new_labels  # Update labels
    }
    
    # Step 4: Convert labels into a feature vector for kernel computation
    label_counts <- table(V(g)$label)
    graph_kernels[[length(graph_kernels) + 1]] <- as.vector(label_counts)
  }
  
  # Step 5: Compute RBF Kernel between graphs
  kernel_matrix <- matrix(0, length(graphs), length(graphs))
  for (i in 1:length(graphs)) {
    for (j in 1:length(graphs)) {
      kernel_matrix[i, j] <- exp(-sum((graph_kernels[[i]] - graph_kernels[[j]])^2))
    }
  }
  
  return(kernel_matrix)
}
