
# re-direct to the desired directory
setwd("...") #for windows working directory
source("source-codes.R")

#### For spectral test
library(mgcv)
library(RMTstat)
source('utils.R')


####################### Input Network ##################
library(igraphdata)
data("karate")
G = karate
G = upgrade_graph(G)

plot(G)
n = gorder(G)

######### Global setting ###########

M = 200
g.kernel= CalculateWLKernel
level = 3
alpha = 0.05

########## Tesiting the fit of an ERMM ###########

test_G = G
C = V(G)$Faction

Q_h = MLE.est(test_G, C)$estimates

p = matrix(Q_h[C,C], length(C), length(C))
diag(p) = 0

results_ERMM = GOF_IRG(test_G, C , p, M, g.kernel, level, alpha)
isTRUE(as.numeric(results_ERMM$LQ_e) >= as.numeric(results_ERMM$obs.test.stat) || as.numeric(results_ERMM$obs.test.stat) >= as.numeric(results_ERMM$UQ_e))
#False means do not reject null

results_ERMM$P_value

########## Tesiting the fit of an ERMM with four groups ###########

load("GOF test for an observed network/Karate_club/four_communities.RData")
source("source-codes.R")

memberships = matrix(0,nrow = n, ncol = 100)
modularities = c()
for (i in 1:100) {
  lc = cluster_louvain(G)
  mem = lc$membership
  memberships[,i] = mem
  modularities[i] = modularity(lc)
}

max(modularities)
id = which(modularities == max(modularities))[1]
mem = memberships[,id]

test_G = G
C = mem
Q_h = MLE.est(test_G, C)$estimates

p = matrix(Q_h[C,C], length(C), length(C))
diag(p) = 0

results_ERMM.L4 = GOF_IRG(test_G, C , p, M, g.kernel, level, alpha)
isTRUE(as.numeric(results_ERMM.L4$LQ_e) >= as.numeric(results_ERMM.L4$obs.test.stat) || as.numeric(results_ERMM.L4$obs.test.stat) >= as.numeric(results_ERMM.L4$UQ_e))
#False means do not reject null

results_ERMM.L4$P_value

########## Tesiting the fit of an DCSBM ###########

test_G = G
C = V(G)$Faction
Y = as_adjacency_matrix(test_G, type = "both")

par = randnet::DCSBM.estimate(Y,C)$Phat
p_h = 1-exp(-par)
diag(p_h) = 0

p = p_h

results_DCSBM = GOF_IRG(test_G, C , p, M, g.kernel, level, alpha)
isTRUE(as.numeric(results_DCSBM$LQ_e) >= as.numeric(results_DCSBM$obs.test.stat) || as.numeric(results_DCSBM$obs.test.stat) >= as.numeric(results_DCSBM$UQ_e))
#False means do not reject null

results_DCSBM$P_value

######Spectral test ##############

test_G = G
A = as_adjacency_matrix(test_G, type = "both")

## group membership provided

set.seed(50)
K0 = 2
test_stat_ST = GoFTestBoot_truecom(A, K0, n.dim = K0, clusters = C, n.b = 50, use.boot = T)$test.stat.truecom #GoFTestBoot(A, K0, n.dim = K0, n.b = 50, use.boot = T)$test.stat.boot[1]
decision_ST = isTRUE(test_stat_ST > qtw(alpha / 2, lower.tail = F))
print(test_stat_ST)
print(decision_ST)
ptw(max(abs(test_stat_ST)), lower.tail=FALSE)

########## Tesiting the fit of an ER ###########

test_G = G
C = rep(1,n)
Q_h = MLE.est(test_G, C)$estimates

p = matrix(Q_h[C,C], length(C), length(C))
diag(p) = 0

results_ER = GOF_IRG(test_G, C , p, M , g.kernel, level , alpha)
isTRUE(as.numeric(results_ER$LQ_e) >= as.numeric(results_ER$obs.test.stat) || as.numeric(results_ER$obs.test.stat) >= as.numeric(results_ER$UQ_e))
#False means do not reject null

results_ER$P_value


