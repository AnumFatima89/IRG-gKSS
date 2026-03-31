
# re-direct to the desired directory
setwd("...") #for windows working directory
source("source-codes.R")

#### For spectral test
library(mgcv)
library(RMTstat)
source('utils.R')

####################### Input Network ##################

data("florentine_m")
G=florentine_m
vertex.attributes(G)
plot(G)

######### Global setting ###########

M = 200 
g.kernel= CalculateWLKernel 
level = 3
alpha = 0.05

########## Tesiting the fit of an ER ###########

test_G = G
n = gorder(test_G)
C = rep(1,n)
Q_h = MLE.est(test_G, C)$estimates

p = matrix(Q_h, length(C), length(C))
diag(p) =0

results_ER = GOF_IRG(test_G, C , p, M , g.kernel, level , alpha)
isTRUE(as.numeric(results_ER$LQ_e) >= as.numeric(results_ER$obs.test.stat) || as.numeric(results_ER$obs.test.stat) >= as.numeric(results_ER$UQ_e))
#False means do not reject null

results_ER$P_value

######Spectral test ##############

test_G = G
A = as_adjacency_matrix(test_G, type = "both")

## group membership provided

set.seed(6)
K0 = 1
test_stat_ST = GoFTestBoot_truecom(A, K0, n.dim = K0, clusters = C, n.b = 50, use.boot = T)$test.stat.truecom 
decision_ST = isTRUE(test_stat_ST > qtw(alpha / 2, lower.tail = F))
print(test_stat_ST)
print(decision_ST)
ptw(max(abs(test_stat_ST)), lower.tail=FALSE)
