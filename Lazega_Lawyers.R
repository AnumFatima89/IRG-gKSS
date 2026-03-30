
##############################
########Methodology###########
##############################

# re-direct to the desired directory
setwd("...") #for windows working directory

load("LazegaLawyers/Lazega_Lawyers.RData")

source("source-codes.R")

#### For spectral test
library(mgcv)
library(RMTstat)
source('utils.R')

####################### Input Network ##################

#G = G_work
#G = G_advice
G = G_friend

# Status: V2
# Gender: V3
# Office: V4
# Practice: V7
# Law school V8

C = vertex_attr$V2
table(C)
n = gorder(G) 

plot(G, vertex.label = "", vertex.color = C)

######### Global setting ###########

M = 200
g.kernel= CalculateWLKernel
level = 3
alpha = 0.05

########## Tesiting the fit of an ERMM ###########

test_G = G
Q_h = MLE.est(test_G, C)$estimates

p = matrix(Q_h[C,C], length(C), length(C))
diag(p) = 0

results_ERMM = GOF_IRG(test_G, C , p, M, g.kernel, level, alpha)
isTRUE(as.numeric(results_ERMM$LQ_e) >= as.numeric(results_ERMM$obs.test.stat) || as.numeric(results_ERMM$obs.test.stat) >= as.numeric(results_ERMM$UQ_e))
#False means do not reject null

########## Tesiting the fit of an DCSBM ###########

test_G = G
Y = as_adjacency_matrix(test_G, type = "both")

par = randnet::DCSBM.estimate(Y,C)$Phat
p_h = 1-exp(-par)
diag(p_h) = 0

p = p_h

results_DCSBM = GOF_IRG(test_G, C , p, M, g.kernel, level, alpha)
isTRUE(as.numeric(results_DCSBM$LQ_e) >= as.numeric(results_DCSBM$obs.test.stat) || as.numeric(results_DCSBM$obs.test.stat) >= as.numeric(results_DCSBM$UQ_e))
#False means do not reject null

######Spectral test ##############

test_G = G
A = as_adjacency_matrix(test_G, type = "both")

## group membership provided

set.seed(6)
K0 = 3
test_stat_ST = GoFTestBoot_truecom(A, K0, n.dim = K0, clusters = C, n.b = 50, use.boot = T)$test.stat.truecom #GoFTestBoot(A, K0, n.dim = K0, n.b = 50, use.boot = T)$test.stat.boot[1]
decision_ST = isTRUE(test_stat_ST > qtw(alpha / 2, lower.tail = F)) #qtw(0.025, beta=1, lower.tail = FALSE, log.p = FALSE))
print(test_stat_ST)
print(decision_ST)
ptw(max(abs(test_stat_ST)), lower.tail=FALSE)

########## Testing the fit of an ER ###########

test_G = G
C = rep(1,n)
Q_h = MLE.est(test_G, C)$estimates


p = matrix(Q_h[C,C], length(C), length(C))
diag(p) = 0

results_ER = GOF_IRG(test_G, rep(1,n) , p, M , g.kernel, level , alpha)
isTRUE(as.numeric(results_ER$LQ_e) >= as.numeric(results_ER$obs.test.stat) || as.numeric(results_ER$obs.test.stat) >= as.numeric(results_ER$UQ_e))
#False means do not reject null

######Spectral test ##############

test_G = G
A = as_adjacency_matrix(test_G, type = "both")

## group membership provided

set.seed(6)
K0 = 3
test_stat_ST = GoFTestBoot_truecom(A, K0, n.dim = K0, clusters = C, n.b = 50, use.boot = T)$test.stat.truecom #GoFTestBoot(A, K0, n.dim = K0, n.b = 50, use.boot = T)$test.stat.boot[1]
decision_ST = isTRUE(test_stat_ST > qtw(alpha / 2, lower.tail = F)) #qtw(0.025, beta=1, lower.tail = FALSE, log.p = FALSE))
print(test_stat_ST)
print(decision_ST)
ptw(max(abs(test_stat_ST)), lower.tail=FALSE)
