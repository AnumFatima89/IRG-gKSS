
# re-direct to the desired directory
setwd("...") #for windows working directory
source("source-codes.R")

#### For spectral test
library(mgcv)
library(RMTstat)
source('utils.R')

####################### Input Network ##################

dolphin <- read.graph(file='https://lipn.univ-paris13.fr/~kanawati/datasets/dolphins.gml', format='gml')
plot(dolphin, vertex.color = V(dolphin)$value, vertex.label="", vertex.label.color = "black", vertex.label.dist=2, edge.width = 1.5, vertex.label.degree= pi/2, vertex.size=15)

G = dolphin
#vertex.attributes(G)
table(V(G)$value)
V(G)$value = replace(V(G)$value,V(G)$value==1,2)
V(G)$value = replace(V(G)$value,V(G)$value==0,1)
table(V(G)$value)

######### Global setting ###########

M = 200
g.kernel= CalculateWLKernel 
level = 3
alpha = 0.05

########## Tesiting the fit of an ERMM ###########

test_G = G
C = V(G)$value
Q_h = MLE.est(test_G, C)$estimates

p = matrix(Q_h[C,C], length(C), length(C))
diag(p) = 0

results_ERMM = GOF_IRG(test_G, C , p, M, g.kernel, level, alpha)
isTRUE(as.numeric(results_ERMM$LQ_e) >= as.numeric(results_ERMM$obs.test.stat) || as.numeric(results_ERMM$obs.test.stat) >= as.numeric(results_ERMM$UQ_e))
#False means do not reject null

results_ERMM$P_value

######Spectral test ##############

test_G = G
A = as_adjacency_matrix(test_G, type = "both")

set.seed(6)
K0 = 2
test_stat_ST = GoFTestBoot_truecom(A, K0, n.dim = K0, clusters = C, n.b = 50, use.boot = T)$test.stat.truecom 
decision_ST = isTRUE(test_stat_ST > qtw(alpha / 2, lower.tail = F))
print(test_stat_ST)
print(decision_ST)
ptw(max(abs(test_stat_ST)), lower.tail=FALSE)

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

##################################################
######### Testing with edge resampling ###########
##################################################

s_size = 30

########## Tesiting the fit of an ER ###########

test_G = G
n = gorder(test_G)
C = rep(1,n)
Q_h = MLE.est(test_G, C)$estimates

p = matrix(Q_h, length(C), length(C))
diag(p) =0

results_ER = GOF_IRG_resamp(test_G, C , p,s_size, M , g.kernel, level , alpha)
isTRUE(as.numeric(results_ER$LQ_e) >= as.numeric(results_ER$obs.test.stat) || as.numeric(results_ER$obs.test.stat) >= as.numeric(results_ER$UQ_e))
#False means do not reject null

results_ER$P_value

########## Tesiting the fit of an ERMM ###########

test_G = G
C = V(G)$value
Q_h = MLE.est(test_G, C)$estimates

p = matrix(Q_h[C,C], length(C), length(C))
diag(p) = 0

results_ERMM = GOF_IRG_resamp(test_G, C , p,s_size, M , g.kernel, level , alpha)
isTRUE(as.numeric(results_ERMM$LQ_e) >= as.numeric(results_ERMM$obs.test.stat) || as.numeric(results_ERMM$obs.test.stat) >= as.numeric(results_ERMM$UQ_e))
#False means do not reject null

results_ERMM$P_value

