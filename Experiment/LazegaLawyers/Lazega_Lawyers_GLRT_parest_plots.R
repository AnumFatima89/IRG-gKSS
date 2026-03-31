
setwd("...") #for windows working directory

load("LazegaLawyers/Lazega_Lawyers.RData")
source("source-codes.R")

################ edge probability matrices #############

G = G_friend
Y = as_adjacency_matrix(G, type = "both")
n = gorder(G) 

C = rep(1 ,n)
Q_h = MLE.est(G, C)$estimates
p = matrix(Q_h[C,C], length(C), length(C))
diag(p) = 0

par = randnet::DCSBM.estimate(Y,C)$Phat
p_h = 1-exp(-par)
diag(p_h) = 0

#### Status ####

C = vertex_attr$V2
Q_h = MLE.est(G, C)$estimates
p2 = matrix(Q_h[C,C], length(C), length(C))
diag(p2) = 0

par = randnet::DCSBM.estimate(Y,C)$Phat
p_h2 = 1-exp(-par)
diag(p_h2) = 0

#### Gender ####

C = vertex_attr$V3
table(C)
Q_h = MLE.est(G, C)$estimates
p3 = matrix(Q_h[C,C], length(C), length(C))
diag(p3) = 0

par = randnet::DCSBM.estimate(Y,C)$Phat
p_h3 = 1-exp(-par)
diag(p_h3) = 0

#### Office ####

C = vertex_attr$V4
Q_h = MLE.est(G, C)$estimates
p4 = matrix(Q_h[C,C], length(C), length(C))
diag(p4) = 0

par = randnet::DCSBM.estimate(Y,C)$Phat
p_h4 = 1-exp(-par)
diag(p_h4) = 0

#### Practice ####

C = vertex_attr$V7
Q_h = MLE.est(G, C)$estimates
p7 = matrix(Q_h[C,C], length(C), length(C))
diag(p7) = 0

par = randnet::DCSBM.estimate(Y,C)$Phat
p_h7 = 1-exp(-par)
diag(p_h7) = 0

#### Law school ####

C = vertex_attr$V8
Q_h = MLE.est(G, C)$estimates
p8 = matrix(Q_h[C,C], length(C), length(C))
diag(p8) = 0

par = randnet::DCSBM.estimate(Y,C)$Phat
p_h8 = 1-exp(-par)
diag(p_h8) = 0

############################## GLRT #################################

edge_loglikelihood <- function(x, p) {
  # Element-wise calculation using the formula
  result <- log((p^x) * ((1 - p)^(1 - x)))
  return(result)
}

############### Actual test #############

## ER against DCSBM
df = (n-1) - 1
## ER against ERMM with two groups
df = 3 - 1
## ER against ERMM with three groups (p4 and p8)
df = 6 - 1
## ERMM against DCSBM with two groups
df = (n-2) + 3 - 3
## ERMM against DCSBM with three groups
df = (n-3) + 6 - 6

alpha = 0.05
test_stat_GLRT = -2*(sum(edge_loglikelihood(Y, p)) - sum(edge_loglikelihood(Y, p_h2)))
decision_GLRT = isTRUE(test_stat_GLRT > qchisq(1 - alpha, df) )

test_stat_GLRT
qchisq(1 - alpha, df)

1 - pchisq(test_stat_GLRT, df) 

####################### Network plots ##################

G1 = G_work
G2 = G_advice

line = -4
cex = 1.5
side = 3
adj=-0.05
par(mfrow=c(1,2), mai = c(0, 1, 0, 0))

plot(G1, vertex.size=10, vertex.label=NA)
mtext("A", side=side, line=line, cex=cex, adj=adj)

plot(G2, vertex.size=10, vertex.label=NA)
mtext("B", side=side, line=line, cex=cex, adj=adj)


G3 = G_friend
l <- layout_nicely(G3)
plot(G3, vertex.size=8, vertex.label=NA, layout = l, margin=0)

line = -4
cex = 1.5
side = 3
adj=-0.05
par(mfrow=c(3,2), mai = c(2, 2, 0, 2))

par(mar = c(2, 2, 0, 2))
plot(G3, vertex.size=8, vertex.label=NA, layout = l, margin=0)
mtext("A", side=side, line=line, cex=cex, adj=adj)
plot(G3, vertex.size= 8, vertex.label=NA, vertex.color = vertex_attr$V2, layout = l, margin=0)
mtext("B", side=side, line=line, cex=cex, adj=adj)

plot(G3, vertex.size= 8, vertex.label=NA, vertex.color = vertex_attr$V3, layout = l, margin=0)
mtext("C", side=side, line=line, cex=cex, adj=adj)
plot(G3, vertex.size= 8, vertex.label=NA, vertex.color = vertex_attr$V4, layout = l, margin=0)
mtext("D", side=side, line=line, cex=cex, adj=adj)

plot(G3, vertex.size= 8, vertex.label=NA, vertex.color = vertex_attr$V7, layout = l, margin=0)
mtext("E", side=side, line=line, cex=cex, adj=adj)
plot(G3, vertex.size= 8, vertex.label=NA, vertex.color = vertex_attr$V8, layout = l, margin=0)
mtext("F", side=side, line=line, cex=cex, adj=adj)

dev.off()

########## Bound on ditance ################

##### ER to ERMM with attributes #####

3*sum(abs(p-p2))/(choose(n,2))

##### ER to DCSBM #####

3*sum(abs(p-p_h))/(choose(n,2))

##### ER to DCSBM with attributes #####

3*sum(abs(p-p_h8))/(choose(n,2))

##### ERMM to DCSBM #####

3*sum(abs(p8-p_h8))/(choose(n,2))



