setwd("...") #for windows working directory
source("source-codes.R")

EL_work = read.table("LazegaLawyers/ELwork.dat")
EL_friend = read.table("LazegaLawyers/ELfriend.dat")
EL_advice = read.table("LazegaLawyers/ELadv.dat")
vertex_attr = read.table("LazegaLawyers/ELattr.dat")

G_work = graph_from_adjacency_matrix(as.matrix(read.table("LazegaLawyers/ELwork.dat")), mode = "max")
G_friend = graph_from_adjacency_matrix(as.matrix(read.table("LazegaLawyers/ELfriend.dat")), mode = "max")
G_advice = graph_from_adjacency_matrix(as.matrix(read.table("LazegaLawyers/ELadv.dat")), mode = "max")

plot(G_work, vertex.label = "", vertex.color = vertex_attr$V2)
plot(G_friend, vertex.label = "", vertex.color = vertex_attr$V2)
plot(G_advice, vertex.label = "", vertex.color = vertex_attr$V2)

