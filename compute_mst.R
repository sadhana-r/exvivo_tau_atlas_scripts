#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(igraph)

adjlist=read.table(args[1]);
#adjlist_1=read.table(args[1]);
#adjlist_2=read.table(args[2]);
#adjlist_3=read.table(args[3]);

N=max(adjlist[,c(1,2)]);

adj=matrix(data=1,nrow=N,ncol=N); #Changed data from NA

#for (i in 1:dim(adjlist)[1]) {
#  adj[adjlist[i,1],adjlist[i,2]] = abs(adjlist[i,5]); #min(adjlist_1[i,3],adjlist_2[i,3],adjlist_3[i,3]);
	#used to be 1 - adjlilst[i,3] when using dsc
#}

for (i in 1:dim(adjlist)[1]) {
  adj[adjlist[i,1],adjlist[i,2]] = 1 - adjlist[i,3]; #min(adjlist_1[i,3],adjlist_2[i,3],adjlist_3[i,3]);
        #used to be 1 - adjlilst[i,3] when using dsc
}

g=graph_from_adjacency_matrix((adj + t(adj)),weighted='weight',diag=F,mode='undirected')

T=minimum.spanning.tree(g, algorithm='prim')

D=distances(g, algorithm='dijkstra')

# Compute the distances from D
#print("Distances from each subject")
#print(colSums(D))
root=which.min(colSums(D));

## If I were to take the second smallest distance as root
#root=order(colSums(D))[2]

# Print the shortest paths
path=shortest_paths(T,from=root)
for (i in 1:N) {
  if(i == root)
    cat(i,"\n")
  else
    cat(rev(path$vpath[[i]]),"\n")
}
