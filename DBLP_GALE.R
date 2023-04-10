library(tidyverse)
library(Matrix)
library(parallel)
library(data.table)
library(permute)
library(dplyr)
library(igraph)

sp.cl <- function(A, k, degree)
{
  # computing the Laplacian
  N <- nrow(A)
  D <- sparseMatrix(i = 1:N, j = 1:N, x = 1/sqrt(degree))
  L <- tcrossprod(crossprod(D, A), D)
  print('Laplace done-eigen started')
  # computing eigen decomposition
  E <- eigen(L, T)$vectors[,1:k]
  #E1 <- t(apply(E,1,function(x){x/sum(x^2)^0.5}))
  #Ei <- eigen(L, T)
  #ord <- sort.list(abs(Ei$values), decreasing = T)[1:k]
  #ord <- ord[sort.list(Ei$values[ord], decreasing = T)]
  #E <- Ei$vectors[,ord]
  print('eigen done.')
  #rss <- rowSums(E^2)
  #print(paste0('rss = ', sum(rss==0)))
  #E1 <- E/sqrt(rss)
  #E1[rss==0,] <- 0
  print('kmeans started.')
  out <- kmeans(E,k,10000000,100)
  print(paste0('kmeans.iter = ', out$iter))
  return(as.integer(out$cluster))
}


spectral.cluster <- function(A, k)
{
  on.exit(gc())
  # find the rows with all 0
  lab <- rep(0, times = nrow(A))
  degree <- colSums(A)
  #tau <- mean(degree)
  index <- (degree == 0)
  lab[index] <- sample(1:k, sum(index), T) # assign the random labels
  l <- (lab==0)
  print(sum(index))
  # perform spectral clustering on the remaining matrix
  lab[l] <- sp.cl(A[l, l], k, degree[l])
  #nde <- degree + tau
  #print(paste0('nde = ', sum(nde==0)))
  #lab <- sp.cl(A,k,nde)
  return(lab)
}

getmode <- function(v) 
{
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

make_basis <- function(l, p) replace(numeric(p), l, TRUE)

mismatch <- function(lab1,lab2) sum(lab1!=lab2)

greedy.match <- function(Z1,Z2, K)  #align Z1 with Z2. Both NxK
{
  M <- as(crossprod(Z1, Z2), 'dMatrix')
  P <- sparseMatrix(i = NULL, j = NULL, dims = c(K,K), x = 1)
  while(max(M) != -1)
  {
    ind <- which(M == max(M), T)[1,]
    P[ind[1],ind[2]] <- 1
    M[ind[1],] <- rep(-1,K)
    M[,ind[2]]  <- rep(-1,K)
  }
  P
}

label.match <- function(lab, fixed)
{
  k <- max(lab, fixed)
  if(identical(lab, fixed)) return(list(lab, sparseMatrix(i = 1:k, j = 1:k, dims = c(k, k))))
  Z1 <- as(t(sapply(lab, make_basis, p = k)), "dgCMatrix")
  Z2 <- as(t(sapply(fixed, make_basis, p = k)), "dgCMatrix")
  P <- greedy.match(Z1, Z2, k)
  list(as.vector(tcrossprod(tcrossprod(Z1, t(P)), rbind(1:k))), P)
}

GALE <- function(A, k, m, s, ncore) # m <- size of each subgraph
  # s <- number of subgraphs
{
  N <- dim(A)[1]
  
  ind <- lapply(1:s, function(i){sample(1:N, m, F)})
  
  ZZ <- mclapply(ind, function(q){
    lab <- spectral.cluster(A[q, q], k)
    sparseMatrix(i = q, j = lab, dims = c(N, k))
  }, mc.cores = ncore)
  
  
  # Forming traversal
  thres <- ceiling(m^2/(2*N))
  
  super <- matrix(rep(0, s^2), s, s)
  super_graph <- sapply(1:s, function(i){
    sapply(1:s, function(j){
      super[i,j] <<- if_else(length(
        intersect(ind[[i]], ind[[j]])) < thres, 0, 1)})
  }) - diag(1, s) 
  
  
  
  x <- as.numeric(dfs(
    graph_from_adjacency_matrix(super_graph,
                                "undirected"),
    root = 1)$order)
  
  J <- length(x)
  unvisit <- rep(1, s)  # Should s be J? No as we're happy as long as all the subgraphs are visited once.
  tau <- rep(1, J)
  
  Zx <- list()
  Zx[[1]] <- ZZ[[x[1]]]
  Zhat <- ZZ[[x[1]]]
  temp <- Zx[[1]]
  temp2 <- sapply(1:N,function(k1){
    as.numeric(k1 %in% ind[[1]])
  })
  for(i in 2:J)
  {
    if(!unvisit[x[i]]) break
    
    unvisit[x[i]] <- 0
    S <- 	intersect(ind[[x[i]]], unlist(ind[x[1:(i-1)]]))
    P <- greedy.match(ZZ[[x[[i]]]][S, ], Zhat[S, ], k)
    Zx[[i]] <- ZZ[[x[[i]]]] %*% P
    
    temp <- temp + Zx[[i]]
    
    temp2 <- temp2 + sapply(1:N, function(k1){
      as.numeric(k1 %in% ind[[i]])
    })
    Zhat <- t(sapply(1:N, function(j){
      if(temp2[j] < tau[i]) 
        return(rep(0,k))
      return(temp[j,]/temp2[j])
    }))
  }
  
  
  final <- apply(Zhat, 1, function(e){
    which(e==max(e))[1]
  })
  
  return(final)
}

################################################################################
################################################################################
## DBLP GALE Code



edges <- read_csv('~/revision/data/dblp-ac-compact.csv')
comm <- read_csv('~/revision/data/author_label.csv')

truth <- comm$a.lab+1

A <- sparseMatrix(i = edges$V1, j = edges$a,
                  dims = c(4057, 4057),
                  symmetric = T)

sim <- 100

parlist <- list(
c(20,300),
c(20,1000)
)

for(i in 1:sim)
for(par in parlist){
tryCatch({
  s <- par[1]
  m <- par[2]
  
  time <- system.time(tmp <- GALE(A = A, k = 4,
            s = s, m = m, ncore = 20))[3]
            
  error <- mean(truth != label.match(tmp, truth)[[1]])
  
  write_csv(data.table(time, error),
            paste0('~/revision/output/data/DBLP/DBLP_GALE_p',
                    s, 'T', m,'.csv'),
            append = T
            )
  gc()},
  error = function(e){cat("ERROR :",conditionMessage(e),"\n")})
}






























