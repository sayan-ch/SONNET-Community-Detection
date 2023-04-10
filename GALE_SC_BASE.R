print("Library read started")
library(parallel)
library(Matrix)
library(dplyr)
library(permute)
library(igraph)
print("library read worked")
## Generating the network

SBM.generate <- function(N, k, p.intra, p.inter, ncore = 1)
{
  on.exit(gc())                                                                 #garbage cleaning in the end       
  m <- sample(1:k, N, T)                                                        #generating community labels
  
  stor <- mclapply(1:(N-1), function(i) {
    tmp <- if_else(m[(i+1):N] == m[i],                                     #for every 
                   rbinom(N-i, 1, p.intra),
                   rbinom(N-i, 1, p.inter))
    
    list(su = sum(tmp), loc = which(tmp == 1))
  }, mc.cores = ncore)
  
  vec <- list('i' = vector(), 'j' = vector())
  count <- 1
  lapply(X = 1:(N-1), FUN = function(i){
    ifelse(stor[[i]]$su > 0,
           {
             count2 <- count + stor[[i]]$su
             vec$i[count:(count2-1)] <<- rep(i, stor[[i]]$su)
             vec$j[count:(count2-1)] <<- i + stor[[i]]$loc
             count <<- count2
           }, 0 )
  })
  
  A <- sparseMatrix(i = vec$i, j = vec$j, dims = c(N,N), symmetric = T)
  return(list('A' = A, 'member' = m))
}


## Spectral clustering

old.sp.cl <- function(A, k, degree)
{
  # computing the Laplacian
  N <- nrow(A)
  D <- sparseMatrix(i = 1:N, j = 1:N, x = 1/sqrt(degree))
  L <- tcrossprod(crossprod(D, A), D)
  
  # computing eigen decomposition
  E <- eigen(L, symmetric = T)$vectors[, 1:k]
  #E1 <- t(apply(E,1,function(x){x/sum(x^2)^0.5}))
  return(as.integer(kmeans(E, k, 100000,10)$cluster))
}


old.spectral.cluster <- function(A, k)
{
  # find the rows with all 0
  lab <- rep(0, times = nrow(A))
  degree <- colSums(A)
  index <- (degree == 0)
  lab[index] <- sample(1:k, sum(index), T) # assign the random labels
  l <- (lab==0)
  
  # perform spectral clustering on the remaining matrix
  lab[l] <- old.sp.cl(A[l, l], k, degree[l])
  return(lab)
}

## Label matching

make_basis <- function(l, p) replace(numeric(p), l, TRUE)

mismatch <- function(lab1,lab2) sum(lab1!=lab2)


getmode <- function(v) 
{
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

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

brute.match <- function(lab, fixed)   #align Z1 with Z2. both Nxk
{
  k <- max(lab, fixed)
  if(identical(lab, fixed)) return(list(lab,
                                        sparseMatrix(i = 1:k, j = 1:k, dims = c(k, k))))
  Z1 <- as(t(sapply(lab, make_basis, p = k)), "dgCMatrix")
  Z2 <- as(t(sapply(fixed, make_basis, p = k)), "dgCMatrix")
  sk <- shuffleSet(k, nset = factorial(k))
  P <- apply(sk, 1, function(x){return(as(x,"pMatrix"))})
  P[[factorial(k)]] <- as(diag(k), "pMatrix")
  l <- lapply(P, function(x){
    return(as.vector(tcrossprod(Z1, crossprod(cbind(1:k), t(x)))))
  })
  
  index <- which.min(sapply(l, function(j) {sum(j != fixed)}))[1]
  list(l[[index]], P[[index]])
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


#####################################################
#####################################################

GALE <- function(A, k, m, s, ncore) # m <- size of each subgraph
  # s <- number of subgraphs
{
  N <- dim(A)[1]
  
  ind <- lapply(1:s, function(i){sample(1:N, m, F)})
  
  ZZ <- mclapply(ind, function(q){
    lab <- old.spectral.cluster(A[q, q], k)
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





























