library(tidyverse)
library(Matrix)
library(parallel)
library(data.table)
library(permute)

sp.cl <- function(A, k, degree)
{
  # computing the Laplacian
  N <- nrow(A)
  D <- sparseMatrix(i = 1:N, j = 1:N, x = 1/sqrt(degree))
  L <- tcrossprod(crossprod(D, A), D)
  
  # computing eigen decomposition
  E <- eigen(L, T)$vectors[,1:k]
  #E1 <- t(apply(E,1,function(x){x/sum(x^2)^0.5}))
  E1 <- E/sqrt(rowSums(E^2))
  out <- kmeans(E1,k,10^8,1000)
  print(out$iter)
  return(as.integer(out$cluster))
}


spectral.cluster <- function(A, k)
{
  # find the rows with all 0
  lab <- rep(0, times = nrow(A))
  degree <- colSums(A)
  index <- (degree == 0)
  lab[index] <- sample(1:k, sum(index), T) # assign the random labels
  l <- (lab==0)
  print(sum(index))
  # perform spectral clustering on the remaining matrix
  lab[l] <- sp.cl(A[l, l], k, degree[l])
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

edges <- read_csv('~/revision/data/dblp-ac-compact.csv')
comm <- read_csv('~/revision/data/author_label.csv')

truth <- comm$a.lab+1

A <- sparseMatrix(i = edges$V1, j = edges$a,
                  dims = c(4057, 4057),
                  symmetric = T)

sim <- 100


for(i in 1:sim){
tryCatch({
  
  time <- system.time(tmp <- spectral.cluster(A, 4))[3]
            
  error <- mean(truth != label.match(tmp, truth)[[1]])
  
  write_csv(data.table(time, error),
            '~/revision/output/data/DBLP/DBLP_whole.csv',
            append = T
            )
  gc()},
  error = function(e){cat("ERROR :",conditionMessage(e),"\n")})
}









