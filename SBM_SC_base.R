library(parallel)
library(Matrix)
library(dplyr)
library(permute)

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
  E <- eigen(L, symmetric = T)$vectors[,1:k]
  #E1 <- t(apply(E,1,function(x){x/sum(x^2)^0.5}))
  return(as.integer(kmeans(E, k, 1000000, 30)$cluster))
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


#############################################################

oldSONNET <- function(A, k, s, O, rep = 0, ncore = 1)
{
  on.exit(gc())
  N <- nrow(A)    # detecting number of nodes
  m <- (N-O)/s 		# Make sure s|N-O
  
  ## SUBGRAPHING STEP for both over and nonover with rep
  
  over <- sample(1:N, O, replace = F) #setting overlap part
  non <- sapply(1:(rep+1), function(x) sample((1:N)[-over], N-O, replace = F))
  # ^^ generating (rep+1) many partitions for nonoverlapping part
  
  output <- sapply(mclapply(1:(s*(rep+1)), function(i){
    r <- ceiling(i / s)   #detecting partition number
    pos <- i %% s + 1   #detecting position of the index
    index <- non[((pos-1)*m+1):(pos*m), r] #obtaining row number
    label <- rep(0, N) # initializing
    
    if(r != 1)  
    {
      label[index] <- old.spectral.cluster(A[index, index], k)
      return(label)
    }else{
      label[c(over,index)] <- old.spectral.cluster(A[c(over, index),
                                                     c(over, index)], k)
      return(label)
    }
  }, mc.cores = ncore), 'c')
  
  # working on first s subgraphs - forming standard labels
  tmp.std <- sapply(2:s, function(j){
    l <- rep(0,N)
    index <- (output[-over, j] != 0)
    t <- label.match(output[over, j], output[over, 1])
    l[over] <- t[[1]]
    l[-over][index] <- as.vector(tcrossprod(crossprod(sapply(output[-over, j][index],
                                                             make_basis, p = k),
                                                      t[[2]]), rbind(1:k)))
    l
  })
  std <- rep(0,N)
  std[over] <- apply(cbind(output[over,1], tmp.std[over,]), 1, getmode)
  std[-over] <- apply(cbind(output[-over,1], tmp.std[-over,]), 1, sum)
  
  if(rep == 0) return(std)
  
  # match the repetitions with standard
  
  ns <- sapply(mclapply((s+1):(s*(rep+1)), function(i){
    index <- (output[, i] != 0)
    lab <- rep(0, N)
    lab[index] <- label.match(output[index, i], std[index])[[1]]
    lab
  }, mc.cores = ncore), 'c')
  
  tmp.ns <- sapply(seq(from = 1, to = (s*rep), by = s), function(j){
    apply(ns[, j:(j+s-1)], 1, sum)
  })
  
  std[-over] <- apply(cbind(std[-over], tmp.ns[-over,]), 1, getmode)
  
  return(std)
}



###########################################################################
###########################################################################


