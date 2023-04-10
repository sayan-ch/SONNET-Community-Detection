library(parallel)
library(Matrix)
library(dplyr)
library(rlist)
library(cluster)
library(RSpectra)


DCBM.fast <- function(N, k, P, Th, ncore = 1)
{
  on.exit(gc())                                                                 #garbage cleaning in the end       
  m <- sample(1:k, N, T)                                                        #generating community labels
  scale <- max(P)*max(Th)^2
  
  stor <- mclapply(1:(N-1), function(i) {
    #tmp <- rbinom(N-i, 1, pmin(P[m[i],m[(i+1):N]]*Th[i]*Th[(i+1):N], 1))
    tmp <- rbinom(N-i, 1, P[m[i],m[(i+1):N]]*Th[i]*Th[(i+1):N]/scale)
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


## Regularized spherical k-median spectral clustering

spher.sp.cl <- function(A, k, degree, tau)
{
  # computing the Laplacian
  N <- nrow(A)
  D <- sparseMatrix(i = 1:N, j = 1:N, x = 1/sqrt(degree + tau))
  L <- tcrossprod(crossprod(D, A), D)
  message('Eigen started')
  # computing eigen decomposition
  E <- eigen(L, symmetric = T)$vectors[,1:k]
  E1 <- t(apply(E, 1, function(x){x/sum(x^2)^0.5}))
  message('Eigen done-pam started')
  i <- which(is.na(E1), T)
  if(nrow(i) == 0) return(as.integer(
    pam(E1, k, diss = F,
        metric = 'manhattan',
        nstart = 1,
        cluster.only = T,
        pamonce = 6)))
  
  i <- unique(i[,1])
  
  out <- vector(length = N)
  out[i] <- sample(1:k, length(i), T)
  
  E2 <- E1[-i,]
  out[-i] <-  as.integer(
    pam(E2, k, diss = F,
        metric = 'manhattan',
        nstart = 1,
        cluster.only = T,
        pamonce = 6))
  message('pam done')
  return(out)
}


spher.spectral.cluster <- function(A, k, tau=-1)
{
  # find the rows with all 0
  lab <- rep(0, times = nrow(A))
  degree <- colSums(A)
  index <- (degree == 0)
  lab[index] <- sample(1:k, sum(index), T) # assign the random labels
  l <- (lab==0)
  
  if(tau < 0) tau <- mean(degree[l])
  # perform spectral clustering on the remaining matrix
  lab[l] <- spher.sp.cl(A[l, l], k, degree[l], tau)
  return(lab)
}


#########################################################

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

label.match <- function(lab, fixed)
{
  k <- max(lab, fixed)
  if(identical(lab, fixed)) return(list(lab, sparseMatrix(i = 1:k, j = 1:k, dims = c(k, k))))
  Z1 <- as(t(sapply(lab, make_basis, p = k)), "dgCMatrix")
  Z2 <- as(t(sapply(fixed, make_basis, p = k)), "dgCMatrix")
  P <- greedy.match(Z1, Z2, k)
  list(as.vector(tcrossprod(tcrossprod(Z1, t(P)), rbind(1:k))), P)
}
#########################################################

snow.sample <- function(A, N, deg, o)
{
  ind <- which(deg == 0)
  if(length(ind) == 0){
    init <- sample(1:N, o, F)
  }  else{
    init <- sample((1:N)[-ind], o, F)
  } 
  
  sapply(init, function(x){
    neighbor <- which(A[x, ] == 1)
    neighbor[sample(length(neighbor), 1, prob = deg[neighbor])]
  })
}

spSONNET <- function(A, k, tau = -1, s, O, rep = 0, ncore = 1)
{
  on.exit(gc())
  N <- nrow(A)    # detecting number of nodes
  m <- (N-O)/s 		# Make sure s|N-O
  #deg <- colSums(A)

  ## SUBGRAPHING STEP for both over and nonover with rep
  
  #over <- snow.sample(A, N, deg, O) #setting overlap part
  over <- sample(1:N, O, F)
  non <- sapply(1:(rep+1), function(x) sample((1:N)[-over], N-O, replace = F))
  # ^^ generating (rep+1) many partitions for nonoverlapping part
  
  output <- sapply(mclapply(1:(s*(rep+1)), function(i){
    r <- ceiling(i / s)   #detecting partition number
    pos <- i %% s + 1   #detecting position of the index
    index <- non[((pos-1)*m+1):(pos*m), r] #obtaining row number
    label <- rep(0, N) # initializing
    
    if(r != 1)  
    {
      label[index] <- spher.spectral.cluster(A[index, index], k, tau)
      return(label)
    }else{
      label[c(over,index)] <- spher.spectral.cluster(A[c(over, index),
                                                       c(over, index)], k, tau)
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


