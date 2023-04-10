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
  print('Laplace done-eigen started')
  # computing eigen decomposition
  #E <- eigen(L, T)$vectors[,1:k]
  #E1 <- t(apply(E,1,function(x){x/sum(x^2)^0.5}))
  Ei <- eigen(L, T)
  #ord <- sort.list(abs(Ei$values), decreasing = T)[1:k]
  #ord <- ord[sort.list(Ei$values[ord], decreasing = T)]
  ord <- c(1:(k/2), (N-k/2+1):N)
  if(k%%2){
  k1 <- floor(k/2)
  if(abs(Ei$values[k1+1]) > abs(Ei$values[N-k1]))
    ord <- c(1:(k1+1), (N-k1+1):N)
  else
    ord <- c(1:k1, (N-k1):N)
    }
  E <- abs(Ei$vectors[,ord])
  print('eigen done.')
  rss <- rowSums(E^2)
  print(paste0('rss = ', sum(rss==0)))
  E1 <- E/sqrt(rss)
  E1[rss==0,] <- 0
  print('kmeans started.')
  out <- kmeans(E1,k,10000000,1000)
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

oldSONNET <- function(A, k, s, O, rep = 0, ncore = 1)
{
  on.exit(gc())
  N <- nrow(A)    # detecting number of nodes
  m <- (N-O)/s 		# Make sure s|N-O
  degree <- colSums(A)
  
  ## SUBGRAPHING STEP for both over and nonover with rep
  
  if(rep == 0)
    over <- sample(1:N, O, replace = F,
                 prob = degree) #setting overlap part
  else
    over <- sample(1:N, O, replace = F)
    
  #over <- snow.sample(A, N, O, degree)
  non <- sapply(1:(rep+1), function(x) sample((1:N)[-over], N-O, replace = F))
  # ^^ generating (rep+1) many partitions for nonoverlapping part
  
  output <- sapply(mclapply(1:(s*(rep+1)), function(i){
    r <- ceiling(i / s)   #detecting partition number
    pos <- i %% s + 1   #detecting position of the index
    index <- non[((pos-1)*m+1):(pos*m), r] #obtaining row number
    label <- rep(0, N) # initializing
    
    if(r != 1)  
    {
      label[index] <- spectral.cluster(A[index, index], k)
      return(label)
    }else{
      label[c(over,index)] <- spectral.cluster(A[c(over, index),
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

################################################################################
################################################################################
## Twitch Code


sim <- 5
for(itera in 1:sim)
for(r in c(0,5)){
edges <- read_csv('~/revision/data/bigadj.csv')
comm <- read_csv('~/revision/data/bigcomm.csv')

A <- sparseMatrix(i = edges$user1, j = edges$user2,
                  dims = c(32407, 32407),
                  symmetric = T)

# s = 2, O = 8425, data.use = 72.62%

print('s=2, o = 8425')
time1 <- system.time({
  outson1 <- oldSONNET(A, 20, 2, 8425, r, 20)
})[3]

error1 <- mean(comm$langnum != label.match(outson1, comm$langnum)[[1]])

write_csv(data.table(s = 2, o = 8425, r = r, time = time1, error = error1),
          paste0('~/revision/output/data/twitch/tw2_s2o8425r',r,'.csv'),
          append = T)
gc()
# s = 8, O = 12879, data.use = 68.23%

print('s=8, o = 12879')
time2 <- system.time({
  outson2 <- oldSONNET(A, 20, 8, 12879, r, 20)
})[3]

error2 <- mean(comm$langnum != label.match(outson2, comm$langnum)[[1]])

write_csv(data.table(s = 8, o = 12879, r = r, time = time2, error = error2),
          paste0('~/revision/output/data/twitch/tw2_s8o12879r',r,'.csv'),
          append = T)

gc()
# s = 5, O = 6002, data.use = 46.89%

print('s=5, o = 6002')
time3 <- system.time({
  outson3 <- oldSONNET(A, 20, 5, 6002, 0, 20)
})[3]

error3 <- mean(comm$langnum != label.match(outson3, comm$langnum)[[1]])

write_csv(data.table(s = 5, o = 6002, r = r, time = time3, error = error3),
          paste0('~/revision/output/data/twitch/tw2_s5o6002r',r,'.csv'),
          append = T)


gc()
# s = 8, O = 6007, data.use = 41.93%

print('s=8, o = 6007')
time4 <- system.time({
  outson4 <- oldSONNET(A, 20, 8, 6007, 0, 20)
})[3]

error4 <- mean(comm$langnum != label.match(outson4, comm$langnum)[[1]])

write_csv(data.table(s = 8, o = 6007, r = r,  time = time4, error = error4),
          paste0('~/revision/output/data/twitch/tw2_s8o6007r',r,'.csv'),
          append = T)

gc()
rm(error1, error2, error3, error4,
time1, time2, time3, time4,
outson1, outson2, outson3, outson4,
A, edges, comm)
}





























