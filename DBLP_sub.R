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

oldSONNET <- function(A, k, s, O, rep = 0, ncore = 1)
{
  on.exit(gc())
  N <- nrow(A)    # detecting number of nodes
  m <- (N-O)/s 		# Make sure s|N-O
  degree <- colSums(A)
  
  ## SUBGRAPHING STEP for both over and nonover with rep
  
  over <- sample(1:N, O, replace = F,
                 prob = degree) #setting overlap part
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



edges <- read_csv('~/revision/data/dblp-ac-compact.csv')
comm <- read_csv('~/revision/data/author_label.csv')

truth <- comm$a.lab+1

A <- sparseMatrix(i = edges$V1, j = edges$a,
                  dims = c(4057, 4057),
                  symmetric = T)

sim <- 100

parlist <- list(
c(59,104),
c(17,113),
c(12,169),
c(5,367),
c(3,532)
)

for(i in 1:sim)
for(r in c(0,2,5))
for(par in parlist){
tryCatch({
  s <- par[1]
  o <- par[2]
  
  time <- system.time(tmp <- oldSONNET(A = A, k = 4,
            s = s, O = o, rep = r, ncore = 20))[3]
            
  error <- mean(truth != label.match(tmp, truth)[[1]])
  
  write_csv(data.table(time, error),
            paste0('~/revision/output/data/DBLP/DBLP_s',
                    s, 'O', o, 'r', r,'.csv'),
            append = T
            )
  gc()},
  error = function(e){cat("ERROR :",conditionMessage(e),"\n")})
}






























