source('~/revision/code/DCBM_RSC_BASE.R')

whole.run <- function(i)
{
  on.exit(gc())
  net <- DCBM.fast(N = 20000, k = 20, 
                   P = matrix(0.03,20,20) + diag(0.07,20),
                   Th = sample(1:100, 20000, T),
                      ncore = 20)
  time <- system.time(
    out <- spher.spectral.cluster(A = net$A, k = 20))[3]
  
  tmp <- c(
    time,
    mean(net$member != 
           label.match(lab = out, fixed = net$member)[[1]])
  )
  
  write.table(data.frame(time = as.numeric(tmp[1]), error = tmp[2]),
            '~/revision/output/sonnet/DCBM/n20k20/TEMP_DCBM_RSC_n20k20_whole.csv',
            col.names = F, append = T, sep = ",")
  
  return(tmp)
}

sim <- 100


output <- sapply(mclapply(1:sim, whole.run, mc.cores = 20),
                 'rbind')

tab <- data.frame('N' = rep(20000, sim),
                  'k' = rep(20, sim),
                  'p.intra' = rep(0.1, sim),
                  'p.inter' = rep(0.03, sim),
                  'method' = rep('RSC (reg) on DCBM', sim),
                  'time.raw' = output[1,],
                  'error.raw' = output[2,])

write.table(tab, '~/revision/output/sonnet/DCBM/n20k20/DCBM_RSC_n20k20_whole.csv',
            col.names = T, sep = ",", append = T)







