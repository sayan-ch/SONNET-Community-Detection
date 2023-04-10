source('~/revision/code/DCBM_RSC_BASE2.R')

whole.run <- function(i)
{
  on.exit(gc())
  message('Generation started')
  net <- DCBM.fast(N = 10000, k = 5, 
                      P = matrix(0.03,5,5)+diag(0.07,5),
                   Th = sample(1:100, 10000, T),
                      ncore = 20)
  message('Generation done - Spher cluster started')
  time <- system.time(
    out <- spher.spectral.cluster(A = net$A, k = 5, tau = -1))[3]
  message('spher cluster done')
  tmp <- c(
    time,
    mean(net$member != 
           label.match(lab = out, fixed = net$member)[[1]])
  )
  
  write.table(data.frame(time = as.numeric(tmp[1]), error = tmp[2]),
            '~/revision/output/sonnet/DCBM/n10k5/TEMP_DCBM_RSC_n10k5_whole.csv',
            col.names = F, append = T, sep = ",")
  
  return(tmp)
}

sim <- 100


output <- sapply(mclapply(1:sim, whole.run, mc.cores = 20),
                 'rbind')

tab <- data.frame('N' = rep(10000, sim),
                  'k' = rep(5, sim),
                  'p.intra' = rep(0.1, sim),
                  'p.inter' = rep(0.03, sim),
                  'method' = rep('RSC (reg) on DCBM', sim),
                  'time.raw' = output[1,],
                  'error.raw' = output[2,])

write.table(tab, '~/revision/output/sonnet/DCBM/n10k5/DCBM_RSC_n10k5_whole.csv',
            col.names = T, sep = ",", append = T)







