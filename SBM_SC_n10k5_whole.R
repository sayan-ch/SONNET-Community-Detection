source('~/revision/code/SBM_SC_base.R')

whole.run <- function(i)
{
  on.exit(gc())
  net <- SBM.generate(N = 10000, k = 5, 
                      p.intra = 0.2, p.inter = 0.05,
                      ncore = 20)
  time <- system.time(
    out <- old.spectral.cluster(A = net$A, k = 5))[3]
  
  tmp <- c(
    time,
    mean(net$member != 
           label.match(lab = out, fixed = net$member)[[1]])
  )
  
  write.table(data.frame(time = as.numeric(tmp[1]), error = tmp[2]),
            '~/revision/output/NEW_SBM_SC_PROC_n10k5_whole.csv',
            col.names = F, append = T, sep = ",")
  
  return(tmp)
}

sim <- 5


output <- sapply(mclapply(1:sim, whole.run, mc.cores = 20),
                 'rbind')

tab <- data.frame('N' = rep(10000, sim),
                  'k' = rep(5, sim),
                  'p.intra' = rep(0.2, sim),
                  'p.inter' = rep(0.05, sim),
                  'method' = rep('SC (unreg) on SBM', sim),
                  'time.raw' = output[1,],
                  'error.raw' = output[2,])

write.table(tab, '~/revision/output/NEW_SBM_SC_PROC_n10k5_whole_Final.csv',
            col.names = T, sep = ",", append = T)







