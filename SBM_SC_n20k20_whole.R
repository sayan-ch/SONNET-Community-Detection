source('~/revision/code/SBM_SC_base.R')

whole.run <- function(i)
{
  on.exit(gc())
  net <- SBM.generate(N = 20000, k = 20, 
                      p.intra = 0.2, p.inter = 0.1,
                      ncore = 20)
  time <- system.time(
    out <- old.spectral.cluster(A = net$A, k = 20))[3]
  
  tmp <- c(
    time,
    mean(net$member != 
           label.match(lab = out, fixed = net$member)[[1]])
  )
  
  write.table(data.frame(time = as.numeric(tmp[1]), error = tmp[2]),
            '~/revision/output/sonnet/n20k20/TEMP_SBM_SC_n20k20_whole.csv',
            col.names = F, append = T, sep = ",")
  
  return(tmp)
}

sim <- 100


output <- sapply(mclapply(1:sim, whole.run, mc.cores = 20),
                 'rbind')

tab <- data.frame('N' = rep(20000, sim),
                  'k' = rep(20, sim),
                  'p.intra' = rep(0.2, sim),
                  'p.inter' = rep(0.1, sim),
                  'method' = rep('SC (unreg) on SBM', sim),
                  'time.raw' = output[1,],
                  'error.raw' = output[2,])

write.table(tab, '~/revision/output/sonnet/n20k20/SBM_SC_n20k20_whole.csv',
            col.names = T, sep = ",", append = T)







