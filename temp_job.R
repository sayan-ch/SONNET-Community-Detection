source('~/revision/code/SBM_SC_base_REG.R')
library(tidyverse)

whole.run <- function(i)
{
  on.exit(gc())
  message('net gen start')
  net <- SBM.generate(N = 3000, k = 5, 
                      p.intra = 0.2, p.inter = 0.05,
                      ncore = 20)
  message('net gen done - reg spec start')
  time <- system.time(
    out <- reg.spectral.cluster(A = net$A, k = 5))[3]
  message('reg spec done - lab mat start')
  
  tmp <- c(
    time,
    mean(net$member != 
           label.match(lab = out, fixed = net$member)[[1]])
  )
  message('lab mat done - write start')
  
  write.table(data.frame(time = as.numeric(tmp[1]), error = tmp[2]),
            '~/revision/TEMP_TRIAL_REG_n3k5_whole.csv',
            col.names = F, append = T, sep = ",")
            
  message('write done')  
  
  return(tmp)
}

sim <- 100


output <- sapply(mclapply(1:sim, whole.run, mc.cores = 20),
                 'rbind')

tab <- data.frame('N' = rep(10000, sim),
                  'k' = rep(5, sim),
                  'p.intra' = rep(0.2, sim),
                  'p.inter' = rep(0.05, sim),
                  'method' = rep('SC (unreg) on SBM', sim),
                  'time.raw' = output[1,],
                  'error.raw' = output[2,])

write.table(tab, '~/revision/output/TRIAL_REG_n3k5_whole_Final.csv',
            col.names = T, sep = ",", append = T)
