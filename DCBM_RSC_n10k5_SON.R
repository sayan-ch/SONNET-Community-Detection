source('~/revision/code/DCBM_RSC_BASE.R')
library(data.table)
library(tidyverse)

f.run <- function(s.run, O.run, r.run)
{
  on.exit(gc())
  net <- DCBM.fast(N = 10000, k = 5, 
                      P = matrix(0.03,5,5)+diag(0.07,5),
                      Th = sample(1:100, 10000, T),
                      ncore = 20)
  time <- system.time(
    out <- spSONNET(A = net$A, k = 5, tau = -1,
                     s = s.run, O = O.run, rep = r.run,
                     ncore = 20))[3]
  
  tmp <- c(
    time,
    mean(net$member != 
           label.match(lab = out, fixed = net$member)[[1]])
  )
  
  write_csv(data.frame(time = as.numeric(tmp[1]), error = tmp[2]),
              paste0(
                '~/revision/output/sonnet/DCBM/n10k5/TEMP_DCBM_RSC_n10k5_SON_s',
                s.run,'O',O.run,'r',r.run,'.csv'),
               append = T)
  
  return(tmp)
}


sim <- 100

parlist <- list(
c(5,905)
)

#for(s in c(3,5,12,17,41))
#  for(O in c(160,276,905))
count <- 1
for(r in c(2, 5))
for(par in parlist){
      s <- par[1]
      O <- par[2]
      output <- sapply(mclapply(1:sim, function(i){
        f.run(s.run = s, O.run = O, r.run = r)
      }, mc.cores = 1),
                       'rbind')
      
      
      tab <- data.table('N' = rep(10000, sim),
                        'k' = rep(5, sim),
                        'p.intra' = rep(0.1, sim),
                        'p.inter' = rep(0.03, sim),
                        'method' = rep('SONNET + RSC (reg) on DCBM', sim),
                        'ncore' = 20,
                        's' = rep(s, sim),
                        'o' = rep(O, sim),
                        'r' = rep(r, sim),
                        'time.raw' = output[1,],
                        'error.raw' = output[2,])
      
      write_csv(tab, paste0(
        '~/revision/output/sonnet/DCBM/n10k5/DCBM_RSC_n10k5_SON_s',s,'O',O,'r',r,'.csv'))
                  
      tab.mean <- data.table('N' = 10000,
                        'k' = 5,
                        'p.intra' = 0.1,
                        'p.inter' = 0.03,
                        'method' = 'SONNET + RSC (reg) on DCBM',
                        'ncore' = 20,
                        's' = s,
                        'o' = O,
                        'r' = r,
                        'time.mean' = round(mean(output[1,]), 4),    
                        'time.se' = round(sd(output[1,])/sqrt(sim), 4),
                        
                        'error_mean' = round(mean(output[2,]), 4),
                        
                        'error.se' = round(sd(output[2,])/sqrt(sim), 4)
                        )
                        
    write_csv(tab.mean, paste0(
        '~/revision/output/sonnet/DCBM/n10k5/MEAN_DCBM_RSC_n10k5_SON_ALL_PARAM.csv'),
        append = T)
        
    count <- count+1
    }
    



