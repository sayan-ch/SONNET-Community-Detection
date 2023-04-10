source('~/revision/code/SBM_SC_base.R')
library(data.table)
library(tidyverse)

f.run <- function(s.run, O.run, r.run)
{
  on.exit(gc())
  net <- SBM.generate(N = 10000, k = 5, 
                      p.intra = 0.2, p.inter = 0.05,
                      ncore = 20)
  time <- system.time(
    out <- oldSONNET(A = net$A, k = 5,
                     s = s.run, O = O.run, rep = r.run,
                     ncore = 20))[3]
  
  tmp <- c(
    time,
    mean(net$member != 
           label.match(lab = out, fixed = net$member)[[1]])
  )
  
  write_csv(data.frame(time = as.numeric(tmp[1]), error = tmp[2]),
              paste0(
                '~/revision/output/sonnet/TEMP_SBM_SC_n10k5_SON_s',
                s,'O',O,'r',r,'.csv'),
               append = T)
  
  return(tmp)
}


sim <- 100

parlist <- list(c(3,1309),
c(5,905),
c(12,400),
c(17,276),
c(41,160)
)

#for(s in c(3,5,12,17,41))
#  for(O in c(160,276,905))
count <- 1
for(par in parlist)
    for(r in c(0, 2, 5))
    {
      s <- par[1]
      O <- par[2]
      output <- sapply(mclapply(1:sim, function(i){
        f.run(s.run = s, O.run = O, r.run = r)
      }, mc.cores = 1),
                       'rbind')
      
      
      tab <- data.table('N' = rep(10000, sim),
                        'k' = rep(5, sim),
                        'p.intra' = rep(0.2, sim),
                        'p.inter' = rep(0.05, sim),
                        'method' = rep('SONNET + SC (unreg) on SBM', sim),
                        'ncore' = 20,
                        's' = rep(s, sim),
                        'o' = rep(O, sim),
                        'r' = rep(r, sim),
                        'time.raw' = output[1,],
                        'error.raw' = output[2,])
      
      write_csv(tab, paste0(
        '~/revision/output/sonnet/SBM_SC_n10k5_SON_s',s,'O',O,'r',r,'.csv'))
                  
      tab.mean <- data.table('N' = 10000,
                        'k' = 5,
                        'p.intra' = 0.2,
                        'p.inter' = 0.05,
                        'method' = 'SONNET + SC (unreg) on SBM',
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
        '~/revision/output/sonnet/MEAN_SBM_SC_n10k5_SON_ALL_PARAM.csv'),
        append = (count != 1))
        
    count <- count+1
    }
    



