source('~/revision/code/SBM_SC_base.R')
library(data.table)
library(tidyverse)

f.run <- function(s.run, O.run, r.run)
{
  on.exit(gc())
  net <- SBM.generate(N = 20000, k = 20, 
                      p.intra = 0.2, p.inter = 0.1,
                      ncore = 1)
  time <- system.time(
    out <- oldSONNET(A = net$A, k = 20,
                     s = s.run, O = O.run, rep = r.run,
                     ncore = 20))[3]
  
  tmp <- c(
    time,
    mean(net$member != 
           label.match(lab = out, fixed = net$member)[[1]])
  )
  
  write_csv(data.frame(time = as.numeric(tmp[1]), error = tmp[2]),
              paste0(
                '~/revision/output/sonnet/n20k20/TEMP_SBM_SC_n20k20_SON_s',
                s.run,'O',O.run,'r',r.run,'.csv'),
               append = T)
  
  return(tmp)
}


sim <- 100

parlist <- list(
c(35,6280)
)


count <- 1
for(par in parlist)
    for(r in c(20, 10))
    {
      sim <- 100
      s <- par[1]
      O <- par[2]
      
      if(s==5 & O == 6605 & (r == 0 | r == 10)) next
      
      if(s==5 & O == 6605 & r == 20) sim <- 36
      
      sim <- 100
      
      output <- sapply(mclapply(1:sim, function(i){
        f.run(s.run = s, O.run = O, r.run = r)
      }, mc.cores = 1),
                       'rbind')
      
      
      tab <- data.table('N' = rep(20000, sim),
                        'k' = rep(20, sim),
                        'p.intra' = rep(0.2, sim),
                        'p.inter' = rep(0.1, sim),
                        'method' = rep('SONNET + SC (unreg) on SBM', sim),
                        'ncore' = 20,
                        's' = rep(s, sim),
                        'o' = rep(O, sim),
                        'r' = rep(r, sim),
                        'time.raw' = output[1,],
                        'error.raw' = output[2,])
      
      write_csv(tab, paste0(
        '~/revision/output/sonnet/n20k20/SBM_SC_n20k20_SON_s',s,'O',O,'r',r,'.csv'))
                  
      tab.mean <- data.table('N' = 20000,
                        'k' = 20,
                        'p.intra' = 0.2,
                        'p.inter' = 0.1,
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
        '~/revision/output/sonnet/n20k20/MEAN_SBM_SC_n20k20_SON_ALL_PARAM.csv'),
        append = T)
        
    count <- count+1
    }
    



