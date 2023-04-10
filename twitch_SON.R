source('~/revision/code/SBM_SC_base.R')

library(tidyverse)
library(data.table)


edge <- read_csv('~/revision/data/twitch_edges.csv')
comm <- read_csv('~/revision/data/twitch_comm.csv')

A <- sparseMatrix(i = edge$user1,
                  j = edge$user2,
                  symmetric = T)

f.run <- function(s.run, O.run, r.run)
{
  on.exit(gc())
  
  time <- system.time(
    out <- oldSONNET(A = A, k = 20,
                     s = s.run, O = O.run, rep = r.run,
                     ncore = 20))[3]
  print(system.time(
  tmp <- c(
    time,
    mean(comm$langnum != 
           label.match(lab = out, fixed = comm$langnum)[[1]])
  )
  ))
  
  write_csv(data.table(time = as.numeric(tmp[1]), error = tmp[2]),
              paste0(
                '~/revision/output/data/TEMP_twitch_SON_s',
                s,'O',O,'r',r,'.csv'),
            append = T)
  
  return(tmp)
}


sim <- 5

parlist <- list(
c(2,8425)
)

for(r in c(0))
for(par in parlist){
      s <- par[1]
      O <- par[2]
      #cut <- par[3]
      
      output <- sapply(mclapply(1:sim, function(i){
        f.run(s.run = s, O.run = O, r.run = r)
      }, mc.cores = 1),
                       'rbind')
      
      
      tab <- data.table('N' = rep(32407, sim),
                        'k' = rep(20, sim),
                        'method' = rep('SONNET + SC + reg', sim),
                        'ncore' = 20,
                        's' = rep(s, sim),
                        'o' = rep(O, sim),
                        'r' = rep(r, sim),
                        'time.raw' = output[1,],
                        'error.raw' = output[2,])
      
      write_csv(tab, paste0(
        '~/revision/output/data/twitch_SON_s',s,'O',O,'r',r,'.csv'),
        append = T)
        
      tab.mean <- data.table('N' = 32407,
                        'k' = 20,
                        'method' = 'SONNET + SC + reg',
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
        '~/revision/output/data/MEAN_twitch_SON_ALL_PARAM.csv'),
        append = (count != 1))
    }







