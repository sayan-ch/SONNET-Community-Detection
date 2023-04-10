source('~/revision/code/SBM_SC_base.R')

f.run <- function(s.run, O.run, r.run)
{
  on.exit(gc())
  net <- SBM.generate(N = 10000, k = 5, 
                      p.intra = 0.2, p.inter = 0.05,
                      ncore = 1)
  time <- system.time(
    out <- oldSONNET(A = net$A, k = 5,
                     s = s.run, O = O.run, rep = r.run,
                     ncore = 20))[3]
  
  tmp <- c(
    time,
    mean(net$member != 
           label.match(lab = out, fixed = net$member)[[1]])
  )
  
  write.table(data.frame(time = as.numeric(tmp[1]), error = tmp[2]),
              paste0(
                '~/revision/output/sonnet/TEMP_SBM_SC_n10k5_SON_s',
                s,'O',O,'r',r,'.csv'),
              col.names = F, append = T, sep = ",")
  
  return(tmp)
}


sim <- 100

for(s in c(10,20,50))
  for(O in c(100, 500, 1000))
    for(r in c(0, 2, 5))
    {
      output <- sapply(mclapply(1:sim, function(i){
        f.run(s.run = s, O.run = O, r.run = r)
      }, mc.cores = 1),
                       'rbind')
      
      
      tab <- data.frame('N' = rep(10000, sim),
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
      
      write.table(tab, paste0(
        '~/revision/output/sonnet/SBM_SC_n10k5_SON_s',s,'O',O,'r',r,'.csv'),
                  col.names = T, sep = ",")
    }







