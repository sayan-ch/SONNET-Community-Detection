source('~/revision/code/GALE_SC_BASE.R')


print("Base File Worked")
f.run <- function(s.run, m.run)
{
  on.exit(gc())
  net <- SBM.generate(N = 10000, k = 5, 
                      p.intra = 0.2, p.inter = 0.05,
                      ncore = 20)
  time <- system.time(
    out <- GALE(A = net$A, k = 5,
                     s = s.run, m = m.run,
                     ncore = 20))[3]
  
  tmp <- c(
    time,
    mean(net$member != 
           label.match(lab = out, fixed = net$member)[[1]])
  )
  
  write.table(data.frame(time = as.numeric(tmp[1]), error = tmp[2]),
              paste0(
                '~/revision/output/gale/TEMP_SBM_SC_n10k5_GALE_t',
                s,'m',m,'.csv'),
              col.names = F, append = T, sep = ",")
  
  return(tmp)
}


sim <- 100

for(s in c(50, 1000))
  for(m in c(300, 500, 1000))
    {
      output <- sapply(mclapply(1:sim, function(i){
        f.run(s.run = s, m.run = m)
      }, mc.cores = 1),
      'rbind')
      
      
      tab <- data.frame('N' = rep(10000, sim),
                        'k' = rep(5, sim),
                        'p.intra' = rep(0.2, sim),
                        'p.inter' = rep(0.05, sim),
                        'method' = rep('GALE + SC (unreg) on SBM', sim),
                        'ncore' = 20,
                        't' = rep(s, sim),
                        'm' = rep(m, sim),
                        'time.raw' = output[1,],
                        'error.raw' = output[2,])
      
      write.table(tab, paste0(
        '~/revision/output/gale/SBM_SC_n10k5_GALE_t',s,'m',m,'.csv'),
        col.names = T, sep = ",", append = T)
    }







