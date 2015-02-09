as.coda.demc <-function(Draws.DEMC){
  require(coda)
  Draws =Draws.DEMC$Draws
  Nchain= Draws.DEMC$Nchain
  Chain= list()
  for (i in seq_len(Nchain)){
    chain.i = seq(from=i, to = ncol(Draws), by=Nchain)
    Chain[[i]]=as.mcmc(t(Draws[,chain.i]))
  }
  mcmc.list(Chain)
}