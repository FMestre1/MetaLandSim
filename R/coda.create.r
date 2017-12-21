coda.create <- function(object,file.name,write.index=TRUE,par.list=list("mupsi1.chain","e.chain","x.chain","b.chain","y.chain","alpha.chain"),niter=101000,nthin=10) {

  p.list.insert <- function(mcmc.list){
    pchain = mcmc.list[[c("p.chain")]]
    lengthm = length(mcmc.list)
    for(i in 1:nrow(pchain)) {
      mcmc.list[[length(mcmc.list)+1]] = pchain[i,]
      names(mcmc.list)[length(mcmc.list)] = paste('p.',i,sep='')
    }
    mcmc.list
  }


  if(nthin!=1) message(paste('NOTE: This function assumes nthin=', nthin, ' has already been applied',sep=''))
  if("p.chain"%in%par.list){
    nyear = nrow(object$p.chain)
    object = p.list.insert(object)
    nvars = length(names(object))
    par.list = c(par.list,names(object)[(nvars-nyear+1):nvars])
  }
  iter.num.short=seq(1,niter,by=nthin)
  value=NULL
  par.list = par.list[par.list!="p.chain"]
  iter.num=rep(iter.num.short,length(par.list))

  for(i in 1:length(par.list)){
    value=c(value,object[[names(object)[names(object)==par.list[i]]]])
  }
  temp=data.frame(iter.num,value)
  write.table(temp,paste(file.name,".txt",sep=""),col.names=FALSE,row.names=FALSE,sep="\t")
  if(write.index==TRUE) {
      row.num=seq(1,length(iter.num),by=length(iter.num.short))
      row.num2=c((row.num[2:length(row.num)]-1),length(iter.num))
      write.table(data.frame(as.character(par.list),row.num,row.num2),file=paste(file.name,"_Index",".txt",sep=""),col.names=FALSE,row.names=FALSE,sep="\t")         }
}