accept.calculate<-function(x,model=c('naive','missing','robust')) {

  if(model=='naive'){
    par.list=list(e.chain=x$e.chain,x.chain=x$x.chain,y.chain=x$y.chain,b.chain=x$b.chain,alpha.chain=x$alpha.chain,deviance.chain=x$deviance.chain)
    acc.list=list(e.chain=NULL,x.chain=NULL,y.chain=NULL,b.chain=NULL,alpha.chain=NULL,deviance.chain=NULL)
      }
  if(model=='missing'){
    par.list=list(muz.missing.chain=x$muz.missing.chain,mupsi1.chain=x$mupsi1.chain,e.chain=x$e.chain,x.chain=x$x.chain, y.chain=x$y.chain,b.chain=x$b.chain,alpha.chain=x$alpha.chain)
    acc.list=list(muz.missing.chain=NULL,mupsi1.chain=NULL,e.chain=NULL,x.chain=NULL, y.chain=NULL,b.chain=NULL,alpha.chain=NULL)
  }
  if(model=='robust'){
    par.list=list(latent.deviance.chain = x$latent.deviance.chain, muz.chain=x$muz.chain,muz.missing.chain=x$muz.missing.chain,p.chain=x$p.chain,mupsi1.chain=x$mupsi1.chain,e.chain=x$e.chain,x.chain=x$x.chain, y.chain=x$y.chain, b.chain=x$b.chain, alpha.chain=x$alpha.chain,latent.deviance.chain=x$latent.deviance.chain)
    acc.list=list(latent.deviance.chain = NULL, muz.chain=NULL,muz.missing.chain=NULL,p.chain=NULL,mupsi1.chain=NULL,e.chain=NULL,x.chain=NULL, y.chain=NULL, b.chain=NULL, alpha.chain=NULL,latent.deviance.chain=NULL)
  }

  for (i in 1:length(par.list)) {
      acc.list[[i]]=acceptance.rate(par.list[[i]])
  }
  return(acc.list)
}
