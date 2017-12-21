log.lik.data.z.p<-function(ntrial.temp,nsucc.temp,z.temp,p.temp)  sum(dbinom(nsucc.temp,ntrial.temp,z.temp*p.temp,log=TRUE),na.rm=TRUE)
