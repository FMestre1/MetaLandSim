\name{simulatedifm}
\alias{simulatedifm}
\docType{data}
\title{
Set of simulated data to use with the IFM parameter estimation functions. The data were generated using the code provided in "details".
}
\description{
This dataset loads several objects:

'sim.area', 'sim.det.20', 'sim.distance', 'z.sim', 'z.sim.20' and 'z.sim.20.fa'.
}
\format{
  \describe{
    \item{\code{sim.area}}{Vector of the areas for each site; here, 100 sites.}
    \item{\code{sim.det.20}}{100 x 10 x 3 array corresponding to nsites x nyears x nvisits. Data simulated with year-specific detection probabilities equal to 0.4,0.6,0.2,0.9,0.3,0.4,0.6,0.2,0.9,0.3.}
    \item{\code{sim.distance}}{nsite x nsite distance matrix.}
    \item{\code{z.sim}}{nyear x nsite occupancy data generated with perfect detection.}
    \item{\code{z.sim.20}}{nyear x nsite occupancy data generated with perfect detection with approximately 20\% of data missing at random.}
    \item{\code{z.sim.20.fa}}{nyear x nsite occupancy data containing false absences, which can be used to explore the bias of ifm.missing.MCMC and ifm.naive.MCMC when there is imperfect detection.}
  }
}
\details{
These datasets were created using the code included in \bold{Examples}.
}
\examples{

\dontrun{

##############################################
# Areas for 100 hundred sites were created from a log normal distribution with mean
# and variance equal to the mean and variance of the area of the sites in the
# Sierra Foothills black rail population.
# Universal Transverse Mercator locations (UTMs) for each site were
# simulated from the mean and variance of the UTM Northing and Easting
#
# Dynamics were simulated for 1000 years. The parameters were chosen such that the metapopulation
# persisted with reasonable turnover. The last 10 years of these data were retained.
# For the detection data sets, we simulated a removal design based on three visits.

#-------------------
# IFM SIMULATE
#------------------
set.seed(123)

#-------
#AREAS
#-------
mean.log.area=-0.75
sd.log.area=1.33
nsite=100
log.sim.area=rnorm(nsite,mean.log.area,sd.log.area)
sim.area = exp(log.sim.area)

#---------
#DISTANCE
#--------
sim.site.x=rnorm(nsite,643930,9000)
sim.site.y=rnorm(nsite,4340500,10500)

sim.site.x.mat.col=matrix(rep(sim.site.x,nsite),ncol=nsite,byrow=TRUE)
sim.site.x.mat.row=matrix(rep(sim.site.x,nsite),ncol=nsite)
sim.site.delta.x=(sim.site.x.mat.col-sim.site.x.mat.row)^2

sim.site.y.mat.col=matrix(rep(sim.site.y,nsite),ncol=nsite,byrow=TRUE)
sim.site.y.mat.row=matrix(rep(sim.site.y,nsite),ncol=nsite)
sim.site.delta.y=(sim.site.y.mat.col-sim.site.y.mat.row)^2

# scale distance so that alpha = 2 is reasonable:
sim.distance=((sim.site.delta.x+sim.site.delta.y)^0.5)/100000
diag(sim.distance)=99

#-------------#
# CREATE SPOM #
#-------------#
spom=function(nsite,nyear.sim,alpha,b,y,e,x) {
	psi=matrix(rep(NA,nsite*nyear.sim),ncol=nyear.sim)
	psi1=rbinom(nsite,1,0.8)
	psi[,1]=psi1
	s.i.temp=exp(-alpha*sim.distance)
	s.i.temp[s.i.temp==1]=0
	e.i=e/sim.area^x
	e.i[e.i>1]=1

	for (i in 2:nyear.sim) {
		s.i=s.i.temp%*%((sim.area^b)*psi[,i-1])
		c.i=s.i^2/(s.i^2+y^2)
		e.i.re=e.i*(1-c.i)
		mu1=psi[,i-1]*(1-e.i.re)+(1-psi[,i-1])*c.i
		psi.temp=rbinom(nsite,1,mu1)
		psi[,i]=psi.temp
		}
	psi
	}

#---------------
# SIMULATE IFM
#---------------

nyear.sim=1000
psi.sim=spom(nsite,nyear.sim,alpha=20,b=0.5,y=7.5,e=0.25,x=0.25)

nyear=10

# Data from this dataset conforms to the assumptions of IFM Naive:
z.sim=psi.sim[,(nyear.sim+1-nyear):nyear.sim]
apply(z.sim,2,mean)


#CREATE DETECTION HISTORY
p=rep(c(0.4,0.6,0.2,0.9,0.3),2)
nrep=3
temp=rep(1,nrep*nsite*nyear)
p.mat=matrix(rep(p,nsite),nrow=nsite,byrow=TRUE)
temp.z.sim=z.sim*p.mat
sim.det=rbinom(temp,1,temp.z.sim)
dim(sim.det)=c(nsite,nyear,nrep)

#ENFORCE REMOVAL DESIGN
for (i in 1:nsite) {
  for(t in 1:nyear) {
    if (sim.det[i,t,1]==1) sim.det[i,t,2]=2
    if (sim.det[i,t,1]==1) sim.det[i,t,3]=2
    if (sim.det[i,t,2]==1) sim.det[i,t,3]=2
  }
}

sim.det[sim.det==2]=NA
sim.det.no.missing.values=sim.det

#RANDOMLY CREATE MISSING DATA; 20%
# Data are missing when a site was never visited in a given year
unif.mat=runif(nyear*nsite)
z.sim.20=z.sim
z.sim.20[unif.mat<0.2]=NA

sim.det.20=sim.det
sim.det.20[rep(is.na(z.sim.20),nrep)]=NA

#CREATE DATASET WITH MISSING VALUES AND FALSE ABSENCES
#THIS IS TO CHECK HOW IFM.NAIVE AND IFM.MISSING
# LEAD TO BIASES
z.sim.20.fa = apply(sim.det.20,c(1,2),sum,na.rm=TRUE)
z.sim.20.fa[unif.mat<0.20]=NA

# z.sim: Perfect detection, one visit per year.
# z.sim.20: Perfect detection, but 20% of sites have missing data.
# z.sim.20.fa: 20% of sites have missing data in each year, sites are visited according to the removal design,
#   then collapsed to a single observation per year equal to one if a detection occurred.
# sim.det.20: 20% of sites have missing data, there is imperfect detection, and the data follow a removal design.
#   The data are arranged in a 3-d array sites x years x visits

save(z.sim,z.sim.20,z.sim.20.fa,sim.det.
20,sim.distance,sim.area,file=paste("SIMULATE_DATA_MDL_CMP",nsite,"_",nyear,
"_",nrep,".RData",sep=""))
}

data(simulatedifm)
ls()

}
\keyword{datasets}
\keyword{ifm}
