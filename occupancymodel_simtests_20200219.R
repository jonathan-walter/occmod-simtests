## Create simulations to test performance of model 

rm(list=ls())

library(rjags)
library(coda)
library(stringr)

## Occupancy only ---------------------------------------------------------------------------------

#simulate a site x reps x spp array of occupancies (0,1) and test whether the model recovers it
#and the underlying parameters

nsites=20
nspp=50
nreps=5

logit<-function(x){
  return(log(x/(1/x)))
}

inv.logit<-function(x){
  return(exp(x)/(1+exp(x)))
}

elevsc<-c(scale(1:nsites))
beta<-rnorm(nspp,0.1,1)
#p.detect<-matrix(runif(nsites*nspp,0,1),nrow=nsites,ncol=nspp) #detection probs by site and species
  # consider treating this as a beta distribution so that there can be many rare spp and a few common
  # ones as is often the case in the real world

p.occur<-matrix(inv.logit(rnorm(n = nsites*nspp
                                , mean= beta*matrix(elevsc #This is a scaled variable but represents change in occurrence probability with elevation
                                              , nrow=nsites #each of these has an elevation attached
                                              , ncol=nspp # I guess this simply recycles the elevation for each species, which has its overall (average) probability of detection given by the beta variable
                                              ))) #sd for this is 1
                ,nsites
                ,nspp)
p.detect<-matrix(inv.logit(runif(nspp,-2,2)),nsites,nspp,byrow=T)
Xtrue<-matrix(NA, nsites,nspp) #initialize the true occupancy array
Xobs<-array(NA, dim=c(nsites,nreps,nspp))
for(site in 1:nsites){ #fill Z with detections
  for(spp in 1:nspp){
    Xtrue[site,spp]<-rbinom(1,1,p.occur[site,spp])
    Xobs[site,,spp]<-rbinom(5,1,p.detect[site,spp]*Xtrue[site,spp])
  }
}
Zobs<-apply(Xobs,c(1,3),function(x){as.numeric(any(x==1))})

#######################Run the model
nburn = 5000
n.iter = 10000
n.chains = 2
thin = 10

###Specify the parameters to be monitored
#sp.params = list("Z","mu.psi","theta","p.fit", "p.fitnew")
sp.params = list("mu.psi")
#Z matrix will store occupancy state
#mu.psi will store occupancy probabilities matrix
#mu.theta will store detection probabilities 
#p.fit and p.fitnew are required for Bayesian p value

sp.params <- as.character(sp.params)

sp.data = list(nspp=nspp, nsite=nsites, nrep=rep(nreps,nsites), X=Xobs, elev=elevsc)

#######################Specify the initial values


# Set 1 - this set of initial values has proven effective, so I'll be using these
sp.inits = function() {
  psi.meanGuess = runif(1,0.001,0.99)
  list(psi.mean=psi.meanGuess, theta.mean=runif(1,0.001,0.99),
       u=rnorm(nspp), v=rnorm(nspp),
       Z = Zobs)
}

# Set 2 - use these for dev2 version
# sp.inits = function() {
#   psi.meanGuess = runif(1,0.001,0.99)
#   list(psi.mean=psi.meanGuess, theta.mean=runif(1,0.001,0.99),
#        u=matrix(rnorm(nsites*nspp),nsites,nspp), 
#        v=matrix(rnorm(nsites*nspp),nsites,nspp),
#        Z = Zobs)
# }

# Set 3 - use these for the dev3 version

ocmod <- jags.model(file = "~/Documents/Research/DATA/BBS/DetectionCorrection/Multisp_model_dev3.txt", inits = sp.inits, data = sp.data, n.chains = n.chains)
update(ocmod, n.iter = nburn)
out <- coda.samples(ocmod, n.iter = n.iter, variable.names = sp.params, thin=thin)
out.mcmc <- as.mcmc(out[[1]])

plot(out)
H<-heidel.diag(out)
print(H)
#G<-gelman.diag(out)

#Z<-out.mcmc[,grepl("Z",colnames(out.mcmc))]
# theta<-out.mcmc[,grepl("theta",colnames(out.mcmc))]
# theta.med<-array(apply(theta,2,mean),dim=c(nsites,nreps,nspp))
# plot(c(p.detect),c(theta.med[,1,]))
# cor(c(p.detect),c(theta.med[,1,]))
# 
mu.psi<-out.mcmc[,grepl("mu.psi",colnames(out.mcmc))]
psi.med<-array(apply(mu.psi,2,median),dim=c(nsites,nreps,nspp))
plot(c(p.occur),c(psi.med[,1,]))
# 
# Z1<-matrix(Z[1,],nrow=nsites,ncol=nspp)
# Zavg<-matrix(colMeans(Z),nrow=nsites,ncol=nspp)
# Zmed<-matrix(apply(Z,2,median),nrow=nsites,ncol=nspp)
# sum(abs(Zobs-Zmed))
# sum(abs(Xtrue-Zmed))

## mu.theta is more bimodally distributed than p.detect, since it's estimated from a binomial dist'n
## Z is very similar to Zobs; the median of Z is perfect

