######### This is code written by JAW 
# to simulate data to test the bird occupancy model developed by Jarzyna and Jetz, we think the one described in
# Jarzyna MA and W Jetz. 2016. Detecting the multiple facets of biodiversity. Trends in Ecology and Evolution, 31(7):527-538
# c.f. https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/e1f60939/
# also take a look at Zipkin

# some questions to answer:

# 1) Did J and J actually do this right?
# 2) for a given set of occupancies and detections, how many sites and reps do you need to get the job done?
# 3) right now, it seems like there might be more variation between spp (which have no covariates) than between sites (which have the elevation gradient). 
  # Will it ever work under these conditions? Does it get harder b/c of this?
# 4) some technical q's: 
  # what is our favorite operational convergence criterion? 
    # R_hat, Gelman-Rubin, Heidler, visual inspection
  # what should trace plots really look like
  # will the convergence problems go away with enough time or are they fundamentally about too few data?
  # How realistic are the simulations? can we add features of realism (i.e. from the BBS data) back in (those that might help model)? 
#load libraries and name a few functions

rm(list=ls())

library(rjags)
library(coda)
library(stringr)

library(R2jags) #I added this one for jags.paralllel which lets you run the mcmc chains in parallel
library(mcmcse)
library(tictoc)

logit<-function(x){
  return(log(x/(1/x)))
}

inv.logit<-function(x){
  return(exp(x)/(1+exp(x)))
}

## Create simulations to test performance of model 

## Occupancy only ---------------------------------------------------------------------------------

#simulate a site x reps x spp array of occupancies (0,1) and test whether the model recovers it
#and the underlying parameters

nsites=20
nspp=50
nreps=50 #In the Jarzyna and Jetz material it sounded like there were actually 50 stops anwyays!
elevsc<-c(scale(1:nsites)) #sets it so that sites are ordered from lowest to higest, and expressed in units of SD so from -1.6 yo 1.6 or so
beta<-rnorm(nspp,0.1,1)/10 #I just divided by 10 since the number of stops instead of sections #This is a scaled variable but represents change in occurrence probability with elevation

#p.detect<-matrix(runif(nsites*nspp,0,1),nrow=nsites,ncol=nspp) #detection probs by site and species
  # consider treating this as a beta distribution so that there can be many rare spp and a few common
  # ones as is often the case in the real world

#This is the actual simulated probability of occurences for each species at each site
p.occur<-matrix(inv.logit #this is used to rescale everything to 0,1
                (rnorm(n = nsites*nspp #I guess rnorm gives random deviates from the species beta *elevation, with sd 1
                          , mean= beta*matrix(elevsc 
                                        , nrow=nsites #each of these has an elevation attached, elevsc
                                        , ncol=nspp # I guess this simply recycles the elevation for each species, 
                                        # which each have their own overall (average) probability of detection given by the beta variable
                                        ))) 
                ,nsites #i think this just means have sites as rows and species as columns
                ,nspp)

p.sum<-apply(p.occur, 2, sum)/20 #average probability of occurrence across sites for each sp.
min(p.sum) #0.3, still high
max(p.sum) #0.6, not that much higher! Nothing has low occurence probability. Hmmm. 

p.detect<-matrix(inv.logit(runif(nspp,-2,2)),nsites,nspp,byrow=T) 

# p.detect #detection does not vary between sites. It ranges from 

min(p.detect)  #0.12 to 
max(p.detect)  #0.87


#There is much more variation between species in detection than there is in occurence, I think. 
# This is probably bad for the model esp. b/c the covariate is on occurence
Xtrue<-matrix(NA, nsites,nspp) #initialize the true occupancy array
Xobs<-array(NA, dim=c(nsites,nreps,nspp))
for(site in 1:nsites){ #fill Z with detections
  for(spp in 1:nspp){
    Xtrue[site,spp]<-rbinom(1,1,p.occur[site,spp]) #occupancy is 1 bernouli trial 
    Xobs[site,,spp]<-rbinom(5,1,p.detect[site,spp]*Xtrue[site,spp]) #detection is 5 bernouli trials times occurence
  }
}

Zobs<-apply(Xobs,c(1,3),function(x){as.numeric(any(x==1))}) #Z is the matrix of presences at the "site" level, not the repeat observation level

#######################Run the model
nburn = 5000
niter = 10000
nchains = 3 #I learned to do at least 3 to assess convergence
thin = 100 # this is to deal with autocorrelation in the chains, can come back to it

###Specify the parameters to be monitored
sp.params = list("Z","mu.psi","mu.theta","p.fit", "p.fitnew")
# sp.params = list("mu.psi")
#Z matrix will store occupancy state
#mu.psi will store occupancy probabilities matrix
#mu.theta will store detection probabilities 
#p.fit and p.fitnew are required for Bayesian p value

sp.params <- as.character(sp.params)

sp.data = list(nspp=nspp, nsite=nsites, nrep=rep(nreps,nsites), X=Xobs, elev=elevsc)

#######################Specify the initial values


# Set 1 - this set of initial values has proven effective, so I'll be using these
#MR: proven effective as in in your own experiments?
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

#I think this will fit the model, I think it will do the chains in series lets see if I can find how I did it in parallel
#trying with R2jags::jags.parallel,b ut this meansa  different models epcificiation
source("Multisp_model_dev3.R")

ocmod <- jags.parallel(data = sp.data
                       , inits = sp.inits
                       , parameters.to.save = sp.params
                       , model.file = jaw_model
                       # have to include object names to export to cluster, 
                       # I determined the membership of this list by trial and error message
                       
                       , export_obj_names = list("nburn", "niter", "nchains", "thin", "Zobs") 
                       , n.chains = nchains
                       , n.iter = niter
                       , n.burnin = nburn
                       , n.thin = thin
) #~/Documents/Research/DATA/BBS/DetectionCorrection/Multisp_model_dev3.txt")

# ocmod <- jags.model(file = "Multisp_model_dev3.txt"
#                     , inits = sp.inits
#                     , data = sp.data
#                     , n.chains = n.chains) #~/Documents/Research/DATA/BBS/DetectionCorrection/Multisp_model_dev3.txt")


traceplot(ocmod, varname="mu.psi") #I think this might not be convergence, seems like it's bucking about wildly.
#I don't think the update is paralllelized.
recompile(ocmod)
tic()
is_upd_parll<-autojags(ocmod)
toc() #this took about 2 minutes and I think it did 3 chains for whatever number of iterations 2 times. 

# looking at traceplots I'm not clear that we have convergence.
traceplot(is_upd_parll, varname="mu.psi") #also look at thetas, 

######
ocmod.mcmc<-as.mcmc(is_upd_parll)
gelman.diag(ocmod.mcmc) #this is taking an impressively long time maybe need more thinning? 
#if I recall correctly this is sort of an ANOVA to compare within-chain and between-chain noisiness.. and the cutoff is somewhere around 1.1... if greater is bad then
# within/between? I could go back and read up on this.
# this failed with much longer chains too. 

# eventually got: Error in chol.default(W) : 
# the leading minor of order 356 is not positive definite
#I think this means too many variables and too few observations, but fitting again. 

#I think this just runs another bunch of iterations, can come back to it. 
# update(ocmod, n.iter = nburn)
#thinning after the fact here... or wait, is coda.samples where the MCMC sampling actually happens?
# out <- coda.samples(ocmod
#                     , n.iter = n.iter
#                     , variable.names = sp.params
#                     , thin=thin)

# out.mcmc <- as.mcmc(out[[1]])

# plot(out)

#
quartz()
plot(is_upd_parll) #Rhats looking close to 1 here!
#Z is confidently at 0 or 1 for each spp (why are there two? look at model file)
# colored points overlapping for mu_psi and theta! that's goodl
dev.off()

H<-heidel.diag(ocmod.mcmc) #looks like it says "failed" a lot more than 5% of everything
print(H)
#this is helpful when going back b/c it informs how many iterations are actually informative, too. 

#with the 50 sites it looks like things pass!

#G<-gelman.diag(out)

Z<-ocmod.mcmc[,grepl("Z",colnames(ocmod.mcmc))]
# theta<-out.mcmc[,grepl("theta",colnames(out.mcmc))]
# theta.med<-array(apply(theta,2,mean),dim=c(nsites,nreps,nspp))
# plot(c(p.detect),c(theta.med[,1,]))
# cor(c(p.detect),c(theta.med[,1,]))
# 
mu.psi<-ocmod.mcmc[,grepl("*psi*",colnames(ocmod.mcmc))]
psi.med<-array(apply(mu.psi,2,median),dim=c(nsites,nreps,nspp))
plot(c(p.occur),c(psi.med[,1,]), xlim=c(0,1), ylim=c(0,1))
# 

gelman.plot(ocmod.mcmc) #all these gelman things have memory issues. Hmmm. Thin more?
######## not sure how to work with my data here!

ocmod.df<-as.data.frame(as.matrix(ocmod.mcmc))

plotPost(ocmod.df$`mu.psi[1,1]`) #this produced a pretty histogram

# Z1<-matrix(Z[1,],nrow=nsites,ncol=nspp)
# Zavg<-matrix(colMeans(Z),nrow=nsites,ncol=nspp)
# Zmed<-matrix(apply(Z,2,median),nrow=nsites,ncol=nspp)
# sum(abs(Zobs-Zmed))
# sum(abs(Xtrue-Zmed))

## mu.theta is more bimodally distributed than p.detect, since it's estimated from a binomial dist'n
## Z is very similar to Zobs; the median of Z is perfect

