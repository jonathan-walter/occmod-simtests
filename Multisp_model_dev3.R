jaw_model<-function() {
    #Prior distributions on the community level occupancy
    #and detection covariates
    psi.mean ~ dunif(0.001,0.99) #vague prior for the hyperparameter of the community-level occupancy covariates
    # does ignoring the boundaries matter here?
    
    a <- log(psi.mean) - log(1-psi.mean) # logit transformation
    
    theta.mean ~ dunif(0.001,0.99) #vague prior for the hyperparameter of the community-level detection covariates
    #theta.mean is the average detection probability?
    
    b <- log(theta.mean) - log(1-theta.mean)
    
    mu.alpha1 ~ dnorm(0, 0.01) #site-level occupancy average
    #need to figure this one out still
    
    tau1 ~ dgamma(10,1) 
    #this is the standard deviation (or precision, find out) for the occupancy distribution
    
    tau2 ~ dgamma(10,1)
    #
    
    tau.alpha1 ~ dgamma(10,1)#Zipkin's original priors #see if we can track down this code. 
    
    rho ~ dunif(-0.99,0.99) # rho is the correlation coefficient... between abundance and occupancy. Seems like counts could be used to verify this. 
    # what is this!
    var.v <- tau2 / (1-(rho^2))
     # ahhh!
    sigma1 <- 1/sqrt(tau1) 
    sigma2 <- 1/sqrt(tau2)
    
    for (i in 1:nspp) {
    #Prior distributions for the occupancy and detection covariates for each species 
        u[i] ~ dnorm(a, tau1) # u[i] is the species occupancy "intercept" (grand mean across all sites?)
         # "a" is the logit tranformation of psi.mean. 
        # tau1 is the variability in occupancy between species, as this is indexed to sp i
        
        mu.v[i] <- b + (rho*sigma2 /sigma1)*(u[i]-a) #subtracting a from u[i] is just centering that variable at 0 
            # This is a complicated bit having to do with the priors, but I think it ultimately has to do with detection? 
            # But why is detection related to u[i] which is occupancy? 
            # what's going on with the tau2 stuff (sigma2, rho) ?
            # why not simply dnorm(b, tau2)? 
        
        #BECAUSE THIS IS A MODEL THAT IGNORES ABUNDANCE AND ALSO IGNORES COUNTS. so Zipkin et al. 2009 J App Ec 
        # induced a correlation between detection and occurence, where the link is basically abundance (occupancy/site abundance correlation)
        #. Because high abundance species are likely to be both easier to detect and more prevalent across the landscape, 
        # we modelled a correlation ðqÞ between occurrence and detection in the model by allowing ui and vi to be jointly distributed such that 
        #... (Dorazio & Royle 2005; Kerry & Royle 2008).
        
        v[i] ~ dnorm(mu.v[i], var.v) #v[i]  (Species-level detection probability)
            # is simply a random deviate from normal with mean mu.v[i], sd var.v, which is given by rho and tau2. 
        
        alpha1[i] ~ dnorm(mu.alpha1, tau.alpha1) # this is beta! confusing that alpha is beta!
    
    #Estimate the occupancy probability (latent Z matrix) for each species #at each point (i.e., route or site)
        for (j in 1:nsite) {
            logit(psi[j,i]) <- u[i] + alpha1[i]*elev[j] #occupancy probability (psi[j,i]) 
            # is given by expected speces occupancy (u[i]) and species-specific elevation response (alpha1[i]*elev[j])
            # the logit transformation is to allow it to be "linearly" related to elevation
            mu.psi[j,i] <- psi[j,i] #i don't yet know why this occurs
            Z[j,i] ~ dbin(psi[j,i], 1)#Z is generally not observed with certainty, instead
           
            #latent Z is the result of 1 bernouli trial from each psi[j,i]
            
            #Estimate the species-specific detection probability, and then the chance that it is observed 
            # at every rep at each point where the species occurs (Z=1)
                for (k in 1:nrep[j]) { 
                logit(theta[j,k,i]) <- v[i]  #we observed data theta[i,j,k] for species i at site j during sampling period k 
                #I guess the idea here is that theta has to be this 3-d array, even though its values are all given by a species-level detection (v[i])
                mu.theta[j,k,i] <- theta[j,k,i]*Z[j,i] #multiply by indicator (1 if occupied)
                X[j,k,i] ~ dbin(mu.theta[j,k,i], 1)  # so mu.theta is the detection parameter
                #X is the 3D array of dependent variable: The detection/non-
                #detection data is defined in a three dimensional
                #array X where the first dimension, j, is the site; the second #dimension, k, is the rep; and the last dimension, 
                # i, is the species.
                Xnew[j,k,i] ~ dbin(mu.theta[j,k,i], 1) 
                #Create second simulated dataset to calculate the Bayesian p-value 
                d[j,k,i] <- abs(X[j,k,i] - mu.theta[j,k,i]) #difference between observation and detection*occupancy (expected)
                dnew[j,k,i] <- abs(Xnew[j,k,i] - mu.theta[j,k,i]) # difference between observation in alternate universe and expected
                d2[j,k,i] <- pow(d[j,k,i],2) #square differences
                dnew2[j,k,i] <- pow(dnew[j,k,i],2)  
                }
            dsum[j,i] <- sum(d2[j,1:nrep[j],i])
            dnewsum[j,i] <- sum(dnew2[j,1:nrep[j],i]) 
            }
        }
    #Calculate the discrepancy measure, which is then defined as the mean(p.fit > p.fitnew)
    p.fit <- sum(dsum[1:nsite,1:nspp]) 
    p.fitnew <- sum(dnewsum[1:nsite,1:nspp]) #this is driven by nrep(j) and the detection probablities... what does it tell us?
    #} } }
    #Estimation of species occupancy (averaged across all the sites)
    for(i in 1:nspp) {
    occ_sp[i] <- sum(Z[1:nsite,i])/nsite 
    }
    #End model specification
    }
    

