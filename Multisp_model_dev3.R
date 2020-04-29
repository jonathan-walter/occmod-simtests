Jarzyna_offset<-function() {
    #Prior distributions on the community level occupancy
    #and detection covariates
    beta.mean ~ dunif(0.001,0.999) #vague prior for the hyperparameter of the community-level occupancy covariates
    # does ignoring the boundaries matter here?
    
    a <- log(beta.mean) - log(1-beta.mean) # logit transformation
    
    u.mean ~ dunif(0.001,0.999) #
    
    q <- log(u.mean) - log(1-u.mean) # logit transformation
    
    p.mean ~ dunif(0.001,0.999) #vague prior for the hyperparameter of the community-level detection covariates
    #theta.mean is the average detection probability?
    
    b <- log(p.mean) - log(1-p.mean)
    
    p.site ~ dunif(0.001,0.999) 
    s <- log(p.site) - log(1-p.site)
    
    # mu.alpha1 ~ dnorm(0, 0.01) #site-level occupancy average
    #need to figure this one out still
    
    #this is the standard deviation (or precision, find out) for the occupancy distribution
    tau.beta ~ dgamma(10,1) #variability in species sensitivity to elevation 
    
    tau.u ~ dgamma(10,1)  # variability in species intercept occupancy (occupancy at average elevation)
    
    tau.p1 ~ dgamma(10,1) #detection offset variability between sites 
    tau.p2 ~ dgamma(10,1) #detection variability between species 
    
    sigma.beta <- 1/sqrt(tau.beta)
    sigma.u <- 1/sqrt(tau.u)
    sigma.p2 <- 1/sqrt(tau.p2)
    sigma.p1 <- 1/sqrt(tau.p1)
    
    for (i in 1:nspp) {
    #Prior distributions for the occupancy and detection covariates for each species 
        beta[i] ~ dnorm(a, tau.beta) # beta[i] is the species sensitivity to elevation
       
        u[i] ~ dnorm(q, tau.u) # u[i] is the species occupancy "intercept" (grand mean across all sites?)
        
        p.sp[i] ~ dnorm(b, tau.p1) #v[i]  (Species-level detection probability)
            # is simply a random deviate from normal with mean mu.v[i], sd var.v, which is given by rho and tau2. 
        
        
    #Estimate the occupancy probability (latent Z matrix) for each species #at each point (i.e., route or site)
        for (j in 1:nsite) {
            logit(psi[j,i]) <- u[i] + beta[i]*elev[j] #occupancy probability (psi[j,i]) 
            # is given by expected species occupancy (u[i]) and species-specific elevation response (beta[i]*elev[j]) #elev[j] comes from data
            # the logit transformation is to allow it to be "linearly" related to elevation
            # mu.psi[j,i] <- psi[j,i] #i don't yet know why this occurs
            Z[j,i] ~ dbin(psi[j,i], 1)#Z is generally not observed with certainty, instead
            p.st[j] ~ dnorm(s, tau.p2)
            #latent Z is the result of 1 bernouli trial from each psi[j,i]
            
            #Estimate the species-specific detection probability, and then the chance that it is observed 
            # at every rep at each point where the species occurs (Z=1)
                for (k in 1:nrep[j]) { 
                logit(theta[j,k,i]) <- p.st[j]+p.sp[j]  # detection theta[i,j,k] given that species species i is poresent at site j during sampling period k 
                mu.theta[j,k,i] <- theta[j,k,i]*Z[j,i] #multiply by indicator (1 if occupied)
                X[j,k,i] ~ dbin(mu.theta[j,k,i], 1)  # 
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
    

