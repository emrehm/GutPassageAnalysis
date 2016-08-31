setwd(getwd()) # replace the getwd() with whatever you want

###
### Packages
###

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


# Load libraries

packages <- c("rjags","boot")
ipak(packages)



###
### Simulate data
###

### Define some parameters

# Sample sizes
n.samples <- 1000
mean.seeds <- 15
n.tree.id <- 20

# "True" survival parameters (that we hope to recapitulate later)
p.surv.fruit <- 0.4
p.surv.gut <- 0.8

# Simulate "true" fruit data
set.seed(1)
true.seeds <- rpois(n.samples, mean.seeds)
hist(true.seeds)

# Simulate "observed" fruit data
tree.id <- rep(1:n.tree.id,each=n.samples/n.tree.id) # This is tree tree.id
p.surv.tree.id <- round(rnorm(n.tree.id,0,0.02),digits=2)
#obs.surv <- rbinom(n = n.samples, size = true.seeds, prob = p.surv.fruit) # This one does not have tree.idual level random effects
obs.surv <- rbinom(n = n.samples, size = true.seeds, prob = p.surv.fruit+rep(p.surv.tree.id,each=n.samples/n.tree.id)) # To make tree.iduals have different survival
est.mean <- mean(true.seeds) # Lets just say that we're getting the mean from the same population
est.mean <- rep(est.mean,length(obs.surv))

# Simulate observed gut-passed data (Nseeds is no longer unknown)
tree.id <- rep(tree.id, 2)
#obs.surv <- c(obs.surv,rbinom(n = n.samples, size = true.seeds, prob = p.surv.gut)) # This one does not have tree.idual level random effects
obs.surv <- c(obs.surv,rbinom(n = n.samples, size = 1, prob = p.surv.gut+rep(p.surv.tree.id,each=n.samples/n.tree.id))) # To make tree.iduals have different survival

whole.fruit <- rep(c(1,0),each=n.samples)
gut.passed <- c(rep(NA, n.samples), rep(0:1, length.out=n.samples))
n.bird <- 10
bird.id <- as.factor(c(rep(NA, n.samples), 
                       rep(letters[1:n.bird], 
                           each=n.bird, 
                           length.out=n.samples)))
# Make half of these not gut passed
bird.id[n.samples + which((n.samples+1):(2*n.samples) %% 2 == T)] <- NA


# Figure out indices for each of the three treatments
wf.indices <- 1:n.samples
mech.indices <- n.samples + which((n.samples+1):(2*n.samples) %% 2 == T)
bird.indices <- n.samples + which((n.samples+1):(2*n.samples) %% 2 == T) + 1
seed.indices <- c(mech.indices,bird.indices)

# Make sure the data simulation makes sense
par(mfrow=c(1,1))
plot(true.seeds, obs.surv[1:n.samples],
     xlim=c(0,range(true.seeds)[2]),
     ylim=c(0,range(true.seeds)[2]),
     main="fruit")
# If there was 100% survival
curve(x*1,add=T)
# This is the mean relationship based on the "true" fruit survival parameter
curve(x*p.surv.fruit,add=T) 


#


###
### Define Model
###

# b0mean - mean survival for fruit
# b0 - the tree individual random effects (contained in b0[tree.id])
# b1 - the effect of gut passage
# b2 - the bird invidual random effects

sink("model.txt")
cat("

model { 
  
    #Likelihood
    
    # First loop over the whole fruits (no bird, yes tree random effects)
    for(i in wf.indices){ 
    
    # Likelihood
    # Process model
    obs.surv[i] ~ dbinom(p[i], n.seed.true[i])
    logit(p[i]) <- b0mean[i] # Deterministic model. Will add other terms.
    b0mean[i] ~ dnorm(b0[tree.id[i]], b0sig)
    
    # Data model
    
    n.seed.true[i] ~ dpois(l.seed) # Probably could just use dnorm
    
    }
    
    
    # Second loop over the mech seeds (no bird, yes tree random effects)
    for(i in mech.indices){ 
    
    # Likelihood
    # Process model
    obs.surv[i] ~ dbinom(p[i], 1) # Could make this multinomial...
    logit(p[i]) <- b0mean[i] + b1[1+gut.passed[i]] # Deterministic model. Will add other terms.
    
    b0mean[i] ~ dnorm(b0[tree.id[i]], b0sig)
    }


    # Third loop over the bird seeds (yes bird, yes tree random effects)
    # Not totally sure if this is best / most efficient. Could be included in the
    # previous forloop or a nested forloop potentially...
    for(i in bird.indices){ 
    
    # Likelihood
    # Process model
    obs.surv[i] ~ dbinom(p[i], 1) # Could make this multinomial...
    logit(p[i]) <- b0mean[i] + b1[1+gut.passed[i]] + b2mean[i] # Deterministic model. Will add other terms.
    
    b0mean[i] ~ dnorm(b0[tree.id[i]], b0sig)
    b2mean[i] ~ dnorm(b2[bird.id[i]], b2sig)
    }
    
    
    
    # Priors
    for(j in 1:max(tree.id)){
    b0[j] ~ dnorm(0,1/2.25)
    }
    b0sig ~ dunif(0,100)
    
    b1[1] ~ dnorm(0,1/2.25)
    b1[2] ~ dnorm(0,1/2.25)
    
    for(j in 1:n.bird){
    b2[j] ~ dnorm(0,1/2.25)
    }
    
    b2sig ~ dunif(0,100)
    
    l.seed ~ dnorm(15,3) # This is where mean seeds per fruit and sd are supplied
    
    
    
    
    }#end of model




    ",fill=TRUE)
sink()



###
### Set up and run model
###

# Make data object for jags.model to use
data <- list(obs.surv=obs.surv, 
             #est.mean=est.mean, # Currently just putting this directly in model...
             gut.passed=gut.passed,
             bird.id=bird.id,
             n.bird=n.bird,
             #n.samples=n.samples,
             wf.indices=wf.indices,
             mech.indices=mech.indices,
             bird.indices=bird.indices,
             #whole.fruit = whole.fruit,
             tree.id = tree.id)

# Make inits object for jags.model to use
# Notice that you have to have high estimates of n.seed.true that's too big (can't be lower than n.obs)
inits <- list(
  list(n.seed.true=rep(30,n.samples),b0=rnorm(n.tree.id,0,2),b1=rnorm(2),l.seed=rpois(n=1,lambda=30),b2=rnorm(n.bird)),#,p=runif(n.tree.id,0,1)),
  list(n.seed.true=rep(30,n.samples),b0=rnorm(n.tree.id,0,2),b1=rnorm(2),l.seed=rpois(n=1,lambda=30),b2=rnorm(n.bird)))#,p=runif(n.tree.id,0,1))


parameters <- c("p","b0") # use this for variable list in coda object


###Compile JAGS model, do updates, obtain coda obect for parameters, and JAGS object for plotting posterior distribution of optimum elevation, Bayesian p value, and probability of occupancy as a function of elevation

M1.model <- jags.model("model.txt",data=data,inits=inits,n.chains=2,n.adapt=1000)

update(M1.model, n.iter=1000)

zj <- jags.samples(M1.model,data=data, inits=inits, n.chains=2,variable.names=c("b0","b1", "b2", "n.seed.true", "p"), n.iter=3000)





###
### Look at model results, see if it recapitulates "true" parameters
###

# Does it provide good estimates of the 'true' seeds?
par(mfrow=c(1,1))
plot(true.seeds ~ summary(zj$n.seed.true,FUN=mean)$stat)
curve(x*1,add=T) # Collapses the distribution slightly...

# There were no true differences among "mech" and "gut" so we expect
# these to be similar.
inv.logit(summary(zj$b1, FUN=mean)$stat)
# They are.



# See if the tree.id random effects terms (y-axis) recapitulate the true parameters (x-axis)
par(mfrow=c(2,1))
plot(inv.logit(summary(zj$b0, FUN=mean)$stat) ~ c(p.surv.tree.id+rep(p.surv.fruit,n.tree.id)),
     xlab="true random effect",
     ylab="model random effect") # tapply(obs.surv,tree.id,mean)
curve(x*1,add=T)
# A different perspective is whether they recaptiulate not the true mean
# seed survival but the observed mean seed survival
plot(inv.logit(summary(zj$b0, FUN=mean)$stat) ~ c(tapply(obs.surv[1:n.samples],tree.id[1:n.samples],mean)/est.mean[1]),
     xlab="observed random effect",
     ylab="model random effect")
curve(x*1,add=T)

mean(summary(zj$p, FUN=mean)$stat) # Overall mean survival - should be mean of 0.8 and 0.4
hist(summary(zj$b2, FUN=mean)$stat) # This is the distribution of bird random effects

inv.logit(mean(summary(zj$b0, FUN=mean)$stat))
inv.logit(mean(summary(zj$b0, FUN=mean)$stat) + mean(summary(zj$b1, FUN=mean)$stat))
inv.logit(mean(summary(zj$b0, FUN=mean)$stat) + mean(summary(zj$b1, FUN=mean)$stat + summary(zj$b2, FUN=mean)$stat))
inv.logit(mean(summary(zj$b0, FUN=mean)$stat))

