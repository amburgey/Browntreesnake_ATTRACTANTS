### SCR analysis for VIS lure vs no lure project
### Data was collected in 2015 where biologists walked transects in CP that either had traps with mice (lure) or no traps (no lure)
### The goal is to understand if detection probability is different between these two scenarios
### In the case of an incipient or suppressed (low density) population, any way that can maximize detection probability of brown treesnakes can be a useful resource for surveyors


rm(list = ls())

library(nimble)

source("PrepData.R")


# Read in capture data and survey data
caps <- read.csv("Captures.csv")[,c("EFFORTID","Date","PITTAG","SVL","TOTAL","WEIGHT","SEX","TRANSECT","LOCATION")]
survs <- read.csv("Surveys.csv")[,c("EFFORTID","Date","TRANSECT","TYPE")]

# Format for analysis
dat <- PrepDat(caps,survs)

# Create trapping grid of CP dimensions (5 ha, 50,000 m2)
locs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
pts <- as.matrix(locs)

## Number of locations (VIS lure/no lure points)
J <- nrow(pts)

## Define state-space of point process. (i.e., where animals live).
## "delta" just adds a fixed buffer to the outer extent of the traps.
## Don't need to estimate state-space since we know it (5 ha enclosed pop)
delta<- 11.874929
Xl<-min(locs[,1]) - delta
Xu<-max(locs[,1]) + delta
Yl<-min(locs[,2]) - delta
Yu<-max(locs[,2]) + delta
# Check area: 
A <- (Xu-Xl)*(Yu-Yl)

# Observations
y <- dat$snks
colnames(y) <- NULL
rownames(y) <- NULL

# Uniquely marked individuals
nind <- nrow(y)

# Active/not active for when transects run (one less day surveyed for TL)
act <- as.matrix(dat$act[,-1])
colnames(act) <- NULL

nActive <- apply(act, 1, sum)  ## number of occasions a camera was active
ActiveOcc <- matrix(NA, J, max(nActive ))
for(j in 1:J){
  ActiveOcc[j,1:nActive[j]] <- which(act[j,]==1)  ## which occasions a camera was active
}

K <- ncol(ActiveOcc)

# Status of surveys (1 = not active, 2 = active but no lure, 3 = active and lure)
stat <- as.matrix(dat$stat[,-1])
colnames(stat) <- NULL

# Number of survey occasions
nocc <- ncol(act)

## Date augmentation
M <- 250
y <-abind(y,array(0,dim=c((M-nrow(y)),ncol(y),nocc)), along = 1)



## Starting values for activity centers
## set inits of AC to mean location of individuals and then to some area within stat space for augmented
set.seed(552020)
sst <- cbind(runif(M,Xl,Xu),runif(M,Yl,Yu))
sumy <- rowSums(y, dims=2)
for (i in 1:nind){
  for(k in 1:nocc){
    sst[i,1] <- mean(pts[sumy[i,]>0,1])
    sst[i,2] <- mean(pts[sumy[i,]>0,2])
  }
}


# MCMC settings
nc = 3; adaptInterval = 500; nb = 10000; ni = 20000+nb; nt = 5



##### NIMBLE model

NimModel <- nimbleCode({

  for(l in 1:2){
    lam0[l]~dunif(0,1) ## Detection model with 1, 2, 3 indicator (based on loop below should skip lamo=1 (inactive))
  }
  sigma ~ dunif(0,50)
  psi ~ dunif(0,1)

  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(Xl,Xu)
    s[i,2] ~ dunif(Yl,Yu)
  
    # d2[i,1:J] <- GetDistance(s=s[i,1:2], pts=pts[1:J,1:2], J=J)
    
    for(j in 1:J) {
      d2[i,j]<- pow(s[i,1]-pts[j,1],2) + pow(s[i,2]-pts[j,2],2)
      
      for(k in 1:nActive[j]){
        p[i,j,ActiveOcc[j,k]] <- z[i]*lam0[STATUS[j,ActiveOcc[j,k]]]*exp(-(d2[i,j])/(2*sigma*sigma))
        y[i,j,ActiveOcc[j,k]] ~ dbern(p[i,j,ActiveOcc[j,k]])
      }
    }
  }
  
  N <- sum(z[1:M])
  D <- N/A
  
})


#### NIMBLE Functions
# GetDistance <- nimbleFunction(
#   run = function(s=double(1), pts=double(2), J=double(0)){
#     returnType(double(1))
#     d2 <- pow(s[1,1]-pts[1:J,1],2) + pow(s[2,2]-pts[1:J,2],2)
#     return(d2)
#   }
# )

## Specifications and data to NIMBLE
constants <- list(pts=pts, M=M, J=J, K=K, Xl=Xl, Xu=Xu, Yl=Yl, Yu=Yu, A=A)
data <- list(y=y, ActiveOcc=ActiveOcc, STATUS=stat-1, nActive=nActive)
## stat-1 since the only functional categories are 2 and 3 and the nActive forloop will skip inactive traps.
inits <- list (sigma=40, z=c(rep(1,nind),rep(0,M-nind)), s=sst, psi=0.5, lam0=c(0.05,0.07))

parameters <- c("sigma","psi","N","D","lam0")

## Compile and run in NIMBLE
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=data, check=FALSE, inits=inits)
conf <- configureMCMC(Rmodel, monitors=parameters, control=list(adaptInterval=adaptInterval), thin=nthin)

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)

## Run multiple chains
out <- runMCMC(Cmcmc, niter=ni, nburnin=nb, nchains=nc, inits=inits, setSeed=FALSE, progressBar=TRUE, samplesAsCodaMCMC=TRUE)
end.time <- Sys.time()
end.time-start.time

## Save
save(out, file="SCRVISlurenolureNIMBLE.RData")

## Summary info
summary(out[,c('sigma','psi','N','D','lam0')])
gelman.diag(out[,c('sigma','psi','N','D','lam0')], multivariate=FALSE)

traceplot(out[,"N"])
traceplot(out[,"sigma"])
traceplot(out[,"lam0"])

## Effective samples per minute
effectiveSize(out[,c('sigma','psi','N','D','lam0')]) / as.numeric(end.time-start.time)


