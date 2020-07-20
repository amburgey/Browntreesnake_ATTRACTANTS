### SCR analysis for VIS lure vs no lure project
### Data was collected in 2015 where biologists walked transects in CP that either had traps with mice (lure) or no traps (no lure)
### The goal is to understand if detection probability is different between these two scenarios
### In the case of an incipient or suppressed (low density) population, any way that can maximize detection probability of brown treesnakes can be a useful resource for surveyors


rm(list = ls())


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

# nolurepts <- X
# lurepts <- X[c(14:26,40:52,66:78,92:104,118:130,144:156,170:182,196:208,222:234,248:260,274:286,300:312,326:338),]

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

# Status of surveys (1 = not active, 2 = active but no lure, 3 = active and lure)
stat <- as.matrix(dat$stat[,-1])

# Number of survey occasions
nocc <- ncol(act)

## Date augmentation
M <- 250
y <-abind(y,array(0,dim=c((M-nrow(y)),ncol(y),nocc)), along = 1)



## Effort
# K <- rowSums(act)

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

## JAGS model
cat(file="SCR0_DataAug.txt","
model {

  for(l in 1:3){
    lam0[l]~dunif(0,1) ## Detection model with 1, 2, 3 indicator
  }
  sigma ~ dunif(0,50)
  psi ~ dunif(0,1)

  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(Xl,Xu)
    s[i,2] ~ dunif(Yl,Yu)
  
      for(j in 1:J) {
        d[i,j,k]<- pow(s[i,1]-pts[j,1],2) + pow(s[i,2]-pts[j,2],2)
        p[i,j,k]<- z[i]*lam0[STATUS[j,k]]*exp(-(d[i,j]*d[i,j])/(2*sigma*sigma))
        y[i,j,k] ~ dbinom(p[i,j,k],act[j,k])
      }
    }
  }
  
  N <- sum(z[])
  D <- N/A
  
}")

# MCMC settings
nc <- 3; nAdapt=100; nb <- 1; ni <- 200+nb; nt <- 1

# data and constants
jags.data <- list (y=y, pts=pts, M=M, J=J, Xl=Xl, Xu=Xu, Yl=Yl, Yu=Yu, A=A, act=act, STATUS=stat, nocc=nocc)

inits <- function(){
  list (sigma=runif(1,1,50), z=c(rep(1,nind),rep(0,M-nind)), s=sst, psi=runif(1), lam0=runif(3,0,0.3))
}

parameters <- c("sigma","psi","N","D","lam0","p")

out <- jagsUI("SCR0_DataAug.txt", data=jags.data, inits=inits, parallel=TRUE,
            n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters)

save(out, file="SCRVISlurenolure.RData")


