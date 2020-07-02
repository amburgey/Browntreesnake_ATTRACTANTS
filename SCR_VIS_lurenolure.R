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
X <- as.matrix(locs)
nolurepts <- X
lurepts <- X[c(14:26,40:52,66:78,92:104,118:130,144:156,170:182,196:208,222:234,248:260,274:286,300:312,326:338),]

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
ynolure <- as.matrix(dat$snksNTL[,-352])
colnames(ynolure) <- NULL
ylure <- as.matrix(dat$snksTL[,-170])
colnames(ylure) <- NULL

## Number of locations (VIS lure/no lure points)
JNL <- nrow(dat$actNTL)
JL <- nrow(dat$actTL)

# Active/not active for when transects run (one less day surveyed for TL)
actNTL <- dat$actNTL[,-c(1,32)]
colnames(actNTL) <- NULL
actTL <- dat$actTL[,-c(1,31)]
colnames(actTL) <- NULL

# Uniquely marked individuals
nindNL <- nrow(ynolure)
nindL <- nrow(ylure)

## Date augmentation
M <- 250
ynolure <-rbind(ynolure,array(0,dim=c((M-nrow(ynolure)),ncol(ynolure))))
ylure <-rbind(ylure,array(0,dim=c((M-nrow(ylure)),ncol(ylure))))

## Effort
KNL <- rowSums(actNTL)
KL <- rowSums(actTL)

## Starting values for activity centers
## set inits of AC to mean location of individuals and then to some area within stat space for augmented
set.seed(552020)
sst <- cbind(runif(M,Xl,Xu),runif(M,Yl,Yu))
for(i in 1:nindNL){
  sst[i,1] <- mean( X[ynolure[i,]>0,1] )
  sst[i,2] <- mean( X[ynolure[i,]>0,2] )
}


## JAGS model
cat(file="SCR0_DataAug.txt","
model {
  lam1 ~ dunif(0,5)  # Lure
  lam2 ~ dunif(0,5)  # No lure
  sigma ~ dunif(0,50)
  psi ~ dunif(0,1)

  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(Xl,Xu)
    s[i,2] ~ dunif(Yl,Yu)
    
    ## NO LURE SURVEYS ##
    for(j in 1:JNL) {
      dNL[i,j]<- pow(s[i,1]-nolurepts[j,1],2) + pow(s[i,2]-nolurepts[j,2],2)
      pNL[i,j]<- z[i]*lam2*exp(-(dNL[i,j])/(2*sigma*sigma))
      ynolure[i,j] ~ dpois(pNL[i,j]*KNL)
    }
    
    ## LURE SURVEYS ##
    for(k in 1:JL){
      dL[i,k]<- pow(pow(s[i,1]-lurepts[k,1],2) + pow(s[i,2]-lurepts[k,2],2),0.5)
      pL[i,k]<- z[i]*lam1*exp(-(dL[i,k]*dL[i,k])/(2*sigma*sigma))
      ylure[i,k] ~ dpois(pL[i,k]*KL)   
    }
  }
  
  N <- sum(z[])
  D <- N/A
  
}")

# MCMC settings
nc <- 3; nAdapt=1000; nb <- 1; ni <- 2000+nb; nt <- 1

# data and constants
jags.data <- list (ynolure=ynolure, ylure=ylure, X=X, M=M, JNL=JNL, JL=JL, Xl=Xl, Xu=Xu, Yl=Yl, Yu=Yu, A=A, KNL=KNL, KL=KL)

inits <- function(){
  list (sigma=runif(1,1,50), z=c(rep(1,nindNL),rep(0,M-nindNL)), s=sst, psi=runif(1), lam1=runif(1,.5,1.5), lam1=runif(1,.5,1.5))
}

parameters <- c("sigma","psi","N","D","lam1","lam2","pNL","pL")

out <- jagsUI("SCR0_DataAug.txt", data=jags.data, inits=inits, parallel=TRUE,
            n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters)

save(out, file="SCRVISlurenolure.RData")

print(out)
traceplot(out$samples[,"lam1"])

# gamma parameters for sigma
a.sigma = (out$mean$sigma)^2/(out$sd$sigma)^2
b.sigma = out$mean$sigma/(out$sd$sigma)^2
N = out$mean$N

plot(rgamma(N, a.sigma, b.sigma), type = "h", xlim=c(0,130))

