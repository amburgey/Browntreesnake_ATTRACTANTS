#### SCR analysis for visual survey - scent trail vs. no scent trail project
## Data was collected in 2016 where biologists walked transects in CP that either had been sprayed with fish fertilizer (scent), been sprayed previously (old scent), or had no fertilizer (no scent)
## The goal is to understand if encounter probability is different between these three scenarios
## In the case of an incipient or suppressed (low density) population, any way that we can maximize encounter rates (and thus detection probability) of brown treesnakes can be a useful resource for surveyors


rm(list = ls())

# Prep data for analysis
source("Data/PrepDataScent.R")

# Read in capture data
caps <- read_csv("Data/CapturesScent.csv")[,c("EFFORTID","Date","PITTAG","SVL","TOTAL","WEIGHT","SEX","TRANSECT","LOCATION")]
# Restrict to 2-month period to meet assumptions of population closure
caps <- caps %>%
  filter(!grepl('Oct', Date)) %>%
  filter(!grepl('Jan', Date))

# Read in survey data
survs <- read_csv("Data/SurveysScent.csv")
survs <- survs %>%
  select(-contains("Oct")) %>%
  select(-contains("Jan"))

# Format for analysis
dat <- PrepDat(caps,survs)

# Create trapping grid of Closed Population (CP) dimensions (5 ha, 50,000 m2)
locs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
pts <- as.matrix(locs)

## Number of locations (no scent/fresh scent/old scent points)
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

# Active/not active for when transects run
act <- as.matrix(dat$act)

# Status of surveys (1 = not active, 2 = active but no scent, 3 = active and fresh scent, 4 = active and old scent)
scent <- as.matrix(dat$scent)

# Number of survey occasions
nocc <- ncol(act)

## Data augmentation
M <- 250
y <-abind(y,array(0,dim=c((M-nrow(y)),ncol(y),nocc)), along = 1)

## Starting values for activity centers
## set inits of AC to mean location of individuals and then to some area within stat space for augmented
set.seed(2252021)
sst <- cbind(runif(M,Xl,Xu),runif(M,Yl,Yu))
sumy <- rowSums(y, dims=2)
for (i in 1:nind){
  for(k in 1:nocc){
    sst[i,1] <- mean(pts[sumy[i,]>0,1])
    sst[i,2] <- mean(pts[sumy[i,]>0,2])
  }
}

## Trick to skip estimating info for inactive transect points
nActive <- apply(act, 1, sum)
ActiveOcc <- matrix(NA, J, max(nActive))
for(j in 1:J){
  ActiveOcc[j,1:nActive[j]] <- which(act[j,]==1)
}

## JAGS model
cat(file="Models/SCR0_DataAugSCENT.txt","
model {

  for(l in 1:3){
    lam0[l]~dunif(0,1) ## Detection model with 1, 2, 3 indicator for treatment
  }
  sigma ~ dunif(0,50)
  psi ~ dunif(0,1)

  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(Xl,Xu)
    s[i,2] ~ dunif(Yl,Yu)
  
    for(j in 1:J) {
      d2[i,j]<- pow(s[i,1]-pts[j,1],2) + pow(s[i,2]-pts[j,2],2)
        
      for(k in 1:nActive[j]){
        p[i,j,ActiveOcc[j,k]] <- z[i]*lam0[STATUS[j,ActiveOcc[j,k]]]*exp(-(d2[i,j])/(2*sigma*sigma))
        y[i,j,ActiveOcc[j,k]] ~ dbern(p[i,j,ActiveOcc[j,k]])
      }
    }
  }
  
  N <- sum(z[])                 ## Abundance
  D <- N/A                      ## Density
  delta1 <- lam0[1] - lam0[2]   ## Difference between no scent and fresh scent enc. rate
  delta2 <- lam0[1] - lam0[3]   ## Difference between no scent and old scent enc. rate
  delta3 <- lam0[3] - lam0[2]   ## Difference between old scent and fresh scent enc. rate
  
}")

# MCMC settings
nc <- 3; nAdapt=1000; nb <- 1; ni <- 2000+nb; nt <- 1

# data and constants
jags.data <- list (y=y, pts=pts, M=M, J=J, Xl=Xl, Xu=Xu, Yl=Yl, Yu=Yu, A=A, STATUS=scent-1, nActive=nActive, ActiveOcc=ActiveOcc)
# scent-1 since the only functional categories are 1, 2 and 3... inactive state becomes 0 this way.
inits <- function(){
  list (sigma=runif(1,40,50), z=c(rep(1,nind),rep(0,M-nind)), s=sst, psi=runif(1), lam0=runif(3,0.05,0.07))
}

parameters <- c("sigma","psi","N","D","lam0","delta1","delta2","delta3")

out <- jagsUI("/Models/SCR0_DataAugSCENT.txt", data=jags.data, inits=inits, parallel=TRUE,
            n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters)

save(out, file="/Results/SCRVISscentnoscent.RData")

