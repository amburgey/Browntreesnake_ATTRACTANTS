### Simulating data to demonstrate how SCR model works for differential encounter probabilities

rm(list = ls())

#### LATENT information (i.e., we don't know this in real life but instead estimate it - however, we can simulate it here to demonstrate our model) ####
# Simulate the true population
# True number of snakes
N <- 109
# Establish parameters for detection, lam0[1] = no trap lure, lam0[2] = trap lure, lam0[3] = trap inactive
sigma <- 32
lam0 <- c(0.003, 0.004, 0)

#### STUDY AREA information ####
# Sampling grid
locs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
# Trap locations (5 transects with lures, 8 transects without, sensu lure project)
pts <- as.matrix(locs)
J <- nrow(pts)
# Define state-space of point process. (i.e., where animals live)
delta<- 11.874929
Xl<-min(locs[,1]) - delta
Xu<-max(locs[,1]) + delta
Yl<-min(locs[,2]) - delta
Yu<-max(locs[,2]) + delta
# Study area
A <- (Xu-Xl)*(Yu-Yl)
# Simulate activity centers of all snakes
sx<-runif(N,Xl,Xu)
sy<-runif(N,Yl,Yu)
S<-cbind(sx,sy) 

#### SAMPLING information ####
# Number of nights trapping
K <- 30
## Matrix of transect points (rows) by surveys (columns),
## Whether a transect point was active (surveyed) is indicated
## by a 1 (active) or 0 (inactive)
act <- matrix(data = sample(c(0,1), size=(J*K), replace=TRUE), nrow = J, ncol = K)
## Matrix of transect points (rows) by surveys (columns) indicating 
## status of surveys (e.g., 1 = not active, 2 = active but no lure, 3 = 
## active and lure)
stat <- matrix(NA, nrow = nrow(act), ncol = ncol(act))
for(i in 1:nrow(act)){
  for(j in 1:ncol(act)){
    stat[i,j] <- ifelse(act[i,j] == 1, 1, sample(2:3, 1))
  }
} 	

#### OBSERVATION information ####
# Function to calculate Euclidean distance between two sets of (x,y) locations
e2dist <- function(A, B)  {
  xdif <- outer(A[, 1], B[, 1], "-")
  ydif <- outer(A[, 2], B[, 2], "-")
  sqrt(xdif^2 + ydif^2)
}
# Distance between each individual AC and each each trap
d2 <- e2dist(S,pts)
# Simulate the encounters of all individuals in traps across the study area
p <- array(NA, dim=c(N,J,K))
y <- array(NA, dim=c(N,J,K))

for(i in 1:N){                  ## for every individual
  for(j in 1:J){                ## for every trap
    for(k in 1:K){              ## for every occasion
      p[i,j,k] <- lam0[stat[j,k]]*exp(-(d2[i,j])/(2*sigma*sigma))
      y[i,j,k] <- rbinom(1,1,p[i,j,k])
    }
  }
}

## NOTE: In this situation, the simulated data has no individual that was never caught (i.e., rowSums[y] > 0 for all individuals). When simulating data, you'd normally need to remove these animals that were never encountered to better represent the data collected.

#### DATA ANALYSIS ####
## Number of observed individuals
nind <- nrow(y)
## Data augmentation (i.e., add a number of uncaptured animals due to unknown abundance)
M <- 250	
y <- abind(y,array(0, dim = c((M - nind), J, K)), along = 1)

## Initial values for snake activity centers
set.seed(552020)
## Generate random state space locations
sst <- cbind(runif(M, Xl, Xu), runif(M, Yl, Yu))
## Find real locations of marked snakes
sumy <- rowSums(y, dims=2)
## Set initial values to mean location of marked individuals
## Set augmented individuals to random location within state space
for (i in 1:nind){
  for (k in 1:K){
    sst[i,1] <- mean(pts[sumy[i,]>0, 1])
    sst[i,2] <- mean(pts[sumy[i,]>0, 2]) 	
  }
}

## Trick to skip estimating info for inactive transect points
nActive <- apply(act, 1, sum)
ActiveOcc <- matrix(NA, J, max(nActive))
for(j in 1:J){
  ActiveOcc[j,1:nActive[j]] <- which(act[j,]==1)
}


