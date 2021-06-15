### Integrating across space to understand how encounter probability scales to probability of detecton on a given night
## Simulate a single snake with an activity center in the center of the study area
## Each treatment is applied uniformly on a given evening (i.e., ever grid cell receives no lure vs. each receives lure)
## Parameters that affect detection (sigma and lam0) are from model results from each study

rm(list = ls())

##### LURE PROJECT #####

load("Results/SCRVISlurenolure.RData")

## Establish parameters for detection, lam0[1] = no trap lure, lam0[2] = trap lure
sigma <- out$sims.list$sigma
lam0 <- out$sims.list$lam0

## Study area grid
locs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
pts <- as.matrix(locs)
J <- nrow(pts)

## Status of all grid cells that evening (1 = no lure, 2 = lure)
STATUS <- c(1,2)

## Define state-space of point process. (i.e., where animals live)
delta<- 11.874929
Xl<-min(locs[,1]) - delta
Xu<-max(locs[,1]) + delta
Yl<-min(locs[,2]) - delta
Yu<-max(locs[,2]) + delta
## Check area: 
A <- (Xu-Xl)*(Yu-Yl)

# Number of nights trapping (just one)
K <- 1

## Snake activity center is in center of grid
S <- as.data.frame(matrix(as.matrix(locs[176,]), nrow = 1, ncol = 2))

## Function to calculate distance between two sets of (x,y) locations
e2dist <- function(A, B)  {
  xdif <- outer(A[, 1], B[, 1], "-")
  ydif <- outer(A[, 2], B[, 2], "-")
  sqrt(xdif^2 + ydif^2)
}

#Distance between each individual AC and each trap
d2 <- e2dist(S,pts)

## Simulate the encounter probabilities of this individual in every trap and repeat for every posterior sample (i.e., 6000 times) in order to calculate mean and uncertainty
p <- array(NA,dim = c(length(sigma),J,length(STATUS)))

set.seed(04042021)

for(l in 1:length(STATUS)){   ## do for a situation where every cell is no lure, lure
  for(j in 1:J){  ## do for every survey location
    for(s in 1:length(sigma)){## do for every element of out$sims.list$sigma
      p[s,j,l] <- lam0[s,STATUS[l]]*exp(-(d2[1,j])/(2*sigma[s]*sigma[s]))
    }
  }
}


## Mean probability of encountering a snake on a grid with no lures
NTLp <- mean(p[,,1])

## Probability of detecting a snake across the area during an evening with no lures
NTLpstar <- vector()
for(s in 1:nrow(p)){
  NTLpstar[s] <- 1-prod(1-p[s,,1])
}
mean(NTLpstar)
quantile(NTLpstar, probs=c(0.025,0.975))

## Mean probability of encountering a snake on a grid with lures
TLp <- mean(p[,,2])

## Probability of detecting a snake across the area during an evening with lures
TLpstar <- vector()
for(s in 1:nrow(p)){
  TLpstar[s] <- 1-prod(1-p[s,,2])
}
mean(TLpstar)
quantile(TLpstar, probs=c(0.025,0.975))



##### SCENT PROJECT #####

rm(list = ls())

load("Results/SCRVISscentnoscent.RData")

## Establish parameters for detection, lam0[1] = no scent, lam0[2] = fresh scent, lam0[3] = old scent
sigma <- out$sims.list$sigma
lam0 <- out$sims.list$lam0

## Study area grid
locs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
pts <- as.matrix(locs)
J <- nrow(pts)

## Status of all grid cells that evening (1 = no scent, 2 = fresh scent, 3 = old scent) 
STATUS <- c(1,2,3)

## Define state-space of point process. (i.e., where animals live)
delta<- 11.874929
Xl<-min(locs[,1]) - delta
Xu<-max(locs[,1]) + delta
Yl<-min(locs[,2]) - delta
Yu<-max(locs[,2]) + delta
## Check area: 
A <- (Xu-Xl)*(Yu-Yl)

# Number of nights trapping (just one)
K <- 1

## Snake activity center is in center of grid
S <- as.data.frame(matrix(as.matrix(locs[176,]), nrow = 1, ncol = 2))

## Function to calculate distance between two sets of (x,y) locations
e2dist <- function(A, B)  {
  xdif <- outer(A[, 1], B[, 1], "-")
  ydif <- outer(A[, 2], B[, 2], "-")
  sqrt(xdif^2 + ydif^2)
}

#Distance between each individual AC and each each trap
d2 <- e2dist(S,pts)

## Simulate the encounter probabilities of this individual in every trap and repeat for every posterior sample (i.e., 6000 times) in order to calculate mean and uncertainty
p <- array(NA,dim = c(length(sigma),J,length(STATUS)))

set.seed(04042021)

for(l in 1:length(STATUS)){   ## do for a situation where every cell is no lure, lure
  for(j in 1:J){  ## do for every survey location
    for(s in 1:length(sigma)){## do for every element of out$sims.list$sigma
      p[s,j,l] <- lam0[s,STATUS[l]]*exp(-(d2[1,j])/(2*sigma[s]*sigma[s]))
    }
  }
}


## Mean probability of encountering a snake on a grid with no scent
NSp <- mean(p[,,1])

## Probability of detecting a snake across the area during an evening with no lures
NSpstar <- vector()
for(s in 1:nrow(p)){
  NSpstar[s] <- 1-prod(1-p[s,,1])
}
mean(NSpstar)
quantile(NSpstar, probs=c(0.025,0.975))

## Mean probability of encountering a snake on a grid with fresh scent
FSp <- mean(p[,,2])

## Probability of detecting a snake across the area during an evening with lures
FSpstar <- vector()
for(s in 1:nrow(p)){
  FSpstar[s] <- 1-prod(1-p[s,,2])
}
mean(FSpstar)
quantile(FSpstar, probs=c(0.025,0.975))

## Mean probability of encountering a snake on a grid with old scent
OSp <- mean(p[,,3])

## Probability of detecting a snake across the area during an evening with lures
OSpstar <- vector()
for(s in 1:nrow(p)){
  OSpstar[s] <- 1-prod(1-p[s,,3])
}
mean(OSpstar)
quantile(OSpstar, probs=c(0.025,0.975))

