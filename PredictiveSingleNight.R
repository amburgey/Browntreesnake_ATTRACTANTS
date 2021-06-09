### Integrating across space to understand how encounter probability on a given night scales to captures across a week

rm(list = ls())

##### LURE PROJECT #####

load("SCRVISlurenolure.RData")

## Establish parameters for detection, lam0[1] = no trap lure, lam0[2] = trap lure
sigma <- out$mean$sigma
lam0 <- c(out$mean$lam0[1], out$mean$lam0[2])

## Grid
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

## Simulate the encounter probability and observations of this individual in every trap and repeat 1000 times
p <- array(NA,dim = c(1,J,length(STATUS)))
y <- array(NA,dim = c(1000,J,length(STATUS)))

set.seed(04042021)

for(l in 1:length(STATUS)){   ## do for a situation where every cell is no lure, lure
  for(s in 1:1000){    ## do 1000 simulations
    for(i in 1:1){
      for(j in 1:J){
        p[i,j,l] <- lam0[STATUS[l]]*exp(-(d2[i,j])/(2*sigma*sigma))
        y[s,j,l] <- rbinom(1,K,p[i,j,l])
      }
    }
  }
}

## Mean probability of encountering a snake on a grid with no lures
NTLp <- mean(p[1,,1])

## Probability of detecting a snake across the area during an evening with no lures
NTLpstar <- 1-prod(1-p[1,,1])

## Mean probability of encountering a snake on a grid with lures
TLp <- mean(p[1,,2])

## Probability of detecting a snake across the area during an evening with lures
TLpstar <- 1-prod(1-p[1,,2])



##### SCENT PROJECT #####

rm(list = ls())

load("SCRVISscentnoscentGroupOwnCat.RData")

## Establish parameters for detection, lam0[1] = no scent, lam0[2] = fresh scent, lam0[3] = old scent
sigma <- out$mean$sigma
lam0 <- c(out$mean$lam0[1], out$mean$lam0[2], out$mean$lam0[3])

## Grid
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

## Simulate the encounter probability and observations of this individual in every trap and repeat 1000 times
p <- array(NA,dim = c(1,J,length(STATUS)))
y <- array(NA,dim = c(1000,J,length(STATUS)))

set.seed(04052021)

for(l in 1:length(STATUS)){   ## do for a situation where every cell is no scent, etc
  for(s in 1:1000){    ## do 100 simulations
    for(i in 1:1){
      for(j in 1:J){
        for(k in 1:K){
          p[i,j,l] <- lam0[STATUS[l]]*exp(-(d2[i,j])/(2*sigma*sigma))
          y[s,j,l] <- rbinom(1,K,p[i,j,l])
        }
      }
    }
  }
}

## Mean probability of capture on a grid with no scent
NSp <- mean(p[1,,1])

## Probability of capturing a snake on grid during an evening with no scent
NSpstar <- 1-prod(1-p[1,,1])

## Mean probability of capture on a grid with fresh scent
FSp <- mean(p[1,,2])

## Probability of capturing a snake on grid during an evening with fresh scent
FSpstar <- 1-prod(1-p[1,,2])

## Mean probability of capture on a grid with old scent
OSp <- mean(p[1,,3])

## Probability of capturing a snake on grid during an evening with old scent
OSpstar <- 1-prod(1-p[1,,3])
