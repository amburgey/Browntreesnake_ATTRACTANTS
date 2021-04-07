### Integrating across space to understand how encounter probability on a given night scales to captures across a week

rm(list = ls())

##### LURE PROJECT #####

load("SCRVISlurenolure.RData")

## Establish parameters for detection, lam0[1] = no trap lure, lam0[2] = trap lure
sigma <- out$mean$sigma
lam0 <- c(out$mean$lam0[1], out$mean$lam0[2])

## Grid
locs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
## Trap locations (5 transects with lures, 8 transects without, sensu lure project)
## BHNTZ traps, AEILPSWAA
# NTL <- locs[c(1:13,53:65,105:117,144:156,196:208,235:247,287:299,339:351),]
# TL <- locs[c(14:26,92:104,170:182,248:260,326:338),]
pts <- as.matrix(locs[c(1:26,53:65,92:117,144:156,170:182,196:208,235:260,287:299,326:351),])
J <- nrow(pts)

# ## When transects are active and inactive (for more than a week)
# act <- matrix(rep(rep(c(1,1,0,0,1,0,0,1,1,0,0,1,0,1,0,1,0,0,1,1,0,0,1,0,0,1,1), each=13),times=7),nrow = 351, ncol = 7)

## Status of transect, constant during the one week of sampling
stat <- rep(c(2,3,2,3,2,2,3,2,2,3,2,3,3), each=13)
STATUS <- stat - 1
# stat <- matrix(rep(rep(c(2,3,1,1,2,1,1,3,2,1,1,2,1,3,1,2,1,1,2,3,1,1,2,1,1,3,3), each=13),times=7),nrow = 351, ncol = 7)

## Define state-space of point process. (i.e., where animals live)
delta<- 11.874929
Xl<-min(locs[,1]) - delta
Xu<-max(locs[,1]) + delta
Yl<-min(locs[,2]) - delta
Yu<-max(locs[,2]) + delta
## Check area: 
A <- (Xu-Xl)*(Yu-Yl)

# Number of nights trapping (one week)
K <- 7

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

## Simulate the encounters of this individual in every trap and repeat 100 times
p <- matrix(NA, nrow=1000,ncol=J)
y <- matrix(NA,nrow=1000,ncol=J)

set.seed(04042021)

for(s in 1:1000){    ## do 1000 simulations
  for(i in 1:1){
    for(j in 1:J){
      p[i,j] <- lam0[STATUS[j]]*exp(-(d2[i,j])/(2*sigma*sigma))
      y[s,j] <- rbinom(1,K,p[i,j])
    }
  }
}

## BHNTZ traps, AEILPSWAA
## Observations on transects without lures
NTLcap <- mean(colSums(as.data.frame(y[,c(1:13,27:39,53:78,92:117,131:143,157:169)])))

## Observations on transects with mouse lures
TLcap <- mean(colSums(as.data.frame(y[,c(14:26,40:52,79:91,118:130,144:156)])))

## CPUE on transects without lures (8 transects of 0.22 km distance for 7 days)
NTLcap/(8*0.22*7)

## CPUE on transects with mouse lures (5 transects of 0.22 km distance for 7 days)
TLcap/(5*0.22*7)



##### SCENT PROJECT #####

rm(list = ls())

load("SCRVISscentnoscentGroupOwnCat.RData")

## Establish parameters for detection
sigma <- out$mean$sigma
lam0 <- c(out$mean$lam0[1], out$mean$lam0[2], out$mean$lam0[3])

## Grid
locs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
## Survey locations (9 transects with scent, 18 transects without, sensu scent project)
## ADGJMPSVY sprayed, BCEFHIKLNOQRTUWXZAA
# NS <- locs[c(14:39,53:78,92:117,131:156,170:195,209:234,248:273,287:312,326:351),]
# S <- locs[c(1:13,40:52,79:91,118:130,157:169,196:208,235:247,274:286,313:325),]
pts <- as.matrix(locs)
J <- nrow(pts)

## Status of transect, 3 spray days and then no spraying during 4 days of sampling
stat <- cbind(rep(rep(c(3,2,2), each=13),times=9),rep(rep(c(3,2,2), each=13),times=9),rep(rep(c(3,2,2), each=13),times=9),rep(rep(c(4,2,2), each=13),times=9))
STATUS <- stat - 1

## Define state-space of point process. (i.e., where animals live)
delta<- 11.874929
Xl<-min(locs[,1]) - delta
Xu<-max(locs[,1]) + delta
Yl<-min(locs[,2]) - delta
Yu<-max(locs[,2]) + delta
## Check area: 
A <- (Xu-Xl)*(Yu-Yl)

# Number of nights trapping ("one week", 4 days and then 3 day break)
K <- 4

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

## Simulate the encounters of this individual in every trap and repeat 100 times
p <- array(NA,dim=c(1000,J,K))
y <- array(NA,dim=c(1000,J,K))

set.seed(04052021)

for(s in 1:1000){    ## do 100 simulations
  for(i in 1:1){
    for(j in 1:J){
      for(k in 1:K){
        p[i,j,k] <- lam0[STATUS[j,k]]*exp(-(d2[i,j])/(2*sigma*sigma))
        y[s,j,k] <- rbinom(1,K,p[i,j,k])
      }
    }
  }
}

## BHNTZ traps, AEILPSWAA
## Observations on transects with no scent
NScap <- mean(colSums(as.data.frame(y[,c(14:39,53:78,92:117,131:156,170:195,209:234,248:273,287:312,326:351),4])))

## Observations on transects with fresh scent
FScap <- mean(colSums(as.data.frame(y[,c(1:13,40:52,79:91,118:130,157:169,196:208,235:247,274:286,313:325),1:3])))

## Observations on transects with old scent
OScap <- mean(colSums(as.data.frame(y[,c(1:13,40:52,79:91,118:130,157:169,196:208,235:247,274:286,313:325),4])))

## CPUE on transects with no scent (18 transects of 0.22 km distance for 4 days)
NScap/(18*0.22*4)

## CPUE on transects with fresh scent (9 transects of 0.22 km distance for 3 days)
FScap/(9*0.22*3)

## CPUE on transects with old scent (9 transects of 0.22 km distance for 1 day)
OScap/(9*0.22*1)
