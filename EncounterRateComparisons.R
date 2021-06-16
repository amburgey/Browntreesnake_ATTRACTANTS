### Results: calculate the number of posterior samples where encounter rate1 > encounter rate 2

library(jagsUI)

## First: Probability that lamba[lure] > lambda[no lure]
load("Results/SCRVISlurenolure.RData")

## Difference between lambda[lure] - lambda[no lure]
temp <- out$sims.list$lam0[,2] - out$sims.list$lam0[,1]
## Number of times where lamba[lure] was larger
larger <- temp[temp > 0]
## Proportion that lure was greater than no lure out of all iterations
length(larger)/length(out$sims.list$lam0[,1])


## Second: Probability that lamba[no scent] > lambda[fresh scent], lambda[no scent] > lambda[old scent],lambda[old scent] > lambda[fresh scent]
load("Results/SCRVISscentnoscent.RData")

## Difference between lambda[no scent] - lambda[fresh scent]
temp <- out$sims.list$lam0[,1] - out$sims.list$lam0[,2]
## Number of times where lamba[no scent] was larger
larger <- temp[temp >= 0]
## Proportion that no scent was greater than fresh scent out of all iterations
length(larger)/length(out$sims.list$lam0[,1])

## Difference between lambda[no scent] - lambda[old scent]
temp <- out$sims.list$lam0[,1] - out$sims.list$lam0[,3]
## Number of times where lamba[no scent] was larger
larger <- temp[temp >= 0]
## Proportion that no scent was greater than old scent out of all iterations
length(larger)/length(out$sims.list$lam0[,1])

## Difference between lambda[old scent] - lambda[fresh scent]
temp <- out$sims.list$lam0[,3] - out$sims.list$lam0[,2]
## Number of times where lamba[old scent] was larger
larger <- temp[temp >= 0]
## Proportion that old scent was greater than fresh scent out of all iterations
length(larger)/length(out$sims.list$lam0[,1])

