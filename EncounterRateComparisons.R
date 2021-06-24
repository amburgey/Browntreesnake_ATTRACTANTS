### Results: calculate the number of posterior samples where encounter rate1 > encounter rate 2

library(jagsUI)

## First: Probability that lamba[lure] > lambda[no lure]
load("Results/SCRVISlurenolure.RData")

## Difference between lambda[lure] - lambda[no lure]
temp1 <- out$sims.list$lam0[,2] - out$sims.list$lam0[,1]
## Number of times where lamba[lure] was larger
larger1 <- temp1[temp1 > 0]
## Proportion that lure was greater than no lure out of all iterations
length(larger1)/length(out$sims.list$lam0[,1])


## Second: Probability that lamba[no scent] > lambda[fresh scent], lambda[no scent] > lambda[old scent],lambda[old scent] > lambda[fresh scent]
load("Results/SCRVISscentnoscent.RData")

## Difference between lambda[no scent] - lambda[fresh scent]
temp2 <- out$sims.list$lam0[,1] - out$sims.list$lam0[,2]
## Number of times where lamba[no scent] was larger
larger2 <- temp2[temp2 >= 0]
## Proportion that no scent was greater than fresh scent out of all iterations
length(larger2)/length(out$sims.list$lam0[,1])

## Difference between lambda[no scent] - lambda[old scent]
temp3 <- out$sims.list$lam0[,1] - out$sims.list$lam0[,3]
## Number of times where lamba[no scent] was larger
larger3 <- temp3[temp3 >= 0]
## Proportion that no scent was greater than old scent out of all iterations
length(larger3)/length(out$sims.list$lam0[,1])

## Difference between lambda[old scent] - lambda[fresh scent]
temp4 <- out$sims.list$lam0[,3] - out$sims.list$lam0[,2]
## Number of times where lamba[old scent] was larger
larger4 <- temp4[temp4 >= 0]
## Proportion that old scent was greater than fresh scent out of all iterations
length(larger4)/length(out$sims.list$lam0[,1])

