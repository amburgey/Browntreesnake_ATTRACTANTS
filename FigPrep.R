## Figure of lure vs. no lure lam0 (baseline encounter rate)

rm(list = ls())

library(ggplot2); library(jagsUI); library(HDInterval)


##### FIGURE ONE #####


cpue <- read.csv("CPUEfromAYAreport.csv")

##### FIGURE TWO #####


load("SCRVISlurenolure.RData")


res <- as.data.frame(cbind(c("Abundance","Encounter Rate w/o lure", "Encounter Rate w/ lure"), c("Abundance", "Encounter Rate", "Encounter Rate"), rbind(out$mean$N, out$mean$lam0[1], out$mean$lam0[2]), rbind(hdi(out$sims.list$N)[1], hdi(out$sims.list$lam0[,1])[1], hdi(out$sims.list$lam0[,2])[1]), rbind(hdi(out$sims.list$N)[2], hdi(out$sims.list$lam0[,1])[2], hdi(out$sims.list$lam0[,2])[2])))

colnames(res) <- c("Parameter", "Type", "MeanEstimate", "Q2.5","Q97.5")

res[,3] <- as.numeric(as.character(res[,3]))
res[,4] <- as.numeric(as.character(res[,4]))
res[,5] <- as.numeric(as.character(res[,5]))

res$Parameter <- factor(res$Parameter, levels = c("Abundance","Encounter Rate w/o lure", "Encounter Rate w/ lure"))

addline_format <- function(x,...){
  gsub('x','\n',x)
}


plot1 <- ggplot(data=res, aes(x=Parameter, y=MeanEstimate, fill=Parameter)) + 
  geom_point(pch = 21, size = 7) + 
  geom_linerange(data=res, aes(ymin=Q2.5, ymax=Q97.5)) + 
  facet_wrap(~ Type, scales = "free") +
  theme(legend.position="none", strip.text.x = element_text(size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14)) + 
  scale_fill_brewer(palette = "Greens") + 
  ylab("Mean Estimate") +
  scale_x_discrete(breaks=unique(res$Parameter), labels=addline_format(c("Abundance", "Encounter Rate x without lure", "Encounter Rate x with lure")))

png(file="Estimates.png",width=8,height=6,units="in",res=300)
plot1
dev.off()
