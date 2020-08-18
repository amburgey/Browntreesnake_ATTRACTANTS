## Figure of lure vs. no lure lam0 (baseline encounter rate)

rm(list = ls())

library(ggplot2); library(jagsUI)

load("SCRVISlurenolure.RData")


res <- as.data.frame(cbind(c("Abundance","Encounter Rate w/o lure", "Encounter Rate w/ lure"), c("Abundance", "Encounter Rate", "Encounter Rate"), rbind(out$mean$N, out$mean$lam0[1], out$mean$lam0[2]), rbind(out$q2.5$N, out$q2.5$lam0[1], out$q2.5$lam0[2]), rbind(out$q97.5$N, out$q97.5$lam0[1], out$q97.5$lam0[2])))

colnames(res) <- c("Parameter", "Type", "MeanEstimate", "Q2.5","Q97.5")

res[,3] <- as.numeric(as.character(res[,3]))
res[,4] <- as.numeric(as.character(res[,4]))
res[,5] <- as.numeric(as.character(res[,5]))


plot1 <- ggplot(data=res, aes(x=Parameter, y=MeanEstimate, fill=Parameter)) + geom_point(pch = 21, size = 7) + 
  geom_linerange(data=res, aes(ymin=Q2.5, ymax=Q97.5)) + facet_wrap(~ Type, scales = "free") +
  theme(legend.position="none") + scale_fill_brewer(palette = "Greens")

png(file="Estimates.png",width=8,height=6,units="in",res=300)
plot1
dev.off()
