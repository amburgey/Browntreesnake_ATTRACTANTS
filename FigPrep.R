## Figure of lure vs. no lure lam0 (baseline encounter rate)

rm(list = ls())

library(ggplot2); library(jagsUI); library(HDInterval)


##### FIGURE ONE #####


cpue <- read.csv("CPUEfromLureTruncated.csv"); names(cpue) <- c("Week","CPUE","CPUE")
cpue <- rbind(as.data.frame(matrix(c("Average",mean(cpue[,2]),mean(cpue[,3])), nrow = 1, ncol = 3, dimnames = list(1, c("Week","CPUE","CPUE")))),cpue)
cpue <- rbind(cpue[,1:2], cpue[,c(1,3)])
cpue$Type <- c(rep("With lure", times=10), rep("Without lure", times=10))
cpue$CPUE <- as.numeric(cpue$CPUE)
cpue$Week <- factor(cpue$Week, levels=c("Average","1","2","3","4","5","6","7","8","9"))
cpue$Symbol <- c("Average",rep("With lure", times=9), "Average", rep("Without lure", times=9))


plot1a <- ggplot(data = cpue, aes(x = Week, y = CPUE, fill = Type, group = Type)) +
  theme_classic() +
  theme(legend.position="bottom", axis.text = element_text(size=12), axis.title = element_text(size=13), legend.text = element_text(size = 13), legend.title = element_blank()) +
  geom_point(aes(shape = Symbol), size = 4) +
  scale_shape_manual(values = c(23,21,21)) +
  scale_fill_manual(values = c("#2E6EA6","#2EA6A2"), breaks = c("With lure", "Without lure"), guide = FALSE) +
  # ylab("Catch-per unit effort (snakes/km)") +
  scale_x_discrete("Week", breaks = c("Average","1","2","3","4","5","6","7","8","9")) +
  scale_y_continuous("Catch-per unit effort (snakes/km)", breaks = c(0,0.5,1,1.5,2,2.5,3), limits = c(0,3)) +
  geom_vline(xintercept = c(1.5, 2.5, 5.5, 8.5, 9.5, 10.5)) +
  annotate("rect",xmin=c(0.5),xmax=c(1.5),ymin=-Inf,ymax=Inf,fill="#999999",alpha=0.2)

# png(file="RawDetections.png",width=8,height=6,units="in",res=300)
# plot1a
# dev.off()

cpuesc <- read.csv("CPUEfromScent.csv"); names(cpuesc) <- c("Week","CPUE","CPUE","CPUE")
cpuesc <- rbind(as.data.frame(matrix(c("Average",mean(cpuesc[,2]),mean(cpuesc[,3]),mean(cpuesc[,4])), nrow = 1, ncol = 4, dimnames = list(1, c("Week","CPUE","CPUE","CPUE")))),cpuesc)
cpuesc <- rbind(cpuesc[,1:2], cpuesc[,c(1,3)], cpuesc[,c(1,4)])
cpuesc$Type <- c(rep("Without scent", times=11), rep("Old scent", times=11), rep("Fresh scent", times=11))
cpuesc$CPUE <- as.numeric(cpuesc$CPUE)
cpuesc$Week <- factor(cpuesc$Week, levels=c("Average","1","2","3","4","5","6","7","8","9","10"))
cpuesc$Symbol <- c("Average",rep("Without scent", times=10), "Average", rep("Old scent", times=10), "Average", rep("Fresh scent", times=10))
## remove 0s so they don't plot
## No old scent in 4th week, no fresh scent in 4th week (0 snakes captured for old scent in week 5)
cpuesc <- cpuesc[-c(16,27),]


plot1b <- ggplot(data = cpuesc, aes(x = Week, y = CPUE, fill = Type, group = Type)) +
  theme_classic() +
  theme(legend.position="bottom", axis.text = element_text(size=12), axis.title = element_text(size=13), legend.text = element_text(size = 13), legend.title = element_blank()) +
  geom_point(aes(shape = Symbol), size = 4) +
  scale_shape_manual(values = c(23,21,21,21)) +
  scale_fill_manual(values = c("#ffc12b","#ee8010","#8e4c09"), breaks = c("Without scent", "Old scent", "Fresh scent"), guide = FALSE) +
  # ylab("Catch-per unit effort (snakes/km)") +
  scale_x_discrete("Week", breaks = c("Average","1","2","3","4","5","6","7","8","9","10")) +
  scale_y_continuous("Catch-per unit effort (snakes/km)", breaks = c(0,0.25,0.5,0.75,1,1.25,1.5), limits = c(0,1.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5)) +
  annotate("rect",xmin=c(0.5),xmax=c(1.5),ymin=-Inf,ymax=Inf,fill="#999999",alpha=0.2) +    geom_text(x=5, y=0.01, label="**")

library(ggpubr)
png(file="RawDetectionsBoth.png",width=8,height=8.5,units="in",res=600)
ggarrange(plot1a, plot1b, nrow = 2, ncol = 1, labels = "AUTO")
dev.off()


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


plot2 <- ggplot(data=res, aes(x=Parameter, y=MeanEstimate, fill=Parameter)) + 
  geom_point(pch = 21, size = 7) + 
  geom_linerange(data=res, aes(ymin=Q2.5, ymax=Q97.5)) + 
  facet_wrap(~ Type, scales = "free") +
  theme(legend.position="none", strip.text.x = element_text(size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14)) + 
  scale_fill_brewer(palette = "Greens") + 
  ylab("Mean Estimate") +
  scale_x_discrete(breaks=unique(res$Parameter), labels=addline_format(c("Abundance", "Encounter Rate x without lure", "Encounter Rate x with lure")))

png(file="EstimatesLures.png",width=8,height=6,units="in",res=300)
plot2
dev.off()



##### FIGURE THREE #####

## No spray (1), Sprayed (2 = grouped sprayed today and yesterday)
load("SCRVISscentenoscentGroupScent.RData")


res <- as.data.frame(cbind(c("Abundance","Encounter Rate w/o scent", "Encounter Rate w/ scent"), c("Abundance", "Encounter Rate", "Encounter Rate"), rbind(out$mean$N, out$mean$lam0[1], out$mean$lam0[2]), rbind(hdi(out$sims.list$N)[1], hdi(out$sims.list$lam0[,1])[1], hdi(out$sims.list$lam0[,2])[1]), rbind(hdi(out$sims.list$N)[2], hdi(out$sims.list$lam0[,1])[2], hdi(out$sims.list$lam0[,2])[2])))

colnames(res) <- c("Parameter", "Type", "MeanEstimate", "Q2.5","Q97.5")

res[,3] <- as.numeric(as.character(res[,3]))
res[,4] <- as.numeric(as.character(res[,4]))
res[,5] <- as.numeric(as.character(res[,5]))

res$Parameter <- factor(res$Parameter, levels = c("Abundance","Encounter Rate w/o scent", "Encounter Rate w/ scent"))

addline_format <- function(x,...){
  gsub('x','\n',x)
}


plot3 <- ggplot(data=res, aes(x=Parameter, y=MeanEstimate, fill=Parameter)) + 
  geom_point(pch = 21, size = 7) + 
  geom_linerange(data=res, aes(ymin=Q2.5, ymax=Q97.5)) + 
  facet_wrap(~ Type, scales = "free") +
  theme(legend.position="none", strip.text.x = element_text(size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(color="black", size=14, face="bold.italic")) + 
  scale_fill_brewer(palette = "Oranges") + 
  ylab("Mean Estimate") +
  scale_x_discrete(breaks=unique(res$Parameter), labels=addline_format(c("Abundance", "Encounter Rate x without scent", "Encounter Rate x with scent"))) +
  ggtitle("With scent category = sprayed today and sprayed yesterday grouped")


## No spray (1 = grouped not sprayed today and sprayed yesterday), Sprayed (2)
load("SCRVISscentnoscentGroupNoScent.RData")


res <- as.data.frame(cbind(c("Abundance","Encounter Rate w/o scent", "Encounter Rate w/ scent"), c("Abundance", "Encounter Rate", "Encounter Rate"), rbind(out$mean$N, out$mean$lam0[1], out$mean$lam0[2]), rbind(hdi(out$sims.list$N)[1], hdi(out$sims.list$lam0[,1])[1], hdi(out$sims.list$lam0[,2])[1]), rbind(hdi(out$sims.list$N)[2], hdi(out$sims.list$lam0[,1])[2], hdi(out$sims.list$lam0[,2])[2])))

colnames(res) <- c("Parameter", "Type", "MeanEstimate", "Q2.5","Q97.5")

res[,3] <- as.numeric(as.character(res[,3]))
res[,4] <- as.numeric(as.character(res[,4]))
res[,5] <- as.numeric(as.character(res[,5]))

res$Parameter <- factor(res$Parameter, levels = c("Abundance","Encounter Rate w/o scent", "Encounter Rate w/ scent"))

addline_format <- function(x,...){
  gsub('x','\n',x)
}


plot4 <- ggplot(data=res, aes(x=Parameter, y=MeanEstimate, fill=Parameter)) + 
  geom_point(pch = 21, size = 7) + 
  geom_linerange(data=res, aes(ymin=Q2.5, ymax=Q97.5)) + 
  facet_wrap(~ Type, scales = "free") +
  theme(legend.position="none", strip.text.x = element_text(size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(color="black", size=14, face="bold.italic")) + 
  scale_fill_brewer(palette = "Oranges") + 
  ylab("Mean Estimate") +
  scale_x_discrete(breaks=unique(res$Parameter), labels=addline_format(c("Abundance", "Encounter Rate x without scent", "Encounter Rate x with scent"))) +
  ggtitle("With scent category = sprayed today, No scent cateogry = not sprayed and sprayed yesterday grouped")


## No spray (1), Sprayed today (2), Sprayed yesterday (3)
load("SCRVISscentnoscentGroupOwnCat.RData")


res <- as.data.frame(cbind(c("Abundance","Encounter Rate w/o scent", "Encounter Rate w/ fresh scent", "Encounter Rate w/ old scent"), c("Abundance", "Encounter Rate", "Encounter Rate", "Encounter Rate"), rbind(out$mean$N, out$mean$lam0[1], out$mean$lam0[2],out$mean$lam0[3]), rbind(hdi(out$sims.list$N)[1], hdi(out$sims.list$lam0[,1])[1], hdi(out$sims.list$lam0[,2])[1], hdi(out$sims.list$lam0[,3])[1]), rbind(hdi(out$sims.list$N)[2], hdi(out$sims.list$lam0[,1])[2], hdi(out$sims.list$lam0[,2])[2], hdi(out$sims.list$lam0[,3])[2])))

colnames(res) <- c("Parameter", "Type", "MeanEstimate", "Q2.5","Q97.5")

res[,3] <- as.numeric(as.character(res[,3]))
res[,4] <- as.numeric(as.character(res[,4]))
res[,5] <- as.numeric(as.character(res[,5]))

res$Parameter <- factor(res$Parameter, levels = c("Abundance","Encounter Rate w/o scent", "Encounter Rate w/ fresh scent", "Encounter Rate w/ old scent"))

addline_format <- function(x,...){
  gsub('x','\n',x)
}


plot5 <- ggplot(data=res, aes(x=Parameter, y=MeanEstimate, fill=Parameter)) + 
  geom_point(pch = 21, size = 7) + 
  geom_linerange(data=res, aes(ymin=Q2.5, ymax=Q97.5)) + 
  facet_wrap(~ Type, scales = "free") +
  theme(legend.position="none", strip.text.x = element_text(size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(color="black", size=14, face="bold.italic")) + 
  scale_fill_brewer(palette = "Oranges") + 
  ylab("Mean Estimate") +
  scale_x_discrete(breaks=unique(res$Parameter), labels=addline_format(c("Abundance", "Encounter Rate x without scent", "Encounter Rate x with fresh scent", "Encounter Rate x with old scent")))  #+
  # ggtitle("No scent, sprayed today, and sprayed yesterday separate")

library(ggpubr)
# png(file="EstimatesScentGroups.png",width=10,height=8,units="in",res=300)
# ggarrange(plot3, plot4, plot5, nrow = 3, ncol = 1)
# dev.off()
# png(file="EstimatesScents.png",width=8,height=6,units="in",res=300)
# plot5
# dev.off()
## All together
png(file="EstimatesLuresScents.png",width=10,height=6,units="in",res=600)
ggarrange(plot2, plot5, nrow = 2, ncol = 1, labels = "AUTO")
dev.off()

