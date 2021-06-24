## Figure details and explanations of content can be found in the published manuscript - see README for more details

rm(list = ls())

library(ggplot2); library(jagsUI); library(HDInterval);library(ggpubr)


##### FIGURE ONE. Catch Per Unit Effort ----

# Catch per unit effort for lure project
cpue <- read.csv("Data/CPUEfromLureTruncated.csv"); names(cpue) <- c("Week","CPUE","CPUE")
# Create a dataframe with details on catch per unit effort by each week
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
  scale_x_discrete("Week", breaks = c("Average","1","2","3","4","5","6","7","8","9")) +
  scale_y_continuous("Catch-per unit effort (snakes/km)", breaks = c(0,0.5,1,1.5,2,2.5,3), limits = c(0,3)) +
  geom_vline(xintercept = c(1.5, 2.5, 5.5, 8.5, 9.5, 10.5)) +
  annotate("rect",xmin=c(0.5),xmax=c(1.5),ymin=-Inf,ymax=Inf,fill="#999999",alpha=0.2)

# Catch per unit effort for scent project
cpuesc <- read.csv("Data/CPUEfromScent.csv"); names(cpuesc) <- c("Week","CPUE","CPUE","CPUE")
# Create a dataframe with details on catch per unit effort by each week
cpuesc <- rbind(as.data.frame(matrix(c("Average",mean(cpuesc[,2]),mean(cpuesc[,3]),mean(cpuesc[,4])), nrow = 1, ncol = 4, dimnames = list(1, c("Week","CPUE","CPUE","CPUE")))),cpuesc)
cpuesc <- rbind(cpuesc[,1:2], cpuesc[,c(1,3)], cpuesc[,c(1,4)])
cpuesc$Type <- c(rep("Without scent", times=11), rep("Old scent", times=11), rep("Fresh scent", times=11))
cpuesc$CPUE <- as.numeric(cpuesc$CPUE)
cpuesc$Week <- factor(cpuesc$Week, levels=c("Average","1","2","3","4","5","6","7","8","9","10"))
cpuesc$Symbol <- c("Average",rep("Without scent", times=10), "Average", rep("Old scent", times=10), "Average", rep("Fresh scent", times=10))
## remove 0s so they don't plot and obscure x axis
## No old scent in 4th week, no fresh scent in 4th week (0 snakes captured for old scent in week 5)
cpuesc <- cpuesc[-c(16,27),]


plot1b <- ggplot(data = cpuesc, aes(x = Week, y = CPUE, fill = Type, group = Type)) +
  theme_classic() +
  theme(legend.position="bottom", axis.text = element_text(size=12), axis.title = element_text(size=13), legend.text = element_text(size = 13), legend.title = element_blank()) +
  geom_point(aes(shape = Symbol), size = 4) +
  scale_shape_manual(values = c(23,21,21,21)) +
  scale_fill_manual(values = c("#ffc12b","#ee8010","#8e4c09"), breaks = c("Without scent", "Old scent", "Fresh scent"), guide = FALSE) +
  scale_x_discrete("Week", breaks = c("Average","1","2","3","4","5","6","7","8","9","10")) +
  scale_y_continuous("Catch-per unit effort (snakes/km)", breaks = c(0,0.25,0.5,0.75,1,1.25,1.5), limits = c(0,1.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5)) +
  annotate("rect",xmin=c(0.5),xmax=c(1.5),ymin=-Inf,ymax=Inf,fill="#999999",alpha=0.2) +    geom_text(x=5, y=0.01, label="**")

png(file="Figures/RawDetectionsBoth.png",width=8,height=8.5,units="in",res=600)
ggarrange(plot1a, plot1b, nrow = 2, ncol = 1, labels = "AUTO")
dev.off()
## Post-processing of figures in image-editor allows for quick coloring in of symbols


##### FIGURE TWO. Comparison of abundances and encounter rates ----

# Estimates for lure project
load("Results/SCRVISlurenolure.RData")

res <- as.data.frame(cbind(c("Abundance","Encounter Rate w/o lure", "Encounter Rate w/ lure"), c("Abundance", "Encounter Rate", "Encounter Rate"), rbind(out$mean$N, out$mean$lam0[1], out$mean$lam0[2]), rbind(hdi(out$sims.list$N)[1], hdi(out$sims.list$lam0[,1])[1], hdi(out$sims.list$lam0[,2])[1]), rbind(hdi(out$sims.list$N)[2], hdi(out$sims.list$lam0[,1])[2], hdi(out$sims.list$lam0[,2])[2])))

colnames(res) <- c("Parameter", "Type", "MeanEstimate", "Q2.5","Q97.5")

res[,3] <- as.numeric(as.character(res[,3]))
res[,4] <- as.numeric(as.character(res[,4]))
res[,5] <- as.numeric(as.character(res[,5]))

res$Parameter <- factor(res$Parameter, levels = c("Abundance","Encounter Rate w/o lure", "Encounter Rate w/ lure"))

# Line breaks in axis labels
addline_format <- function(x,...){
  gsub('x','\n',x)
}

plot2a <- ggplot(data=res, aes(x=Parameter, y=MeanEstimate, fill=Parameter)) + 
  geom_point(pch = 21, size = 7) + 
  geom_linerange(data=res, aes(ymin=Q2.5, ymax=Q97.5)) + 
  facet_wrap(~ Type, scales = "free", dir = "v") +
  theme(legend.position="none", strip.text.x = element_text(size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14)) + 
  scale_fill_brewer(palette = "Greens") + 
  ylab(c("Mean Estimate")) +
  scale_x_discrete(breaks=unique(res$Parameter), labels=addline_format(c(" ", "Encounter Rate x without lure", "Encounter Rate x with lure")))


# Estimates for scent project
load("Results/SCRVISscentnoscent.RData")

res <- as.data.frame(cbind(c("Abundance","Encounter Rate w/o scent", "Encounter Rate w/ fresh scent", "Encounter Rate w/ old scent"), c("Abundance", "Encounter Rate", "Encounter Rate", "Encounter Rate"), rbind(out$mean$N, out$mean$lam0[1], out$mean$lam0[2],out$mean$lam0[3]), rbind(hdi(out$sims.list$N)[1], hdi(out$sims.list$lam0[,1])[1], hdi(out$sims.list$lam0[,2])[1], hdi(out$sims.list$lam0[,3])[1]), rbind(hdi(out$sims.list$N)[2], hdi(out$sims.list$lam0[,1])[2], hdi(out$sims.list$lam0[,2])[2], hdi(out$sims.list$lam0[,3])[2])))

colnames(res) <- c("Parameter", "Type", "MeanEstimate", "Q2.5","Q97.5")

res[,3] <- as.numeric(as.character(res[,3]))
res[,4] <- as.numeric(as.character(res[,4]))
res[,5] <- as.numeric(as.character(res[,5]))

res$Parameter <- factor(res$Parameter, levels = c("Abundance","Encounter Rate w/o scent", "Encounter Rate w/ fresh scent", "Encounter Rate w/ old scent"))

plot2b <- ggplot(data=res, aes(x=Parameter, y=MeanEstimate, fill=Parameter)) + 
  geom_point(pch = 21, size = 7) + 
  geom_linerange(data=res, aes(ymin=Q2.5, ymax=Q97.5)) + 
  facet_wrap(~ Type, scales = "free", dir = "v") +
  theme(legend.position="none", strip.text.x = element_text(size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(color="black", size=14, face="bold.italic")) + 
  scale_fill_brewer(palette = "Oranges") + 
  ylab("Mean Estimate") +
  scale_x_discrete(breaks=unique(res$Parameter), labels=addline_format(c(" ", "Encounter Rate x without scent", "Encounter Rate x with fresh scent", "Encounter Rate x with old scent")))

png(file="Figures/EstimatesLuresScents.png",width=10,height=6,units="in",res=600)
ggarrange(plot2a, plot2b, nrow = 1, ncol = 2, labels = "AUTO")
dev.off()
