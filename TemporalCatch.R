### Figure for number of animals caught each day

rm(list = ls())

library(ggplot2); library(dplyr); library(tidyverse); library(reshape2)


### Lure Trials
caps <- read.csv("Captures.csv")[,c("EFFORTID","Date","PITTAG","SVL","TOTAL","WEIGHT","SEX","TRANSECT","LOCATION")]
caps <- caps %>%
  filter(!grepl('Apr', Date))

survs <- read.csv("Surveys.csv")[,c("EFFORTID","Date","TRANSECT","TYPE")]
survs <- survs %>%
  filter(!grepl('Apr', Date))

##### CAPTURE DATA #####
## Remove transects that aren't "core"
caps <- caps[!(caps$TRANSECT=="SWE" | caps$TRANSECT=="NEE"),]
caps <- droplevels(caps)
caps$Date2 <- as.Date(as.character(caps$Date), format = "%d-%b-%y")
## Add transect ID
caps$TRANID <- paste(caps$TRANSECT,caps$LOCATION, sep="")
## One EFFORTID seems to be incorrect for a survey - changed EFFORTID 12219 for X5 on 2015-03-03 to EFFORTID 12229
caps$EFFORTID[which(caps$EFFORTID == "12219" & caps$Date2 == "2015-03-03" & caps$TRANID == "X5")] <- "12229"


##### SURVEY DATA #####
survs$Date2 <- as.Date(as.character(survs$Date), format = "%d-%b-%y")
## Expand survey records to have a record per transect point (1-13) rather than just for overall transect
survpts <- survs[rep(seq_len(nrow(survs)), each = 13), ]
survpts$TRANID <- paste(survpts$TRANSECT, rep(1:13, times = (dim(survpts)[1]/13)), sep = "")

## Create vector to use for sorting
siteord <- c(unique(survpts$TRANID)[1:13], unique(survpts$TRANID)[27:52], unique(survpts$TRANID)[313:325], unique(survpts$TRANID)[53:182], unique(survpts$TRANID)[326:338], unique(survpts$TRANID)[183:247], unique(survpts$TRANID)[339:351], unique(survpts$TRANID)[248:312], unique(survpts$TRANID)[14:26])


##### RESHAPE FOR ANALYSIS #####
## Add in zero captures where survey done and no animal found
## Keep in mind that transects with lures are those with cameras (so missing info on TL/NTL can be inferred)
all <- merge(caps, survpts, by=c("EFFORTID","Date","Date2","TRANSECT","TRANID"), all = TRUE)

## Remove duplicate surveys for AA, one with capture and one with no capture causes unnecessary rows
all2 <- unique(all[,c(1:12)])

## FLRX, BHNZ, DJPV, FLRX are camera transects that are TL
## T transect is TL when BHNZ is being run (23-March to 30-March) but NTL otherwise
## Remove two snakes found on forest edge and concrete barrier (i.e., not on transect)
all3 <- all2[!(!is.na(all2$PITTAG) & is.na(all2$LOCATION)),]
## all3 has a weird number of rows that doesn't match surveypts because sometimes more than one snake found at a point
## Some snakes found on random other transects while searching (i.e., wouldn't show up as being surveyed in the survey file) > removed
all3 <- all3[!(all3$EFFORTID == "12096" & all3$TRANID == "Q8"),]
all3 <- all3[!(all3$EFFORTID == "12144" & all3$TRANID == "T12"),]
all3 <- all3[!(all3$EFFORTID == "12191" & all3$TRANID == "W7"),]
all3 <- all3[!(all3$EFFORTID == "12352" & all3$TRANID == "Z5"),]
all3 <- all3[!(all3$EFFORTID == "12335" & all3$TRANID == "B13"),]
all3 <- all3[!(all3$EFFORTID == "12279" & all3$TRANID == "C1"),]
all3<- all3[!(all3$EFFORTID == "12318" & all3$TRANID == "T3"),]
all3 <- all3[!(all3$EFFORTID == "12319" & all3$TRANID == "H6"),]

## Add TYPE == NTL to all with missing survey info, found in MASTER-ARCHIVAL-MARKREL-ALL file
for(i in 1:nrow(all3)){
  if(!is.na(all3[i,6]) & is.na(all3[i,12])){
    all3[i,12] <- "NTL"
  }
}
all3$TRANID <- factor(all3$TRANID,levels=siteord)
all3 <- all3[order(all3$TRANID),]
all3$Snk <- 1

## Subset to captured snakes
snks <- subset(all3, !is.na(PITTAG))
snkssum <- aggregate(snks$Snk, by=list(snks$Date2, snks$TYPE), sum)
snkssum <- snkssum[order(as.Date(snkssum$Group.1, format="%Y-%m-%d")),]
colnames(snkssum) <- c("Date","Type","Count")
## No snakes caught at one instance so add zero
snkssum <- snkssum %>% add_row(Date = as.Date("2015-03-29"), Type = "TL", Count = 0, .before = 48)

## Add weeks and same trap location
snkssum$Week <- rep(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,9,9), each = 2)
snkssum$Location <- rep(c(1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,4,5,5), each = 2)


## Plot snakes by date
plot1 <- ggplot(snkssum, aes(x=Date, y=Count, color=Type, shape=Type)) + 
  geom_point(size=4) +
  scale_shape_manual(values=c(16, 17)) +
  scale_color_manual(values=c("#2E6EA6","#2EA6A2")) +
  theme(axis.text.x = element_text(angle = 90)) +  #, panel.background = element_rect(fill="white", colour = "white")
  scale_x_date(breaks = "3 days") +
  geom_vline(xintercept = c(as.numeric(as.Date("2015-02-06")), as.numeric(as.Date("2015-02-27")), as.numeric(as.Date("2015-03-23")), as.numeric(as.Date("2015-03-27"))))



### Scent Trials
# Read in capture data
caps <- read_csv("CapturesScent.csv")[,c("EFFORTID","Date","PITTAG","SVL","TOTAL","WEIGHT","SEX","TRANSECT","LOCATION")]
# Restrict to 2-month period to meet assumptions of population closure
caps <- caps %>%
  filter(!grepl('Oct', Date)) %>%
  filter(!grepl('Jan', Date))

# Read in survey data
survs <- read_csv("ScentNoScentTrail.csv")
survs <- survs %>%
  select(-contains("Oct")) %>%
  select(-contains("Jan"))

## Remove transects that aren't "core"
caps <- caps[!(caps$TRANSECT=="SWE" | caps$TRANSECT=="NEE"| caps$TRANSECT=="PR" | caps$LOCATION=="0" | caps$LOCATION=="14"),]
caps <- droplevels(caps)
caps$Date2 <- as.Date(as.character(caps$Date), format = "%d-%b-%y")
## Add transect ID
caps$TRANID <- paste(caps$TRANSECT,caps$LOCATION, sep="")

##### SURVEY DATA #####
survs <- survs[!(survs$TRANID=="SWE" | survs$TRANID=="NEE"),]
## Expand survey records to have a record per transect point (1-13) rather than just for overall transect
survpts <- survs[rep(seq_len(nrow(survs)), each = 13), ]
survpts$TRANID <- paste(survpts$TRANID, rep(1:13, times = 27), sep = "")

## Create vector to use for sorting
siteord <- survpts[,1]
colnames(siteord) <- c("TRANID")

## Add 0 or 1 to indicate captured snake at survey
caps$TRANID <- factor(caps$TRANID, levels = as.vector(as.matrix(siteord)))
caps$Act <- 1
caps <- caps[order(caps$TRANID),]
caps$Snk <- 1
## Reshape survey info to type in row
survpts <- melt(survpts)
survpts$value <- ifelse(survpts$value == 1, "NS",
                  ifelse(survpts$value == 2, "FS",
                         ifelse(survpts$value == 3, "OS", -9999)))
survpts$Date2 <- as.Date(survpts$variable, format ="%d-%B-%y")
colnames(survpts)[3] <- c("TYPE")
caps <- merge(caps, survpts[,-2], by=c("TRANID","Date2"))

## Subset to captured snakes
snks2 <- subset(caps, !is.na(PITTAG))
snkssum2 <- aggregate(snks2$Snk, by=list(snks2$Date2, snks2$TYPE), sum)
snkssum2 <- snkssum2[order(as.Date(snkssum2$Group.1, format="%Y-%m-%d")),]
colnames(snkssum2) <- c("Date","Type","Count")


## Plot snakes by date
plot2 <- ggplot(snkssum2, aes(x=Date, y=Count, color=Type, shape=Type)) + 
  geom_point(size=4) +
  scale_shape_manual(values=c(16, 17, 18)) +
  scale_color_manual(values=c("#ffc12b","#ee8010","#8e4c09")) +
  theme(axis.text.x = element_text(angle = 90)) +  #, panel.background = element_rect(fill="white", colour = "white")
  scale_x_date(breaks = "3 days") +
  geom_vline(xintercept = c(as.numeric(as.Date("2016-11-05")), as.numeric(as.Date("2016-11-12")), as.numeric(as.Date("2016-11-16")), as.numeric(as.Date("2016-11-25")), as.numeric(as.Date("2016-12-03")), as.numeric(as.Date("2016-12-10")), as.numeric(as.Date("2016-12-17")), as.numeric(as.Date("2016-12-24"))))

