### Preparing data for SCR analysis for VIS lure vs no lure project
### Data was collected in 2015 where biologists walked transects in CP that either had traps with mice (lure) or no traps (no lure)
### The goal is to understand if detection probability is different between these two scenarios
### In the case of an incipient or suppressed (low density) population, any way that can maximize detection probability of brown treesnakes can be a useful resource for surveyors

rm(list = ls())

library(dplyr)

##### CAPTURE DATA #####
caps <- read.csv("Captures.csv")[,c("EFFORTID","Date","PITTAG","SVL","TOTAL","WEIGHT","SEX","TRANSECT","LOCATION")]
## Remove transects that aren't "core"
caps <- caps[!(caps$TRANSECT=="SWE" | caps$TRANSECT=="NEE"),]
caps <- droplevels(caps)
caps$Date2 <- as.Date(as.character(caps$Date), format = "%d-%b-%y")
## Add transect ID
caps$TRANID <- paste(caps$TRANSECT,caps$LOCATION, sep="")


##### SURVEY DATA #####
survs <- read.csv("Surveys.csv")[,c("EFFORTID","Date","TRANSECT","BI","TYPE")]
survs$Date2 <- as.Date(as.character(survs$Date), format = "%d-%b-%y")
## Expand survey records to have a record per transect point (1-13) rather than just for overall transect
survpts <- survs[rep(seq_len(nrow(survs)), each = 13), ]
survpts$TRANID <- paste(survpts$TRANSECT, rep(1:13, times = 708), sep = "")


##### RESHAPE FOR ANALYSIS #####
## Add in zero captures where survey done and no animal found
## Keep in mind that transects with lures are those with cameras (so missing info on TL/NTL can be inferred)
all <- merge(caps, survpts, by=c("EFFORTID","Date","Date2","TRANSECT","TRANID"), all = TRUE)

## Remove duplicate surveys for AA, one with capture and one with no capture causes unnecessary rows
all2 <- unique(all[,c(1:11,13)])

## FLRX, BHNZ, DJPV, FLRX are camera transects that are TL
## T transect is TL when BHNZ is being run (23-March to 30-March) but NTL otherwise
## Remove two snakes found on forest edge and concrete barrier (i.e., not on transect)
all3 <- all2[!(!is.na(all2$PITTAG) & is.na(all2$LOCATION)),]
## Add TYPE == NTL to all with missing survey info, found in MASTER-ARCHIVAL-MARKREL-ALL file
for(i in 1:nrow(all3)){
  if(!is.na(all3[i,5]) & is.na(all3[i,12])){
    all3[i,12] <- "NTL"
  }
}

##### Create Transect grid cell by Date dataframe to indicate active/inactive (walked/not walked) #####



## Add 0 or 1 to indicate captured snake at survey
all3$Cap <- ifelse(is.na(all3$PITTAG), 0, 1)

##### Create PITTAG by Date dataframe #####
test <- dcast(data = all3, formula = PITTAG ~ Date2, length)


all3 <- all3[order(all3$Date2, all3$TRANSECT, all3$TRANID),]
survpts <- survpts[order(survpts$Date2, survpts$TRANSECT, survpts$TRANID),]

comp <- matrix(NA)
for (i in 1:nrow(all3)){
  if(all3[i,1] == survpts[i,1] & all3[i,2] == survpts[i,2] & all3[i,3] == survpts[i,6] & all3[i,4] == survpts[i,3] & all3[i,5] == survpts[i,7] & all3[i,12] == survpts[i,5]){
    comp[i] <- TRUE
  } else {
  comp[i] <- FALSE
  }
}



comparedf(test1,test2)

