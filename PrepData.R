### Preparing data for SCR analysis for VIS lure vs no lure project
### Data was collected in 2015 where biologists walked transects in CP that either had traps with mice (lure) or no traps (no lure)
### The goal is to understand if detection probability is different between these two scenarios
### In the case of an incipient or suppressed (low density) population, any way that can maximize detection probability of brown treesnakes can be a useful resource for surveyors

library(dplyr)
library(reshape2)
library(tidyverse)
library(secr)
library(jagsUI)
library(data.table)


PrepDat <- function(caps,survs){

  ##### CAPTURE DATA #####
  ## Remove transects that aren't "core"
  caps <- caps[!(caps$TRANSECT=="SWE" | caps$TRANSECT=="NEE"),]
  caps <- droplevels(caps)
  caps$Date2 <- as.Date(as.character(caps$Date), format = "%d-%b-%y")
  ## Add transect ID
  caps$TRANID <- paste(caps$TRANSECT,caps$LOCATION, sep="")
  
  
  ##### SURVEY DATA #####
  survs$Date2 <- as.Date(as.character(survs$Date), format = "%d-%b-%y")
  ## Expand survey records to have a record per transect point (1-13) rather than just for overall transect
  survpts <- survs[rep(seq_len(nrow(survs)), each = 13), ]
  survpts$TRANID <- paste(survpts$TRANSECT, rep(1:13, times = 708), sep = "")
  
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
  ## all3 has a weird number of rows that doesn't match surveypts because sometimes more than one snake found at a point and also some snakes found on random other transects while searching (i.e., wouldn't show up as being surveyed in the survey file)
  ## Add TYPE == NTL to all with missing survey info, found in MASTER-ARCHIVAL-MARKREL-ALL file
  for(i in 1:nrow(all3)){
    if(!is.na(all3[i,6]) & is.na(all3[i,12])){
      all3[i,12] <- "NTL"
    }
  }
  
  
  ##### Split into TL and NTL datasets to create capture and active dataframes #####
  ntl <- subset(all3, TYPE == "NTL")
  tl <- subset(all3, TYPE == "TL")
  ## out of interest, quick sum of number of snakes caught on each kind of survey
  # ct1 <- ntl %>% summarise(NTLcount = sum(!is.na(PITTAG)))
  # sur1 <- ntl %>% summarise(NTLsurvs = sum(!is.na(Date2)))
  # ct2 <- tl %>% summarise(TLcount = sum(!is.na(PITTAG)))
  # sur2 <- tl %>% summarise(TLsurvs = sum(!is.na(Date2)))
  # ratioNTL <- ct1/sur1
  # ratioTL <- ct2/sur2  # Raw data looks slightly (very) higher rate of picking up snakes on VIS lure
  
  ##### Create Transect grid cell by Date dataframe to indicate active/inactive (walked/not walked) #####
  actNTL <- ntl[,c("Date2","TRANID")]
  actNTL$Active <- 1
  actNTL <- dcast(actNTL, TRANID ~ Date2)
  actNTL <- actNTL %>% mutate_if(is.numeric, ~1 * (. > 0))  
  ## Sort here so everything matching below
  actNTL$test <- factor(actNTL$TRANID,levels=siteord)
  actNTL <- actNTL[order(actNTL$test),]
  
  actTL <- tl[,c("Date2","TRANID")]
  actTL$Active <- 1
  actTL <- dcast(actTL, TRANID ~ Date2)
  actTL <- actTL %>% mutate_if(is.numeric, ~1 * (. > 0))
  ## Sort here so everything matching below
  actTL$test <- factor(actTL$TRANID,levels=siteord)
  actTL <- actTL[order(actTL$test),]
  
  
  ##### Create PITTAG by Date dataframe #####
  ## Add 0 or 1 to indicate captured snake at survey
  ## Missing 028347041 from original captures file, never seen during this period of time aside from once when not on transect (removed above)
  ## Sum of captures equals capture raw data - 2 snakes removed above not on transects
  ntl$Cap <- ifelse(is.na(ntl$PITTAG), 0, 1)
  NL <- ntl[,12]
  snksNTL <- dcast(data = ntl, formula = PITTAG ~ TRANID)[-95,]
  snksNTL <- setcolorder(snksNTL, siteord)
  
  tl$Cap <- ifelse(is.na(tl$PITTAG), 0, 1)
  L <- tl[,12]
  snksTL <- dcast(data = tl, formula = PITTAG ~ TRANID)[-63,]
  snksTL <- setcolorder(snksTL, actTL$TRANID)
  
  
  prepdat <- list(actNTL = actNTL, snksNTL = snksNTL, actTL = actTL, snksTL = snksTL)
  
  return(prepdat)

}

