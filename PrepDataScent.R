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
library(abind)


PrepDat <- function(caps,survs){

  ##### CAPTURE DATA #####
  ## Remove transects that aren't "core"
  caps <- caps[!(caps$TRANSECT=="SWE" | caps$TRANSECT=="NEE"| caps$TRANSECT=="PR" | caps$LOCATION=="0"),]
  caps <- droplevels(caps)
  caps$Date2 <- as.Date(as.character(caps$Date), format = "%d-%b-%y")
  ## Add transect ID
  caps$TRANID <- paste(caps$TRANSECT,caps$LOCATION, sep="")
  
  
  ##### SURVEY DATA #####
  survs <- survs[!(survs$TRANID=="SWE" | survs$TRANID=="NEE"),]
  ## Expand survey records to have a record per transect point (1-13) rather than just for overall transect
  survpts <- survs[rep(seq_len(nrow(survs)), each = 13), ]
  survpts$TRANID <- paste(survpts$TRANID, rep(1:13, times = 27), sep = "")
  rownames(survpts) <- NULL
  colnames(survpts) <- NULL
  
  ## Create matrix of active/inactive traps
  act <- ifelse(as.matrix(survpts[,-1]) > 0, 1, 0)
  
  ## Create matrix of spray vs unsprayed transects (1 = inactive, 2 = active and unsprayed, 3 = active and sprayed)
  scent <- ifelse(as.matrix(survpts[,-1]) >= 2, 3, 
                  ifelse(as.matrix(survpts[,-1]) == 1, 2,
                  ifelse(as.matrix(survpts[,-1]) == 0, 1, 999)))
  
  ## Create vector to use for sorting
  siteord <- survpts[,1]
  
  
  ##### RESHAPE FOR ANALYSIS #####
  ## Add in zero captures where survey done and no animal found
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
  
  ##### Split into TL and NTL datasets to create capture and active dataframes #####
  ## out of interest, quick sum of number of snakes caught on each kind of survey
  # ct1 <- ntl %>% summarise(NTLcount = sum(!is.na(PITTAG)))
  # sur1 <- ntl %>% summarise(NTLsurvs = sum(!is.na(Date2)))
  # ct2 <- tl %>% summarise(TLcount = sum(!is.na(PITTAG)))
  # sur2 <- tl %>% summarise(TLsurvs = sum(!is.na(Date2)))
  # ratioNTL <- ct1/sur1
  # ratioTL <- ct2/sur2  # Raw data looks slightly (very) higher rate of picking up snakes on VIS lure
  
  ##### Create Transect grid cell by Date dataframe to indicate active/inactive (walked/not walked) #####
  act <- all3[,c("Date2","TRANID")]
  act$Active <- 1
  act <- reshape2::dcast(act, TRANID ~ Date2, fun.aggregate = length, value.var = "Active")
  act <- act %>% mutate_if(is.numeric, ~1 * (. > 0))
  
  ##### Create TRAP by Date dataframe #####
  ## STATUS of trap on each date
  ## Inactive = 1
  ## Active with no lure = 2
  ## Active with lure = 3
  suball <- all3[,c(3,5,12)]
  suball$TYPE <- ifelse(suball$TYPE == "NTL", 2, 3)
  stat <- reshape2::dcast(suball, TRANID ~ Date2, value.var = "TYPE", fun.aggregate = unique, fill = 1)
  
  ##### Create PITTAG by Date dataframe #####
  ## Add 0 or 1 to indicate captured snake at survey
  ## Missing 028347041 from original captures file, never seen during this period of time aside from once when not on transect (removed above)
  ## Sum of captures equals capture raw data - 2 snakes removed above not on transects
  snks <- reshape2::acast(data = all3, formula = PITTAG ~ TRANID ~ Date2, fun.aggregate = length, value.var = "TYPE")[-110,,]
  
  prepdat <- list(act = act, snks = snks, stat = stat)
  
  return(prepdat)

}

