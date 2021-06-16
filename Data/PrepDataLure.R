### Preparing data for SCR analysis for visual survey - lure vs no lure project
### Data was collected in 2015 where biologists walked transects in CP that either had traps with mice (lure) or no traps (no lure)
### The goal is to understand if encounter probability is different between these two scenarios
### In the case of an incipient or suppressed (low density) population, any way that can maximize encounter rates (and thus detection probability) of brown treesnakes can be a useful resource for surveyors

library(dplyr)
library(reshape2)
library(tidyverse)
library(secr)
library(jagsUI)
library(data.table)
library(abind)


PrepDat <- function(caps,survs){

  ##### CAPTURE DATA #####
  ## Remove transects that aren't "core" - i.e., not part of the main survey grid
  caps <- caps[!(caps$TRANSECT=="SWE" | caps$TRANSECT=="NEE"),]
  caps <- droplevels(caps)
  ## Make sure recognized as date
  caps$Date2 <- as.Date(as.character(caps$Date), format = "%d-%b-%y")
  ## Add transect ID
  caps$TRANID <- paste(caps$TRANSECT,caps$LOCATION, sep="")
  ## One EFFORTID seems to be incorrect for a survey - changed EFFORTID 12219 for X5 on 2015-03-03 to EFFORTID 12229
  caps$EFFORTID[which(caps$EFFORTID == "12219" & caps$Date2 == "2015-03-03" & caps$TRANID == "X5")] <- "12229"
  
  
  ##### SURVEY DATA #####
  ## Remove transects that aren't "core" - i.e., not part of the main survey grid
  survs <- survs[!(survs$TRANSECT=="SWE" | survs$TRANSECT=="NEE"),]
  ## Make sure recognized as date
  survs$Date2 <- as.Date(as.character(survs$Date), format = "%d-%b-%y")
  ## Expand survey records to have a record per transect point (1-13) rather than just for overall transect
  survpts <- survs[rep(seq_len(nrow(survs)), each = 13), ]
  survpts$TRANID <- paste(survpts$TRANSECT, rep(1:13, times = (dim(survpts)[1]/13)), sep = "")
  
  ## Create vector to use for sorting
  siteord <- c(unique(survpts$TRANID)[1:13], unique(survpts$TRANID)[27:52], unique(survpts$TRANID)[313:325], unique(survpts$TRANID)[53:182], unique(survpts$TRANID)[326:338], unique(survpts$TRANID)[183:247], unique(survpts$TRANID)[339:351], unique(survpts$TRANID)[248:312], unique(survpts$TRANID)[14:26])
  
  
  ##### RESHAPE FOR ANALYSIS #####
  ## Add in zero captures where survey done and no animal found
  all <- merge(caps, survpts, by=c("EFFORTID","Date","Date2","TRANSECT","TRANID"), all = TRUE)
  
  ## Remove duplicate surveys for AA, one with capture and one with no capture causes unnecessary rows
  all2 <- unique(all[,c(1:12)])
  
  ## FLRX, BHNZ, DJPV, FLRX are transects that have mouse lures
  ## T transect has a lure when BHNZ is being run (23-March to 30-March) but no lure otherwise
  ## Remove two snakes found on forest edge and concrete barrier (i.e., not on transect)
  all3 <- all2[!(!is.na(all2$PITTAG) & is.na(all2$LOCATION)),]
  ## Rows of all3 don't match surveypts because sometimes more than one snake found at a point
  ## Removed incidental snakes found on other transects while searching (i.e., wouldn't show up as being surveyed in the survey file and difficult to quantify effort here)
  all3 <- all3[!(all3$EFFORTID == "12096" & all3$TRANID == "Q8"),]
  all3 <- all3[!(all3$EFFORTID == "12144" & all3$TRANID == "T12"),]
  all3 <- all3[!(all3$EFFORTID == "12191" & all3$TRANID == "W7"),]
  all3 <- all3[!(all3$EFFORTID == "12352" & all3$TRANID == "Z5"),]
  all3 <- all3[!(all3$EFFORTID == "12335" & all3$TRANID == "B13"),]
  all3 <- all3[!(all3$EFFORTID == "12279" & all3$TRANID == "C1"),]
  all3<- all3[!(all3$EFFORTID == "12318" & all3$TRANID == "T3"),]
  all3 <- all3[!(all3$EFFORTID == "12319" & all3$TRANID == "H6"),]
  
  ## Add TYPE == no lure (NTL) to all with missing survey info (based on survey datasheets)
  for(i in 1:nrow(all3)){
    if(!is.na(all3[i,6]) & is.na(all3[i,12])){
      all3[i,12] <- "NTL"
    }
  }
  ## Specify Transect is factor and order
  all3$TRANID <- factor(all3$TRANID,levels=siteord)
  all3 <- all3[order(all3$TRANID),]
  
  ##### Create Transect grid cell by Date dataframe to indicate active/inactive (walked/not walked) #####
  act <- all3[,c("Date2","TRANID")]
  act$Active <- 1
  act <- reshape2::dcast(act, TRANID ~ Date2, fun.aggregate = length, value.var = "Active")
  ## Change to 1 (active) or 0 (inactive)
  act <- act %>% mutate_if(is.numeric, ~1 * (. > 0))
  
  ##### Create TRAP by Date dataframe #####
  ## STATUS of trap on each date
  ## Inactive = 1
  ## Active with no lure = 2
  ## Active with lure = 3
  ## Subset to relevant columns
  suball <- all3[,c("Date2","TRANID","TYPE")]
  ## Give numeric indicator
  suball$TYPE <- ifelse(suball$TYPE == "NTL", 2, 3)
  stat <- reshape2::dcast(suball, TRANID ~ Date2, value.var = "TYPE", fun.aggregate = unique, fill = 1)
  
  ##### Create PITTAG by Date dataframe #####
  ## Add 0 or 1 to indicate captured snake at survey
  ## Missing 028347041 from original captures file, never seen during this period of time aside from once when not on transect (removed above)
  ## Sum of captures equals capture raw data - 2 snakes removed above not on transects
  snks <- reshape2::acast(data = all3, formula = PITTAG ~ TRANID ~ Date2, fun.aggregate = length, value.var = "TYPE")
  snks <- snks[-dim(snks)[1],,]
  
  prepdat <- list(act = act, snks = snks, stat = stat)
  
  return(prepdat)

}

