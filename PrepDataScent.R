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
  dates <- colnames(survpts)[-1]
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
  colnames(siteord) <- c("TRANID")
  
  
  ##### RESHAPE FOR ANALYSIS #####
  ## Check that no dates are missing from capture data (i.e., an animal was never found during a survey)
  if(length(dates) != length(unique(caps$Date2)))
     stop("Mismatch in number of surveys in cap and surv data")
  
  ##### Create PITTAG by Date dataframe #####
  ## Add in TRANID where no snake ever captured
  sched <- cbind(as.data.frame(rep(as.matrix(siteord), each=length(dates))), as.data.frame(rep(dates, times = 27)))
  colnames(sched) <- c("TRANID","Date2")
  sched$Date2 <- as.Date(sched$Date2, format = "%d-%b-%y")
  sched$TRANID <- as.factor(sched$TRANID)
  sched$Act <- 0
  
  ## Find names of missing columns
  Missing <- setdiff(as.vector(as.matrix(siteord)),unique(caps$TRANID))
  
  ## Add 0 or 1 to indicate captured snake at survey
  caps$TRANID <- factor(caps$TRANID, levels = as.vector(as.matrix(siteord)))
  caps$Act <- 1
  caps <- merge(caps, sched, by = c("TRANID","Date2","Act"), all = TRUE)
  caps <- caps[order(caps$TRANID),]
  snks <- reshape2::acast(data = caps, formula = PITTAG ~ TRANID ~ Date2, fun.aggregate = length, value.var = "EFFORTID")[-97,,]
  ## Because added negative data, need to remove row with PITTAG == NA
  
  prepdat <- list(act = act, snks = snks, scent = scent)
  
  return(prepdat)

}

