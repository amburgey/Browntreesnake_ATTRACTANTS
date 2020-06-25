### SCR analysis for VIS lure vs no lure project
### Data was collected in 2015 where biologists walked transects in CP that either had traps with mice (lure) or no traps (no lure)
### The goal is to understand if detection probability is different between these two scenarios
### In the case of an incipient or suppressed (low density) population, any way that can maximize detection probability of brown treesnakes can be a useful resource for surveyors


rm(list = ls())

library(dplyr)
library(reshape2)
library(tidyverse)

source("PrepData.R")


# Read in capture data and survey data
caps <- read.csv("Captures.csv")[,c("EFFORTID","Date","PITTAG","SVL","TOTAL","WEIGHT","SEX","TRANSECT","LOCATION")]
survs <- read.csv("Surveys.csv")[,c("EFFORTID","Date","TRANSECT","TYPE")]

# Format for analysis
dat <- PrepDat(caps,survs)


