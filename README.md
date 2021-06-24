In the case of an incipient (aka, low density) population of brown treesnakes (*Boiga irregularis*), any way that we can maximize the encounter rate and thus detection probability of brown treesnakes across the landscape can be a useful resource for surveyors. These two projects were our attempts to understand if hotspots of increased detection probability could be created using attractants to facilitate a rapid response to a brown treesnake invasion. Both of these studies were conducted in the Closed Population (CP) on the island of Gu&aring;han.

## Table of Contents
Data folder - contains capture, survey, and catch per unit effort information for both studies in addition to scripts to prepare data for analysis

Scripts folder - contains scripts to run analyses for each study

Models folder - contains model text file for each study

Results folder - contains saved results from the analysis of each study

EncounterRateComparison - script to run post-analysis comparison

PredictiveSingleNight - script to calculate detection probabilities from encounter rates

FigPrep - script for preparing figures

Figures folder - output from FigPrep

## Required Packages and Versions Used
dplyr_1.4.4

reshape2_1.4.4

tidyverse_1.3.0

secr_4.3.0

jagsUI_1.5.1

data.table_1.13.0

abind_1.4.5

ggplot2_3.3.2

ggpubr_0.4.0

HDInterval_0.2.2

## Attractant One: visual survey (VIS) lure/no lure

Data were collected in 2015 where biologists walked transects in CP that either had no traps (no lure) or traps with mice (lure)

The goal was to understand if encounter rate and detection probability were different between these two scenarios

We use an indicator variable (STATUS) to denote when 1 = a transect was not surveyed, 2 = a transect was surveyed but no lures were used, and 3 = a transect was surveyed and lures were used 

## Attractant Two: visual survey (VIS) scent/no scent

Data were collected in 2016 where biologists walked transects in CP that either had no scent applied (no scent), were freshly sprayed earlier that evening with fish fertilizer (fresh scent), or had been sprayed 24-hours previously with that scent (old scent)

The goal was to understand if encounter rate and detection probability were different between these three scenarios

We use an indicator variable (STATUS) to denote when 0 = a transect was not surveyed, 2 = a transect was surveyed but no scent was sprayed, 3 = a transect was surveyed and fresh scent was applied, and 4 = a transect was surveyed and scent had been applied the day before

We fit two models in a Spatially Explicit Capture-Recapture framework. Details of this work can be found in the published journal article on this topic.

## Details of article TBD

## How to use this repository

Start in Scripts to run data prep and analyses; Files to calculate results are EncounterRateComparisons and PredictiveSingleNight; FigPrep creates figures of results
