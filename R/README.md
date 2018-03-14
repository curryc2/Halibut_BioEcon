# Halibut Bioconomic Model Functions
### Contains functions for used in simulation model.

***

Script Name                                | Description
-------------------------------------------|-----------------------------------------
fisheryFootPrint_plus.R                    | Series of functions to calculate selectivity, maturation schedules, weight@age, length@age, ect. Original versions of these functions may be found on Steve's [github](https://github.com/seastateinc/fisheryFootprint).
Halibut_Plot_Fxns.R                        | Series of functions to plot different biological and fishery parameters by sex and age.
beverton-holt-recruit.R                    | Beverton-Holt recruitment function parameterized with steepness.
ricker-recruit.R                           | Ricker recruitment function parameterizied with steepness.
C-to-F.R                                   | Finds fishing mortality rate that would produce a specified catch for each fishing sector, given abundance and selectivity.
F-to-C.R                                   | Calculates catch by fishing sector for a specified fishing mortality rate, given abundance and selectivity.
read-update-params.R                       | Function that reads in parameter values from Halibut Model Inputs.xlsx, and updates life history and selectivity parameters in the halibut object. 
extract-params.R                           | Function to extract and rename parameters from the Halibut object for use in the bio-fishery simulation framework.
create-sim-objects.R                       | Function to create data objects for simulation.              
calc-init-age-prop.R                       | Function to calculate initial sex-specific proportions of biomass at age, based on unfished equilibrium
calc-SBpR.R                                | Function to calculate spawning biomass per recruit for harvest control rule calculations.

***
### Harvest Control Rule Options
Functions for alternative harvest control rules used for setting fishing mortality rate **F** or **biomass/abundance** level.

* Contained in [R/HCR](https://github.com/curryc2/Halibut_BioEcon/tree/master/R/HCR)

Script Name                                | Control Rule
-------------------------------------------|-----------------------------------------
HCR-linear.R                               | Sloping linear control rule with specified floor, ceiling, and transition zone.
HCR-threshold.R                            | Threshold harvest control rule, with specified floor and ceiling.

***
### Observation Error Functions
Functions for adding observation error to observed quantities.

* Contained in [R/Obs_Error](https://github.com/curryc2/Halibut_BioEcon/tree/master/R/Obs_Error)

Credit:
=============
### Dr. Steve Martell
Many of the functions contained in fisheryFootPrint_plus.R used were originally written by Steve Martell.
The original code may be found on [github](https://github.com/seastateinc/fisheryFootprint).
