# COVID-19-ON-model

**Note:** This repository is currently under construction.

Simulation and model-fitting code for manuscript by Fair, Karatayev, Anand, and Bauch (2021). 

All data needed to run scripts are included in the "InputFiles" folder, with the exception of the parameter set inputs for covidHierV46_2Github_basesim.R, as the large number of files makes this infeasible. However, such parameter sets can be generated using covidHierV46_2Github_parmfit.R. Expected runtime varies between scripts, but an individual simulation should run within 1 minute.

**Scripts**

* covidHierV46_2Github_parmfit.R runs fitting procedure for obtaining parameter sets, outputs files containing the best parameter set as determined by that run of the fitting algorithms.

* covidHierV46_2Github_basesim.R evaluates output from covidHierV46_2Github_parmfit.R to identify parameter sets resulting in the best fit to data, runs simulations using those parameter sets, and outputs files containing the results of that identification and simulation.

* covidHierV46_2Github_firstcounterfactual.R runs simulations for counterfactuals covering the first time period (10 March 2020 to 15 August 2020).

* covidHierV46_2Github_secondcounterfactual.R runs simulations for counterfactuals covering the second time period (12 June 2020 to 28 February 2021).

* covidHierV46_2Github_indexp.R runs simulations for our model extension to consider individual-level NPI adherence during first time period (10 March 2020 to 15 August 2020). 


**Input Files**

* simsFITPARMS_v46v2phuNUDGE_V2.rds contains best-fit parameter sets (generated from covidHierV46_2Github_parmfit.R and covidHierV46_2Github_basesim.R), used for all counterfactual simulations

* cmat_data_weighted.rds contains contact matrices; data adapted from https://github.com/kieshaprem/synthetic-contact-matrices

* PHUtoREGIONlinker_numeric.rds contains a key for linking Ontario's census divisions (regions) to the Public Health Units (PHU) they belong to.

* Ondat_newKage_to2021-03-02.rds contains data on new reported COVID-19 cases in Ontario, grouped by age-class (data from https://covid-19.ontario.ca/data).

* ONdat_casesbyPHU_to2021-03-01.rs contains data on new reported COVID-19 cases in Ontario, grouped by PHU (data from https://covid-19.ontario.ca/data).

* mobilitydat-REAL_to2021-02-27_v4.rds contains data on mobility trends for retail and recreational locations (adapted from https://www.google.com/covid19/mobility/).

* covidHierData.rds contains the following objects (note that some are artefacts of a previous version of the model retained within the file to preseve the dimensions of certain data structures - these are denoted with an asterisk):

  * Distd50*: the proportion of people voluntarily distancing by day 21 of the outbreak (data from https://news.gallup.com/opinion/gallup/298310/americans-step-social-distancing-even-further.aspx)

  * testRat: estimate of the ratio of true infections to positive tests (from https://doi.org/10.1101/2020.03.14.20036178)

  * Msave: Travel matrix, where entries are numbers of commuters between a home and workplace county. Data from Statistics canada, 2016 census, catalogue no. 98-400-x2016391.

  * storeSatIncTrm.65*: Pre-calculated values of baseline transmission probability (across values of c and xi) which leads to 65% of population being infected after 1 year without mitigation

  * caseCtBin*: the daily number of reported cases, with smaller counties grouped by population density (data from https://www.publichealthontario.ca/en/data-and-analysis/infectious-disease/covid-19-data-surveillance/covid-19-data-tool)

  * testT3.1*: function giving the (smoothed) testing intensity on any given day t from the day the 50th case positive case was reported (i.e., t=0 if fewer than 50 cumulative cases reported). These values increase from 0 to 1following the (smoothed) relative increase in number of tests taken each day based on data from https://www.ontario.ca/page/2019-novel-coronavirussection-0

  * CountyPopDensities*: population densities in each county (data from https://www12.statcan.gc.ca/health-sante/82-228/details/page.cfm?Lang=E&Tab=1&Geo1=HR&Code1=3570&Geo2=PR&Code2=35&Data=Rate&SearchText=York%20Regional%20Health%20Unit&SearchType=Contains&SearchPR=01&B1=All&Custom=&B2=All&B3=All)

**System & hardware requirements**

 * Windows 10 Pro Version 2004

 * Sufficient RAM to support storage of data structures during simulations

**Installation guides**

All software should install within a few minutes on a standard computer, the versions listed here are those the scripts have been tested on.

 * R Version 4.0.3 https://www.r-project.org/

 * R Studio Version 1.2.5019 (IDE for R) https://rstudio.com/ 
