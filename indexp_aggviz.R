### DATA VIZ
library(tidyverse); library(paletteer); library(khroma); library(colorBlindness); library(ggdark); library(Rfast); library(ggthemes); library(scales); library(Hmisc); library(ggnewscale);
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

`%notin%` <- Negate(`%in%`);

setwd("C:/Users/jogad/OneDrive - University of Waterloo/S20/CovidModelling/OntarioProj/COVID-19-Hierarchy-model-master/SharcnetSims/")

DATA=readRDS("covidHierData.rds"); #See readme file for description of contents
pops=colSums(DATA$Msave)

if (exists('rw.weekly')==FALSE) #Skips running this data collating step if already done (marked by existence of final object generated in this step)
{
  # DATA=readRDS("covidHierData.rds"); #See readme file for description of contents
  # pops=colSums(DATA$Msave)
  
  #Read in empircal data for comparison
  rw.ind.0<-read.csv("conposcovidloc_mar32021.csv", na.strings=c(""," ","NA"))
  rw.tots<-read.csv("covidtesting_mar32021.csv")
  
  #Drop all dates in "grey" period where there may be reporting lags, etc
  rw.ind<-rw.ind.0[rw.ind.0$Case_Reported_Date<="2021-02-28",]
  rw.tots<-rw.tots[rw.tots$Reported.Date<="2021-02-28",]
  
  #Fill in missing days for data
  rw.tots<- rw.tots %>%
    mutate(Date = as.Date(Reported.Date)) %>%
    complete(Date = seq.Date(min(Date), max(Date), by="day"))
  
  # Aggregate age class data
  agetagger=function(x) {
    if (is.na(x)==TRUE) {return(NA)}
    if(x == "<20") {return(1)}
    if(x %in% c("20s", "30s")) {return(2)}
    if (x %in% c("40s", "50s")) {return(3)}
    if (x %in% c("60s", "70s")) {return(4)}
    if (x %in% c("80s", "90s", "90+")) {return(5)}
    if (x == "UNKNOWN") {return("UNKNOWN")}}
  rw.ind$ageclass <- NA
  
  rw.ind$ageclass<-sapply(rw.ind$Age_Group, function(x) agetagger(x))
  
  ### Expand ind data to cover all possible dates and age classes
  rw.ind$Date<-as.Date(rw.ind$Case_Reported_Date)
  temp50<- rw.ind %>% group_by(ageclass, Date) %>% dplyr::summarize (newcases=n(), deaths=sum(Outcome1=="Fatal"))
  temp50 <- temp50[order(temp50$Date),]
  temp50$totalcases<-cumsum(temp50$newcases)
  df.dates<-data.frame(Date=(seq(as.Date("2020-01-01"), as.Date("2021-03-31"), by = "days")))
  df.50<-data.frame(Date=(seq(as.Date("2020-01-01"), as.Date("2021-03-31"), by = "days")), ts50=1:length(seq(as.Date("2020-01-01"), as.Date("2021-03-31"), by = "days")) - length(seq(as.Date("2020-01-01"),temp50$Date[min(which(temp50$totalcases>=50))], by="days")))
  rw.ind.dates<-merge(rw.ind, df.dates, by="Date", all=TRUE)
  template.full<- rw.ind.dates %>% tidyr::expand(ageclass, Date, Reporting_PHU)
  rw.ind.full<-merge(template.full, rw.ind.dates, by=c("ageclass", "Date", "Reporting_PHU"), all=T)
  rw.ind.full<-merge(rw.ind.full, df.50, by="Date")
  
  template.full.OB<- rw.ind.dates[rw.ind.dates$Outbreak_Related=="Yes",] %>% tidyr::expand(ageclass, Date, Reporting_PHU)
  rw.ind.OB<-merge(template.full.OB, rw.ind.dates[rw.ind.dates$Outbreak_Related=="Yes",], by=c("ageclass", "Date", "Reporting_PHU"), all=T)
  rw.ind.OB<-merge(rw.ind.OB, df.50, by="Date")
  
  ##Aggregate individual case reporting data - deaths not accurate, taken out
  rw.agg<- rw.ind.full %>% group_by(ageclass, Date,ts50) %>% dplyr::summarize (newcases=sum(!is.na(Case_Reported_Date))) #, deaths=sum(Outcome1=="Fatal", na.rm=T))
  rw.agg<-rw.agg[order(rw.agg$Date),]
  
  rw.agg.OB<- rw.ind.OB %>% group_by(ageclass, Date,ts50) %>% dplyr::summarize (newcases.OB=sum(!is.na(Case_Reported_Date))) #, deaths.OB=sum(Outcome1=="Fatal", na.rm=T))
  rw.agg.OB<-rw.agg.OB[order(rw.agg.OB$Date),]
  
  # rw.agg.region<- rw.ind.full %>% group_by(ageclass, Date, Reporting_PHU, ts50) %>% dplyr::summarize (newcases=sum(!is.na(Case_Reported_Date))) #,  deaths=sum(Outcome1=="Fatal", na.rm=T))
  # rw.agg.region<-rw.agg.region[order(rw.agg.region$Date),]
  
  rw.tots$new.K<-rw.tots$Total.Cases-lag(rw.tots$Total.Cases)
  rw.tots<-merge(rw.tots, df.50, by="Date", all.x=TRUE)
  
  #Create death summary across age classes based on ind data - NOTE THAT THE DATES DO NOT CORRESPOND TO DATE OF DEATH
  # rw.agg.sum <- rw.agg %>% group_by(Date, ts50) %>% dplyr::summarize(new.deaths.rw.ind=sum(deaths))
  #rw.agg.sum$total.deaths.rw.ind<-cumsum(rw.agg.sum$new.deaths.rw.ind)
  #rw.agg.region.sum <- rw.agg.region %>% group_by(Date, ts50, Reporting_PHU) %>% dplyr::summarize(new.deaths.rw.ind=sum(deaths), new.cases.rw.ind=sum(newcases))
  #rw.agg.region.sum$total.deaths.rw.ind<-cumsum(rw.agg.region.sum$new.deaths.rw.ind)
  
  #Write a subset of rw.agg for fitting
  rw.agg.full<-merge(rw.agg, df.dates, by="Date", all.x=TRUE, all.y=TRUE)
  rw.agg.fit<-rw.agg.full[order(rw.agg.full$ageclass), colnames(rw.agg.full) %notin% "deaths"] %>%
    tidyr::pivot_wider(names_from = ageclass, values_from = newcases,  names_prefix="newK.")
  rw.agg.fit<-rw.agg.fit[order(rw.agg.fit$ts50),]
  rw.agg.fit<-rw.agg.fit[rw.agg.fit$Date<=max(as.Date(rw.ind.0$Case_Reported_Date),na.rm=TRUE),] #Drop any dates from the fitting that go past the final day of data we have from PHO
  rw.agg.fit<-rw.agg.fit[rowSums(is.na(rw.agg.fit)) != ncol(rw.agg.fit),] #Get rid of rows that are all NAs
  rw.agg.fit$total.K<-cumsum(rowSums(rw.agg.fit[,3:9], na.rm=TRUE))
  rw.agg.fit$new.K<-rowSums(rw.agg.fit[,3:9], na.rm=TRUE)
  
  rw.agg.OB.full<-merge(rw.agg.OB, df.50, by="Date", all.x=TRUE, all.y=TRUE)
  colnames(rw.agg.OB.full)[colnames(rw.agg.OB.full)=="ts50.y"]<-"ts50"
  rw.agg.OB.fit<-rw.agg.OB.full[order(rw.agg.OB.full$ageclass), colnames(rw.agg.OB.full) %notin% "deaths.OB"] %>%
    tidyr::pivot_wider(names_from = ageclass, values_from = newcases.OB,  names_prefix="newK.OB.")
  rw.agg.OB.fit<-rw.agg.OB.fit[order(rw.agg.OB.fit$ts50),]
  rw.agg.OB.fit<-rw.agg.OB.fit[rw.agg.OB.fit$Date<=max(as.Date(rw.ind.0$Case_Reported_Date),na.rm=TRUE),] #Drop any dates from the fitting that go past the final day of data we have from PHO
  rw.agg.OB.fit<-rw.agg.OB.fit[rowSums(is.na(rw.agg.OB.fit)) != ncol(rw.agg.OB.fit),] #Get rid of rows that are all NAs
  rw.agg.OB.fit$total.K.OB<-cumsum(rowSums(rw.agg.OB.fit[,4:10], na.rm=TRUE))
  rw.agg.OB.fit$new.K.OB<-rowSums(rw.agg.OB.fit[,4:10], na.rm=TRUE)
  
  rw.agg.fit<-merge(rw.agg.fit, rw.agg.OB.fit[,c("ts50", "total.K.OB", "new.K.OB")], by="ts50")
  rw.agg<-merge(rw.agg, rw.agg.OB, by=c("ageclass", "Date", "ts50"), all.x=T, all.y=T)
  rw.agg$newcases.OB[is.na(rw.agg$newcases.OB)]<-0
  #rw.agg$deaths.OB[is.na(rw.agg$deaths.OB)]<-0
  rw.agg<-rw.agg[order(rw.agg$Date),]
  
  ### Updating with new death outcome data by age
  rw.hosp<-readRDS("RW_hospbyage_to2021-03-02.rds")
  
  #Should note down the citations for these - written in old notebook
  #Corrected the lags ot account for the fact that we're now not adding +5 for the 5 dats from exposure to onset
  #Instead we are subtracting the avg of 5 days from onset to case becoming known, for a total of -10 off previous values for lags
  # for "leaves" stay as-is because these are times to recover (presumably from time of admission?)
  hosplag<-2; #12; 
  hospleave<- 1/10; #1/15; #1/18;
  iculag<-4; #14; 
  iculeave<- 1/13.25; #1/13; #1/16;
  deathlag<-6; #16;
  
  rw.hosp.fit<-rw.tots[, c("ts50", "Deaths", "Number.of.patients.hospitalized.with.COVID.19", "Number.of.patients.in.ICU.with.COVID.19", "Number.of.patients.in.ICU.on.a.ventilator.with.COVID.19" )]
  rw.hosp.fit$newDeaths.rw.tot<-rw.hosp.fit$Deaths-lag(rw.hosp.fit$Deaths) #Assuming death counts reported in total are accurate to that day, new deaths should be a simple lag1
  rw.hosp.fit$newDeaths.rw.tot[min(which(rw.hosp.fit$Deaths>=0))]<-rw.hosp.fit$Deaths[min(which(rw.hosp.fit$Deaths>=0))]
  #rw.hosp.fit3<-merge(rw.hosp.fit,rw.agg.sum, by="ts50", all=T)
  #rw.hosp.fit3$newDeaths.rw.ind.lagged<-lag(rw.hosp.fit3$total.deaths.rw.ind,deathlag)-lag(lag(rw.hosp.fit3$total.deaths.rw.ind,deathlag))
  
  #Knock out dates in rw data we're unceratain of
  rw.agg<-rw.agg[rw.agg$Date<=as.Date("2021-02-28"),]
  rw.agg.fit<-rw.agg.fit[rw.agg.fit$Date<=as.Date("2021-02-28"),]
  #rw.agg.region<-rw.agg.region[rw.agg.region$Date<=as.Date("2021-01-01"),]
  
  ###### Rest of new CFR stuff
  
  rw.tots$week<-as.Date(cut(rw.tots$Date,
                            breaks = "week"))
  rw.tots$deaths[is.na(rw.tots$Deaths)]<-0
  rw.tots$newdeaths<-rw.tots$Deaths - lag(rw.tots$Deaths)
  rw.weekly <- rw.tots %>% group_by(week) %>% dplyr::summarise(newcases.weekly.actual=sum(new.K), newdeaths.weekly.actual=sum(newdeaths),
                                                               mean.icu=mean(Number.of.patients.in.ICU.with.COVID.19 + Number.of.patients.in.ICU.on.a.ventilator.with.COVID.19))
  rw.weekly$cases.l2.actual=lag(rw.weekly$newcases.weekly.actual,2)
  rw.weekly$weekly.cfr.actual<-rw.weekly$newdeaths.weekly.actual/rw.weekly$cases.l2.actual
  rw.weekly$binned.icu<-cut(rw.weekly$mean.icu, breaks=c(0,150,300,750))
  rw.weekly$wave<-cut(rw.weekly$week, breaks=c(as.Date("2020-01-20"), as.Date("2020-08-10"),as.Date("2021-02-22")), right=T)
  rw.weekly<-rw.weekly[is.na(rw.weekly$wave)==FALSE,]
  
  ######
  
}

# 
# path = "C:/Users/jogad/OneDrive - University of Waterloo/S20/CovidModelling/OntarioProj/COVID-19-Hierarchy-model-master/SharcnetSims/ncincr_experiment_V46"
# fitparms.file.names <- dir(path, pattern ="simsFITPARMS_suscRAND_ncincr_v46v2phuNUDGE_susceffect0.50_", full.names=TRUE)
# simprovince.file.names <- dir(path, pattern ="simsPROVINCE_suscRAND_ncincr_v46v2phuNUDGE_susceffect0.50_", full.names=TRUE)
# 
# if(exists('df.model.sum.comb')==FALSE || (exists('df.model.sum.comb')==TRUE & dim(df.model.sum.comb)[1]!=(length(fitparms.file.names)*10))) {
# 
#   for(fpick in 1:length(fitparms.file.names))
#   {
#     #Read in simulation data
#     #Select how many of the top fits you want
#     topnum<-10
# 
#     fitparms00=readRDS(fitparms.file.names[fpick]);
#     simdata00=readRDS(simprovince.file.names[fpick]);
# 
#     if(fitparms.file.names[fpick]=="C:/Users/jogad/OneDrive - University of Waterloo/S20/CovidModelling/OntarioProj/COVID-19-Hierarchy-model-master/SharcnetSims/ncincr_experiment_V46/simsFITPARMS_suscRAND_ncincr_v46v2phuNUDGE_susceffect0.50_0.99_1.00.rds") {
#       fitparms00=fitparms00[fitparms00$ncperc==1.00,]
#       simdata00=simdata00[simdata00$ncperc==1.00,]
#     }
# 
#     bestcutoff<-if(length(unique(simdata00$LL))>(topnum+1)) {Rfast::nth(unique(simdata00$LL), topnum+1, descending = F)} else {1e8}
#     fitparms0 <- fitparms00[fitparms00$LL<bestcutoff,]
#     simdata0 <- simdata00[simdata00$LL<bestcutoff,]
# 
#     #Get rid of any dup parm sets (ID by LL)
#     checker<-fitparms0[!duplicated(fitparms0[ , "LL" ] ), ]
#     fitparms<- fitparms0[fitparms0$fit %in% checker$fit,]
#     simdata<- simdata0[simdata0$fit %in% checker$fit,]
# 
#     # ### Tags
#     # #Define state space for compliant individuals
#     # Susc=c("S1", "S2", "S3", "S4", "S5"); E=c("E1", "E2","E3", "E4", "E5"); R=c("R1", "R2", "R3", "R4","R5");
#     # Da=c("A1","Ak1","SA1","SAk1", "A2","Ak2","SA2","SAk2", "A3","Ak3","SA3","SAk3", "A4","Ak4","SA4","SAk4", "A5","Ak5","SA5","SAk5")
#     # Di=c("I1","Ik1","SI1","SIk1", "I2","Ik2","SI2","SIk2", "I3","Ik3","SI3","SIk3", "I4","Ik4","SI4","SIk4", "I5","Ik5","SI5","SIk5")
#     # Tn=c("A1","SA1","I1","SI1","A2","SA2","I2","SI2","A3","SA3","I3","SI3","A4","SA4","I4","SI4","A5","SA5","I5","SI5")
#     # Tk=c("Ak1","SAk1","Ik1","SIk1","Ak2","SAk2","Ik2","SIk2","Ak3","SAk3","Ik3","SIk3","Ak4","SAk4","Ik4","SIk4","Ak5","SAk5","Ik5","SIk5")
#     # All1<-c("S1", "E1", "R1", "A1","Ak1","SA1","SAk1", "I1","Ik1","SI1","SIk1"); All2<-c("S2", "E2", "R2", "A2","Ak2","SA2","SAk2", "I2","Ik2","SI2","SIk2"); All3<-c("S3", "E3", "R3", "A3","Ak3","SA3","SAk3", "I3","Ik3","SI3","SIk3"); All4<-c("S4", "E4", "R4", "A4","Ak4","SA4","SAk4", "I4","Ik4","SI4","SIk4"); All5<-c("S5", "E5", "R5", "A5","Ak5","SA5","SAk5", "I5","Ik5","SI5","SIk5");
#     # sick1<-c("E1", "A1","Ak1","SA1","SAk1", "I1","Ik1","SI1","SIk1"); sick2<-c("E2", "A2","Ak2","SA2","SAk2", "I2","Ik2","SI2","SIk2"); sick3<-c("E3", "A3","Ak3","SA3","SAk3", "I3","Ik3","SI3","SIk3"); sick4<-c("E4", "A4","Ak4","SA4","SAk4", "I4","Ik4","SI4","SIk4"); sick5<-c("E5", "A5","Ak5","SA5","SAk5", "I5","Ik5","SI5","SIk5");
#     # symp1<-c("I1","Ik1","SI1","SIk1"); symp2<-c("I2","Ik2","SI2","SIk2"); symp3<-c("I3","Ik3","SI3","SIk3"); symp4<-c("I4","Ik4","SI4","SIk4"); symp5<-c("I5","Ik5","SI5","SIk5");
#     # #Define state space for non-compliant individuals
#     # Suscnc=c("S1nc", "S2nc", "S3nc", "S4nc", "S5nc"); Enc=c("E1nc", "E2nc","E3nc", "E4nc", "E5nc"); Rnc=c("R1nc", "R2nc", "R3nc", "R4nc","R5nc");
#     # Danc=c("A1nc","Ak1nc","SA1nc","SAk1nc", "A2nc","Ak2nc","SA2nc","SAk2nc", "A3nc","Ak3nc","SA3nc","SAk3nc", "A4nc","Ak4nc","SA4nc","SAk4nc", "A5nc","Ak5nc","SA5nc","SAk5nc")
#     # Dinc=c("I1nc","Ik1nc","SI1nc","SIk1nc", "I2nc","Ik2nc","SI2nc","SIk2nc", "I3nc","Ik3nc","SI3nc","SIk3nc", "I4nc","Ik4nc","SI4nc","SIk4nc", "I5nc","Ik5nc","SI5nc","SIk5nc")
#     # Tnnc=c("A1nc","SA1nc","I1nc","SI1nc","A2nc","SA2nc","I2nc","SI2nc","A3nc","SA3nc","I3nc","SI3nc","A4nc","SA4nc","I4nc","SI4nc","A5nc","SA5nc","I5nc","SI5nc")
#     # Tknc=c("Ak1nc","SAk1nc","Ik1nc","SIk1nc","Ak2nc","SAk2nc","Ik2nc","SIk2nc","Ak3nc","SAk3nc","Ik3nc","SIk3nc","Ak4nc","SAk4nc","Ik4nc","SIk4nc","Ak5nc","SAk5nc","Ik5nc","SIk5nc")
#     # All1nc<-c("S1nc", "E1nc", "R1nc", "A1nc","Ak1nc","SA1nc","SAk1nc", "I1nc","Ik1nc","SI1nc","SIk1nc"); All2nc<-c("S2nc", "E2nc", "R2nc", "A2nc","Ak2nc","SA2nc","SAk2nc", "I2nc","Ik2nc","SI2nc","SIk2nc"); All3nc<-c("S3nc", "E3nc", "R3nc", "A3nc","Ak3nc","SA3nc","SAk3nc", "I3nc","Ik3nc","SI3nc","SIk3nc"); All4nc<-c("S4nc", "E4nc", "R4nc", "A4nc","Ak4nc","SA4nc","SAk4nc", "I4nc","Ik4nc","SI4nc","SIk4nc"); All5nc<-c("S5nc", "E5nc", "R5nc", "A5nc","Ak5nc","SA5nc","SAk5nc", "I5nc","Ik5nc","SI5nc","SIk5nc");
#     # sick1nc<-c("E1nc", "A1nc","Ak1nc","SA1nc","SAk1nc", "I1nc","Ik1nc","SI1nc","SIk1nc"); sick2nc<-c("E2nc", "A2nc","Ak2nc","SA2nc","SAk2nc", "I2nc","Ik2nc","SI2nc","SIk2nc"); sick3nc<-c("E3nc", "A3nc","Ak3nc","SA3nc","SAk3nc", "I3nc","Ik3nc","SI3nc","SIk3nc"); sick4nc<-c("E4nc", "A4nc","Ak4nc","SA4nc","SAk4nc", "I4nc","Ik4nc","SI4nc","SIk4nc"); sick5nc<-c("E5nc", "A5nc","Ak5nc","SA5nc","SAk5nc", "I5nc","Ik5nc","SI5nc","SIk5nc");
#     # symp1nc<-c("I1nc","Ik1nc","SI1nc","SIk1nc"); symp2nc<-c("I2nc","Ik2nc","SI2nc","SIk2nc"); symp3nc<-c("I3nc","Ik3nc","SI3nc","SIk3nc"); symp4nc<-c("I4nc","Ik4nc","SI4nc","SIk4nc"); symp5nc<-c("I5nc","Ik5nc","SI5nc","SIk5nc");
#     # Snm=c(Susc,E,Da,Di,R,"Nt","Cw", "Cs", "C", "Nt1","Nt2","Nt3","Nt4","Nt5", "VD"); Snmnc=c(Suscnc,Enc,Danc,Dinc,Rnc,"Ntnc", "Nt1nc","Nt2nc","Nt3nc","Nt4nc","Nt5nc");
#     # Snm.comb=c(Snm,Snmnc, "probcc1", "probcnc1", "probncc1", "probncnc1", "probcc2", "probcnc2", "probncc2", "probncnc2", "probcc3", "probcnc3", "probncc3", "probncnc3", "probcc4", "probcnc4", "probncc4", "probncnc4", "probcc5", "probcnc5", "probncc5", "probncnc5");
# 
#     ### Calcualted probabilties of hospital outcomes using data from "COVID-19 Epidemeological Summary Severity Report"
#     hosp<-as.data.frame(rbind(c(4+11, 1, 0)/601, c(69+115, 11+21, 2+6)/5251, c(224+444, 65+148, 16+61)/6693, c(511+570, 178+140, 141+309)/4659, c(571+259, 60+9, 676+614)/4702))
#     colnames(hosp)<-c("hospital", "icu", "death.old")
#     #hosp$death<-rw.hosp$rw.deaths/rw.hosp$rw.known
#     df.cfr<-readRDS("ONdat_severeoutcomes_to2021-08-31.rds")
# 
#     ### create a data frame with lagged hospital outcomes
#     df.hosp.0<-merge(simdata, df.cfr, by="ts50", all.x=T)
#     #df.hosp.region.0<-merge(simregions, df.cfr, by="ts50", all.x=T)
# 
#     #Create columns for hospital outcomes to fill
#     df.hosp.0$cfr1[is.na(df.hosp.0$cfr1)] <-0;
#     df.hosp.0$cfr2[is.na(df.hosp.0$cfr2)] <-0;
#     df.hosp.0$cfr3[is.na(df.hosp.0$cfr3)] <-0;
#     df.hosp.0$cfr4[is.na(df.hosp.0$cfr4)] <-0;
#     df.hosp.0$cfr5[is.na(df.hosp.0$cfr5)] <-0;
# 
#     df.hosp.0$currenthospratio.on.1[is.na(df.hosp.0$currenthospratio.on.1)] <-0;
#     df.hosp.0$currenthospratio.on.2[is.na(df.hosp.0$currenthospratio.on.2)] <-0;
#     df.hosp.0$currenthospratio.on.3[is.na(df.hosp.0$currenthospratio.on.3)] <-0;
#     df.hosp.0$currenthospratio.on.4[is.na(df.hosp.0$currenthospratio.on.4)] <-0;
#     df.hosp.0$currenthospratio.on.5[is.na(df.hosp.0$currenthospratio.on.5)] <-0;
# 
#     df.hosp.0$currenticuratio.on.1[is.na(df.hosp.0$currenticuratio.on.1)] <-0;
#     df.hosp.0$currenticuratio.on.2[is.na(df.hosp.0$currenticuratio.on.2)] <-0;
#     df.hosp.0$currenticuratio.on.3[is.na(df.hosp.0$currenticuratio.on.3)] <-0;
#     df.hosp.0$currenticuratio.on.4[is.na(df.hosp.0$currenticuratio.on.4)] <-0;
#     df.hosp.0$currenticuratio.on.5[is.na(df.hosp.0$currenticuratio.on.5)] <-0;
# 
#     #Create columns for hospital outcomes to fill
#     df.hosp.0$hospK.1<-0; df.hosp.0$hospK.2<-0; df.hosp.0$hospK.3<-0; df.hosp.0$hospK.4<-0; df.hosp.0$hospK.5<-0;
#     df.hosp.0$icuK.1<-0; df.hosp.0$icuK.2<-0; df.hosp.0$icuK.3<-0; df.hosp.0$icuK.4<-0; df.hosp.0$icuK.5<-0;
#     df.hosp.0$new.deathK.1<-0; df.hosp.0$new.deathK.2<-0; df.hosp.0$new.deathK.3<-0; df.hosp.0$new.deathK.4<-0; df.hosp.0$new.deathK.5<-0;
# 
#     # df.hosp.region.0$hosp.1<-0; df.hosp.region.0$hosp.2<-0; df.hosp.region.0$hosp.3<-0; df.hosp.region.0$hosp.4<-0; df.hosp.region.0$hosp.5<-0;
#     # df.hosp.region.0$icu.1<-0; df.hosp.region.0$icu.2<-0; df.hosp.region.0$icu.3<-0; df.hosp.region.0$icu.4<-0; df.hosp.region.0$icu.5<-0;
#     # df.hosp.region.0$new.death.1<-0; df.hosp.region.0$new.death.2<-0; df.hosp.region.0$new.death.3<-0; df.hosp.region.0$new.death.4<-0; df.hosp.region.0$new.death.5<-0;
#     # df.hosp.region.0$hospK.1<-0; df.hosp.region.0$hospK.2<-0; df.hosp.region.0$hospK.3<-0; df.hosp.region.0$hospK.4<-0; df.hosp.region.0$hospK.5<-0;
#     # df.hosp.region.0$icuK.1<-0; df.hosp.region.0$icuK.2<-0; df.hosp.region.0$icuK.3<-0; df.hosp.region.0$icuK.4<-0; df.hosp.region.0$icuK.5<-0;
#     # df.hosp.region.0$new.deathK.1<-0; df.hosp.region.0$new.deathK.2<-0; df.hosp.region.0$new.deathK.3<-0; df.hosp.region.0$new.deathK.4<-0; df.hosp.region.0$new.deathK.5<-0;
# 
#     ### Hospital outcomes generated using binom probabiltiies
#     nsamp<-1 #number of samples to take for hospital outcomes
#     set.seed(1e3); #set seed for rbinom
#     # ### Selection individual regions to calculate outcomes for
#     # regionselect<-NULL
# 
#     addon<-0
# 
#     for (j in 1:length(unique(df.hosp.0$fit))) {
#       for (k in 1:length(unique(df.hosp.0$run))) {
#         for (m in 1:length(unique(df.hosp.0$ncperc))) {
#           for (n in 1:length(unique(df.hosp.0$susceffect))) {
# 
#             df.hosp<-df.hosp.0[(df.hosp.0$fit==unique(df.hosp.0$fit)[j] & df.hosp.0$run==unique(df.hosp.0$run)[k] & df.hosp.0$ncperc==unique(df.hosp.0$ncperc)[m] & df.hosp.0$susceffect==unique(df.hosp.0$susceffect)[n]), ]
# 
#             # tmp.1<-matrix(rep(df.hosp$hosp.1,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             # tmp.2<-matrix(rep(df.hosp$hosp.2,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             # tmp.3<-matrix(rep(df.hosp$hosp.3,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             # tmp.4<-matrix(rep(df.hosp$hosp.4,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             # tmp.5<-matrix(rep(df.hosp$hosp.5,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             #
#             # for (i in (13+addon):length(df.hosp$ts50))  {
#             #   tmp.1[i,] <- tmp.1[i-1,] + rbinom(nsamp, lag(df.hosp$new.sick.1,hosplag)[i], hosp[1,1]) - rbinom(nsamp, tmp.1[i-1,], hospleave)
#             #   tmp.2[i,] <- tmp.2[i-1,] + rbinom(nsamp, lag(df.hosp$new.sick.2,hosplag)[i], hosp[2,1]) - rbinom(nsamp, tmp.2[i-1,], hospleave)
#             #   tmp.3[i,] <- tmp.3[i-1,] + rbinom(nsamp, lag(df.hosp$new.sick.3,hosplag)[i], hosp[3,1]) - rbinom(nsamp, tmp.3[i-1,], hospleave)
#             #   tmp.4[i,] <- tmp.4[i-1,] + rbinom(nsamp, lag(df.hosp$new.sick.4,hosplag)[i], hosp[4,1]) - rbinom(nsamp, tmp.4[i-1,], hospleave)
#             #   tmp.5[i,] <- tmp.5[i-1,] + rbinom(nsamp, lag(df.hosp$new.sick.5,hosplag)[i], hosp[5,1]) - rbinom(nsamp, tmp.5[i-1,], hospleave)
#             # }
#             #
#             # df.hosp$hosp.1<-rowMeans(tmp.1)
#             # df.hosp$hosp.2<-rowMeans(tmp.2)
#             # df.hosp$hosp.3<-rowMeans(tmp.3)
#             # df.hosp$hosp.4<-rowMeans(tmp.4)
#             # df.hosp$hosp.5<-rowMeans(tmp.5)
#             #
#             tmp.1<-matrix(rep(df.hosp$hospK.1,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             tmp.2<-matrix(rep(df.hosp$hospK.2,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             tmp.3<-matrix(rep(df.hosp$hospK.3,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             tmp.4<-matrix(rep(df.hosp$hospK.4,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             tmp.5<-matrix(rep(df.hosp$hospK.5,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
# 
#             for (i in (13+addon):length(df.hosp$ts50))  {
# 
#               tmp.1[i,] <- tmp.1[i-1,] + rbinom(nsamp, lag(df.hosp$new.K.1.c + df.hosp$new.K.1.nc,hosplag)[i], df.hosp$currenthospratio.on.1[i]) - rbinom(nsamp, tmp.1[i-1,], hospleave)
#               tmp.2[i,] <- tmp.2[i-1,] + rbinom(nsamp, lag(df.hosp$new.K.2.c + df.hosp$new.K.2.nc,hosplag)[i], df.hosp$currenthospratio.on.2[i]) - rbinom(nsamp, tmp.2[i-1,], hospleave)
#               tmp.3[i,] <- tmp.3[i-1,] + rbinom(nsamp, lag(df.hosp$new.K.3.c + df.hosp$new.K.3.nc,hosplag)[i], df.hosp$currenthospratio.on.3[i]) - rbinom(nsamp, tmp.3[i-1,], hospleave)
#               tmp.4[i,] <- tmp.4[i-1,] + rbinom(nsamp, lag(df.hosp$new.K.4.c + df.hosp$new.K.4.nc,hosplag)[i], df.hosp$currenthospratio.on.4[i]) - rbinom(nsamp, tmp.4[i-1,], hospleave)
#               tmp.5[i,] <- tmp.5[i-1,] + rbinom(nsamp, lag(df.hosp$new.K.5.c + df.hosp$new.K.5.nc,hosplag)[i], df.hosp$currenthospratio.on.5[i]) - rbinom(nsamp, tmp.5[i-1,], hospleave)
#             }
# 
#             df.hosp$hospK.1<-rowMeans(tmp.1)
#             df.hosp$hospK.2<-rowMeans(tmp.2)
#             df.hosp$hospK.3<-rowMeans(tmp.3)
#             df.hosp$hospK.4<-rowMeans(tmp.4)
#             df.hosp$hospK.5<-rowMeans(tmp.5)
# 
#             # tmp.1<-matrix(rep(df.hosp$icu.1,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             # tmp.2<-matrix(rep(df.hosp$icu.2,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             # tmp.3<-matrix(rep(df.hosp$icu.3,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             # tmp.4<-matrix(rep(df.hosp$icu.4,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             # tmp.5<-matrix(rep(df.hosp$icu.5,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             #
#             # for (i in (15+addon):length(df.hosp$ts50))  {
#             #
#             #   tmp.1[i,] <- tmp.1[i-1,] + rbinom(nsamp, lag(df.hosp$new.sick.1,iculag)[i], df.hosp$currenticuratio.on.1[i]) - rbinom(nsamp, tmp.1[i-1,], iculeave)
#             #   tmp.2[i,] <- tmp.2[i-1,] + rbinom(nsamp, lag(df.hosp$new.sick.2,iculag)[i], df.hosp$currenticuratio.on.2[i]) - rbinom(nsamp, tmp.2[i-1,], iculeave)
#             #   tmp.3[i,] <- tmp.3[i-1,] + rbinom(nsamp, lag(df.hosp$new.sick.3,iculag)[i], df.hosp$currenticuratio.on.3[i]) - rbinom(nsamp, tmp.3[i-1,], iculeave)
#             #   tmp.4[i,] <- tmp.4[i-1,] + rbinom(nsamp, lag(df.hosp$new.sick.4,iculag)[i], df.hosp$currenticuratio.on.4[i]) - rbinom(nsamp, tmp.4[i-1,], iculeave)
#             #   tmp.5[i,] <- tmp.5[i-1,] + rbinom(nsamp, lag(df.hosp$new.sick.5,iculag)[i], df.hosp$currenticuratio.on.5[i]) - rbinom(nsamp, tmp.5[i-1,], iculeave)
#             #
#             #   # df.hosp$icu.1[i] <- df.hosp$icu.1[i-1] + rbinom(1, lag(df.hosp$new.sick.1,iculag)[i], hosp[1,2]) - rbinom(1, df.hosp$icu.1[i-1], iculeave)
#             #   # df.hosp$icu.2[i] <- df.hosp$icu.2[i-1] + rbinom(1, lag(df.hosp$new.sick.2,iculag)[i], hosp[2,2]) - rbinom(1, df.hosp$icu.2[i-1], iculeave)
#             #   # df.hosp$icu.3[i] <- df.hosp$icu.3[i-1] + rbinom(1, lag(df.hosp$new.sick.3,iculag)[i], hosp[3,2]) - rbinom(1, df.hosp$icu.3[i-1], iculeave)
#             #   # df.hosp$icu.4[i] <- df.hosp$icu.4[i-1] + rbinom(1, lag(df.hosp$new.sick.4,iculag)[i], hosp[4,2]) - rbinom(1, df.hosp$icu.4[i-1], iculeave)
#             #   # df.hosp$icu.5[i] <- df.hosp$icu.5[i-1] + rbinom(1, lag(df.hosp$new.sick.5,iculag)[i], hosp[5,2]) - rbinom(1, df.hosp$icu.5[i-1], iculeave)
#             # }
#             #
#             # df.hosp$icu.1<-rowMeans(tmp.1)
#             # df.hosp$icu.2<-rowMeans(tmp.2)
#             # df.hosp$icu.3<-rowMeans(tmp.3)
#             # df.hosp$icu.4<-rowMeans(tmp.4)
#             # df.hosp$icu.5<-rowMeans(tmp.5)
# 
#             tmp.1<-matrix(rep(df.hosp$icuK.1,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             tmp.2<-matrix(rep(df.hosp$icuK.2,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             tmp.3<-matrix(rep(df.hosp$icuK.3,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             tmp.4<-matrix(rep(df.hosp$icuK.4,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             tmp.5<-matrix(rep(df.hosp$icuK.5,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
# 
#             for (i in (15+addon):length(df.hosp$ts50))  {
# 
#               tmp.1[i,] <- tmp.1[i-1,] + rbinom(nsamp, lag(df.hosp$new.K.1.c + df.hosp$new.K.1.nc,iculag)[i], df.hosp$currenticuratio.on.1[i]) - rbinom(nsamp, tmp.1[i-1,], iculeave)
#               tmp.2[i,] <- tmp.2[i-1,] + rbinom(nsamp, lag(df.hosp$new.K.2.c + df.hosp$new.K.2.nc,iculag)[i], df.hosp$currenticuratio.on.2[i]) - rbinom(nsamp, tmp.2[i-1,], iculeave)
#               tmp.3[i,] <- tmp.3[i-1,] + rbinom(nsamp, lag(df.hosp$new.K.3.c + df.hosp$new.K.3.nc,iculag)[i], df.hosp$currenticuratio.on.3[i]) - rbinom(nsamp, tmp.3[i-1,], iculeave)
#               tmp.4[i,] <- tmp.4[i-1,] + rbinom(nsamp, lag(df.hosp$new.K.4.c + df.hosp$new.K.4.nc,iculag)[i], df.hosp$currenticuratio.on.4[i]) - rbinom(nsamp, tmp.4[i-1,], iculeave)
#               tmp.5[i,] <- tmp.5[i-1,] + rbinom(nsamp, lag(df.hosp$new.K.5.c + df.hosp$new.K.5.nc,iculag)[i], df.hosp$currenticuratio.on.5[i]) - rbinom(nsamp, tmp.5[i-1,], iculeave)
#             }
# 
#             df.hosp$icuK.1<-rowMeans(tmp.1)
#             df.hosp$icuK.2<-rowMeans(tmp.2)
#             df.hosp$icuK.3<-rowMeans(tmp.3)
#             df.hosp$icuK.4<-rowMeans(tmp.4)
#             df.hosp$icuK.5<-rowMeans(tmp.5)
# 
#             # tmp.1<-matrix(rep(df.hosp$new.death.1,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             # tmp.2<-matrix(rep(df.hosp$new.death.2,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             # tmp.3<-matrix(rep(df.hosp$new.death.3,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             # tmp.4<-matrix(rep(df.hosp$new.death.4,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             # tmp.5<-matrix(rep(df.hosp$new.death.5,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             #
#             # for (i in (17+addon):length(df.hosp$ts50))  {
#             #
#             #   tmp.1[i,] <- rbinom(nsamp, lag(df.hosp$new.sick.1,deathlag)[i], df.hosp$cfr1[i])
#             #   tmp.2[i,] <- rbinom(nsamp, lag(df.hosp$new.sick.2,deathlag)[i], df.hosp$cfr2[i])
#             #   tmp.3[i,] <- rbinom(nsamp, lag(df.hosp$new.sick.3,deathlag)[i], df.hosp$cfr3[i])
#             #   tmp.4[i,] <- rbinom(nsamp, lag(df.hosp$new.sick.4,deathlag)[i], df.hosp$cfr4[i])
#             #   tmp.5[i,] <- rbinom(nsamp, lag(df.hosp$new.sick.5,deathlag)[i], df.hosp$cfr5[i])
#             # }
#             #
#             # df.hosp$new.death.1<-rowMeans(tmp.1)
#             # df.hosp$new.death.2<-rowMeans(tmp.2)
#             # df.hosp$new.death.3<-rowMeans(tmp.3)
#             # df.hosp$new.death.4<-rowMeans(tmp.4)
#             # df.hosp$new.death.5<-rowMeans(tmp.5)
# 
#             tmp.1<-matrix(rep(df.hosp$new.deathK.1,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             tmp.2<-matrix(rep(df.hosp$new.deathK.2,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             tmp.3<-matrix(rep(df.hosp$new.deathK.3,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             tmp.4<-matrix(rep(df.hosp$new.deathK.4,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
#             tmp.5<-matrix(rep(df.hosp$new.deathK.5,nsamp), nrow=nrow(df.hosp), ncol=nsamp, byrow=F)
# 
#             for (i in (17+addon):length(df.hosp$ts50))  {
# 
#               tmp.1[i,] <- rbinom(nsamp, lag(df.hosp$new.K.1.c + df.hosp$new.K.1.nc,deathlag)[i], lag(df.hosp$cfr1)[i])
#               tmp.2[i,] <- rbinom(nsamp, lag(df.hosp$new.K.2.c + df.hosp$new.K.2.nc,deathlag)[i], lag(df.hosp$cfr2)[i])
#               tmp.3[i,] <- rbinom(nsamp, lag(df.hosp$new.K.3.c + df.hosp$new.K.3.nc,deathlag)[i], lag(df.hosp$cfr3)[i])
#               tmp.4[i,] <- rbinom(nsamp, lag(df.hosp$new.K.4.c + df.hosp$new.K.4.nc,deathlag)[i], lag(df.hosp$cfr4)[i])
#               tmp.5[i,] <- rbinom(nsamp, lag(df.hosp$new.K.5.c + df.hosp$new.K.5.nc,deathlag)[i], lag(df.hosp$cfr5)[i])
#             }
# 
#             df.hosp$new.deathK.1<-rowMeans(tmp.1)
#             df.hosp$new.deathK.2<-rowMeans(tmp.2)
#             df.hosp$new.deathK.3<-rowMeans(tmp.3)
#             df.hosp$new.deathK.4<-rowMeans(tmp.4)
#             df.hosp$new.deathK.5<-rowMeans(tmp.5)
# 
#             if (j==1 & k==1 && m==1 && n==1) {
#               df.hosp.new<-df.hosp
#             } else {
#               df.hosp.new<-rbind(df.hosp.new, df.hosp)
#             }
#           }
#         }
#       }
#       #print(j)
#     }
# 
#     # Pivot data for easier hospital outcome plotting
#     df.hosp.new.1 <- df.hosp.new[, colnames(df.hosp.new) %in% c("fit", "run", "ncperc", "susceffect", "ts50", "new.sick.1.c", "new.sick.2.c", "new.sick.3.c", "new.sick.4.c", "new.sick.5.c",
#                                                                 "new.sick.1.nc", "new.sick.2.nc", "new.sick.3.nc", "new.sick.4.nc", "new.sick.5.nc")] %>%
#       tidyr::pivot_longer(cols = starts_with("new.sick"), names_to = "ageclass", values_to = "value",  names_prefix="new.sick.")
#     df.hosp.new.1$type<-"new.sick"
# 
#     # df.hosp.new.2 <- df.hosp.new[, colnames(df.hosp.new) %in% c("fit", "run", "ncperc", "susceffect", "ts50", "hosp.1", "hosp.2", "hosp.3", "hosp.4", "hosp.5")] %>%
#     #   tidyr::pivot_longer(cols = starts_with("hosp"), names_to = "ageclass", values_to = "value", names_prefix="hosp.")
#     # df.hosp.new.2$type<-"hosp"
#     #
#     # df.hosp.new.3 <- df.hosp.new[, colnames(df.hosp.new) %in% c("fit", "run", "ncperc", "susceffect", "ts50", "icu.1", "icu.2", "icu.3", "icu.4", "icu.5")] %>%
#     #   tidyr::pivot_longer(cols = starts_with("icu"), names_to = "ageclass", values_to = "value", names_prefix="icu.")
#     # df.hosp.new.3$type<-"icu"
# 
#     # df.hosp.new.4 <- df.hosp.new[, colnames(df.hosp.new) %in% c("fit", "run", "ncperc", "susceffect", "ts50", "new.death.1", "new.death.2", "new.death.3", "new.death.4", "new.death.5")] %>%
#     #   tidyr::pivot_longer(cols = starts_with("new.death"), names_to = "ageclass", values_to = "value", names_prefix="new.death.")
#     # df.hosp.new.4$type<-"new.death"
# 
#     df.hosp.new.5 <- df.hosp.new[, colnames(df.hosp.new) %in% c("fit", "run", "ncperc", "susceffect", "ts50", "hospK.1", "hospK.2", "hospK.3", "hospK.4", "hospK.5")] %>%
#       tidyr::pivot_longer(cols = starts_with("hospK"), names_to = "ageclass", values_to = "value", names_prefix="hospK.")
#     df.hosp.new.5$type<-"hospK"
# 
#     df.hosp.new.6 <- df.hosp.new[, colnames(df.hosp.new) %in% c("fit", "run", "ncperc", "susceffect", "ts50", "icuK.1", "icuK.2", "icuK.3", "icuK.4", "icuK.5")] %>%
#       tidyr::pivot_longer(cols = starts_with("icuK"), names_to = "ageclass", values_to = "value", names_prefix="icuK.")
#     df.hosp.new.6$type<-"icuK"
# 
#     df.hosp.new.7 <- df.hosp.new[, colnames(df.hosp.new) %in% c("fit", "run", "ncperc", "susceffect", "ts50", "new.deathK.1", "new.deathK.2", "new.deathK.3", "new.deathK.4", "new.deathK.5")] %>%
#       tidyr::pivot_longer(cols = starts_with("new.deathK"), names_to = "ageclass", values_to = "value", names_prefix="new.deathK.")
#     df.hosp.new.7$type<-"new.deathK"
# 
#     df.hosp.new.8 <- df.hosp.new[, colnames(df.hosp.new) %in% c("fit", "run", "ncperc", "susceffect", "ts50", "new.K.1.c", "new.K.2.c", "new.K.3.c", "new.K.4.c", "new.K.5.c",
#                                                                 "new.K.1.nc", "new.K.2.nc", "new.K.3.nc", "new.K.4.nc", "new.K.5.nc")] %>%
#       tidyr::pivot_longer(cols = starts_with("new.K"), names_to = "ageclass", values_to = "value",  names_prefix="new.K.")
#     df.hosp.new.8$type<-"new.K"
# 
#     # df.sum.hosp<-rbind(df.hosp.new.1, df.hosp.new.2, df.hosp.new.3, df.hosp.new.4, df.hosp.new.5, df.hosp.new.6, df.hosp.new.7, df.hosp.new.8)
#     df.sum.hosp<-rbind(df.hosp.new.1, df.hosp.new.5, df.hosp.new.6, df.hosp.new.7, df.hosp.new.8)
#     rm(list = ls()[grepl("df.hosp.new.", ls())])
# 
#     ### Add some calendar dates for ease of reading
#     df.dates.big<-data.frame(Date=(seq(as.Date("2019-06-01"), as.Date("2021-01-31"), by = "days")), ts50=1:length(seq(as.Date("2019-06-01"), as.Date("2021-01-31"), by = "days")) - length(seq(as.Date("2019-06-01"),rw.agg$Date[min(which(cumsum(rw.agg$newcases)>=50))], by="days")))
#     simdata<-merge(simdata, df.dates.big, by="ts50")
#     df.sum.hosp<-merge(df.sum.hosp, df.dates.big, by="ts50")
#     df.hosp.new <- merge(df.hosp.new, df.dates.big, by="ts50")
#     df.hosp.new <- df.hosp.new %>% mutate(#new.death=new.death.1+new.death.2+new.death.3+new.death.4+new.death.5,
#       new.deathK=new.deathK.1+new.deathK.2+new.deathK.3+new.deathK.4+new.deathK.5,
#       # hosp=hosp.1+hosp.2+hosp.3+hosp.4+hosp.5,
#       hospK=hospK.1+hospK.2+hospK.3+hospK.4+hospK.5,
#       # icu=icu.1+icu.2+icu.3+icu.4+icu.5,
#       icuK=icuK.1+icuK.2+icuK.3+icuK.4+icuK.5
#     )
#     # model.hosp.fit<-df.hosp.new[, c("ts50", "fit", "run", "new.death", "new.deathK", "hosp", "hospK", "icu", "icuK", "ncperc", "susceffect")]
#     model.hosp.fit<-df.hosp.new[, c("ts50", "fit", "run", "new.deathK",  "hospK", "icuK", "ncperc", "susceffect")]
#     df.fit<-merge(model.hosp.fit, rw.hosp.fit, by="ts50")
#     df.fit.plotting<-merge(model.hosp.fit, rw.hosp.fit, by="ts50", all.x=T)
#     df.fit.plotting<-merge(df.fit.plotting[, colnames(df.fit.plotting)!="Date"], df.dates.big, by="ts50")
#     df.fit<- df.fit %>% drop_na(Deaths,Number.of.patients.hospitalized.with.COVID.19,Number.of.patients.in.ICU.with.COVID.19,Number.of.patients.in.ICU.on.a.ventilator.with.COVID.19)
#     df.fit$weights<-0
#     df.fit$weights[df.fit$Date>="2020-08-01"]<-1
#     df.fit$weights[df.fit$Date<"2020-08-01"]<-1
# 
#     if("Date" %notin% colnames(rw.hosp.fit)) {rw.hosp.fit<-merge(rw.hosp.fit, df.dates.big, by="ts50", all.x=T, all.y=F)}
# 
#     comp.limselect<-c(as.Date("2020-03-10"), as.Date("2020-08-15"))
#     limselect<-comp.limselect; x.iculabel<-as.Date("2020-03-10"); dbreak<-"1 month";
# 
#     df.sum.hosp$cleanageclass<-sapply(df.sum.hosp$ageclass, function(x) unlist(strsplit(x, split='[.]'))[1])
#     df.sum.hosp$cperc<-1-df.sum.hosp$ncperc
# 
#     df.sum.hosp<-df.sum.hosp[(df.sum.hosp$Date<=limselect[2] & df.sum.hosp$Date>=limselect[1]),]
# 
#     df.model.sum.1 <- df.sum.hosp %>% group_by(run,fit,ncperc,susceffect) %>% dplyr::summarise(#model.deaths=sum(value[type=="new.death"], na.rm=TRUE),
#       model.deathsK=sum(value[type=="new.deathK"], na.rm=TRUE),
#       model.sick=sum(value[type=="new.sick"], na.rm=TRUE),
#       model.known=sum(value[type=="new.K"], na.rm=TRUE))
# 
# 
#     df.model.sum.1$cperc<-1-df.model.sum.1$ncperc
# 
#     # df.model.sum.2 <- df.model.sum.1 %>% group_by(ncperc,susceffect) %>% dplyr::summarize(mean.model.deaths=mean(model.deaths), sd.model.deaths=sd(model.deaths), min.model.deaths=min(model.deaths), max.model.deaths=max(model.deaths),
#     #                                                                                       mean.model.deathsK=mean(model.deathsK), sd.model.deathsK=sd(model.deathsK), min.model.deathsK=min(model.deathsK), max.model.deathsK=max(model.deathsK),
#     #                                                                                       mean.model.known=mean(model.known) ,  sd.model.known=sd(model.known), min.model.known=min(model.known), max.model.known=max(model.known),
#     #                                                                                       mean.model.infected=mean(model.sick), sd.model.infected=sd(model.sick), min.model.infected=min(model.sick), max.model.infected=max(model.sick))
#     #
#     # df.model.sum.2$cperc<-1-df.model.sum.2$ncperc
# 
#     ########### new CFR stuff
#     df.parms<-readRDS("parmfits_cfrcase_fcnform1_v2.rds")
#     df.parms$nu<-df.parms$nu/1e2
#     fcnparms.w1<-df.parms[(df.parms$fcnform==1 & df.parms$lag==2 & df.parms$wave==1 & df.parms$SSE==min(df.parms$SSE[df.parms$fcnform==1 & df.parms$lag==2 & df.parms$wave==1])), c(1,2,3)]
#     fcnparms.w2<-df.parms[(df.parms$fcnform==1 & df.parms$lag==2 & df.parms$wave==2 & df.parms$SSE==min(df.parms$SSE[df.parms$fcnform==1 & df.parms$lag==2 & df.parms$wave==2])), c(1,2,3)]
#     fcnparms.w3<-df.parms[(df.parms$fcnform==1 & df.parms$lag==2 & df.parms$wave==3 & df.parms$SSE==min(df.parms$SSE[df.parms$fcnform==1 & df.parms$lag==2 & df.parms$wave==3])), c(1,2,3)]
# 
#     cfrfcn = function(cases2lag, week) {
#       parms=data.frame(week=week);
#       parms$v1<-ifelse(parms$week<="2020-05-25", fcnparms.w1[[1]], ifelse(parms$week<="2020-08-10",fcnparms.w2[[1]], fcnparms.w3[[1]]));
#       parms$v2<-ifelse(parms$week<="2020-05-25", fcnparms.w1[[2]], ifelse(parms$week<="2020-08-10",fcnparms.w2[[2]], fcnparms.w3[[2]]));
#       parms$v3<-ifelse(parms$week<="2020-05-25", fcnparms.w1[[3]], ifelse(parms$week<="2020-08-10",fcnparms.w2[[3]], fcnparms.w3[[3]]));
#       x=parms$v1 + parms$v2*(1-exp(-parms$v3*cases2lag));  return(x);}
# 
#     simdata$week<-as.Date(cut(simdata$Date,
#                               breaks = "week"))
#     simdeath <- simdata %>% group_by(run,fit,ncperc,susceffect,week) %>% dplyr::summarise(newcases.weekly=sum(new.K.c + new.K.nc))
#     simdeath <- simdeath %>% group_by(run,fit,ncperc,susceffect) %>% mutate(cases.l2=lag(newcases.weekly,2), weekly.deaths.worst=cfrfcn(cases.l2, week)*cases.l2) %>% ungroup()
#     simdeath<-merge(simdeath, rw.weekly, by="week")
#     simdeath$weekly.deaths.best<-simdeath$weekly.cfr.actual*simdeath$cases.l2
# 
#     ###########
# 
#     if (fpick==1) {
#       df.model.sum.comb<-df.model.sum.1
#       simdeath.comb<-simdeath
#     } else {
#       df.model.sum.comb<-rbind(df.model.sum.comb, df.model.sum.1)
#       simdeath.comb<-rbind(simdeath.comb, simdeath)
#     }
# 
#     if (any(round(df.sum.hosp$cperc,2) %in% round(seq(0,1,0.1),2))) {
#       print("saving near est.")
#       if (exists('df.sum.hosp.save')==FALSE) {
#         df.sum.hosp.save<-df.sum.hosp[round(df.sum.hosp$cperc,2) %in% round(seq(0,1,0.1),2),]
#       } else {
#         df.sum.hosp.save<-rbind(df.sum.hosp.save,df.sum.hosp[round(df.sum.hosp$cperc,2) %in% round(seq(0,1,0.1),2),])
#       }
# 
#       print(unique(df.sum.hosp.save$cperc))
#     }
# 
#     print(fpick/length(fitparms.file.names))
#   }
# 
# }
# 
# 
# #Clear from workspace all objects that aren't needed for plotting
# rm(list=setdiff(ls(), c("rw.agg","rw.agg.fit",  "fitparms", "pops", "%notin%","df.model.sum.comb",
#                         "hosplag", "hospleave", "iculag", "iculeave", "deathlag", "limselect",
#                         "df.sum.hosp.save", "rw.hosp.fit", "simdeath.comb", "rw.weekly", "rw.tots")))
# 
# df.sum.rw<-data.frame(ageclass=c(1,2,3,4,5, NA), rw.known=c(sum(rw.agg$newcases[(rw.agg$Date<=limselect[2] & rw.agg$Date>=limselect[1] & rw.agg$ageclass==1)], na.rm = TRUE),sum(rw.agg$newcases[(rw.agg$Date<=limselect[2] & rw.agg$Date>=limselect[1] & rw.agg$ageclass==2)], na.rm = TRUE),sum(rw.agg$newcases[(rw.agg$Date<=limselect[2] & rw.agg$Date>=limselect[1] & rw.agg$ageclass==3)], na.rm = TRUE),sum(rw.agg$newcases[(rw.agg$Date<=limselect[2] & rw.agg$Date>=limselect[1] & rw.agg$ageclass==4)], na.rm = TRUE),sum(rw.agg$newcases[(rw.agg$Date<=limselect[2] & rw.agg$Date>=limselect[1] & rw.agg$ageclass==5)], na.rm = TRUE), 0),
#                       rw.deaths=c(0,0,0,0,0,rw.hosp.fit$Deaths[(rw.hosp.fit$Date==limselect[2] & is.na(rw.hosp.fit$Date)==FALSE)]))
# 
# #df.totals.rw <- df.sum.rw %>% summarize(total.known = sum(rw.known), total.deaths=sum(rw.deaths))
# 
# saveRDS(df.sum.hosp.save, file = sprintf("dfIND_sum_hosp_save_to%s_v2.rds", limselect[2]))
# saveRDS(df.sum.rw, file = sprintf("dfIND_sum_rw_to%s_v2.rds", limselect[2]))
# saveRDS(df.model.sum.comb, file = sprintf("dfIND_model_sum_comb_to%s_v2.rds", limselect[2]))
# saveRDS(simdeath.comb, file = sprintf("dfIND_model_simdeath_to%s_v2.rds", limselect[2]))

## Viz
library(RColorBrewer); library(patchwork); library(Rfast);

rm(list=setdiff(ls(), c("rw.agg","rw.agg.fit",  "fitparms", "pops", "%notin%","df.model.sum.comb", "hosplag", "hospleave", "iculag", "iculeave", "deathlag", "limselect", "df.sum.hosp.save", "rw.hosp.fit", "rw.weekly", "rw.tots")))
limselect<-c(as.Date("2020-03-10"), as.Date("2020-08-15")); dbreak<-"1 month";

df.sum.hosp.save = readRDS("dfIND_sum_hosp_save_to2020-08-15_v2.rds")
df.sum.rw = readRDS("dfIND_sum_rw_to2020-08-15_v2.rds");
df.model.sum.comb = readRDS("dfIND_model_sum_comb_to2020-08-15_v2.rds");
simdeath.comb = readRDS("dfIND_model_simdeath_to2020-08-15_v2.rds");
simdeath.comb<-simdeath.comb[(simdeath.comb$week>=limselect[1] & simdeath.comb$week<=limselect[2]),]
simdeath.comb$cperc<-1-simdeath.comb$ncperc
simdeath.comb.sum0 <- simdeath.comb %>% group_by(run,fit,ncperc,susceffect,cperc) %>% dplyr::summarise(model.deathsK.worst=sum(weekly.deaths.worst, na.rm=TRUE),
                                                                                                       model.deathsK.best=sum(weekly.deaths.best, na.rm=TRUE))
simdeath.comb.sum <- simdeath.comb.sum0 %>% group_by(run,fit) %>% arrange(fit,run,cperc) %>% mutate(incr.deaths.best=(model.deathsK.best) - lag(model.deathsK.best),
                                                                                                    lead.incr.deaths.best=lead(model.deathsK.best) - (model.deathsK.best),
                                                                                                    red.deathsK.best=model.deathsK.best-model.deathsK.best[round(cperc,2)==0],
                                                                                                    incr.deaths.worst=(model.deathsK.worst) - lag(model.deathsK.worst),
                                                                                                    lead.incr.deaths.worst=lead(model.deathsK.worst) - (model.deathsK.worst),
                                                                                                    red.deathsK.worst=model.deathsK.worst-model.deathsK.worst[round(cperc,2)==0],
                                                                                                    incr.cperc=(cperc) - lag(cperc), lead.incr.cperc=lead(cperc),
                                                                                                    lead.incr.cperc=lead(cperc) - (cperc))


df.model.sum.comb<-df.model.sum.comb %>% group_by(fit, run) %>% arrange(fit, run,cperc) %>% mutate(incr.known=(model.known) - lag(model.known), incr.infected=(model.sick) - lag(model.sick),
                                                                                                   incr.deaths=(model.deathsK) - lag(model.deathsK),
                                                                                                   lead.incr.known=lead(model.known) - (model.known), lead.incr.infected=lead(model.sick) - (model.sick),
                                                                                                   lead.incr.deaths=lead(model.deathsK) - (model.deathsK),
                                                                                                   incr.cperc=(cperc) - lag(cperc), lead.incr.cperc=lead(cperc) - (cperc),
                                                                                                   red.known=model.known-model.known[round(cperc,2)==0],
                                                                                                   red.infected=model.sick-model.sick[round(cperc,2)==0],
                                                                                                   red.deathsK=model.deathsK-model.deathsK[round(cperc,2)==0])

rw.weekly<-rw.weekly[(rw.weekly$week>=limselect[1] & rw.weekly$week<=limselect[2]),]
# rw.weekly.sum<- rw.weekly %>% dplyr::summarise(rw.deaths=sum(newdeaths.weekly, na.rm=T))

# df.model.sum.comb.2<- df.model.sum.comb %>% group_by(fit, cperc) %>% dplyr::summarize(mean.deaths=mean(model.deaths), sd.deaths=sd(model.deaths), min.deaths=min(model.deaths), max.deaths=max(model.deaths),
#                                                                                 mean.deathsK=mean(model.deathsK), sd.deathsK=sd(model.deathsK), min.deathsK=min(model.deathsK), max.deathsK=max(model.deathsK),
#                                                                                 mean.known=mean(model.known),  sd.known=sd(model.known), min.known=min(model.known), max.known=max(model.known),
#                                                                                 mean.infected=mean(model.sick), sd.infected=sd(model.sick), min.infected=min(model.sick), max.infected=max(model.sick),
#                                                                                 mean.incr.infected=mean(incr.infected), sd.incr.infected=sd(incr.infected), min.incr.infected=min(incr.infected), max.incr.infected=max(incr.infected),
#                                                                                 mean.incr.known=mean(incr.known), sd.incr.known=sd(incr.known), min.incr.known=min(incr.known), max.incr.known=max(incr.known),
#                                                                                 mean.incr.deaths=mean(incr.deaths), sd.incr.deaths=sd(incr.deaths), min.incr.deaths=min(incr.deaths), max.incr.deaths=max(incr.deaths),
#                                                                                 mean.incr.cperc=mean(incr.cperc), sd.incr.cperc=sd(incr.cperc), min.incr.cperc=min(incr.cperc), max.incr.cperc=max(incr.cperc))



df.hosp.ts <- df.sum.hosp.save[(df.sum.hosp.save$Date>=limselect[1] & df.sum.hosp.save$Date<=limselect[2]),] %>% group_by(fit, run, cperc, ts50) %>% dplyr::summarise(#model.deaths=sum(value[type=="new.death"], na.rm=TRUE),
  model.hospK=sum(value[type=="hospK"], na.rm=TRUE),
  model.icuK=sum(value[type=="icuK"], na.rm=TRUE),
  model.deathsK=sum(value[type=="new.deathK"], na.rm=TRUE),
  model.sick=sum(value[type=="new.sick"], na.rm=TRUE),
  model.known=sum(value[type=="new.K"], na.rm=TRUE))

rm("df.sum.hosp.save"); #Get rid, due to large file size


df.dates.big<-data.frame(Date=(seq(as.Date("2019-06-01"), as.Date("2021-01-31"), by = "days")), ts50=1:length(seq(as.Date("2019-06-01"), as.Date("2021-01-31"), by = "days")) - length(seq(as.Date("2019-06-01"),rw.agg$Date[min(which(cumsum(rw.agg$newcases)>=50))], by="days")))
if("Date" %notin% colnames(df.hosp.ts)) {df.hosp.ts<-merge(df.hosp.ts, df.dates.big, by="ts50", all.x=T, all.y=F)}
if("Date" %notin% colnames(rw.hosp.fit)) {rw.hosp.fit<-merge(rw.hosp.fit, df.dates.big, by="ts50", all.x=T, all.y=F)}

#Knock out dates in rw data we're unceratain of
rw.agg<-rw.agg[rw.agg$Date<=as.Date("2021-01-01"),]
rw.hosp.fit<-rw.hosp.fit[rw.hosp.fit$Date<=as.Date("2021-01-01"),]

custom.pal<- c("black",(brewer.pal(4,"RdBu")))

p.known.ts<- ggplot() +
  geom_line(data=df.hosp.ts[(df.hosp.ts$cperc>=0.4 & df.hosp.ts$cperc<=0.7),], aes(x=Date, y=model.known, group=interaction(cperc, fit, run), colour=as.factor(100*(cperc))), alpha=0.1) +
  stat_summary( data=df.hosp.ts[(df.hosp.ts$cperc>=0.4 & df.hosp.ts$cperc<=0.7),], aes(x=Date, y=model.known, group=interaction(cperc), colour=as.factor(100*(cperc))), fun = median, geom="line", size=0.75) +
  scale_color_manual("Scenario", values=custom.pal[-1], guide = guide_legend(order = 2), labels=c("40% of individuals adhering to NPIs", "50% of individuals adhering to NPIs", "60% of individuals adhering to NPIs", "70% of individuals adhering to NPIs")) +
  new_scale_color() +
  #stat_summary( data=df.hosp.ts[df.hosp.ts$cperc>=0.5,], aes(x=Date, y=model.known, group=interaction(cperc), colour=as.factor(100*(1-cperc))), fun = mean, geom="line", linetype="dashed") +
  geom_line(data=rw.agg.fit, aes(x=Date, y=new.K, colour=as.factor("Actual") ), size=0.75) +
  scale_color_manual("", values=custom.pal[1], guide = guide_legend(order = 1)) +
  #geom_line(data=rw.tots, aes(x=Date, y=new.K), colour="black", linetype="dashed") +
  xlab("Date") +
  ylab("New confirmed cases") +
  # scale_x_date(expand = c(0,0), date_breaks = dbreak, limits=c(min(df.hosp.ts$Date, rw.agg.fit$Date[rw.agg.fit$new.K>=0]), limselect[2]), date_labels = "%b '%y") +
  scale_x_date(expand = c(0,0), date_breaks = dbreak, limits=c(limselect[1], limselect[2]), date_labels = "%b '%y") +
  #scale_y_continuous(expand = c(0.01,0), limits=c(0,max(df.hosp.ts$model.known[(df.hosp.ts$Date<=max(limselect) & df.hosp.ts$cperc>=0.6)], rw.agg.fit$new.K[rw.agg.fit$Date<=max(limselect)], na.rm=T))) +
  scale_y_continuous(labels=comma, expand=c(0.01,0)) +
  coord_cartesian(ylim=c(0,2500))+
  theme_light()

# p.hosp.ts<- ggplot() +
#   #geom_line(data=df.hosp.ts, aes(x=Date, y=model.deathsK, group=interaction(cperc, fit, run), colour=as.factor(1-cperc)), alpha=0.1) +
#   stat_summary( data=df.hosp.ts[df.hosp.ts$cperc>=0.5,], aes(x=Date, y=model.hospK, group=interaction(cperc), colour=as.factor(100*(1-cperc))), fun = median, geom="line") +
#   stat_summary( data=df.hosp.ts[df.hosp.ts$cperc>=0.5,], aes(x=Date, y=model.hospK, group=interaction(cperc), colour=as.factor(100*(1-cperc))), fun = mean, geom="line", linetype="dashed") +
#   geom_line(data=rw.hosp.fit, aes(x=Date, y=Number.of.patients.hospitalized.with.COVID.19), colour="black") +
#   xlab("Date") +
#   ylab("Individuals hospitalized") +
#   # scale_x_date(expand = c(0,0), date_breaks = dbreak, limits=c(min(df.hosp.ts$Date, rw.hosp.fit$Date[rw.hosp.fit$Number.of.patients.hospitalized.with.COVID.19>=0], na.rm=T), limselect[2]), date_labels = "%b '%y") +
#   scale_x_date(expand = c(0,0), date_breaks = dbreak, limits=c(limselect[1], limselect[2]), date_labels = "%b '%y") +
#   #scale_y_continuous(expand = c(0.01,0), limits=c(0,max(df.hosp.ts$model.deathsK[(df.hosp.ts$Date<=max(limselect) & df.hosp.ts$cperc>=0.6)], rw.hosp.fit$newDeaths.rw.tot[rw.hosp.fit$Date<=max(limselect)], na.rm=T))) +
#   scale_color_brewer("% of population not distancing", palette="RdBu") +
#   theme_light()
# 
# p.icu.ts<- ggplot() +
#   #geom_line(data=df.hosp.ts, aes(x=Date, y=model.deathsK, group=interaction(cperc, fit, run), colour=as.factor(1-cperc)), alpha=0.1) +
#   stat_summary( data=df.hosp.ts[df.hosp.ts$cperc>=0.5,], aes(x=Date, y=model.icuK, group=interaction(cperc), colour=as.factor(100*(1-cperc))), fun = median, geom="line") +
#   stat_summary( data=df.hosp.ts[df.hosp.ts$cperc>=0.5,], aes(x=Date, y=model.icuK, group=interaction(cperc), colour=as.factor(100*(1-cperc))), fun = mean, geom="line", linetype="dashed") +
#   geom_line(data=rw.hosp.fit, aes(x=Date, y=Number.of.patients.in.ICU.on.a.ventilator.with.COVID.19+Number.of.patients.in.ICU.with.COVID.19), colour="black") +
#   xlab("Date") +
#   ylab("Individuals in ICU") +
#   # scale_x_date(expand = c(0,0), date_breaks = dbreak, limits=c(min(df.hosp.ts$Date, rw.hosp.fit$Date[(rw.hosp.fit$Number.of.patients.in.ICU.on.a.ventilator.with.COVID.19 + rw.hosp.fit$Number.of.patients.in.ICU.with.COVID.19) >=0], na.rm=T), limselect[2]), date_labels = "%b '%y") +
#   scale_x_date(expand = c(0,0), date_breaks = dbreak, limits=c(limselect[1], limselect[2]), date_labels = "%b '%y") +
#   #scale_y_continuous(expand = c(0.01,0), limits=c(0,max(df.hosp.ts$model.deathsK[(df.hosp.ts$Date<=max(limselect) & df.hosp.ts$cperc>=0.6)], rw.hosp.fit$newDeaths.rw.tot[rw.hosp.fit$Date<=max(limselect)], na.rm=T))) +
#   scale_color_brewer("% of population not distancing", palette="RdBu") +
#   theme_light()

p.death.ts<- ggplot() +
  geom_line(data=df.hosp.ts[(df.hosp.ts$cperc>=0.4 & df.hosp.ts$cperc<=0.7),], aes(x=Date, y=model.deathsK, group=interaction(cperc, fit, run), colour=as.factor(100*(cperc))), alpha=0.1) +
  #stat_summary( data=df.hosp.ts[df.hosp.ts$cperc>=0.5,], aes(x=Date, y=model.deathsK, group=interaction(cperc), colour=as.factor(100*(1-cperc))), fun = median, geom="line") +
  stat_summary( data=df.hosp.ts[(df.hosp.ts$cperc>=0.4 & df.hosp.ts$cperc<=0.7),], aes(x=Date, y=model.deathsK, group=interaction(cperc), colour=as.factor(100*(cperc))), fun = median, geom="line", size=0.75) +
  scale_color_manual("Scenario", values=custom.pal[-1], guide = guide_legend(order = 2), labels=c("40% of individuals adhering to NPIs", "50% of individuals adhering to NPIs", "60% of individuals adhering to NPIs", "70% of individuals adhering to NPIs")) +
  new_scale_color() +
  geom_line(data=rw.hosp.fit, aes(x=Date, y=newDeaths.rw.tot, colour="Actual"), size=0.75) +
  scale_color_manual("", values=custom.pal[1], guide = guide_legend(order = 1)) +
  xlab("Date") +
  ylab("Deaths, CFR=f(t)") +
  #scale_x_date(expand = c(0,0), date_breaks = dbreak, limits=c(min(df.hosp.ts$Date, rw.hosp.fit$Date[rw.hosp.fit$newDeaths.rw.tot>=0], na.rm=T), limselect[2]), date_labels = "%b '%y") +
  scale_x_date(expand = c(0,0), date_breaks = dbreak, limits=c(limselect[1], limselect[2]), date_labels = "%b '%y") +
  #scale_y_continuous(expand = c(0.01,0), limits=c(0,max(df.hosp.ts$model.deathsK[(df.hosp.ts$Date<=max(limselect) & df.hosp.ts$cperc>=0.6)], rw.hosp.fit$newDeaths.rw.tot[rw.hosp.fit$Date<=max(limselect)], na.rm=T))) +
  scale_y_continuous(labels=comma, expand=c(0.01,0)) +
  #scale_color_manual("Individuals adhering to NPIs (% of population)", values=custom.pal) +
  theme_light()+
  ggtitle("Best-case scenario")

p.death.weekly.ts.best<-ggplot() +
  geom_line(data=simdeath.comb[(round(simdeath.comb$cperc,2) %in% c(0.4,0.5,0.6,0.7)),], aes(x=week, y=weekly.deaths.best, group=interaction(run,fit,cperc), colour=as.factor(100*(cperc))), alpha=0.1) +
  #stat_summary(data=simdeath.comb[round(simdeath.comb$cperc,2) %in% seq(0.5,1,0.1),], aes(x=week, y=weekly.deaths, group=interaction(cperc), colour=as.factor(100*(1-cperc))), fun=median,geom="line") +
  stat_summary(data=simdeath.comb[(round(simdeath.comb$cperc,2) %in% c(0.4,0.5,0.6,0.7)),], aes(x=week, y=weekly.deaths.best, group=interaction(cperc), colour=as.factor(100*(cperc))), fun=median,geom="line", size=0.75) +
  scale_color_manual("Scenario", values=custom.pal[-1], guide = guide_legend(order = 2), labels=c("40% of individuals adhering to NPIs", "50% of individuals adhering to NPIs", "60% of individuals adhering to NPIs", "70% of individuals adhering to NPIs")) +
  new_scale_color() +
  geom_line(data=rw.weekly, aes(x=week, y=newdeaths.weekly.actual, colour="Actual"), size=0.75) +
  scale_color_manual("", values=custom.pal[1], guide = guide_legend(order = 1)) +
  xlab("Date") +
  ylab("Weekly deaths") +
  scale_x_date(expand = c(0,0), limits=limselect, date_breaks = dbreak, date_labels = "%b '%y") +
  scale_y_continuous(labels=comma, expand=c(0.01,0), limits=c(0,3000)) +
  #scale_color_manual("Individuals adhering to NPIs (% of population)", values=custom.pal) +
  theme_light()+
  ggtitle("Best-case scenario")


p.death.weekly.ts.worst<-ggplot() +
  geom_line(data=simdeath.comb[(round(simdeath.comb$cperc,2) %in% c(0.4,0.5,0.6,0.7)),], aes(x=week, y=weekly.deaths.worst, group=interaction(run,fit,cperc), colour=as.factor(100*(cperc))), alpha=0.1) +
  #stat_summary(data=simdeath.comb[round(simdeath.comb$cperc,2) %in% seq(0.5,1,0.1),], aes(x=week, y=weekly.deaths, group=interaction(cperc), colour=as.factor(100*(1-cperc))), fun=median,geom="line") +
  stat_summary(data=simdeath.comb[(round(simdeath.comb$cperc,2) %in% c(0.4,0.5,0.6,0.7)),], aes(x=week, y=weekly.deaths.worst, group=interaction(cperc), colour=as.factor(100*(cperc))), fun=median,geom="line", size=0.75) +
  scale_color_manual("Scenario", values=custom.pal[-1], guide = guide_legend(order = 2), labels=c("40% of individuals adhering to NPIs", "50% of individuals adhering to NPIs", "60% of individuals adhering to NPIs", "70% of individuals adhering to NPIs")) +
  new_scale_color() +
  geom_line(data=rw.weekly, aes(x=week, y=newdeaths.weekly.actual, colour="Actual"), size=0.75) +
  scale_color_manual("", values=custom.pal[1], guide = guide_legend(order = 1)) +
  xlab("Date") +
  ylab("Weekly deaths, CFR=f(cases)") +
  scale_x_date(expand = c(0,0), limits=limselect, date_breaks = dbreak, date_labels = "%b '%y") +
  scale_y_continuous(labels=comma, expand=c(0.01,0), limits=c(0,3000)) +
  #scale_color_manual("Individuals adhering to NPIs (% of population)", values=custom.pal) +
  theme_light()+
  ggtitle("Worst-case scenario")




###THIS IS WHERE YOU GOT TO!!!###

psize<-0.025
lsize<-0.75

p.infected<- ggplot(df.model.sum.comb, aes(x=100*cperc, y=model.sick)) +
  geom_point(colour="#0072BB", alpha=0.1) +
  #stat_summary(aes(group=fit), fun =median, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
  stat_summary( fun = median, geom="line", colour="#0072BB", size=lsize) +
  #stat_summary( fun = mean, geom="line", colour="#0072BB", linetype="dashed", size=lsize) +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Total infections (by Aug 15, 2020)") +
  #scale_colour_brewer(palette="Set3") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0.01,0)) +
  theme_light()

p.known<- ggplot(df.model.sum.comb, aes(x=100*cperc, y=model.known)) +
  geom_point(colour="#0072BB", alpha=0.1) +
  #stat_summary(aes(group=fit), fun =median, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
  stat_summary( fun = median, geom="line", colour="black", size=lsize) +
  #stat_summary( fun = mean, geom="line", colour="black", linetype="dashed", size=lsize) +
  geom_hline(data=df.sum.rw, aes(yintercept=sum(rw.known)), colour="black", linetype="dotted", size=lsize) +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Total confirmed cases") +
  #scale_colour_brewer(palette="Set3") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0.01,0)) +
  theme_light()

p.death<- ggplot(df.model.sum.comb, aes(x=100*cperc, y=model.deathsK)) +
  geom_point(colour="#0072BB", alpha=0.1) +
  #stat_summary(aes(group=fit), fun =median, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
  stat_summary( fun = median, geom="line", colour="black", size=lsize) +
  # stat_summary( fun = mean, geom="line", colour="black", linetype="dashed", size=lsize) +
  geom_hline(data=df.sum.rw, aes(yintercept=sum(rw.deaths)), colour="black", linetype="dotted", size=lsize) +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Total deaths") +
  #scale_colour_brewer(palette="Set3") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0.01,0), limits=c(0,230000)) +
  theme_light()+
  ggtitle("Best-case scenario")

p.death.weekly.best<- ggplot(simdeath.comb.sum, aes(x=100*cperc, y=model.deathsK.best)) +
  geom_point(colour="#0072BB", alpha=0.1) +
  #stat_summary(aes(group=fit), fun =median, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
  stat_summary( fun = median, geom="line", colour="black", size=lsize) +
  #stat_summary( fun = mean, geom="line", colour="black", linetype="dashed", size=lsize) +
  geom_hline(data=rw.weekly, aes(yintercept=sum(newdeaths.weekly.actual, na.rm=T)), colour="black", linetype="dotted", size=lsize) +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Total deaths") +
  #scale_colour_brewer(palette="Set3") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0.01,0), limits=c(0,230000)) +
  theme_light()+
  ggtitle("Best-case scenario")

p.death.weekly.worst<- ggplot(simdeath.comb.sum, aes(x=100*cperc, y=model.deathsK.worst)) +
  geom_point(colour="#0072BB", alpha=0.1) +
  #stat_summary(aes(group=fit), fun =median, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
  stat_summary( fun = median, geom="line", colour="black", size=lsize) +
  #stat_summary( fun = mean, geom="line", colour="black", linetype="dashed", size=lsize) +
  geom_hline(data=rw.weekly, aes(yintercept=sum(newdeaths.weekly.actual, na.rm=T)), colour="black", linetype="dotted", size=lsize) +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Total deaths, CFR=f(cases)") +
  #scale_colour_brewer(palette="Set3") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0.01,0), limits=c(0,230000)) +
  theme_light()+
  ggtitle("Worst-case scenario")


p.death.saved<- ggplot(df.model.sum.comb, aes(x=100*cperc, y=-red.deathsK, colour=fit)) +
  #stat_summary(aes(group=fit), fun =mean, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
  stat_summary( fun = median, geom="line", colour="black", size=lsize) +
  #stat_summary( fun = mean, geom="line", colour="black", linetype="dashed", size=lsize) +
  geom_vline(aes(xintercept=45), colour="black", linetype="dotted", size=lsize) +
  geom_vline(aes(xintercept=60), colour="black", linetype="dotted", size=lsize) +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Total deaths prevented (by Aug 15, 2020)") +
  #scale_colour_brewer(palette="Set3") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0.01,0)) +
  theme_light()+
  ggtitle("Best-case scenario")

p.death.weekly.saved.best<- ggplot(simdeath.comb.sum, aes(x=100*cperc, y=-red.deathsK.best, colour=fit)) +
  #stat_summary(aes(group=fit), fun =mean, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
  stat_summary( fun = median, geom="line", colour="black", size=lsize) +
  #stat_summary( fun = mean, geom="line", colour="black", linetype="dashed", size=lsize) +
  geom_vline(aes(xintercept=45), colour="black", linetype="dotted", size=lsize) +
  geom_vline(aes(xintercept=60), colour="black", linetype="dotted", size=lsize) +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Total deaths prevented") +
  #scale_colour_brewer(palette="Set3") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0.01,0)) +
  theme_light() +
  ggtitle("Worst-case scenario")

p.death.weekly.saved.worst<- ggplot(simdeath.comb.sum, aes(x=100*cperc, y=-red.deathsK.worst, colour=fit)) +
  #stat_summary(aes(group=fit), fun =mean, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
  stat_summary( fun = median, geom="line", colour="black", size=lsize) +
  #stat_summary( fun = mean, geom="line", colour="black", linetype="dashed", size=lsize) +
  geom_vline(aes(xintercept=45), colour="black", linetype="dotted", size=lsize) +
  geom_vline(aes(xintercept=60), colour="black", linetype="dotted", size=lsize) +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Total deaths prevented, CFR=f(cases)") +
  #scale_colour_brewer(palette="Set3") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0.01,0)) +
  theme_light() +
  ggtitle("Worst-case scenario")

p.infected.saved<- ggplot(df.model.sum.comb, aes(x=100*cperc, y=-red.infected, colour=fit)) +
  #stat_summary(aes(group=fit), fun =mean, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
  stat_summary( fun = median, geom="line", colour="black", size=lsize) +
  # stat_summary( fun = mean, geom="line", colour="black", linetype="dashed", size=lsize) +
  geom_vline(aes(xintercept=45), colour="black", linetype="dotted", size=lsize) +
  geom_vline(aes(xintercept=60), colour="black", linetype="dotted", size=lsize) +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Total infections prevented (by Aug 15, 2020)") +
  #geom_vline(aes(xintercept=87), colour="black", linetype="dotted", size=lsize) +
  #scale_colour_brewer(palette="Set3") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0.01,0)) +
  theme_light()+
  ggtitle("Best-case scenario")

#Data frame to make our rectangles
d1=data.frame(x1=c(55,55), x2=c(59,59), y1=c(0,0), y2=c(70,5000))

p.infected.tosave<- ggplot() +
  geom_rect(data=d1[d1$y2==70,], aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey", alpha=0.3) +
  #stat_summary(aes(group=fit), fun =mean, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
  stat_summary(data=df.model.sum.comb, aes(x=100*(cperc), y=-round(sum(pops)*lead.incr.cperc)/lead.incr.infected, colour=fit),  fun = median, geom="point", colour="#0072BB") +
  #stat_summary( fun = mean, geom="line", colour="black", linetype="dashed") +
  xlab("Individuals adhering to NPIs (% of population)") +
  # ylab(paste("Median # of additional individuals distancing \nto prevent 1 infection (to Aug 15, 2020)", sep='')) +
  ylab("Median # of additional individuals adhering to NPIs to prevent 1 infection") +
  # geom_vline(aes(xintercept=45), colour="black", linetype="dotted", size=lsize) +
  # geom_vline(aes(xintercept=60), colour="black", linetype="dotted", size=lsize) +
  #scale_colour_brewer(palette="Set3") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0,0)) +
  #coord_cartesian(expand=FALSE, xlim=c(0,100), ylim=c(0,25)) +
  theme_light()

p.death.tosave<- ggplot() +
  geom_rect(data=d1[d1$y2==5000,], aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey", alpha=0.3) +
  # stat_summary(aes(group=fit), fun =mean, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
  stat_summary(data=df.model.sum.comb, aes(x=100*(cperc), y=-round(sum(pops)*lead.incr.cperc)/lead.incr.deaths, colour=fit), fun = median, geom="point", colour="#0072BB") +
  #stat_summary( fun = mean, geom="line", colour="black", linetype="dashed") +
  xlab("Individuals adhering to NPIs (% of population)") +
  # ylab(paste("Median # of additional individuals distancing \nto prevent 1 death (to Aug 15, 2020)", sep='')) +
  ylab("Median # of additional individuals adhering to NPIs to prevent 1 death") +
  # geom_vline(aes(xintercept=45), colour="black", linetype="dotted", size=lsize) +
  # geom_vline(aes(xintercept=60), colour="black", linetype="dotted", size=lsize) +
  #scale_colour_brewer(palette="Set3") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0,0)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,750)) +
  #coord_cartesian(expand=FALSE, xlim=c(0,100), ylim=c(0,750)) +
  coord_cartesian(ylim=c(0,5000)) +
  theme_light()+
  ggtitle("Best-case scenario")

p.death.weekly.tosave.best<- ggplot() +
  geom_rect(data=d1[d1$y2==5000,], aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey", alpha=0.3) +
  # stat_summary(aes(group=fit), fun =mean, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
  stat_summary(data=simdeath.comb.sum, aes(x=100*(cperc), y=-round(sum(pops)*lead.incr.cperc)/lead.incr.deaths.best, colour=fit), fun = median, geom="point", colour="#0072BB") +
  #stat_summary( fun = mean, geom="line", colour="black", linetype="dashed") +
  xlab("Individuals adhering to NPIs (% of population)") +
  # ylab(paste("Median # of additional individuals distancing \nto prevent 1 death (to Aug 15, 2020)", sep='')) +
  ylab("Median # of additional individuals adhering to NPIs to prevent 1 death") +
  # geom_vline(aes(xintercept=45), colour="black", linetype="dotted", size=lsize) +
  # geom_vline(aes(xintercept=60), colour="black", linetype="dotted", size=lsize) +
  #scale_colour_brewer(palette="Set3") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0,0)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,750)) +
  coord_cartesian(ylim=c(0,5000)) +
  theme_light()+
  ggtitle("Best-case scenario")

p.death.weekly.tosave.worst<- ggplot() +
  geom_rect(data=d1[d1$y2==5000,], aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey", alpha=0.3) +
  # stat_summary(aes(group=fit), fun =mean, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
  stat_summary(data=simdeath.comb.sum, aes(x=100*(cperc), y=-round(sum(pops)*lead.incr.cperc)/lead.incr.deaths.worst, colour=fit), fun = median, geom="point", colour="#0072BB") +
  #stat_summary( fun = mean, geom="line", colour="black", linetype="dashed") +
  xlab("Individuals adhering to NPIs (% of population)") +
  # ylab(paste("Median # of additional individuals distancing \nto prevent 1 death (to Aug 15, 2020)", sep='')) +
  ylab("Median # of additional individuals adhering to NPIs to prevent 1 death, CFR=f(cases)") +
  # geom_vline(aes(xintercept=45), colour="black", linetype="dotted", size=lsize) +
  # geom_vline(aes(xintercept=60), colour="black", linetype="dotted", size=lsize) +
  #scale_colour_brewer(palette="Set3") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0,0)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,750)) +
  coord_cartesian(ylim=c(0,5000)) +
  theme_light()+
  ggtitle("Worst-case scenario")

# p.comb<-p.infected +   p.death +  p.known + p.infected.saved +  p.death.saved + plot_spacer() +  p.infected.tosave + p.death.tosave +
#   plot_layout(guides = 'collect', nrow=3) +
#   plot_annotation(tag_levels = 'a', title=sprintf("Importance of voluntary distancing by susceptible: %f", unique(df.model.sum.comb$susceffect))) &
#   theme(legend.position = 'bottom',  legend.margin = margin(0, 0, 0, -5, "pt"), plot.margin = margin(0, 7, 2, 0, "pt"))


p.comb.save<-p.infected.tosave + p.death.weekly.tosave.best + p.death.weekly.tosave.worst +
  plot_layout(guides = 'collect', nrow=1) +
  plot_annotation(tag_levels = 'a') & #, title=sprintf("Importance of voluntary distancing by susceptible: %f", unique(df.model.sum.comb$susceffect))) &
  theme(legend.position = 'bottom',  legend.margin = margin(0, 0, 0, -5, "pt"), plot.margin = margin(0, 7, 2, 0, "pt"))


png("Plot_IndividualExperiment_main_v2.png", width=28, height=16, units="cm", res=500)
print(p.comb.save)
dev.off()


p.comb.totals<-p.known + p.death.weekly.best + p.death.weekly.worst +
  plot_layout(guides = 'collect', nrow=1) +
  plot_annotation(tag_levels = 'a') & #, title=sprintf("Importance of voluntary distancing by susceptible: %f", unique(df.model.sum.comb$susceffect))) &
  theme(legend.position = 'bottom',  legend.margin = margin(0, 0, 0, -5, "pt"), plot.margin = margin(0, 7, 2, 0, "pt"))


png("Plot_IndividualExperiment_si1_v2.png", width=29.25, height=12.5, units="cm", res=500)
print(p.comb.totals)
dev.off()

p.comb.ts<-p.known.ts + guide_area()  + p.death.weekly.ts.best + theme(legend.position="none") + p.death.weekly.ts.worst  + theme(legend.position="none") + #p.hosp.ts + p.icu.ts +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a') &
  theme(legend.margin = margin(0, 0, 0, -5, "pt"), plot.margin = margin(0, 7, 2, 0, "pt"))

png("Plot_IndividualExperiment_si2_v2.png", width=22.5, height=22.5, units="cm", res=500)
print(p.comb.ts)
dev.off()


hel=lo

df.model.sum.comb.final <- df.model.sum.comb %>% mutate(case.tosave= -round(sum(pops)*lead.incr.cperc)/lead.incr.infected)

df.model.sum.comb.final.sum<-df.model.sum.comb.final %>% group_by(cperc) %>% dplyr::summarise(med.case.tosave = median(case.tosave), med.case.total=median(model.known))


simdeath.comb.sum.final <- simdeath.comb.sum %>% mutate(
  death.tosave.best=-round(sum(pops)*lead.incr.cperc)/lead.incr.deaths.best, 
  death.tosave.worst=-round(sum(pops)*lead.incr.cperc)/lead.incr.deaths.worst)

simdeath.comb.sum.final.sum<-simdeath.comb.sum.final %>% group_by(cperc) %>% dplyr::summarise(
  med.death.tosave.best=median(death.tosave.best),
  med.death.total.best=median(model.deathsK.best),
  med.death.tosave.worst=median(death.tosave.worst),
  med.death.total.worst=median(model.deathsK.worst))

# p.infected.incr<- ggplot(df.model.sum.comb, aes(x=100*cperc, y=incr.infected, colour=fit)) +
#   stat_summary(aes(group=fit), fun =mean, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
#   stat_summary( fun = median, geom="line", colour="black") +
#   stat_summary( fun = mean, geom="line", colour="black", linetype="dashed") +
#   xlab("Individuals voluntarily distancing (% of population)") +
#   ylab(expression(Delta~"Total infections (by Aug 15, 2020)")) +
#   geom_vline(aes(xintercept=75), colour="black", linetype="dotted") +
#   #scale_colour_brewer(palette="Set3") +
#   guides(colour=F) +
#   scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
#   theme_light()
#
# p.known.incr<- ggplot(df.model.sum.comb, aes(x=100*cperc, y=incr.known, colour=fit)) +
#   stat_summary(aes(group=fit), fun =mean, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
#   stat_summary( fun = median, geom="line", colour="black") +
#   stat_summary( fun = mean, geom="line", colour="black", linetype="dashed") +
#   xlab("Individuals voluntarily distancing (% of population)") +
#   ylab(expression(Delta~"Confirmed cases (by Aug 15, 2020)")) +
#   geom_vline(aes(xintercept=75), colour="black", linetype="dotted") +
#   #scale_colour_brewer(palette="Set3") +
#   guides(colour=F) +
#   scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
#   theme_light()
#
# p.death.incr<- ggplot(df.model.sum.comb, aes(x=100*cperc, y=incr.deaths, colour=fit)) +
#   stat_summary(aes(group=fit), fun =mean, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
#   stat_summary( fun = median, geom="line", colour="black") +
#   stat_summary( fun = mean, geom="line", colour="black", linetype="dashed") +
#   xlab("Individuals voluntarily distancing (% of population)") +
#   ylab(expression(Delta~"Deaths (by Aug 15, 2020)")) +
#   geom_vline(aes(xintercept=75), colour="black", linetype="dotted") +
#   #scale_colour_brewer(palette="Set3") +
#   guides(colour=F) +
#   scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
#   theme_light()
#
#
# p.infected.marginal<- ggplot(df.model.sum.comb, aes(x=100*cperc, y=incr.infected/round(sum(pops)*incr.cperc), colour=fit)) +
#   stat_summary(aes(group=fit), fun =mean, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
#   stat_summary( fun = median, geom="line", colour="black") +
#   stat_summary( fun = mean, geom="line", colour="black", linetype="dashed") +
#   xlab("Individuals voluntarily distancing (% of population)") +
#   ylab("Est. marginal infections (by Aug 15, 2020)") +
#   geom_vline(aes(xintercept=75), colour="black", linetype="dotted") +
#   #scale_colour_brewer(palette="Set3") +
#   guides(colour=F) +
#   scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
#   theme_light()
#
# p.known.marginal<- ggplot(df.model.sum.comb, aes(x=100*cperc, y=incr.known/round(sum(pops)*incr.cperc), colour=fit)) +
#   stat_summary(aes(group=fit), fun =mean, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
#   stat_summary( fun = median, geom="line", colour="black") +
#   stat_summary( fun = mean, geom="line", colour="black", linetype="dashed") +
#   xlab("Individuals voluntarily distancing (% of population)") +
#   ylab("Est. marginal confirmed cases (by Aug 15, 2020)") +
#   geom_vline(aes(xintercept=75), colour="black", linetype="dotted") +
#   #scale_colour_brewer(palette="Set3") +
#   guides(colour=F) +
#   scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
#   theme_light()
#
# p.death.marginal<- ggplot(df.model.sum.comb, aes(x=100*cperc, y=incr.deaths/round(sum(pops)*incr.cperc), colour=fit)) +
#   stat_summary(aes(group=fit), fun =mean, fun.min = function(x) min(x),  fun.max = function(x) max(x), geom = "pointrange", size=psize) +
#   stat_summary( fun = median, geom="line", colour="black") +
#   stat_summary( fun = mean, geom="line", colour="black", linetype="dashed") +
#   xlab("Individuals voluntarily distancing (% of population)") +
#   ylab("Est. marginal deaths (by Aug 15, 2020)") +
#   geom_vline(aes(xintercept=75), colour="black", linetype="dotted") +
#   #scale_colour_brewer(palette="Set3") +
#   guides(colour=F) +
#   scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
#   theme_light()

# p.comb<-wrap_plots(p.infected, p.known, p.death, p.infected.incr, p.known.incr, p.death.incr,  p.infected.marginal,  p.known.marginal,  p.death.marginal, nrow=3) +
#   plot_layout(guides = 'collect') +
#   plot_annotation(tag_levels = 'a', title=sprintf("Importance of voluntary distancing by susceptible: %f", unique(df.model.sum.comb$susceffect))) &
#   theme(legend.position = 'bottom',  legend.margin = margin(0, 0, 0, -5, "pt"), plot.margin = margin(0, 7, 2, 0, "pt"))


# p.comb.save<-wrap_plots(p.infected.save, p.death.save, nrow=1) +
#   plot_layout(guides = 'collect') +
#   plot_annotation(tag_levels = 'a', title=sprintf("Importance of voluntary distancing by susceptible: %f", unique(df.model.sum.comb$susceffect))) &
#   theme(legend.position = 'bottom',  legend.margin = margin(0, 0, 0, -5, "pt"), plot.margin = margin(0, 7, 2, 0, "pt"))
#
# png("Plot_ncincrSAVES_V37.png", width=25, height=17.5, units="cm", res=500)
# print(p.comb.save)
# dev.off()


# p.infected.single<- ggplot(df.model.sum.comb[(df.model.sum.comb$run==1 & df.model.sum.comb$fit==unique(df.model.sum.comb$fit)[1]),], aes(x=100*cperc, y=model.sick, colour=as.factor(run))) +
#   geom_line(alpha=1) +
#   # stat_summary(fun = median,geom="line", colour="black") +
#   # stat_summary(fun = mean, geom="line", colour="black", linetype="dashed") +
#   xlab("Individuals voluntarily distancing (% of population)") +
#   ylab("Total infections (by Aug 15, 2020)") +
#   guides(colour=F) +
#   scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
#   theme_light()
#
# p.incr.infected.single<- ggplot(df.model.sum.comb[(df.model.sum.comb$run==1 & df.model.sum.comb$fit==unique(df.model.sum.comb$fit)[1]),], aes(x=100*cperc, y=incr.infected, colour=as.factor(run))) +
#   geom_point(alpha=1) +
#   # stat_summary(fun = median,geom="line", colour="black") +
#   # stat_summary(fun = mean,geom="line", colour="black", linetype="dashed") +
#   xlab("Individuals voluntarily distancing (% of population)") +
#   ylab(expression(Delta~"Total infections (by Aug 15, 2020)")) +
#   guides(colour=F) +
#   scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
#   theme_light()
#
# p.comb.single<-wrap_plots(p.infected.single,p.incr.infected.single, nrow=1) +
#   plot_layout(guides = 'collect') +
#   plot_annotation(tag_levels = 'a', title=sprintf("Importance of voluntary distancing by susceptible: %f", unique(df.model.sum.comb$susceffect))) &
#   theme(legend.position = 'bottom',  legend.margin = margin(0, 0, 0, -5, "pt"), plot.margin = margin(0, 7, 2, 0, "pt"))
#
# png("Plot_ncincr_single_constseed.png", width=25, height=10, units="cm", res=500)
# print(p.comb.single)
# dev.off()
