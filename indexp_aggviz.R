#     Script condenses and visualizes output from covidHierV46_2Github_indexp.R
#     Copyright (C) 2021  Kathyrn R Fair
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

### DATA VIZ
library(tidyverse); library(Rfast); library(scales); library(Hmisc); library(ggnewscale); library(RColorBrewer); library(patchwork);
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

`%notin%` <- Negate(`%in%`);

### Set whether you need to condense simulation output ("yes"), or whether this has already been achieved ("no")
condense<-"no"

DATA=readRDS("covidHierData.rds"); #See readme file for description of contents
pops=colSums(DATA$Msave)

if (exists('rw.weekly')==FALSE) #Skips running this data collating step if already done (marked by existence of final object generated in this step)
{
  
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
  
  rw.agg<- rw.ind.full %>% group_by(ageclass, Date,ts50) %>% dplyr::summarize (newcases=sum(!is.na(Case_Reported_Date)))
  rw.agg<-rw.agg[order(rw.agg$Date),]
  
  rw.agg.OB<- rw.ind.OB %>% group_by(ageclass, Date,ts50) %>% dplyr::summarize (newcases.OB=sum(!is.na(Case_Reported_Date)))
  rw.agg.OB<-rw.agg.OB[order(rw.agg.OB$Date),]

  rw.tots$new.K<-rw.tots$Total.Cases-lag(rw.tots$Total.Cases)
  rw.tots<-merge(rw.tots, df.50, by="Date", all.x=TRUE)
  
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
  rw.agg<-rw.agg[order(rw.agg$Date),]
  
  rw.hosp.fit<-rw.tots[, c("Date", "ts50", "Deaths", "Number.of.patients.hospitalized.with.COVID.19", "Number.of.patients.in.ICU.with.COVID.19", "Number.of.patients.in.ICU.on.a.ventilator.with.COVID.19" )]
  rw.hosp.fit$newDeaths.rw.tot<-rw.hosp.fit$Deaths-lag(rw.hosp.fit$Deaths) 
  rw.hosp.fit$newDeaths.rw.tot[min(which(rw.hosp.fit$Deaths>=0))]<-rw.hosp.fit$Deaths[min(which(rw.hosp.fit$Deaths>=0))]
  
  #Drop all dates in "grey" period where there may be reporting lags, etc
  rw.agg<-rw.agg[rw.agg$Date<=as.Date("2021-02-28"),]
  rw.agg.fit<-rw.agg.fit[rw.agg.fit$Date<=as.Date("2021-02-28"),]
  
  ###### First portion of CFR calculations
  
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

if (condense=="yes")
{
  ### Condenses output from indexp simulations
path = getwd()
fitparms.file.names <- dir(path, pattern ="simsFITPARMS_suscRAND_ncincr_v46v2phuNUDGE_susceffect0.50_", full.names=TRUE)
simprovince.file.names <- dir(path, pattern ="simsPROVINCE_suscRAND_ncincr_v46v2phuNUDGE_susceffect0.50_", full.names=TRUE)

if(exists('df.model.sum.comb')==FALSE || (exists('df.model.sum.comb')==TRUE & dim(df.model.sum.comb)[1]!=(length(fitparms.file.names)*10))) {

  for(fpick in 1:length(fitparms.file.names))
  {
    #Read in simulation data
    #Select how many of the top fits you want (i.e. topnum<-10 indicates you want the 10 best parameter sets)
    topnum<-10

    fitparms00=readRDS(fitparms.file.names[fpick]);
    simdata00=readRDS(simprovince.file.names[fpick]);
    
    bestcutoff<-if(length(unique(simdata00$LL))>(topnum+1)) {Rfast::nth(unique(simdata00$LL), topnum+1, descending = F)} else {1e8}
    fitparms0 <- fitparms00[fitparms00$LL<bestcutoff,]
    simdata0 <- simdata00[simdata00$LL<bestcutoff,]

    #Get rid of any dup parm sets (ID by LL)
    checker<-fitparms0[!duplicated(fitparms0[ , "LL" ] ), ]
    fitparms<- fitparms0[fitparms0$fit %in% checker$fit,]
    simdata<- simdata0[simdata0$fit %in% checker$fit,]

    ### Add some calendar dates for ease of reading
    df.dates.big<-data.frame(Date=(seq(as.Date("2019-06-01"), as.Date("2021-01-31"), by = "days")), ts50=1:length(seq(as.Date("2019-06-01"), as.Date("2021-01-31"), by = "days")) - length(seq(as.Date("2019-06-01"),rw.agg$Date[min(which(cumsum(rw.agg$newcases)>=50))], by="days")))
    simdata<-merge(simdata, df.dates.big, by="ts50")

    limselect<-c(as.Date("2020-03-10"), as.Date("2020-08-15"))
    dbreak<-"1 month";
    simdata$cperc<-1-simdata$ncperc
    simdata<-simdata[(simdata$Date<=limselect[2] & simdata$Date>=limselect[1]),]
    simdata$new.K<-simdata$new.K.c + simdata$new.K.nc
    simdata$new.sick<-simdata$new.sick.c + simdata$new.sick.nc

    df.model.sum <- simdata %>% group_by(run,fit,cperc,susceffect) %>% dplyr::summarise(
      model.sick=sum(new.sick, na.rm=TRUE),
      model.known=sum(new.K, na.rm=TRUE))

    ########### Second part of CFR calculations
    df.parms<-readRDS("parmfits_cfrcase_fcnform1_v2.rds")
    df.parms$nu<-df.parms$nu/1e2
    fcnparms.w1<-df.parms[(df.parms$fcnform==1 & df.parms$lag==2 & df.parms$wave==1 & df.parms$SSE==min(df.parms$SSE[df.parms$fcnform==1 & df.parms$lag==2 & df.parms$wave==1])), c(1,2,3)]
    fcnparms.w2<-df.parms[(df.parms$fcnform==1 & df.parms$lag==2 & df.parms$wave==2 & df.parms$SSE==min(df.parms$SSE[df.parms$fcnform==1 & df.parms$lag==2 & df.parms$wave==2])), c(1,2,3)]
    fcnparms.w3<-df.parms[(df.parms$fcnform==1 & df.parms$lag==2 & df.parms$wave==3 & df.parms$SSE==min(df.parms$SSE[df.parms$fcnform==1 & df.parms$lag==2 & df.parms$wave==3])), c(1,2,3)]

    cfrfcn = function(cases2lag, week) {
      parms=data.frame(week=week);
      parms$v1<-ifelse(parms$week<="2020-05-25", fcnparms.w1[[1]], ifelse(parms$week<="2020-08-10",fcnparms.w2[[1]], fcnparms.w3[[1]]));
      parms$v2<-ifelse(parms$week<="2020-05-25", fcnparms.w1[[2]], ifelse(parms$week<="2020-08-10",fcnparms.w2[[2]], fcnparms.w3[[2]]));
      parms$v3<-ifelse(parms$week<="2020-05-25", fcnparms.w1[[3]], ifelse(parms$week<="2020-08-10",fcnparms.w2[[3]], fcnparms.w3[[3]]));
      x=parms$v1 + parms$v2*(1-exp(-parms$v3*cases2lag));  return(x);}

    simdata$week<-as.Date(cut(simdata$Date,
                              breaks = "week"))
    simdeath <- simdata %>% group_by(run,fit,ncperc,susceffect,week) %>% dplyr::summarise(newcases.weekly=sum(new.K.c + new.K.nc))
    simdeath <- simdeath %>% group_by(run,fit,ncperc,susceffect) %>% mutate(cases.l2=lag(newcases.weekly,2), weekly.deaths.worst=cfrfcn(cases.l2, week)*cases.l2) %>% ungroup()
    simdeath<-merge(simdeath, rw.weekly, by="week")
    simdeath$weekly.deaths.best<-simdeath$weekly.cfr.actual*simdeath$cases.l2

    ###########

    if (fpick==1) {
      df.model.sum.comb<-df.model.sum
      simdeath.comb<-simdeath
    } else {
      df.model.sum.comb<-rbind(df.model.sum.comb, df.model.sum)
      simdeath.comb<-rbind(simdeath.comb, simdeath)
    }

    print(fpick/length(fitparms.file.names))
  }

}


df.sum.rw<-data.frame(ageclass=c(1,2,3,4,5, NA), rw.known=c(sum(rw.agg$newcases[(rw.agg$Date<=limselect[2] & rw.agg$Date>=limselect[1] & rw.agg$ageclass==1)], na.rm = TRUE),sum(rw.agg$newcases[(rw.agg$Date<=limselect[2] & rw.agg$Date>=limselect[1] & rw.agg$ageclass==2)], na.rm = TRUE),sum(rw.agg$newcases[(rw.agg$Date<=limselect[2] & rw.agg$Date>=limselect[1] & rw.agg$ageclass==3)], na.rm = TRUE),sum(rw.agg$newcases[(rw.agg$Date<=limselect[2] & rw.agg$Date>=limselect[1] & rw.agg$ageclass==4)], na.rm = TRUE),sum(rw.agg$newcases[(rw.agg$Date<=limselect[2] & rw.agg$Date>=limselect[1] & rw.agg$ageclass==5)], na.rm = TRUE), 0),
                      rw.deaths=c(0,0,0,0,0,rw.hosp.fit$Deaths[(rw.hosp.fit$Date==limselect[2] & is.na(rw.hosp.fit$Date)==FALSE)]))

saveRDS(df.sum.rw, file = sprintf("dfIND_sum_rw_to%s_v2.rds", limselect[2]))
saveRDS(df.model.sum.comb, file = sprintf("dfIND_model_sum_comb_to%s_v2.rds", limselect[2]))
saveRDS(simdeath.comb, file = sprintf("dfIND_model_simdeath_to%s_v2.rds", limselect[2]))
}

### Viz
limselect<-c(as.Date("2020-03-10"), as.Date("2020-08-15")); dbreak<-"1 month";

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
                                                                                                   lead.incr.known=lead(model.known) - (model.known), lead.incr.infected=lead(model.sick) - (model.sick),
                                                                                                   incr.cperc=(cperc) - lag(cperc), lead.incr.cperc=lead(cperc) - (cperc),
                                                                                                   red.known=model.known-model.known[round(cperc,2)==0],
                                                                                                   red.infected=model.sick-model.sick[round(cperc,2)==0])

rw.weekly<-rw.weekly[(rw.weekly$week>=limselect[1] & rw.weekly$week<=limselect[2]),]


df.dates.big<-data.frame(Date=(seq(as.Date("2019-06-01"), as.Date("2021-01-31"), by = "days")), ts50=1:length(seq(as.Date("2019-06-01"), as.Date("2021-01-31"), by = "days")) - length(seq(as.Date("2019-06-01"),rw.agg$Date[min(which(cumsum(rw.agg$newcases)>=50))], by="days")))
if("Date" %notin% colnames(rw.hosp.fit)) {rw.hosp.fit<-merge(rw.hosp.fit, df.dates.big, by="ts50", all.x=T, all.y=F)}

#Drop all dates in "grey" period where there may be reporting lags, etc
rw.agg<-rw.agg[rw.agg$Date<=as.Date("2021-02-28"),]
rw.hosp.fit<-rw.hosp.fit[rw.hosp.fit$Date<=as.Date("2021-02-28"),]

custom.pal<- c("black",(brewer.pal(4,"RdBu")))

psize<-0.025
lsize<-0.75

p.known<- ggplot(df.model.sum.comb, aes(x=100*cperc, y=model.known)) +
  geom_point(colour="#0072BB", alpha=0.1) +
  stat_summary( fun = median, geom="line", colour="black", size=lsize) +
  geom_hline(data=df.sum.rw, aes(yintercept=sum(rw.known)), colour="black", linetype="dotted", size=lsize) +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Total confirmed cases") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0.01,0)) +
  theme_light()

p.death.weekly.best<- ggplot(simdeath.comb.sum, aes(x=100*cperc, y=model.deathsK.best)) +
  geom_point(colour="#0072BB", alpha=0.1) +
  stat_summary( fun = median, geom="line", colour="black", size=lsize) +
  geom_hline(data=rw.weekly, aes(yintercept=sum(newdeaths.weekly.actual, na.rm=T)), colour="black", linetype="dotted", size=lsize) +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Total deaths") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0.01,0), limits=c(0,230000)) +
  theme_light()+
  ggtitle("Best-case scenario")

p.death.weekly.worst<- ggplot(simdeath.comb.sum, aes(x=100*cperc, y=model.deathsK.worst)) +
  geom_point(colour="#0072BB", alpha=0.1) +
  stat_summary( fun = median, geom="line", colour="black", size=lsize) +
  geom_hline(data=rw.weekly, aes(yintercept=sum(newdeaths.weekly.actual, na.rm=T)), colour="black", linetype="dotted", size=lsize) +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Total deaths, CFR=f(cases)") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0.01,0), limits=c(0,230000)) +
  theme_light()+
  ggtitle("Worst-case scenario")

### Calcs to check values for shaded area
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

#Data frame to make our shaded areas
d1=data.frame(x1=c(55,55), x2=c(59,59), y1=c(0,0), y2=c(70,5000))

p.infected.tosave<- ggplot() +
  geom_rect(data=d1[d1$y2==70,], aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey", alpha=0.3) +
  stat_summary(data=df.model.sum.comb, aes(x=100*(cperc), y=-round(sum(pops)*lead.incr.cperc)/lead.incr.infected, colour=fit),  fun = median, geom="point", colour="#0072BB") +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Median # of additional individuals adhering to NPIs to prevent 1 infection") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0,0)) +
  theme_light()

p.death.tosave<- ggplot() +
  geom_rect(data=d1[d1$y2==5000,], aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey", alpha=0.3) +
  stat_summary(data=df.model.sum.comb, aes(x=100*(cperc), y=-round(sum(pops)*lead.incr.cperc)/lead.incr.deaths, colour=fit), fun = median, geom="point", colour="#0072BB") +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Median # of additional individuals adhering to NPIs to prevent 1 death") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0,0)) +
  coord_cartesian(ylim=c(0,5000)) +
  theme_light()+
  ggtitle("Best-case scenario")

p.death.weekly.tosave.best<- ggplot() +
  geom_rect(data=d1[d1$y2==5000,], aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey", alpha=0.3) +
  stat_summary(data=simdeath.comb.sum, aes(x=100*(cperc), y=-round(sum(pops)*lead.incr.cperc)/lead.incr.deaths.best, colour=fit), fun = median, geom="point", colour="#0072BB") +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Median # of additional individuals adhering to NPIs to prevent 1 death") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0,0)) +
  coord_cartesian(ylim=c(0,5000)) +
  theme_light()+
  ggtitle("Best-case scenario")

p.death.weekly.tosave.worst<- ggplot() +
  geom_rect(data=d1[d1$y2==5000,], aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey", alpha=0.3) +
  stat_summary(data=simdeath.comb.sum, aes(x=100*(cperc), y=-round(sum(pops)*lead.incr.cperc)/lead.incr.deaths.worst, colour=fit), fun = median, geom="point", colour="#0072BB") +
  xlab("Individuals adhering to NPIs (% of population)") +
  ylab("Median # of additional individuals adhering to NPIs to prevent 1 death, CFR=f(cases)") +
  guides(colour=F) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,100)) +
  scale_y_continuous(labels=comma, expand=c(0,0)) +
  coord_cartesian(ylim=c(0,5000)) +
  theme_light()+
  ggtitle("Worst-case scenario")


p.comb.save<-p.infected.tosave + p.death.weekly.tosave.best + p.death.weekly.tosave.worst +
  plot_layout(guides = 'collect', nrow=1) +
  plot_annotation(tag_levels = 'a') &
  theme(legend.position = 'bottom',  legend.margin = margin(0, 0, 0, -5, "pt"), plot.margin = margin(0, 7, 2, 0, "pt"))


png("Plot_IndividualExperiment_main_v2.png", width=28, height=16, units="cm", res=500)
print(p.comb.save)
dev.off()


p.comb.totals<-p.known + p.death.weekly.best + p.death.weekly.worst +
  plot_layout(guides = 'collect', nrow=1) +
  plot_annotation(tag_levels = 'a') &
  theme(legend.position = 'bottom',  legend.margin = margin(0, 0, 0, -5, "pt"), plot.margin = margin(0, 7, 2, 0, "pt"))


png("Plot_IndividualExperiment_si_v2.png", width=29.25, height=12.5, units="cm", res=500)
print(p.comb.totals)
dev.off()
