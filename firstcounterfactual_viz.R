#     Script visualizes output from covidHierV46_2Github_firstcounterfactual.R
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
library(tidyverse); library(scales); library(ggdark); library(Hmisc); library(ggnewscale);
`%notin%` <- Negate(`%in%`);

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

#Select set of sims to visualize
code.list<-c("v46v2phuNUDGE_bothopen_restricted_vdON","v46v2phuNUDGE_bothopen_restricted_vdOFF", "v46v2phuNUDGE_neitheropen_restricted_vdON","v46v2phuNUDGE_neitheropen_restricted_vdOFF")

if (exists('mobilitydat')==FALSE) #Skips running this step if has already been done (marked by existence of final object generated in this step)
{
  DATA=readRDS("covidHierData.rds"); #See readme file for description of contents
  pops=colSums(DATA$Msave)
  
  #Read in empircal data for comparison
  rw.ind.0<-read.csv("conposcovidloc_mar32021_simple.csv", na.strings=c(""," ","NA"))
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
  for (i in 1:length(rw.ind$Age_Group))  {rw.ind$ageclass[i] <- agetagger(rw.ind$Age_Group[i])}
  
  ### Expand ind data to cover all possible dates and age classes
  rw.ind$Date<-as.Date(rw.ind$Case_Reported_Date)
  temp50<- rw.ind %>% group_by(ageclass, Date) %>% dplyr::summarize (newcases=n(), deaths=sum(Outcome1=="Fatal"))
  temp50 <- temp50[order(temp50$Date),]
  temp50$totalcases<-cumsum(temp50$newcases)
  df.dates<-data.frame(Date=(seq(as.Date("2020-01-01"), as.Date("2021-03-31"), by = "days")))
  df.50<-data.frame(Date=(seq(as.Date("2020-01-01"), as.Date("2021-03-31"), by = "days")), ts50=1:length(seq(as.Date("2020-01-01"), as.Date("2021-03-31"), by = "days")) - length(seq(as.Date("2020-01-01"),temp50$Date[min(which(temp50$totalcases>=50))], by="days")))
  rw.ind.dates<-merge(rw.ind, df.dates, by="Date", all=TRUE)
  template.full<- rw.ind.dates %>% expand(ageclass, Date, Reporting_PHU)
  rw.ind.full<-merge(template.full, rw.ind.dates, by=c("ageclass", "Date", "Reporting_PHU"), all=T)
  rw.ind.full<-merge(rw.ind.full, df.50, by="Date")
  
  template.full.OB<- rw.ind.dates[rw.ind.dates$Outbreak_Related=="Yes",] %>% expand(ageclass, Date, Reporting_PHU)
  rw.ind.OB<-merge(template.full.OB, rw.ind.dates[rw.ind.dates$Outbreak_Related=="Yes",], by=c("ageclass", "Date", "Reporting_PHU"), all=T)
  rw.ind.OB<-merge(rw.ind.OB, df.50, by="Date")
  
  rw.agg<- rw.ind.full %>% group_by(ageclass, Date,ts50) %>% dplyr::summarize (newcases=sum(!is.na(Case_Reported_Date)))
  rw.agg<-rw.agg[order(rw.agg$Date),]
  
  rw.agg.OB<- rw.ind.OB %>% group_by(ageclass, Date,ts50) %>% dplyr::summarize (newcases.OB=sum(!is.na(Case_Reported_Date)))
  rw.agg.OB<-rw.agg.OB[order(rw.agg.OB$Date),]
  
  rw.agg.region<- rw.ind.full %>% group_by(ageclass, Date, Reporting_PHU, ts50) %>% dplyr::summarize (newcases=sum(!is.na(Case_Reported_Date))) 
  rw.agg.region<-rw.agg.region[order(rw.agg.region$Date),]
  
  rw.tots$new.K<-rw.tots$Total.Cases-lag(rw.tots$Total.Cases)
  rw.tots<-merge(rw.tots, df.50, by="Date", all.x=TRUE)
 
  #Write a subset of rw.agg for fitting
  rw.agg.full<-merge(rw.agg, df.dates, by="Date", all.x=TRUE, all.y=TRUE)
  rw.agg.fit<-rw.agg.full[order(rw.agg.full$ageclass), colnames(rw.agg.full) %notin% "deaths"] %>%
    pivot_wider(names_from = ageclass, values_from = newcases,  names_prefix="newK.")
  rw.agg.fit<-rw.agg.fit[order(rw.agg.fit$ts50),]
  rw.agg.fit<-rw.agg.fit[rw.agg.fit$Date<=max(as.Date(rw.ind.0$Case_Reported_Date),na.rm=TRUE),] #Drop any dates from the fitting that go past the final day of data we have from PHO
  rw.agg.fit<-rw.agg.fit[rowSums(is.na(rw.agg.fit)) != ncol(rw.agg.fit),] #Get rid of rows that are all NAs
  rw.agg.fit$total.K<-cumsum(rowSums(rw.agg.fit[,3:9], na.rm=TRUE))
  rw.agg.fit$new.K<-rowSums(rw.agg.fit[,3:9], na.rm=TRUE)
  
  rw.agg.OB.full<-merge(rw.agg.OB, df.50, by="Date", all.x=TRUE, all.y=TRUE)
  colnames(rw.agg.OB.full)[colnames(rw.agg.OB.full)=="ts50.y"]<-"ts50"
  rw.agg.OB.fit<-rw.agg.OB.full[order(rw.agg.OB.full$ageclass), colnames(rw.agg.OB.full) %notin% "deaths.OB"] %>%
    pivot_wider(names_from = ageclass, values_from = newcases.OB,  names_prefix="newK.OB.")
  rw.agg.OB.fit<-rw.agg.OB.fit[order(rw.agg.OB.fit$ts50),]
  rw.agg.OB.fit<-rw.agg.OB.fit[rw.agg.OB.fit$Date<=max(as.Date(rw.ind.0$Case_Reported_Date),na.rm=TRUE),] #Drop any dates from the fitting that go past the final day of data we have from PHO
  rw.agg.OB.fit<-rw.agg.OB.fit[rowSums(is.na(rw.agg.OB.fit)) != ncol(rw.agg.OB.fit),] #Get rid of rows that are all NAs
  rw.agg.OB.fit$total.K.OB<-cumsum(rowSums(rw.agg.OB.fit[,4:10], na.rm=TRUE))
  rw.agg.OB.fit$new.K.OB<-rowSums(rw.agg.OB.fit[,4:10], na.rm=TRUE)
  
  rw.agg.fit<-merge(rw.agg.fit, rw.agg.OB.fit[,c("ts50", "total.K.OB", "new.K.OB")], by="ts50")
  rw.agg<-merge(rw.agg, rw.agg.OB, by=c("ageclass", "Date", "ts50"), all.x=T, all.y=T)
  rw.agg$newcases.OB[is.na(rw.agg$newcases.OB)]<-0
  rw.agg<-rw.agg[order(rw.agg$Date),]
  
  #Read in other info
  mobilitydat<-readRDS("mobilitydat_REAL_to2021-02-27_V4.rds") #read in mobility data
}


if (exists('calccheck')==FALSE) #Skips running this data collating/hosp outcome generating step if we're looking at the same codeeversion and this step has already been done (marked by existence of final object generated in this step)
{
  
  for(codeselect in 1:length(code.list))  {
    
    ###Set code version
    codeversion<-code.list[codeselect]

    #Read in simulation data
    #Select how many of the top fits you want (i.e. topnum<-10 indicates you want the 10 best parameter sets)
    topnum<-10
    
    fitparms00=readRDS(sprintf("simsFITPARMS_all_COUNTERFACTUAL_%s_V2.rds", codeversion));
    simdata00=readRDS(sprintf("simsPROVINCE_all_COUNTERFACTUAL_%s_V2.rds", codeversion));
    simdata00 <- simdata00 %>% rename(meanVD=V66, sdVD=V67)

    bestcutoff<-if(length(unique(simdata00$LL))>(topnum+1)) {Rfast::nth(unique(simdata00$LL), topnum+1, descending = F)} else {1e8}
    fitparms0 <- fitparms00[fitparms00$LL<bestcutoff,]
    simdata0 <- simdata00[simdata00$LL<bestcutoff,]
    
    #Get rid of any dup parm sets (ID by LL)
    checker<-fitparms0[!duplicated(fitparms0[ , "LL" ] ), ]
    fitparms<- fitparms0[fitparms0$fit %in% checker$fit,]
    simdata<- simdata0[simdata0$fit %in% checker$fit,]
   
    ### Add some calendar dates for ease of reading
    df.dates.big<-data.frame(Date=(seq(as.Date("2019-06-01"), as.Date("2021-03-31"), by = "days")), ts50=1:length(seq(as.Date("2019-06-01"), as.Date("2021-03-31"), by = "days")) - length(seq(as.Date("2019-06-01"),rw.agg$Date[min(which(cumsum(rw.agg$newcases)>=50))], by="days")))
    simdata<-merge(simdata, df.dates.big, by="ts50")
    
    ########### CFR-based death calculations
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
    
    rw.tots$week<-as.Date(cut(rw.tots$Date,
                              breaks = "week")) 
    rw.tots$deaths[is.na(rw.tots$Deaths)]<-0
    rw.tots$newdeaths<-rw.tots$Deaths - lag(rw.tots$Deaths)
    rw.weekly <- rw.tots %>% group_by(week) %>% dplyr::summarise(newcases.weekly.actual=sum(new.K), newdeaths.weekly.actual=sum(newdeaths),
                                                                 mean.icu.actual=mean(Number.of.patients.in.ICU.with.COVID.19 + Number.of.patients.in.ICU.on.a.ventilator.with.COVID.19))
    rw.weekly$cases.l2.actual=lag(rw.weekly$newcases.weekly.actual,2)
    rw.weekly$weekly.cfr.actual<-rw.weekly$newdeaths.weekly.actual/rw.weekly$cases.l2.actual
    rw.weekly$binned.icu<-cut(rw.weekly$mean.icu.actual, breaks=c(0,150,300,750))
    rw.weekly$wave<-cut(rw.weekly$week, breaks=c(as.Date("2020-01-20"), as.Date("2020-08-10"),as.Date("2021-02-22")), right=T)
    rw.weekly<-rw.weekly[is.na(rw.weekly$wave)==FALSE,]
    
    simdata$week<-as.Date(cut(simdata$Date,
                              breaks = "week")) 
    simdeath <- simdata %>% group_by(run,fit,week,cftype,reopeningtype,vdtype) %>% dplyr::summarise(newcases.weekly=sum(new.K))
    simdeath <- simdeath %>% group_by(run,fit,cftype,reopeningtype,vdtype) %>% mutate(cases.l2=lag(newcases.weekly,2), weekly.deaths.worst=cfrfcn(cases.l2, week)*cases.l2) %>% ungroup()
    simdeath<-merge(simdeath, rw.weekly, by="week")
    simdeath$weekly.deaths.best<-simdeath$weekly.cfr.actual*simdeath$cases.l2
    ###########
    
    if (codeselect==1) {
      simdata.comb<-simdata
      simdeath.comb<-simdeath
    } else {
      simdata.comb<-rbind(simdata, simdata.comb)
      simdeath.comb<-rbind(simdeath.comb, simdeath)
    }
    
    print(codeselect)
    
    #To check if we've switched versions before allowing shortcut
    if(length(unique(interaction(simdata.comb$cftype,simdata.comb$reopeningtype)))==length(code.list)) {calccheck<-"good"}
    
  }
  
  rw.hosp.fit<-merge(rw.hosp.fit,df.dates.big, by="ts50", all.x=T, all.y=F)
  
  #Knock out dates in rw data we're unceratain of
  rw.agg<-rw.agg[rw.agg$Date<=as.Date("2021-02-21"),]
  rw.agg.fit<-rw.agg.fit[rw.agg.fit$Date<=as.Date("2021-02-21"),]
  rw.agg.region<-rw.agg.region[rw.agg.region$Date<=as.Date("2021-02-21"),]
}

### Viz
library(RColorBrewer); library(patchwork); library(Rfast);

pal.select<-c( "black", rev(brewer.pal(4,"RdBu")))

limselect<-c(as.Date("2020-03-10"), as.Date("2020-08-15")); #Set span of dates to plot
dbreak<-"1 month";
rw.agg.scenario<-rw.agg[(rw.agg$Date >= limselect[1] & rw.agg$Date <= limselect[2]),]

#Set transparency of lines from individual sims
alphaval<-0.15; rw.line<-0.5; summary.line<-0.5;

#Summarize provincial outcomes

simdata.comb <- simdata.comb %>% mutate(tag=factor(interaction(cftype,reopeningtype,vdtype), levels = c("rw", "neitheropen.restricted.vdON", "bothopen.restricted.vdON", "neitheropen.restricted.vdOFF","bothopen.restricted.vdOFF")), ordered=T)

p.cases.1<-ggplot() +
  geom_line(data=simdata.comb, aes(x=Date, y=new.K, group=interaction(run,fit,cftype,reopeningtype,vdtype), colour=tag), alpha=alphaval, key_glyph = "rect") +
  stat_summary(data=simdata.comb, fun=median,geom="line",aes(x=Date, y=new.K, group=interaction(cftype,reopeningtype,vdtype), colour=tag), size=summary.line, key_glyph = "rect")+
  geom_vline(aes(xintercept=as.Date("2020-10-09"))) +
  scale_colour_manual("Scenario", values=c("bothopen.restricted.vdON"=pal.select[3], "neitheropen.restricted.vdON"=pal.select[2], "bothopen.restricted.vdOFF"=pal.select[5], "neitheropen.restricted.vdOFF"=pal.select[4], "rw"=pal.select[1]),
                      labels=c("bothopen.restricted.vdON"="NPI adherence without closures", "neitheropen.restricted.vdON"="NPI adherence with closures",
                               "bothopen.restricted.vdOFF"="No NPI adherence or closures", "neitheropen.restricted.vdOFF"="Closures without NPI adherence",
                               "rw"="Actual"), limits=c("neitheropen.restricted.vdON", "bothopen.restricted.vdON", "neitheropen.restricted.vdOFF","bothopen.restricted.vdOFF") , guide = guide_legend(order = 2)) +
  new_scale_color() +
  scale_colour_manual("", values=c("bothopen.restricted.vdON"=pal.select[3], "neitheropen.restricted.vdON"=pal.select[2], "bothopen.restricted.vdOFF"=pal.select[5], "neitheropen.restricted.vdOFF"=pal.select[4], "rw"=pal.select[1]),
                      labels=c("bothopen.restricted.vdON"="NPI adherence without closures", "neitheropen.restricted.vdON"="NPI adherence with closures",
                               "bothopen.restricted.vdOFF"="No NPI adherence or closures", "neitheropen.restricted.vdOFF"="Closures without NPI adherence",
                               "rw"="Actual"), limits=c("rw") , guide = guide_legend(order = 1) ) +
  geom_line(data=rw.tots, aes(x=Date, y=new.K, colour=factor("rw", levels = c("rw"), ordered=T)), key_glyph = "rect") +
  xlab("Date") +
  ylab("New confirmed cases") +
  scale_x_date(expand = c(0,0), limits=limselect, date_breaks = dbreak, date_labels = "%b '%y") +
  scale_y_continuous(labels=comma, expand = c(0.01,0)) +
  theme_light()

p.cases.2<-ggplot() +
  geom_line(data=simdata.comb[interaction(simdata.comb$cftype,simdata.comb$reopeningtype,simdata.comb$vdtype)=="neitheropen.restricted.vdON",], aes(x=Date, y=new.K, group=interaction(run,fit,cftype,reopeningtype,vdtype), colour=interaction(cftype,reopeningtype,vdtype)), alpha=alphaval, key_glyph = "rect") +
  stat_summary(data=simdata.comb[interaction(simdata.comb$cftype,simdata.comb$reopeningtype,simdata.comb$vdtype)=="neitheropen.restricted.vdON",], fun=median,geom="line",aes(x=Date, y=new.K, group=interaction(cftype,reopeningtype,vdtype), colour=interaction(cftype,reopeningtype,vdtype)), size=summary.line, key_glyph = "rect")+
  geom_vline(aes(xintercept=as.Date("2020-10-09"))) +
  geom_line(data=rw.tots, aes(x=Date, y=new.K, colour="rw")) +
  xlab("Date") +
  ylab("New confirmed cases") +
  scale_x_date(expand = c(0,0), limits=limselect, date_breaks = dbreak, date_labels = "%b '%y") +
  scale_y_continuous(labels=comma, expand = c(0.01,0), limits=c(0,875)) + 
  scale_colour_manual("Source", values=c("bothopen.restricted.vdON"=pal.select[3], "neitheropen.restricted.vdON"=pal.select[2], "bothopen.restricted.vdOFF"=pal.select[5], "neitheropen.restricted.vdOFF"=pal.select[4], "rw"=pal.select[1]),
                      labels=c("bothopen.restricted.vdON"="NPI adherence without closures", "neitheropen.restricted.vdON"="NPI adherence with closures", 
                               "bothopen.restricted.vdOFF"="No NPI adherence or closures", "neitheropen.restricted.vdOFF"="Closures without NPI adherence",
                               "rw"="Actual")) +
  theme_light() + guides(colour="none") + theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0, 5, 2, 0, "pt"))

p.cases<- p.cases.1 + inset_element(p.cases.2, 0.575, 0.49, 0.995, 0.995 , ignore_tag = TRUE)

## Bar charts
simdeath.sum.0 <- simdeath.comb[(simdeath.comb$week>=limselect[1] & simdeath.comb$week<=limselect[2]),] %>% group_by(run,fit,cftype,reopeningtype,vdtype) %>% summarise(model.deathsweeklyK.best=sum(weekly.deaths.best, na.rm=TRUE),
                                                                                                                                                                        model.deathsweeklyK.worst=sum(weekly.deaths.worst, na.rm=TRUE))

simdeath.sum.ci <- simdeath.sum.0 %>% group_by(cftype,reopeningtype,vdtype) %>% dplyr::summarize(
                                                                                                meanci.deathsweeklyK.best=smean.cl.boot(model.deathsweeklyK.best, conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)[1],
                                                                                                lowci.deathsweeklyK.best=smean.cl.boot(model.deathsweeklyK.best, conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)[2],
                                                                                                highci.deathsweeklyK.best=smean.cl.boot(model.deathsweeklyK.best, conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)[3] ,
                                                                                                meanci.deathsweeklyK.worst=smean.cl.boot(model.deathsweeklyK.worst, conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)[1],
                                                                                                lowci.deathsweeklyK.worst=smean.cl.boot(model.deathsweeklyK.worst, conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)[2],
                                                                                                highci.deathsweeklyK.worst=smean.cl.boot(model.deathsweeklyK.worst, conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)[3])

simdeath.sum.1 <- simdeath.sum.0 %>% group_by(cftype,reopeningtype,vdtype) %>% dplyr::summarize(mean.model.deathsweeklyK.best=mean(model.deathsweeklyK.best),
                                                                                                min.model.deathsweeklyK.best=min(model.deathsweeklyK.best),
                                                                                                max.model.deathsweeklyK.best=max(model.deathsweeklyK.best),
                                                                                                mean.model.deathsweeklyK.worst=mean(model.deathsweeklyK.worst),
                                                                                                min.model.deathsweeklyK.worst=min(model.deathsweeklyK.worst),
                                                                                                max.model.deathsweeklyK.worst=max(model.deathsweeklyK.worst))

rw.weekly.sum <- rw.weekly[(rw.weekly$week>=limselect[1] & rw.weekly$week<=limselect[2]),]  %>% dplyr::summarize(mean.model.deathsweeklyK.best = sum(newdeaths.weekly.actual, na.rm=T), min.model.deathsweeklyK.best=NA, max.model.deathsweeklyK.best=NA)
rw.weekly.sum$mean.model.deathsweeklyK.worst<-rw.weekly.sum$mean.model.deathsweeklyK.best
rw.weekly.sum$min.model.deathsweeklyK.worst<-rw.weekly.sum$min.model.deathsweeklyK.best
rw.weekly.sum$max.model.deathsweeklyK.worst<-rw.weekly.sum$max.model.deathsweeklyK.best

rw.weekly.sum$cftype=rw.weekly.sum$reopeningtype=rw.weekly.sum$vdtype<-"rw"

simdeath.sum.2<-rbind(simdeath.sum.1, rw.weekly.sum)
simdeath.sum.2$groups <- factor(interaction(simdeath.sum.2$cftype,simdeath.sum.2$reopeningtype,simdeath.sum.2$vdtype),
                                levels = c("rw.rw.rw", "neitheropen.restricted.vdON", "bothopen.restricted.vdON", "neitheropen.restricted.vdOFF","bothopen.restricted.vdOFF"))

h.deathweekly.best<- ggplot() +
  geom_bar(data =simdeath.sum.2 , aes(x=groups, y= mean.model.deathsweeklyK.best, fill=groups), stat = "identity", position = position_dodge()) +
  geom_errorbar(data=simdeath.sum.2 , aes(x=groups, ymin=min.model.deathsweeklyK.best, ymax=max.model.deathsweeklyK.best, group=groups),width=.5,
                position=position_dodge(.9), colour= "black") +
  xlab("Source") +
  ylab("Deaths") +
  scale_x_discrete(labels=NULL) +
  scale_y_continuous(labels = comma, expand=c(0,0.01), limits=c(0,617000)) +
  scale_fill_manual("Source", values=c("bothopen.restricted.vdON"=pal.select[3], "neitheropen.restricted.vdON"=pal.select[2], "bothopen.restricted.vdOFF"=pal.select[5], "neitheropen.restricted.vdOFF"=pal.select[4], "rw.rw.rw"=pal.select[1]),
                    labels=c("bothopen.restricted.vdON"="NPI adherence without closures", "neitheropen.restricted.vdON"="NPI adherence with closures",
                             "bothopen.restricted.vdOFF"="No NPI adherence or closures", "neitheropen.restricted.vdOFF"="Closures without NPI adherence",
                             "rw.rw.rw"="Actual")) +
  theme_light() +
  ggtitle("Best-case scenario")

h.deathweekly.worst<- ggplot() +
  geom_bar(data =simdeath.sum.2 , aes(x=groups, y= mean.model.deathsweeklyK.worst, fill=groups), stat = "identity", position = position_dodge()) +
  geom_errorbar(data=simdeath.sum.2 , aes(x=groups, ymin=min.model.deathsweeklyK.worst, ymax=max.model.deathsweeklyK.worst, group=groups),width=.5,
                position=position_dodge(.9), colour= "black") +
  xlab("Source") +
  ylab("Deaths, CFR=f(cases)") +
  scale_x_discrete(labels=NULL) +
  scale_y_continuous(labels = comma, expand=c(0,0.01), limits=c(0,617000)) +
  scale_fill_manual("Source", values=c("bothopen.restricted.vdON"=pal.select[3], "neitheropen.restricted.vdON"=pal.select[2], "bothopen.restricted.vdOFF"=pal.select[5], "neitheropen.restricted.vdOFF"=pal.select[4], "rw.rw.rw"=pal.select[1]),
                    labels=c("bothopen.restricted.vdON"="NPI adherence without closures", "neitheropen.restricted.vdON"="NPI adherence with closures",
                             "bothopen.restricted.vdOFF"="No NPI adherence or closures", "neitheropen.restricted.vdOFF"="Closures without NPI adherence",
                             "rw.rw.rw"="Actual")) +
  theme_light() +
  ggtitle("Worst-case scenario")

p.comb.main<- (p.cases + guide_area() + plot_layout(guides="collect",widths=c(0.7,0.3))) / (h.deathweekly.best + guides(fill="none")  + h.deathweekly.worst + guides(fill="none") )   +
  plot_layout(heights=c(0.65,0.35)) +
  plot_annotation(tag_levels = 'a') & theme(legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0, "pt"),
                                            plot.margin = margin(0, 5, 2, 0, "pt"))

png("CFplot_wave1_main_V2.png", width=20, height=20, units="cm", res=500)
print(p.comb.main)
dev.off()

