#     Script visualizes output from covidHierV46_2Github_basesim.R
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


library(tidyverse); library(Hmisc);  library(scales); library(RColorBrewer); library(patchwork); library(Rfast);
`%notin%` <- Negate(`%in%`);
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)
#Select set of simulations to visualize
codeversion<-c("v46v2phuNUDGE")

DATA=readRDS("covidHierData.rds"); #See readme file for description of contents
pops=colSums(DATA$Msave)

#Read in empircal data for comparison
rw.ind.0<-read.csv("conposcovidloc_mar32021.csv", na.strings=c(""," ","NA"))
rw.tots<-read.csv("covidtesting_mar32021.csv")

#Drop all dates where totals are impacted by reporting lags
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

##Aggregate reporting data
rw.agg<- rw.ind.full %>% group_by(ageclass, Date,ts50) %>% dplyr::summarize (newcases=sum(!is.na(Case_Reported_Date)))
rw.agg<-rw.agg[order(rw.agg$Date),]

rw.agg.OB<- rw.ind.OB %>% group_by(ageclass, Date,ts50) %>% dplyr::summarize (newcases.OB=sum(!is.na(Case_Reported_Date)))
rw.agg.OB<-rw.agg.OB[order(rw.agg.OB$Date),]

rw.agg.region<- rw.ind.full %>% group_by(ageclass, Date, Reporting_PHU, ts50) %>% dplyr::summarize (newcases=sum(!is.na(Case_Reported_Date)))
rw.agg.region<-rw.agg.region[order(rw.agg.region$Date),]


rw.tots$new.K<-rw.tots$Total.Cases-lag(rw.tots$Total.Cases)
rw.tots<-merge(rw.tots, df.50, by="Date", all.x=TRUE)

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

rw.hosp.fit<-rw.tots[, c("ts50", "Deaths", "Number.of.patients.hospitalized.with.COVID.19", "Number.of.patients.in.ICU.with.COVID.19", "Number.of.patients.in.ICU.on.a.ventilator.with.COVID.19" )]
rw.hosp.fit$newDeaths.rw.tot<-rw.hosp.fit$Deaths-lag(rw.hosp.fit$Deaths)
rw.hosp.fit$newDeaths.rw.tot[min(which(rw.hosp.fit$Deaths>=0))]<-rw.hosp.fit$Deaths[min(which(rw.hosp.fit$Deaths>=0))]

#Read in simulation data
#Select how many of the top fits you want (i.e. topnum<-10 indicates you want the 10 best parameter sets)
topnum<-10

fitparms00=readRDS(sprintf("simsFITPARMS_%s_V2.rds", codeversion));
simdata00=readRDS(sprintf("simsPROVINCE_%s_V2.rds", codeversion));
simdata00 <- simdata00 %>% dplyr::rename(meanVD=V66, sdVD=V67)

simregions00=readRDS(sprintf("simsREGIONS_%s_V2.rds", codeversion));
bestcutoff<-if(length(unique(simdata00$LL))>(topnum+1)) {Rfast::nth(unique(simdata00$LL), topnum+1, descending = F)} else {1e8}
fitparms0 <- fitparms00[fitparms00$LL<bestcutoff,]
simdata0 <- simdata00[simdata00$LL<bestcutoff,]
simregions0<- simregions00[simregions00$LL<bestcutoff,]

#Get rid of any dup parm sets (ID by LL)
checker<-fitparms0[!duplicated(fitparms0[ , "LL" ] ), ]
fitparms<- fitparms0[fitparms0$fit %in% checker$fit,]
simdata<- simdata0[simdata0$fit %in% checker$fit,]
simregions<- simregions0[simregions0$fit %in% checker$fit,]

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
rw.weekly <- rw.tots %>% group_by(week) %>% dplyr::summarise(newcases.weekly=sum(new.K), newdeaths.weekly=sum(newdeaths), 
                                                             mean.icu=mean(Number.of.patients.in.ICU.with.COVID.19 + Number.of.patients.in.ICU.on.a.ventilator.with.COVID.19))

rw.weekly$binned.icu<-cut(rw.weekly$mean.icu, breaks=c(0,150,300,750))
rw.weekly$wave<-cut(rw.weekly$week, breaks=c(as.Date("2020-01-20"), as.Date("2020-08-10"),as.Date("2021-02-22")), right=T)
rw.weekly<-rw.weekly[is.na(rw.weekly$wave)==FALSE,]

simdata$week<-as.Date(cut(simdata$Date,
                          breaks = "week")) 
simdeath <- simdata %>% group_by(run,fit,week) %>% dplyr::summarise(newcases.weekly=sum(new.K))
simdeath <- simdeath %>% group_by(run,fit) %>% mutate(cases.l2=lag(newcases.weekly,2), weekly.deaths=cfrfcn(cases.l2, week)*cases.l2) %>% ungroup()

########################


### Viz
rw.hosp.fit<-merge(rw.hosp.fit,df.dates.big, by="ts50", all.x=T, all.y=F)

#Drop all dates where totals are impacted by reporting lags
rw.agg<-rw.agg[rw.agg$Date<=as.Date("2021-02-28"),]
rw.agg.fit<-rw.agg.fit[rw.agg.fit$Date<=as.Date("2021-02-28"),]
rw.agg.region<-rw.agg.region[rw.agg.region$Date<=as.Date("2021-02-28"),]
rw.hosp.fit<-rw.hosp.fit[rw.hosp.fit$Date<=as.Date("2021-02-28"),]

pal.select<-c("#0072BB", "#E30022", "darkorange", "black")

limselect<-c(as.Date("2020-03-01"), as.Date("2021-02-28"));
dbreak<-"1 month";
rw.agg.scenario<-rw.agg[(rw.agg$Date >= limselect[1] & rw.agg$Date <= limselect[2]),]
rw.hosp.fit.scenario<-rw.hosp.fit[(rw.hosp.fit$Date >= limselect[1] & rw.hosp.fit$Date <= limselect[2]),]

alphaval<-0.15; rw.line<-0.5; summary.line<-1;

#Age-class plot

agelabels=c("0-19", "20-39", "40-59", "60-79", "80+")
agealpha<-0.4*alphaval;


p.case.list<-lapply(1:5,function(i)
{
  ggplot() +
    geom_line(data=simdata, aes_string(x="Date", y=sprintf("new.K.%i",i), group="interaction(run,fit)", colour=shQuote("modelK")), alpha=agealpha, key_glyph = "rect") +
    stat_summary(data=simdata, aes_string(x="Date", y=sprintf("new.K.%i",i), colour=shQuote("modelK")), fun=median,geom="line", size=rw.line, key_glyph = "rect")+
    geom_line(data=rw.agg[rw.agg$ageclass==i,], aes(x=Date, y=newcases, colour="rw"), size=rw.line, key_glyph = "rect") +
    xlab("Date") +
    ylab(sprintf("New confirmed cases, age %s", agelabels[i])) +
    scale_x_date(expand = c(0,0), limits=limselect, date_breaks = dbreak, date_labels = "%b '%y") +
    scale_y_continuous(labels=comma, expand = c(0.01,0), limits=c(0,max(simdata[simdata$Date<=max(limselect), sprintf("new.K.%i",i)], rw.agg$newcases[(rw.agg$ageclass==i  & rw.agg$Date<=max(limselect))], na.rm=T))) +
    scale_colour_manual("", values=c("modelT"=pal.select[2], "modelK"=pal.select[1], "optfit"=pal.select[3],  "rw"=pal.select[4]), labels=c("modelT"="Infections (model)", "modelK"="Model", "optfit"="Fit based on model",  "rw"="Actual")) +
    theme_light()
})


p.hcomb.cases<-wrap_plots(
  p.case.list[[1]]+ guides(colour = "none"),
  p.case.list[[2]]+ guides(colour = "none"),
  p.case.list[[3]]+ guides(colour = "none"),
  p.case.list[[4]]+ guides(colour = "none"),
  p.case.list[[5]],
  nrow=2) + guide_area() +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a') & theme(legend.margin = margin(0, 0, 0, -5, "pt"), plot.margin = margin(0, 5, 2, 0, "pt"),axis.text.x = element_text(angle = 45, hjust = 1))

png(sprintf("Plot_casesbyage_%s_V2.png", codeversion), width=25, height=16, units="cm", res=500)
print(p.hcomb.cases)
dev.off()


### PHU level plots

region.model.all<-merge(simregions, df.dates.big, by="ts50")
regionlink<-read.csv("simplePHUtoREGIONlinker.csv")
simregions.linked<-merge(region.model.all,regionlink[, colnames(regionlink) != "pop2016"], by.x="region", by.y="ï..regionid")
simregions.linked.summary <- simregions.linked %>% group_by(Date,allcase.phu,run,fit) %>% dplyr::summarize(newcases.allage=sum(new.K))

df.phu<-read.csv("Allcasetrendsdata_mar32021.csv")
df.phu$Dateclean<-as.Date(df.phu$Date , format = "%d-%m-%Y")
df.phu<-df.phu[df.phu$Public.Health.Unit!="Ontario",]

phu.labels<-c("Algoma", "Brant County", "Chatham-Kent", "Hamilton", "Durham", "Eastern Ontario", "Grey Bruce", "Haldimand-Norfolk", paste("Haliburton, Kawartha, \nPine Ridge", sep=""), "Halton", "Hastings Prince Edward",
              "Huron Perth", paste("Kingston, Frontenac, \nLennox & Addington", sep=""), "Lambton", "Leeds, Grenville & Lanark", "Middlesex-London", "Niagara", "North Bay Parry Sound", "Northwestern", "Ottawa",
              "Peel", "Peterborough", "Porcupine", "Sudbury & Districts", "Waterloo", "Renfrew County", "Simcoe Muskoka", "Southwestern", "Thunder Bay", "Timiskaming", "Toronto", "Wellington-Dufferin-Guelph",
              "Windsor-Essex", "York")


p.regiontotal.list<-lapply(1:length(unique(df.phu$Public.Health.Unit)), function(i)
{
  ggplot() +
    geom_line(data=simregions.linked.summary[simregions.linked.summary$allcase.phu==unique(df.phu$Public.Health.Unit)[i],], aes(x=Date, y=newcases.allage, group=interaction(run,fit), colour="modelK"), alpha=0.1, key_glyph = "rect") +
    stat_summary(data=simregions.linked.summary[simregions.linked.summary$allcase.phu==unique(df.phu$Public.Health.Unit)[i],], aes(x=Date, y=newcases.allage, colour="modelK"), fun=mean,geom="line", size=rw.line, key_glyph = "rect")+
    geom_line(data=df.phu[df.phu$Public.Health.Unit==unique(df.phu$Public.Health.Unit)[i],], aes(x=Dateclean, y=Cases.by.reported.date, colour="rw"), size=rw.line, key_glyph = "rect") +
    xlab("Date") +
    ylab("New confirmed cases") +
    scale_colour_manual("", values=c("modelT"=pal.select[2], "modelK"=pal.select[1], "optfit"=pal.select[3],  "rw"=pal.select[4]), labels=c("modelT"="Infections (model)", "modelK"="Model", "optfit"="Fit based on model",  "rw"="Actual")) +
    scale_x_date(expand = c(0,0), limits=limselect, date_breaks = dbreak, date_labels = "%b '%y") +
    scale_y_continuous(labels=comma, expand = c(0.01,0))+
    # ggtitle(str_trim(unique(df.phu$Public.Health.Unit)[i], side = c("right"))) +       
    ggtitle(phu.labels[i]) +
    theme_light()
})

p.regions1<-wrap_plots(p.regiontotal.list[1:17]) + guide_area() +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels="a") & theme(legend.box="vertical", plot.margin = margin(0, 5, 2, 0, "pt"),axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=8))

png(sprintf("Plot_phu1_%s_V2.png", codeversion), width=26, height=25, units="cm", res=500)
print(p.regions1)
dev.off()

p.regions2<-wrap_plots(p.regiontotal.list[18:34]) +  guide_area() + 
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels="a") & theme(legend.box="vertical", plot.margin = margin(0, 5, 2, 0, "pt"),axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=8))

png(sprintf("Plot_phu2_%s_V2.png", codeversion), width=26, height=25, units="cm", res=500)
print(p.regions2)
dev.off()
