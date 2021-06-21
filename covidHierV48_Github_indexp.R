#     covidHierV48_Github_indexp.R runs simulations for model extension using individual-level differences in NPI adherence
#     Copyright (C) 2021  Kathyrn R Fair, Vadim A Karatayev
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

library(abind); library(plyr); library(dplyr); library(data.table);

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

DATA=readRDS("covidHierData.rds"); #See readme file for description of contents
cmat=readRDS("cmat_data_weighted.rds");
regionid=readRDS("PHUtoREGIONlinker_numeric.rds");
agedat=readRDS("ONdat_newKage_to2021-03-02.rds"); ### Contains new age specific known case data for fitting
agedat$tot<-cumsum(rowSums(agedat[,3:7]))
regiondat=readRDS("ONdat_casesbyPHU_to2021-03-01.rds"); ### Contains new region specific known case data for fitting
regiondat$tot<-cumsum(rowSums(regiondat[,2:35]))

omgdat<-readRDS("mobilitydat_REAL_to2021-02-27_V4.rds") #read in mobility data
#Drop all dates where totals are impacted by reporting lags
agedat<-agedat[agedat$Date<="2021-02-28",]
regiondat<-regiondat[regiondat$Date<="2021-02-28",]

#Set code version
codeversion<-"v48"

### Set counterfactual type i.e. schools remain open ("schoolopen"), workplaces remain open ("workopen"), both remain open ("bothopen") or both are shut ("neitheropen", corresponds to what actually occured in the province)
# Set reopening type; with ("restricted") or without ("unrestricted") NPIs in schools/workplaces
cftype<-"neitheropen"; reopeningtype<-"restricted"

#Read in data on best parameter sets from fitting
#Select how many of the top fits you want (i.e. topnum<-10 indicates you want the 10 best parameter sets)
topnum<-10

fitparms00=readRDS(sprintf("simsFITPARMS_%s.rds", codeversion));
bestcutoff<-if(length(unique(fitparms00$LL))>(topnum+1)) {dplyr::nth(sort(unique(fitparms00$LL)), topnum+1)} else {1e8}
fitparms0 <- fitparms00[fitparms00$LL<bestcutoff,]
#Get rid of any dup parm sets (ID by LL)
checker<-fitparms0[!duplicated(fitparms0[ , "LL" ] ), ]
fitparms<- fitparms0[fitparms0$fit %in% checker$fit,]

tstart<-Sys.time()

#args takes commandline argument for which 2% band for % of population who are not adhering NPIs to consider 
#can expand size of band but output data frame size may become prohibitively large, depending on user's RAM
#can be replaced by args <- "value in [0.01,1]" to run w/o command line (e.g. args<-50 will run simulations with 49,50% of population non-compliant)
args <- commandArgs(trailingOnly = TRUE)
print(c((as.numeric(args[1])/100 - 0.02), (as.numeric(args[1])/100)))

goodfits<-1;
for (fit in 1:topnum) {
  
  parms=as.matrix(fitparms[fitparms$fit==unique(fitparms$fit)[fit],1:88]); # use one of our "best" parameter sets
  
  #Some parameter name differences from paper:
  #initial testing tauI_t0=cvTl, final testing tauI_tf=tauI, tauA=tauEAf*tauI
  # a1,a2,a3,a4,a5 are age-specific susceptibility modifiers (gamma_i, i=1,2,3,4,5 in the paper)
  # boost is L0 (impact of stay at home orders on NPI adherence)
  #Mfact scales how much travel happens compared to reality; Mfact=1 is just movement from commuting data
  #epsP allows individual NPI adherence efficacy to differ from that of closures (epsP=eps in paper)
  #Tg and Tl are the gobal and local closure thresholds, n0 is fraction initially infected
  #Msave is the travel matrix. Here Msave entries are numbers of commuters; msim3 converts them to proportions
pops=colSums(DATA$Msave)

###Adjust to 2020 Q4 pop estimate (from https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1710000901), as region pop estimates are from 2016 census
pop2020<-14733119
popratio<-pop2020/sum(pops)

popadj2020=round(popratio*pops)

# Specification of closure/reopening events
inCstart=50; #all events occur at a distance from the day total cases >=50 (March 10th)
inClen_w1=79; #Duration of initial provincial workplace closure (March 25-June 11, Phase 2 commenced June 12)
inClen_s1=178; #Duration of initial provincial school closure (Mar 14 - Sept 7, Schools reopen Sept 8)
inClen_w2=46; #Duration of second provincial workplace closure, to first part of staggered reopening (Dec 26th 2020 - Feb 9th 2021)
inClen_s2=51; #Duration of second provincial school closure, to first part of staggered reopening (Dec 21st 2020 - Feb 9th 2021)

workstaggergap1<-6; #Gap for first 3 reopening stages (Feb 10/16/22)
workstaggergap2<-20; #Gap between Feb 16th and march 8th
schoolstaggergap<-8; #Gap between 2 stages of return to in-person learning (Feb 8/16)
inClen_w2=rep(46+workstaggergap1,49); #days between Dec 26th lockdown and Feb 16th main reopening day (46 days is to the 10th)
inClen_s2=rep(49,49); #days between Dec 21st first day of christmas holiday, feb 8th first reopening day

#Add in mods for workplaces in Hastings Prince Edward (Regions 8,9), Kingston, Frontenac and Lennox & Addington (Regions 6,7), and Renfrew County (Region 38) which reopen 6 days earlier on Feb 10
inClen_w2[c(6,7,8,9,38)]<-(inClen_w2[c(6,7,8,9,38)]-workstaggergap1)
#Add in mod for workplaces in York (14) reopening 6 days later on Feb 22
inClen_w2[c(14)]<-(inClen_w2[c(14)]+workstaggergap1)
#Add in mods for workplaces in Toronto (15), Peel (16), North Bay - Parry Sound (39,40) reopening 20 days later on March 8th
inClen_w2[c(15,16,39,40)]<-(inClen_w2[c(15,16,39,40)]+workstaggergap2)
#Add in mods for schools in Toronto (15), Peel (16), and York (14) reopening 8 days later on feb 16th
inClen_s2[c(14,15,16)]<-(inClen_s2[c(14,15,16)]+schoolstaggergap)

workclosuregap1=14 #Gap between day total cases >= 50 (Mar 10th) and day workplaces closed (March 25th)
schoolclosuregap1=3 #Gap between day total cases >= 50 (Mar 10th) and day schools closed for March Break (March 14th)
workclosuregap2=290 #Gap between day total cases >= 50 (Mar 10th) and day workplaces re-closed (Dec 26th)
schoolclosuregap2=285 #Gap between day total cases >= 50 (Mar 10th) and day schools re-closed (Dec 21st)

omggap=144 #Gap between day total cases >= 50 (Mar 10th) and day we switch to wave 2 omega value (Aug 1st)
Resurge_w=50 #Duration of additional closures 
Resurge_s=50
NP=ncol(parms)-nrow(parms)

#Dummy start values for closures s.t. we can do counterfactuals for first wave
inCstart_w=50;
inCstart_s=50;
if (cftype %in% c("workopen", "bothopen")) {inCstart_w=1e8;} #Set so workplaces never close

years<- 1 
agepopdist<-c(3019645,3475990,3855065,2505545,592260)
ncperc=0 #Set number of compliant susceptibles
susceffect=0.5 #Reduction in VD efficacy if susceptible is not distancing

#Defining model state space. Tn, Tk, Da, and Di are all untested, tested, asymptomatics, and infecteds respectively.
#Nt tracks cumulative # positive cases (including those recovered)
#C tracks # days since last closure, but in output msim converts all positive C entries into eps
#VD indicates level of NPI adherence from individuals
#Sick indicates all individuals who are in any of E, P, A, I (i.e. all those either exposed or currently infectious)
#Define state space for compliant individuals
Susc=c("S1", "S2", "S3", "S4", "S5"); E=c("E1", "E2","E3", "E4", "E5"); R=c("R1", "R2", "R3", "R4","R5");
Da=c("A1","Ak1","SA1","SAk1", "A2","Ak2","SA2","SAk2", "A3","Ak3","SA3","SAk3", "A4","Ak4","SA4","SAk4", "A5","Ak5","SA5","SAk5")
Di=c("I1","Ik1","SI1","SIk1", "I2","Ik2","SI2","SIk2", "I3","Ik3","SI3","SIk3", "I4","Ik4","SI4","SIk4", "I5","Ik5","SI5","SIk5")
Tn=c("A1","SA1","I1","SI1","A2","SA2","I2","SI2","A3","SA3","I3","SI3","A4","SA4","I4","SI4","A5","SA5","I5","SI5")
Tk=c("Ak1","SAk1","Ik1","SIk1","Ak2","SAk2","Ik2","SIk2","Ak3","SAk3","Ik3","SIk3","Ak4","SAk4","Ik4","SIk4","Ak5","SAk5","Ik5","SIk5")
All1<-c("S1", "E1", "R1", "A1","Ak1","SA1","SAk1", "I1","Ik1","SI1","SIk1"); All2<-c("S2", "E2", "R2", "A2","Ak2","SA2","SAk2", "I2","Ik2","SI2","SIk2"); All3<-c("S3", "E3", "R3", "A3","Ak3","SA3","SAk3", "I3","Ik3","SI3","SIk3"); All4<-c("S4", "E4", "R4", "A4","Ak4","SA4","SAk4", "I4","Ik4","SI4","SIk4"); All5<-c("S5", "E5", "R5", "A5","Ak5","SA5","SAk5", "I5","Ik5","SI5","SIk5");
sick1<-c("E1", "A1","Ak1","SA1","SAk1", "I1","Ik1","SI1","SIk1"); sick2<-c("E2", "A2","Ak2","SA2","SAk2", "I2","Ik2","SI2","SIk2"); sick3<-c("E3", "A3","Ak3","SA3","SAk3", "I3","Ik3","SI3","SIk3"); sick4<-c("E4", "A4","Ak4","SA4","SAk4", "I4","Ik4","SI4","SIk4"); sick5<-c("E5", "A5","Ak5","SA5","SAk5", "I5","Ik5","SI5","SIk5");
symp1<-c("I1","Ik1","SI1","SIk1"); symp2<-c("I2","Ik2","SI2","SIk2"); symp3<-c("I3","Ik3","SI3","SIk3"); symp4<-c("I4","Ik4","SI4","SIk4"); symp5<-c("I5","Ik5","SI5","SIk5");
#Define state space for non-compliant individuals
Suscnc=c("S1nc", "S2nc", "S3nc", "S4nc", "S5nc"); Enc=c("E1nc", "E2nc","E3nc", "E4nc", "E5nc"); Rnc=c("R1nc", "R2nc", "R3nc", "R4nc","R5nc");
Danc=c("A1nc","Ak1nc","SA1nc","SAk1nc", "A2nc","Ak2nc","SA2nc","SAk2nc", "A3nc","Ak3nc","SA3nc","SAk3nc", "A4nc","Ak4nc","SA4nc","SAk4nc", "A5nc","Ak5nc","SA5nc","SAk5nc")
Dinc=c("I1nc","Ik1nc","SI1nc","SIk1nc", "I2nc","Ik2nc","SI2nc","SIk2nc", "I3nc","Ik3nc","SI3nc","SIk3nc", "I4nc","Ik4nc","SI4nc","SIk4nc", "I5nc","Ik5nc","SI5nc","SIk5nc")
Tnnc=c("A1nc","SA1nc","I1nc","SI1nc","A2nc","SA2nc","I2nc","SI2nc","A3nc","SA3nc","I3nc","SI3nc","A4nc","SA4nc","I4nc","SI4nc","A5nc","SA5nc","I5nc","SI5nc")
Tknc=c("Ak1nc","SAk1nc","Ik1nc","SIk1nc","Ak2nc","SAk2nc","Ik2nc","SIk2nc","Ak3nc","SAk3nc","Ik3nc","SIk3nc","Ak4nc","SAk4nc","Ik4nc","SIk4nc","Ak5nc","SAk5nc","Ik5nc","SIk5nc")
All1nc<-c("S1nc", "E1nc", "R1nc", "A1nc","Ak1nc","SA1nc","SAk1nc", "I1nc","Ik1nc","SI1nc","SIk1nc"); All2nc<-c("S2nc", "E2nc", "R2nc", "A2nc","Ak2nc","SA2nc","SAk2nc", "I2nc","Ik2nc","SI2nc","SIk2nc"); All3nc<-c("S3nc", "E3nc", "R3nc", "A3nc","Ak3nc","SA3nc","SAk3nc", "I3nc","Ik3nc","SI3nc","SIk3nc"); All4nc<-c("S4nc", "E4nc", "R4nc", "A4nc","Ak4nc","SA4nc","SAk4nc", "I4nc","Ik4nc","SI4nc","SIk4nc"); All5nc<-c("S5nc", "E5nc", "R5nc", "A5nc","Ak5nc","SA5nc","SAk5nc", "I5nc","Ik5nc","SI5nc","SIk5nc");
sick1nc<-c("E1nc", "A1nc","Ak1nc","SA1nc","SAk1nc", "I1nc","Ik1nc","SI1nc","SIk1nc"); sick2nc<-c("E2nc", "A2nc","Ak2nc","SA2nc","SAk2nc", "I2nc","Ik2nc","SI2nc","SIk2nc"); sick3nc<-c("E3nc", "A3nc","Ak3nc","SA3nc","SAk3nc", "I3nc","Ik3nc","SI3nc","SIk3nc"); sick4nc<-c("E4nc", "A4nc","Ak4nc","SA4nc","SAk4nc", "I4nc","Ik4nc","SI4nc","SIk4nc"); sick5nc<-c("E5nc", "A5nc","Ak5nc","SA5nc","SAk5nc", "I5nc","Ik5nc","SI5nc","SIk5nc");
symp1nc<-c("I1nc","Ik1nc","SI1nc","SIk1nc"); symp2nc<-c("I2nc","Ik2nc","SI2nc","SIk2nc"); symp3nc<-c("I3nc","Ik3nc","SI3nc","SIk3nc"); symp4nc<-c("I4nc","Ik4nc","SI4nc","SIk4nc"); symp5nc<-c("I5nc","Ik5nc","SI5nc","SIk5nc");
#Combine to create overall state space
Snm=c(Susc,E,Da,Di,R,"Nt","Cw", "Cs", "C", "Nt1","Nt2","Nt3","Nt4","Nt5", "VD"); Snmnc=c(Suscnc,Enc,Danc,Dinc,Rnc,"Ntnc", "Nt1nc","Nt2nc","Nt3nc","Nt4nc","Nt5nc"); 
Snm.comb=c(Snm,Snmnc, "probcc1", "probcnc1", "probncc1", "probncnc1", "probcc2", "probcnc2", "probncc2", "probncnc2", "probcc3", "probcnc3", "probncc3", "probncnc3", "probcc4", "probcnc4", "probncc4", "probncnc4", "probcc5", "probcnc5", "probncc5", "probncnc5");

#these functions handle all the state transitions. Input x is a matrix where 1st half of columns are source states and 2nd half are destination states
gpTransB=function(x,Prs,seed,nc=ncol(x)){ 
  xvp=cbind(as.vector(x[,1:(nc/2)]),as.vector(Prs))
  if(max(xvp[,1])==0) return(x); nz=(xvp[,1]*xvp[,2])>0; xvp[!nz,1]=0; 
  set.seed(seed); xvp[nz,1]=apply(matrix(xvp[nz,],ncol=2),1,function(y) rbinom(1,y[1],y[2])); 
  return(x+matrix(c(-xvp[,1],xvp[,1]),ncol=nc))
}

#gpTrans is a simplified version where one transition probability applies to all states. 
#If recovery=TRUE have individuals from 1st columns of x all transitioning into the last column
gpTrans=function(x,Pr,seed,Recovery=FALSE, nc=ncol(x)){
  xv=as.vector(x[,1:c(nc/2,nc-1)[1+Recovery]])
  if(max(xv)==0) return(x); set.seed(seed); xv[xv>0]=sapply(xv[xv>0], function(y) rbinom(1,y,Pr));
  if(Recovery){ Tr=matrix(xv,ncol=nc-1); return(x+cbind(-Tr,rowSums(Tr))); }; return(x+matrix(c(-xv,xv),ncol=nc));
}

#Transition probabilities for travel. Multinomial faster than many binomials. rs=TRUE returns just total # people going to each province
reshfl2=function(x,M,seed,rs=TRUE,L=ncol(M),xnm=diag(x)){
  set.seed(seed); if(max(x)>0) xnm[,x>0]=apply(matrix(rbind(x,M)[,x>0],nrow=L+1), 2, function(y) rmultinom(1,y[1],y[-1]));
  if(rs) return(rowSums(xnm)); return(xnm); 
}

#Modifier of travel matrix as people sick and/or tested positive less likely to travel by a proportion pstay
Mstay=function(M,pstay,Mod=M*(-(diag(ncol(M))-1)*(1-pstay))) Mod+diag(1-colSums(Mod))

meansd=function(x,wtR=1+x*0,wts=wtR/sum(wtR)){ xm=weighted.mean(x,wts); return(c(xm,sqrt(sum(wts*(x-xm)^2)))); }

#Main function handling all within-day transitions over state space S
m3iter=function(S,parms,timeinfo,seed,Ns=parms[,"N"]){

  ## Ramp up testing
  # For epi states P, A
  testPA1<- if(timeinfo[2]<0) {parms[,"taumax"]*0;} else {parms[1,"kappa"]*parms[,"taumax"]*(1-exp(-parms[1,"psi1"]*timeinfo[2]));} 
  testPA2<- if(timeinfo[2]<0) {parms[,"taumax"]*0;} else {parms[1,"kappa"]*parms[,"taumax"]*(1-exp(-parms[1,"psi2"]*timeinfo[2]));} 
  testPA3<- if(timeinfo[2]<0) {parms[,"taumax"]*0;} else {parms[1,"kappa"]*parms[,"taumax"]*(1-exp(-parms[1,"psi3"]*timeinfo[2]));} 
  testPA4<- if(timeinfo[2]<0) {parms[,"taumax"]*0;} else {parms[1,"kappa"]*parms[,"taumax"]*(1-exp(-parms[1,"psi4"]*timeinfo[2]));} 
  testPA5<- if(timeinfo[2]<0) {parms[,"taumax"]*0;} else {parms[1,"kappa"]*parms[,"taumax"]*(1-exp(-parms[1,"psi5"]*timeinfo[2]));} 
  # For epi state I
  testI1<- if(timeinfo[2]<0) {parms[,"tau0"];} else {parms[,"taumax"] - (parms[,"taumax"] - parms[,"tau0"])*exp(-parms[1,"psi1"]*timeinfo[2]);}
  testI2<- if(timeinfo[2]<0) {parms[,"tau0"];} else {parms[,"taumax"] - (parms[,"taumax"] - parms[,"tau0"])*exp(-parms[1,"psi2"]*timeinfo[2]);}
  testI3<- if(timeinfo[2]<0) {parms[,"tau0"];} else {parms[,"taumax"] - (parms[,"taumax"] - parms[,"tau0"])*exp(-parms[1,"psi3"]*timeinfo[2]);}
  testI4<- if(timeinfo[2]<0) {parms[,"tau0"];} else {parms[,"taumax"] - (parms[,"taumax"] - parms[,"tau0"])*exp(-parms[1,"psi4"]*timeinfo[2]);}
  testI5<- if(timeinfo[2]<0) {parms[,"tau0"];} else {parms[,"taumax"] - (parms[,"taumax"] - parms[,"tau0"])*exp(-parms[1,"psi5"]*timeinfo[2]);}

  if(testPA1 > testI1 || testPA2 > testI2 || testPA3 > testI3 || testPA4 > testI4 || testPA5 > testI5) {print("problem with testing rates"); invisible(readline(prompt="Press [enter] to continue"));}
  
  #Implement testing (test results coming back from previous day)
  #Compliant individuals
  S0k=S[,Tk]; S0k_1=S[,Tk[1:4]]; S0k_2=S[,Tk[5:8]]; S0k_3=S[,Tk[9:12]]; S0k_4=S[,Tk[13:16]]; S0k_5=S[,Tk[17:20]];
  S0n_1=S[,Tn[1:4]]; S0n_2=S[,Tn[5:8]]; S0n_3=S[,Tn[9:12]]; S0n_4=S[,Tn[13:16]]; S0n_5=S[,Tn[17:20]];
  S[,c(Tn,Tk)]=gpTransB(S[,c(Tn,Tk)],	as.vector(cbind(testPA1,testPA1,testI1,testI1,testPA2,testPA2,testI2,testI2,testPA3,testPA3,testI3,testI3,testPA4,testPA4,testI4,testI4,testPA5,testPA5,testI5,testI5)),seed+40); 
  #Non-compliant individuals
  S0knc=S[,Tknc]; S0knc_1=S[,Tknc[1:4]]; S0knc_2=S[,Tknc[5:8]]; S0knc_3=S[,Tknc[9:12]]; S0knc_4=S[,Tknc[13:16]]; S0knc_5=S[,Tknc[17:20]];
  S0nnc_1=S[,Tnnc[1:4]]; S0nnc_2=S[,Tnnc[5:8]]; S0nnc_3=S[,Tnnc[9:12]]; S0nnc_4=S[,Tnnc[13:16]]; S0nnc_5=S[,Tnnc[17:20]];
  S[,c(Tnnc,Tknc)]=gpTransB(S[,c(Tnnc,Tknc)],	as.vector(cbind(testPA1,testPA1,testI1,testI1,testPA2,testPA2,testI2,testI2,testPA3,testPA3,testI3,testI3,testPA4,testPA4,testI4,testI4,testPA5,testPA5,testI5,testI5)),seed+52); 
  
  #calculate pos, the vector of local active case prevalence
  #Compliant individuals
  S[,"Nt"]=S[,"Nt"]+rowSums(S[,Tk]-S0k);
  S[,"Nt1"]=S[,"Nt1"]+rowSums(S[,Tk[1:4]]-S0k_1);
  S[,"Nt2"]=S[,"Nt2"]+rowSums(S[,Tk[5:8]]-S0k_2);
  S[,"Nt3"]=S[,"Nt3"]+rowSums(S[,Tk[9:12]]-S0k_3);
  S[,"Nt4"]=S[,"Nt4"]+rowSums(S[,Tk[13:16]]-S0k_4);
  S[,"Nt5"]=S[,"Nt5"]+rowSums(S[,Tk[17:20]]-S0k_5);
  #Non-compliant individuals
  S[,"Ntnc"]=S[,"Ntnc"]+rowSums(S[,Tknc]-S0knc);
  S[,"Nt1nc"]=S[,"Nt1nc"]+rowSums(S[,Tknc[1:4]]-S0knc_1);
  S[,"Nt2nc"]=S[,"Nt2nc"]+rowSums(S[,Tknc[5:8]]-S0knc_2);
  S[,"Nt3nc"]=S[,"Nt3nc"]+rowSums(S[,Tknc[9:12]]-S0knc_3);
  S[,"Nt4nc"]=S[,"Nt4nc"]+rowSums(S[,Tknc[13:16]]-S0knc_4);
  S[,"Nt5nc"]=S[,"Nt5nc"]+rowSums(S[,Tknc[17:20]]-S0knc_5);
  
  #Calculate total active cases across compliant and non-compliant individuals
  pos=rowSums(S[,c(Tk, Tknc)])/Ns; posglobal=(sum(S[,c(Tk,Tknc)])/sum(Ns));

  #Disease progression: zeta is the fraction of people who never show symptoms (modeled implicitly)
  zeta=0.2;
  ######## Compliant individuals
  ### Symptomatic removals by age class
  S[,c(Di[1:4],R[1])]=gpTrans(S[,c(Di[1:4],R[1])],parms[1,"rho"],seed,TRUE); 
  S[,c(Di[5:8],R[2])]=gpTrans(S[,c(Di[5:8],R[2])],parms[1,"rho"],seed+1,TRUE); 
  S[,c(Di[9:12],R[3])]=gpTrans(S[,c(Di[9:12],R[3])],parms[1,"rho"],seed+2,TRUE); 
  S[,c(Di[13:16],R[4])]=gpTrans(S[,c(Di[13:16],R[4])],parms[1,"rho"],seed+3,TRUE); 
  S[,c(Di[17:20],R[5])]=gpTrans(S[,c(Di[17:20],R[5])],parms[1,"rho"],seed+4,TRUE); 
  ### Asymptomatic removals by age class
  S[,c(Da[1:4],R[1])]=gpTrans(S[,c(Da[1:4],R[1])],zeta*prod(parms[1,c("sig","rho")]),seed,TRUE); 
  S[,c(Da[5:8],R[2])]=gpTrans(S[,c(Da[5:8],R[2])],zeta*prod(parms[1,c("sig","rho")]),seed+1,TRUE); 
  S[,c(Da[9:12],R[3])]=gpTrans(S[,c(Da[9:12],R[3])],zeta*prod(parms[1,c("sig","rho")]),seed+2,TRUE); 
  S[,c(Da[13:16],R[4])]=gpTrans(S[,c(Da[13:16],R[4])],zeta*prod(parms[1,c("sig","rho")]),seed+3,TRUE); 
  S[,c(Da[17:20],R[5])]=gpTrans(S[,c(Da[17:20],R[5])],zeta*prod(parms[1,c("sig","rho")]),seed+4,TRUE); 
  ### Transition from pre-symptomatic to symptomatic by age class
  S[,c(Da[1:4],Di[1:4])]=gpTrans(S[,c(Da[1:4],Di[1:4])],parms[1,"sig"]*(1-zeta),seed+5); 
  S[,c(Da[5:8],Di[5:8])]=gpTrans(S[,c(Da[5:8],Di[5:8])],parms[1,"sig"]*(1-zeta),seed+6); 
  S[,c(Da[9:12],Di[9:12])]=gpTrans(S[,c(Da[9:12],Di[9:12])],parms[1,"sig"]*(1-zeta),seed+7); 
  S[,c(Da[13:16],Di[13:16])]=gpTrans(S[,c(Da[13:16],Di[13:16])],parms[1,"sig"]*(1-zeta),seed+8); 
  S[,c(Da[17:20],Di[17:20])]=gpTrans(S[,c(Da[17:20],Di[17:20])],parms[1,"sig"]*(1-zeta),seed+9); 
  ### Shift from exposed to asymptomatic by age class
  S[,c("E1","A1")]=gpTrans(S[,c("E1","A1")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+10); 
  S[,c("E2","A2")]=gpTrans(S[,c("E2","A2")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+11); 
  S[,c("E3","A3")]=gpTrans(S[,c("E3","A3")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+12); 
  S[,c("E4","A4")]=gpTrans(S[,c("E4","A4")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+13); 
  S[,c("E5","A5")]=gpTrans(S[,c("E5","A5")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+14); 
  ### Shift from exposed to superspreader asymptomatic by age class
  S[,c("E1","SA1")]=gpTrans(S[,c("E1","SA1")],parms[1,"alpha"]*parms[1,"s"],seed+15);
  S[,c("E2","SA2")]=gpTrans(S[,c("E2","SA2")],parms[1,"alpha"]*parms[1,"s"],seed+16);
  S[,c("E3","SA3")]=gpTrans(S[,c("E3","SA3")],parms[1,"alpha"]*parms[1,"s"],seed+17);
  S[,c("E4","SA4")]=gpTrans(S[,c("E4","SA4")],parms[1,"alpha"]*parms[1,"s"],seed+18);
  S[,c("E5","SA5")]=gpTrans(S[,c("E5","SA5")],parms[1,"alpha"]*parms[1,"s"],seed+19);
  
  ######## Non-Compliant individuals
  ### Symptomatic removals by age class
  S[,c(Dinc[1:4],Rnc[1])]=gpTrans(S[,c(Dinc[1:4],Rnc[1])],parms[1,"rho"],seed+53,TRUE); 
  S[,c(Dinc[5:8],Rnc[2])]=gpTrans(S[,c(Dinc[5:8],Rnc[2])],parms[1,"rho"],seed+54,TRUE); 
  S[,c(Dinc[9:12],Rnc[3])]=gpTrans(S[,c(Dinc[9:12],Rnc[3])],parms[1,"rho"],seed+55,TRUE); 
  S[,c(Dinc[13:16],Rnc[4])]=gpTrans(S[,c(Dinc[13:16],Rnc[4])],parms[1,"rho"],seed+56,TRUE); 
  S[,c(Dinc[17:20],Rnc[5])]=gpTrans(S[,c(Dinc[17:20],Rnc[5])],parms[1,"rho"],seed+57,TRUE); 
  ### Asymptomatic removals by age class
  S[,c(Danc[1:4],Rnc[1])]=gpTrans(S[,c(Danc[1:4],Rnc[1])],zeta*prod(parms[1,c("sig","rho")]),seed+53,TRUE); 
  S[,c(Danc[5:8],Rnc[2])]=gpTrans(S[,c(Danc[5:8],Rnc[2])],zeta*prod(parms[1,c("sig","rho")]),seed+54,TRUE); 
  S[,c(Danc[9:12],Rnc[3])]=gpTrans(S[,c(Danc[9:12],Rnc[3])],zeta*prod(parms[1,c("sig","rho")]),seed+55,TRUE); 
  S[,c(Danc[13:16],Rnc[4])]=gpTrans(S[,c(Danc[13:16],Rnc[4])],zeta*prod(parms[1,c("sig","rho")]),seed+56,TRUE); 
  S[,c(Danc[17:20],Rnc[5])]=gpTrans(S[,c(Danc[17:20],Rnc[5])],zeta*prod(parms[1,c("sig","rho")]),seed+57,TRUE); 
  ### Transition from pre-symptomatic to symptomatic by age class
  S[,c(Danc[1:4],Dinc[1:4])]=gpTrans(S[,c(Danc[1:4],Dinc[1:4])],parms[1,"sig"]*(1-zeta),seed+58); 
  S[,c(Danc[5:8],Dinc[5:8])]=gpTrans(S[,c(Danc[5:8],Dinc[5:8])],parms[1,"sig"]*(1-zeta),seed+59); 
  S[,c(Danc[9:12],Dinc[9:12])]=gpTrans(S[,c(Danc[9:12],Dinc[9:12])],parms[1,"sig"]*(1-zeta),seed+60); 
  S[,c(Danc[13:16],Dinc[13:16])]=gpTrans(S[,c(Danc[13:16],Dinc[13:16])],parms[1,"sig"]*(1-zeta),seed+61); 
  S[,c(Danc[17:20],Dinc[17:20])]=gpTrans(S[,c(Danc[17:20],Dinc[17:20])],parms[1,"sig"]*(1-zeta),seed+62); 
  ### Shift from exposed to asymptomatic by age class
  S[,c("E1nc","A1nc")]=gpTrans(S[,c("E1nc","A1nc")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+63); 
  S[,c("E2nc","A2nc")]=gpTrans(S[,c("E2nc","A2nc")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+64); 
  S[,c("E3nc","A3nc")]=gpTrans(S[,c("E3nc","A3nc")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+65); 
  S[,c("E4nc","A4nc")]=gpTrans(S[,c("E4nc","A4nc")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+66); 
  S[,c("E5nc","A5nc")]=gpTrans(S[,c("E5nc","A5nc")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+67); 
  ### Shift from exposed to superspreader asymptomatic by age class
  S[,c("E1nc","SA1nc")]=gpTrans(S[,c("E1nc","SA1nc")],parms[1,"alpha"]*parms[1,"s"],seed+68);
  S[,c("E2nc","SA2nc")]=gpTrans(S[,c("E2nc","SA2nc")],parms[1,"alpha"]*parms[1,"s"],seed+19);
  S[,c("E3nc","SA3nc")]=gpTrans(S[,c("E3nc","SA3nc")],parms[1,"alpha"]*parms[1,"s"],seed+70);
  S[,c("E4nc","SA4nc")]=gpTrans(S[,c("E4nc","SA4nc")],parms[1,"alpha"]*parms[1,"s"],seed+71);
  S[,c("E5nc","SA5nc")]=gpTrans(S[,c("E5nc","SA5nc")],parms[1,"alpha"]*parms[1,"s"],seed+72);
  
  # Conditional statements for opening/closing workplaces and schools
  if(cftype %in% c("workopen", "neitheropen")){
    if(timeinfo[3]==-1 || timeinfo[2]<schoolclosuregap1) closed_s = 0*parms[,"eps_s"]; #We have not yet reached >=50 cases (prior to Mar 10) or >=50 cases reached, schools not yet closed (Mar 10 - Mar 13)
    if(timeinfo[3]>=schoolclosuregap1 && timeinfo[3]<schoolclosuregap1+inClen_s1) closed_s = parms[,"eps_s"]; #Initial school closure started (March 14-Sept 7)
    if(timeinfo[3]>=schoolclosuregap1+inClen_s1 && timeinfo[3]< schoolclosuregap2  && reopeningtype=="restricted") closed_s =parms[1,"eps_s"]*parms[,"reopen_s"]; #Schools reopen with covid regs in place
    if(timeinfo[3]>=schoolclosuregap1+inClen_s1 && timeinfo[3]< schoolclosuregap2  && reopeningtype=="unrestricted") closed_s =0*parms[,"eps_s"]; #Schools reopen WITHOUT covid regs in place
    
    if (timeinfo[3]>=schoolclosuregap2) {
      if (reopeningtype=="restricted") {closed_s <- ifelse(timeinfo[3]<schoolclosuregap2+inClen_s2, parms[,"eps_s"], parms[1,"eps_s"]*parms[,"reopen_s"]);} #Christmas holiday to Boxing day lockdown to stay at home, with staggered reopening
      if (reopeningtype=="unrestricted") {closed_s <- ifelse(timeinfo[3]<schoolclosuregap2+inClen_s2, parms[,"eps_s"], 0*parms[,"reopen_s"]);} #Christmas holiday to Boxing day lockdown to stay at home, with staggered reopening
    }
  }
  
  if(cftype %in% c("schoolopen", "bothopen")){
    if(timeinfo[3]==-1 || timeinfo[2]<schoolbreakgap) closed_s = 0*parms[,"eps_s"]; #We have not yet reached >=50 cases (prior to Mar 10) or >=50 cases reached, schools not yet closed (Mar 10 - June 30)
    if(timeinfo[3]>=schoolbreakgap && timeinfo[3]<schoolbreakgap+inBlen_s) closed_s = parms[,"eps_s"]; #Initial summer school break started (July 1st-Sept 8th)
    if(timeinfo[3]>=schoolbreakgap+inBlen_s && timeinfo[3]< schoolclosuregap2 && reopeningtype=="restricted") closed_s =parms[1,"eps_s"]*parms[,"reopen_s"]; #Schools reopen with covid regs in place
    if(timeinfo[3]>=schoolbreakgap+inBlen_s && timeinfo[3]< schoolclosuregap2 && reopeningtype=="unrestricted") closed_s =0*parms[,"eps_s"]; #Schools reopen WITHOUT covid regs in place
    
    if (timeinfo[3]>=schoolclosuregap2) {
      if (reopeningtype=="restricted") {closed_s <- ifelse(timeinfo[3]<schoolclosuregap2+inClen_s2, parms[,"eps_s"], parms[1,"eps_s"]*parms[,"reopen_s"]);} #Christmas holiday to Boxing day lockdown to stay at home, with staggered reopening
      if (reopeningtype=="unrestricted") {closed_s <- ifelse(timeinfo[3]<schoolclosuregap2+inClen_s2, parms[,"eps_s"], 0*parms[,"reopen_s"]);} #Christmas holiday to Boxing day lockdown to stay at home, with staggered reopening
    }
    
  }

  if(timeinfo[4]==-1 || timeinfo[2]<workclosuregap1) closed_w = 0*parms[,"eps_w"]; #We have not yet reached >=50 cases (prior to Mar 10) or >=50 cases reached, workplaces not yet closed (Mar 10 - Mar 24)
  if(timeinfo[4]>=workclosuregap1 && timeinfo[4]<workclosuregap1+inClen_w1) {closed_w = parms[,"eps_w"];}#Initial workplace closure started (March 25-June 11th)
  if(timeinfo[4]>=workclosuregap1+inClen_w1 && timeinfo[4]< workclosuregap2 && reopeningtype=="restricted") closed_w =parms[1,"eps_w"]*parms[,"reopen_w"]; #Workplaces reopen with covid regs in place
  if(timeinfo[4]>=workclosuregap1+inClen_w1 && timeinfo[4]< workclosuregap2 && reopeningtype=="unrestricted") closed_w = 0*parms[,"eps_w"]; #Workplaces reopen WITHOUT covid regs in place
  
  if (timeinfo[4]>=workclosuregap2) {
    if (reopeningtype=="restricted") {closed_w <- ifelse(timeinfo[4]<workclosuregap2+inClen_w2, parms[,"eps_w"], parms[1,"eps_w"]*parms[,"reopen_w"]);} #Boxing day lockdown to stay at home, with staggered reopening
    if (reopeningtype=="unrestricted") {closed_w <- ifelse(timeinfo[4]<workclosuregap2+inClen_w2, parms[,"eps_w"], 0*parms[,"reopen_w"]);} #Boxing day lockdown to stay at home, with staggered reopening
  }
  
  S[,"Cs"] = closed_s;
  S[,"Cw"] = closed_w;

  #Make modified travel matrices for people feeling sick and/or tested positive, then implement travel
  M=Mc=parms[,-(1:NP)]; Mc=M[1:nrow(M),]*(1-closed_w)*(1-closed_s); diag(Mc)=diag(Mc)+1-colSums(Mc);
  ### Revamping McA to include age specific travel rates, first of each pair is for old/young, second for middle
  McA=abind(Mstay(Mc,parms[1,"atravel"]), Mc, Mstay(Mc,1-(1-parms[1,"atravel"])*(1-parms[1,"eta"])), Mstay(Mc,parms[1,"eta"]), Mstay(Mc,1-(1-parms[1,"atravel"])*(1-parms[1,"r"])), Mstay(Mc,parms[1,"r"]), Mstay(Mc,1-(1-parms[1,"eta"])*(1-parms[1,"r"])*(1-parms[1,"atravel"])), Mstay(Mc,1-(1-parms[1,"eta"])*(1-parms[1,"r"])), along=3);

  #Implement travel
  ### Do Ss for however many age classes you have, added atravel for age-class specific travel
  ##Ss for compliant susceptibles
  Ss=abind(reshfl2(S[,"S1"],Mstay(Mc,parms[1,"atravel"]),seed+20,FALSE),reshfl2(S[,"S2"],Mc,seed+21,FALSE),reshfl2(S[,"S3"],Mc,seed+22,FALSE),reshfl2(S[,"S4"],Mstay(Mc,parms[1,"atravel"]),seed+23,FALSE),reshfl2(S[,"S5"],Mstay(Mc,parms[1,"atravel"]),seed+24,FALSE),along=3);
  ##Ss for non-compliant susceptibles
  Ssnc=abind(reshfl2(S[,"S1nc"],Mstay(Mc,parms[1,"atravel"]),seed+73,FALSE),reshfl2(S[,"S2nc"],Mc,seed+74,FALSE),reshfl2(S[,"S3nc"],Mc,seed+75,FALSE),reshfl2(S[,"S4nc"],Mstay(Mc,parms[1,"atravel"]),seed+76,FALSE),reshfl2(S[,"S5nc"],Mstay(Mc,parms[1,"atravel"]),seed+77,FALSE),along=3);
  #Shuffle all other compliant compartments
  Rearr=apply(rbind(seed+(25:34), c(c(1,2,2,1,1),c(1,2,2,1,1),c(1,3,1,3,2,4,2,4,2,4,2,4,1,3,1,3,1,3,1,3),c(5,7,5,7,6,8,6,8,6,8,6,8,5,7,5,7,5,7,5,7)), S[,c(R,E,Da,Di)]), 2, function(x) reshfl2(x[-(1:2)],McA[,,x[2]],x[1]))
  #Shuffle all other non-compliant compartments
  Rearrnc=apply(rbind(seed+(78:87), c(c(1,2,2,1,1),c(1,2,2,1,1),c(1,3,1,3,2,4,2,4,2,4,2,4,1,3,1,3,1,3,1,3),c(5,7,5,7,6,8,6,8,6,8,6,8,5,7,5,7,5,7,5,7)), S[,c(Rnc,Enc,Danc,Dinc)]), 2, function(x) reshfl2(x[-(1:2)],McA[,,x[2]],x[1]))
  
  #Age specific contacts, rows are age classes, columns are their contact age classes
  ageSpecifics_w<-cmat[,,1] ## age specific contacts for work
  ageSpecifics_s<-cmat[,,2] ## age specific contacts for school
  ageSpecifics_h<-cmat[,,3] ## age specific contacts for home
  ageSpecifics_o<-cmat[,,4] ## age specific contacts for other
  

  #Seasonal forcing component
  scomp<-1+parms[,"B"]*cos((2*pi/365)*(timeinfo[1] + parms[,"phi"]))

  S[, "VD"]<-0
  
  #Infection probability for S compliant/infectious compliant interaction
  Infect1cc = scomp*parms[,"a1"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[1,]) + (1-closed_s)%*%t(ageSpecifics_s[1,]) + (1-parms[,"eps_h"])%*%t(ageSpecifics_h[1,]) + (1-parms[,"eps_o"])%*%t(ageSpecifics_o[1,]))
  Infect2cc = scomp*parms[,"a2"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[2,]) + (1-closed_s)%*%t(ageSpecifics_s[2,]) + (1-parms[,"eps_h"])%*%t(ageSpecifics_h[2,]) + (1-parms[,"eps_o"])%*%t(ageSpecifics_o[2,]))
  Infect3cc = scomp*parms[,"a3"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[3,]) + (1-closed_s)%*%t(ageSpecifics_s[3,]) + (1-parms[,"eps_h"])%*%t(ageSpecifics_h[3,]) + (1-parms[,"eps_o"])%*%t(ageSpecifics_o[3,]))
  Infect4cc = scomp*parms[,"a4"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[4,]) + (1-closed_s)%*%t(ageSpecifics_s[4,]) + (1-parms[,"eps_h"])%*%t(ageSpecifics_h[4,]) + (1-parms[,"eps_o"])%*%t(ageSpecifics_o[4,]))
  Infect5cc = scomp*parms[,"a5"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[5,]) + (1-closed_s)%*%t(ageSpecifics_s[5,]) + (1-parms[,"eps_h"])%*%t(ageSpecifics_h[5,]) + (1-parms[,"eps_o"])%*%t(ageSpecifics_o[5,]))
  #Infection probability for S non-compliant/infectious compliant interaction
  Infect1ncc = scomp*parms[,"a1"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[1,]) + (1-closed_s)%*%t(ageSpecifics_s[1,]) + (1-((1-susceffect)*parms[,"eps_h"]))%*%t(ageSpecifics_h[1,]) + (1-((1-susceffect)*parms[,"eps_o"]))%*%t(ageSpecifics_o[1,]))
  Infect2ncc = scomp*parms[,"a2"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[2,]) + (1-closed_s)%*%t(ageSpecifics_s[2,]) + (1-((1-susceffect)*parms[,"eps_h"]))%*%t(ageSpecifics_h[2,]) + (1-((1-susceffect)*parms[,"eps_o"]))%*%t(ageSpecifics_o[2,]))
  Infect3ncc = scomp*parms[,"a3"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[3,]) + (1-closed_s)%*%t(ageSpecifics_s[3,]) + (1-((1-susceffect)*parms[,"eps_h"]))%*%t(ageSpecifics_h[3,]) + (1-((1-susceffect)*parms[,"eps_o"]))%*%t(ageSpecifics_o[3,]))
  Infect4ncc = scomp*parms[,"a4"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[4,]) + (1-closed_s)%*%t(ageSpecifics_s[4,]) + (1-((1-susceffect)*parms[,"eps_h"]))%*%t(ageSpecifics_h[4,]) + (1-((1-susceffect)*parms[,"eps_o"]))%*%t(ageSpecifics_o[4,]))
  Infect5ncc = scomp*parms[,"a5"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[5,]) + (1-closed_s)%*%t(ageSpecifics_s[5,]) + (1-((1-susceffect)*parms[,"eps_h"]))%*%t(ageSpecifics_h[5,]) + (1-((1-susceffect)*parms[,"eps_o"]))%*%t(ageSpecifics_o[5,]))
  #Infection probability for S compliant/infectious non-compliant interaction
  Infect1cnc = scomp*parms[,"a1"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[1,]) + (1-closed_s)%*%t(ageSpecifics_s[1,])+ (1-(susceffect*parms[,"eps_h"]))%*%t(ageSpecifics_h[1,]) + (1-(susceffect*parms[,"eps_o"]))%*%t(ageSpecifics_o[1,]))
  Infect2cnc = scomp*parms[,"a2"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[2,]) + (1-closed_s)%*%t(ageSpecifics_s[2,])+ (1-(susceffect*parms[,"eps_h"]))%*%t(ageSpecifics_h[2,]) + (1-(susceffect*parms[,"eps_o"]))%*%t(ageSpecifics_o[2,]))
  Infect3cnc = scomp*parms[,"a3"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[3,]) + (1-closed_s)%*%t(ageSpecifics_s[3,])+ (1-(susceffect*parms[,"eps_h"]))%*%t(ageSpecifics_h[3,]) + (1-(susceffect*parms[,"eps_o"]))%*%t(ageSpecifics_o[3,]))
  Infect4cnc = scomp*parms[,"a4"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[4,]) + (1-closed_s)%*%t(ageSpecifics_s[4,])+ (1-(susceffect*parms[,"eps_h"]))%*%t(ageSpecifics_h[4,]) + (1-(susceffect*parms[,"eps_o"]))%*%t(ageSpecifics_o[4,]))
  Infect5cnc = scomp*parms[,"a5"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[5,]) + (1-closed_s)%*%t(ageSpecifics_s[5,])+ (1-(susceffect*parms[,"eps_h"]))%*%t(ageSpecifics_h[5,]) + (1-(susceffect*parms[,"eps_o"]))%*%t(ageSpecifics_o[5,]))
  #Infection probability for S non-compliant/infectious non-compliant interaction
  Infect1ncnc = scomp*parms[,"a1"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[1,]) + (1-closed_s)%*%t(ageSpecifics_s[1,]) + (1-(0*parms[,"eps_h"]))%*%t(ageSpecifics_h[1,]) + (1-(0*parms[,"eps_o"]))%*%t(ageSpecifics_o[1,]))
  Infect2ncnc = scomp*parms[,"a2"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[2,]) + (1-closed_s)%*%t(ageSpecifics_s[2,]) + (1-(0*parms[,"eps_h"]))%*%t(ageSpecifics_h[2,]) + (1-(0*parms[,"eps_o"]))%*%t(ageSpecifics_o[2,]))
  Infect3ncnc = scomp*parms[,"a3"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[3,]) + (1-closed_s)%*%t(ageSpecifics_s[3,]) + (1-(0*parms[,"eps_h"]))%*%t(ageSpecifics_h[3,]) + (1-(0*parms[,"eps_o"]))%*%t(ageSpecifics_o[3,]))
  Infect4ncnc = scomp*parms[,"a4"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[4,]) + (1-closed_s)%*%t(ageSpecifics_s[4,]) + (1-(0*parms[,"eps_h"]))%*%t(ageSpecifics_h[4,]) + (1-(0*parms[,"eps_o"]))%*%t(ageSpecifics_o[4,]))
  Infect5ncnc = scomp*parms[,"a5"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[5,]) + (1-closed_s)%*%t(ageSpecifics_s[5,]) + (1-(0*parms[,"eps_h"]))%*%t(ageSpecifics_h[5,]) + (1-(0*parms[,"eps_o"]))%*%t(ageSpecifics_o[5,]))
  
  #Class-specific infection modifiers
  modK=1-parms[1,"eta"]; modS=1/parms[1,"s"]-1; 
  
  #Introduce age specific modA 3-D arrays for not superspreader or known, not supersrpeader but known, superspreader not known, superspreader and known
  modA1cc=abind(Infect1cc,modK*Infect1cc,modS*Infect1cc,modS*modK*Infect1cc,along=3);
  modA2cc=abind(Infect2cc,modK*Infect2cc,modS*Infect2cc,modS*modK*Infect2cc,along=3);
  modA3cc=abind(Infect3cc,modK*Infect3cc,modS*Infect3cc,modS*modK*Infect3cc,along=3);
  modA4cc=abind(Infect4cc,modK*Infect4cc,modS*Infect4cc,modS*modK*Infect4cc,along=3);
  modA5cc=abind(Infect5cc,modK*Infect5cc,modS*Infect5cc,modS*modK*Infect5cc,along=3);
  
  modA1ncc=abind(Infect1ncc,modK*Infect1ncc,modS*Infect1ncc,modS*modK*Infect1ncc,along=3);
  modA2ncc=abind(Infect2ncc,modK*Infect2ncc,modS*Infect2ncc,modS*modK*Infect2ncc,along=3);
  modA3ncc=abind(Infect3ncc,modK*Infect3ncc,modS*Infect3ncc,modS*modK*Infect3ncc,along=3);
  modA4ncc=abind(Infect4ncc,modK*Infect4ncc,modS*Infect4ncc,modS*modK*Infect4ncc,along=3);
  modA5ncc=abind(Infect5ncc,modK*Infect5ncc,modS*Infect5ncc,modS*modK*Infect5ncc,along=3);
  
  modA1cnc=abind(Infect1cnc,modK*Infect1cnc,modS*Infect1cnc,modS*modK*Infect1cnc,along=3);
  modA2cnc=abind(Infect2cnc,modK*Infect2cnc,modS*Infect2cnc,modS*modK*Infect2cnc,along=3);
  modA3cnc=abind(Infect3cnc,modK*Infect3cnc,modS*Infect3cnc,modS*modK*Infect3cnc,along=3);
  modA4cnc=abind(Infect4cnc,modK*Infect4cnc,modS*Infect4cnc,modS*modK*Infect4cnc,along=3);
  modA5cnc=abind(Infect5cnc,modK*Infect5cnc,modS*Infect5cnc,modS*modK*Infect5cnc,along=3);
  
  modA1ncnc=abind(Infect1ncnc,modK*Infect1ncnc,modS*Infect1ncnc,modS*modK*Infect1ncnc,along=3);
  modA2ncnc=abind(Infect2ncnc,modK*Infect2ncnc,modS*Infect2ncnc,modS*modK*Infect2ncnc,along=3);
  modA3ncnc=abind(Infect3ncnc,modK*Infect3ncnc,modS*Infect3ncnc,modS*modK*Infect3ncnc,along=3);
  modA4ncnc=abind(Infect4ncnc,modK*Infect4ncnc,modS*Infect4ncnc,modS*modK*Infect4ncnc,along=3);
  modA5ncnc=abind(Infect5ncnc,modK*Infect5ncnc,modS*Infect5ncnc,modS*modK*Infect5ncnc,along=3);
 
  #Flatten into 2-D arrays where order is based on age class to match Rearr which is ordered as (all pre/asympt by age class, all sympt by age class)
  modAbind1cc<-cbind(modA1cc[,1,],modA1cc[,2,],modA1cc[,3,],modA1cc[,4,],modA1cc[,5,]);
  modAbind2cc<-cbind(modA2cc[,1,],modA2cc[,2,],modA2cc[,3,],modA2cc[,4,],modA2cc[,5,]);
  modAbind3cc<-cbind(modA3cc[,1,],modA3cc[,2,],modA3cc[,3,],modA3cc[,4,],modA3cc[,5,]);
  modAbind4cc<-cbind(modA4cc[,1,],modA4cc[,2,],modA4cc[,3,],modA4cc[,4,],modA4cc[,5,]);
  modAbind5cc<-cbind(modA5cc[,1,],modA5cc[,2,],modA5cc[,3,],modA5cc[,4,],modA5cc[,5,]);
  
  modAbind1ncc<-cbind(modA1ncc[,1,],modA1ncc[,2,],modA1ncc[,3,],modA1ncc[,4,],modA1ncc[,5,]);
  modAbind2ncc<-cbind(modA2ncc[,1,],modA2ncc[,2,],modA2ncc[,3,],modA2ncc[,4,],modA2ncc[,5,]);
  modAbind3ncc<-cbind(modA3ncc[,1,],modA3ncc[,2,],modA3ncc[,3,],modA3ncc[,4,],modA3ncc[,5,]);
  modAbind4ncc<-cbind(modA4ncc[,1,],modA4ncc[,2,],modA4ncc[,3,],modA4ncc[,4,],modA4ncc[,5,]);
  modAbind5ncc<-cbind(modA5ncc[,1,],modA5ncc[,2,],modA5ncc[,3,],modA5ncc[,4,],modA5ncc[,5,]);
  
  modAbind1cnc<-cbind(modA1cnc[,1,],modA1cnc[,2,],modA1cnc[,3,],modA1cnc[,4,],modA1cnc[,5,]);
  modAbind2cnc<-cbind(modA2cnc[,1,],modA2cnc[,2,],modA2cnc[,3,],modA2cnc[,4,],modA2cnc[,5,]);
  modAbind3cnc<-cbind(modA3cnc[,1,],modA3cnc[,2,],modA3cnc[,3,],modA3cnc[,4,],modA3cnc[,5,]);
  modAbind4cnc<-cbind(modA4cnc[,1,],modA4cnc[,2,],modA4cnc[,3,],modA4cnc[,4,],modA4cnc[,5,]);
  modAbind5cnc<-cbind(modA5cnc[,1,],modA5cnc[,2,],modA5cnc[,3,],modA5cnc[,4,],modA5cnc[,5,]);
  
  modAbind1ncnc<-cbind(modA1ncnc[,1,],modA1ncnc[,2,],modA1ncnc[,3,],modA1ncnc[,4,],modA1ncnc[,5,]);
  modAbind2ncnc<-cbind(modA2ncnc[,1,],modA2ncnc[,2,],modA2ncnc[,3,],modA2ncnc[,4,],modA2ncnc[,5,]);
  modAbind3ncnc<-cbind(modA3ncnc[,1,],modA3ncnc[,2,],modA3ncnc[,3,],modA3ncnc[,4,],modA3ncnc[,5,]);
  modAbind4ncnc<-cbind(modA4ncnc[,1,],modA4ncnc[,2,],modA4ncnc[,3,],modA4ncnc[,4,],modA4ncnc[,5,]);
  modAbind5ncnc<-cbind(modA5ncnc[,1,],modA5ncnc[,2,],modA5ncnc[,3,],modA5ncnc[,4,],modA5ncnc[,5,]);
  
  noprobcc1<-apply((1-cbind(modAbind1cc, 2*modAbind1cc))^Rearr[,-(1:10)], 1, prod);
  noprobcnc1<-apply((1-cbind(modAbind1cnc, 2*modAbind1cnc))^Rearrnc[,-(1:10)], 1, prod);
  noprobncc1<-apply((1-cbind(modAbind1ncc, 2*modAbind1ncc))^Rearr[,-(1:10)], 1, prod);
  noprobncnc1<-apply((1-cbind(modAbind1ncnc, 2*modAbind1ncnc))^Rearrnc[,-(1:10)], 1, prod);
  
  noprobcc2<-apply((1-cbind(modAbind2cc, 2*modAbind2cc))^Rearr[,-(1:10)], 1, prod);
  noprobcnc2<-apply((1-cbind(modAbind2cnc, 2*modAbind2cnc))^Rearrnc[,-(1:10)], 1, prod);
  noprobncc2<-apply((1-cbind(modAbind2ncc, 2*modAbind2ncc))^Rearr[,-(1:10)], 1, prod);
  noprobncnc2<-apply((1-cbind(modAbind2ncnc, 2*modAbind2ncnc))^Rearrnc[,-(1:10)], 1, prod);
  
  noprobcc3<-apply((1-cbind(modAbind3cc, 2*modAbind3cc))^Rearr[,-(1:10)], 1, prod);
  noprobcnc3<-apply((1-cbind(modAbind3cnc, 2*modAbind3cnc))^Rearrnc[,-(1:10)], 1, prod);
  noprobncc3<-apply((1-cbind(modAbind3ncc, 2*modAbind3ncc))^Rearr[,-(1:10)], 1, prod);
  noprobncnc3<-apply((1-cbind(modAbind3ncnc, 2*modAbind3ncnc))^Rearrnc[,-(1:10)], 1, prod);
  
  noprobcc4<-apply((1-cbind(modAbind4cc, 2*modAbind4cc))^Rearr[,-(1:10)], 1, prod);
  noprobcnc4<-apply((1-cbind(modAbind4cnc, 2*modAbind4cnc))^Rearrnc[,-(1:10)], 1, prod);
  noprobncc4<-apply((1-cbind(modAbind4ncc, 2*modAbind4ncc))^Rearr[,-(1:10)], 1, prod);
  noprobncnc4<-apply((1-cbind(modAbind4ncnc, 2*modAbind4ncnc))^Rearrnc[,-(1:10)], 1, prod);
  
  noprobcc5<-apply((1-cbind(modAbind5cc, 2*modAbind5cc))^Rearr[,-(1:10)], 1, prod);
  noprobcnc5<-apply((1-cbind(modAbind5cnc, 2*modAbind5cnc))^Rearrnc[,-(1:10)], 1, prod);
  noprobncc5<-apply((1-cbind(modAbind5ncc, 2*modAbind5ncc))^Rearr[,-(1:10)], 1, prod);
  noprobncnc5<-apply((1-cbind(modAbind5ncnc, 2*modAbind5ncnc))^Rearrnc[,-(1:10)], 1, prod);
  
  S[,"probcc1"]<-1-noprobcc1; S[,"probcnc1"]<-1-noprobcnc1; S[,"probncc1"]<-1-noprobncc1; S[,"probncnc1"]<-1-noprobncnc1;
  S[,"probcc2"]<-1-noprobcc2; S[,"probcnc2"]<-1-noprobcnc2; S[,"probncc2"]<-1-noprobncc2; S[,"probncnc2"]<-1-noprobncnc2;
  S[,"probcc3"]<-1-noprobcc3; S[,"probcnc3"]<-1-noprobcnc3; S[,"probncc3"]<-1-noprobncc3; S[,"probncnc3"]<-1-noprobncnc3;
  S[,"probcc4"]<-1-noprobcc4; S[,"probcnc4"]<-1-noprobcnc4; S[,"probncc4"]<-1-noprobncc4; S[,"probncnc4"]<-1-noprobncnc4;
  S[,"probcc5"]<-1-noprobcc5; S[,"probcnc5"]<-1-noprobcnc5; S[,"probncc5"]<-1-noprobncc5; S[,"probncnc5"]<-1-noprobncnc5;
  
  Infects1c=1 - (noprobcc1*noprobcnc1);
  Infects2c=1 - (noprobcc2*noprobcnc2);
  Infects3c=1 - (noprobcc3*noprobcnc3);
  Infects4c=1 - (noprobcc4*noprobcnc4);
  Infects5c=1 - (noprobcc5*noprobcnc5);
  
  Infects1nc=1 - (noprobncc1*noprobncnc1);
  Infects2nc=1 - (noprobncc2*noprobncnc2);
  Infects3nc=1 - (noprobncc3*noprobncnc3);
  Infects4nc=1 - (noprobncc4*noprobncnc4);
  Infects5nc=1 - (noprobncc5*noprobncnc5);

  #Implement infection and move susceptibles and newly exposeds back to home county 
  #For each Ss[,,i], entires in each column are individuals from the same region, with the rows showing if/where they travelled (so row sums give current total people in a region including visitors)
  #Exposure for compliant individuals
  S[,c("S1","E1")]=cbind(0,S[,"E1"]) + colSums(gpTransB(cbind(Ss[,,1],0*Ss[,,1]),round(Infects1c,10),seed+42))
  S[,c("S2","E2")]=cbind(0,S[,"E2"]) + colSums(gpTransB(cbind(Ss[,,2],0*Ss[,,2]),round(Infects2c,10),seed+43))
  S[,c("S3","E3")]=cbind(0,S[,"E3"]) + colSums(gpTransB(cbind(Ss[,,3],0*Ss[,,3]),round(Infects3c,10),seed+44))
  S[,c("S4","E4")]=cbind(0,S[,"E4"]) + colSums(gpTransB(cbind(Ss[,,4],0*Ss[,,4]),round(Infects4c,10),seed+45))
  S[,c("S5","E5")]=cbind(0,S[,"E5"]) + colSums(gpTransB(cbind(Ss[,,5],0*Ss[,,5]),round(Infects5c,10),seed+46))

  #Exposure for non-compliant individuals
  S[,c("S1nc","E1nc")]=cbind(0,S[,"E1nc"]) + colSums(gpTransB(cbind(Ssnc[,,1],0*Ssnc[,,1]),round(Infects1nc,10),seed+93))
  S[,c("S2nc","E2nc")]=cbind(0,S[,"E2nc"]) + colSums(gpTransB(cbind(Ssnc[,,2],0*Ssnc[,,2]),round(Infects2nc,10),seed+94))
  S[,c("S3nc","E3nc")]=cbind(0,S[,"E3nc"]) + colSums(gpTransB(cbind(Ssnc[,,3],0*Ssnc[,,3]),round(Infects3nc,10),seed+95))
  S[,c("S4nc","E4nc")]=cbind(0,S[,"E4nc"]) + colSums(gpTransB(cbind(Ssnc[,,4],0*Ssnc[,,4]),round(Infects4nc,10),seed+96))
  S[,c("S5nc","E5nc")]=cbind(0,S[,"E5nc"]) + colSums(gpTransB(cbind(Ssnc[,,5],0*Ssnc[,,5]),round(Infects5nc,10),seed+97))
  
  if(sum(S[,c(Susc,E,Da,Di,R,Suscnc,Enc,Danc,Dinc,Rnc)])!=sum(parms[,"N"]))  hel=lo
  return(S)
}; FUN=m3iter

#Implement initial closure and changes in testing probability
closeinappl=function(parms,TS,tm=dim(TS)[3],delayInit=10){
  
  Nta=colSums(t(t(TS[,"Nt",] + TS[,"Ntnc",]))); 
  
  if( (max(Nta)<inCstart)) parms[,"eps_h"]=parms[,"eps_o"]=0; #No NPI adherence before start of pandemic (as inCstart is March 10th when total cases>=50) 
  
  return(parms);
}

#Pulls info on what timestep it is
gettime=function(TS,tm=dim(TS)[3]){
  
  tsNt=colSums(t(t(TS[,"Nt",] + TS[,"Ntnc",]))); 
  if(max(tsNt)>=inCstart) {tepi<-length(tsNt[tsNt>=inCstart]) -1} else {tepi<--1} #Calc days into epidemic (since >=50 cases)
  if(max(tsNt)>=inCstart_s) {tstart_s<-length(tsNt[tsNt>=inCstart_s]) -1} else {tstart_s<--1} #Calc days past trigger date for school closures (==tepi for no counterfactuals)
  if(max(tsNt)>=inCstart_w) {tstart_w<-length(tsNt[tsNt>=inCstart_w]) -1} else {tstart_w<--1} #Calc days past trigger date for work closures (==tepi for no counterfactuals)
  
  timeinfo=c(tm, tepi, tstart_s, tstart_w)
  return(timeinfo);
}

#Function to implement simulations. InitInf sets which stages the initially infected people are in, InitInf=c(6:10, 71:75) indicates all initial infections correspond to individuals who are exposed, but split between different age classes and compliants/non-compliants
msim3=function(FUN,parms,Trun=365,seed0=11,plotgive=TRUE,InitInf=c(6:10,71:75)){
  #assign initial infections
  L=nrow(parms); nI=InitInf;
  set.seed(seed0);
  #Distribute population across age classes
  NsA=sapply(popadj2020, function(x) rmultinom(1,x,(agepopdist/13448505)));  # pop age demographics, StatCan Census 2016
  #Distribute population between compliants and non-compliants
  NsAnc=apply(NsA, c(1,2), function(x) rbinom(1,x,ncperc)); NsAc=NsA-NsAnc;

  #Distirbute Intial infections across age classes
  ### For compliants
  set.seed(seed0+1);
  Infs1c=NI0_1c=apply(cbind(NsAc[1,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[1]),x[1],prob=x[2]/length(nI[1])));
  set.seed(seed0+2);
  Infs2c=NI0_2c=apply(cbind(NsAc[2,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[2]),x[1],prob=x[2]/length(nI[2])));
  set.seed(seed0+3);
  Infs3c=NI0_3c=apply(cbind(NsAc[3,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[3]),x[1],prob=x[2]/length(nI[3])));
  set.seed(seed0+4);
  Infs4c=NI0_4c=apply(cbind(NsAc[4,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[4]),x[1],prob=x[2]/length(nI[4])));
  set.seed(seed0+5);
  Infs5c=NI0_5c=apply(cbind(NsAc[5,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[5]),x[1],prob=x[2]/length(nI[5])));
  ### For non-compliants
  set.seed(seed0+7);
  Infs1nc=NI0_1nc=apply(cbind(NsAnc[1,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[6]),x[1],prob=x[2]/length(nI[6])));
  set.seed(seed0+8);
  Infs2nc=NI0_2nc=apply(cbind(NsAnc[2,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[7]),x[1],prob=x[2]/length(nI[7])));
  set.seed(seed0+9);
  Infs3nc=NI0_3nc=apply(cbind(NsAnc[3,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[8]),x[1],prob=x[2]/length(nI[8])));
  set.seed(seed0+10);
  Infs4nc=NI0_4nc=apply(cbind(NsAnc[4,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[9]),x[1],prob=x[2]/length(nI[9])));
  set.seed(seed0+11);
  Infs5nc=NI0_5nc=apply(cbind(NsAnc[5,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[10]),x[1],prob=x[2]/length(nI[10])));

  #Create object to store simulation results
  TS=array(0,dim=c(L,length(Snm.comb),1)); 
  #Seed infections amongst compliant individuals
  TS[,c(1,InitInf[1]),1]=cbind(NsAc[1,]-NI0_1c,Infs1c); TS[,c(2,InitInf[2]),1]=cbind(NsAc[2,]-NI0_2c,Infs2c); TS[,c(3,InitInf[3]),1]=cbind(NsAc[3,]-NI0_3c,Infs3c); TS[,c(4,InitInf[4]),1]=cbind(NsAc[4,]-NI0_4c,Infs4c);	TS[,c(5,InitInf[5]),1]=cbind(NsAc[5,]-NI0_5c,Infs5c);
  #Seed infectious amongst non-compliant individuals
  TS[,c(66,InitInf[6]),1]=cbind(NsAnc[1,]-NI0_1nc,Infs1nc); TS[,c(67,InitInf[7]),1]=cbind(NsAnc[2,]-NI0_2nc,Infs2nc); TS[,c(68,InitInf[8]),1]=cbind(NsAnc[3,]-NI0_3nc,Infs3nc); TS[,c(69,InitInf[9]),1]=cbind(NsAnc[4,]-NI0_4nc,Infs4nc);	TS[,c(70,InitInf[10]),1]=cbind(NsAnc[5,]-NI0_5nc,Infs5nc);
  
  TS[,57:59,1]=0; #Set Cw=Cs=C=0
  colnames(TS)=Snm.comb; 
  #Modify travel probability as needed
  Mn=round(parms[1,"Mfact"]*parms[,-(1:NP)]*(1-diag(L))); diag(Mn)=parms[,"N"]-colSums(Mn); parms[,-(1:NP)]=Mn%*%diag(1/colSums(Mn));
  
  #Implement simulation
  set.seed(seed0+6); Seeds=rnorm(Trun,1e6,1e6);
  for(i in 2:Trun) TS=abind(TS,FUN(TS[,,i-1],closeinappl(parms,TS),gettime(TS),seed=Seeds[i]),along=3); TS[,"C",]<-apply(TS[,c("Cw","Cs"),],c(1,3),sum); #TS[,"Cw",]=parms[1,"eps_w"]*(TS[,"Cw",]>0); TS[,"Cs",]=parms[1,"eps_s"]*(TS[,"Cs",]>0);
  

  #Different levels of aggregation in model output. plotgive=3 or 3.5 give shortest output form (tracking only infections and costs) in integer format to reduce output size
  if(plotgive=="TS") return(TS); if(plotgive==TRUE) return(Resagg(TS,parms,plotgive==1));
  if(plotgive%in%c(3,3.5)){ TS[,"C",]=TS[,"C",]*matrix(parms[,"N"],nrow(parms),Trun); TSn=apply(TS,c(2,3),sum); out=matrix(as.integer(TSn),nrow(TSn),ncol(TSn)); if(plotgive==3.5) return(out[c(1:5,59),]); return(rbind(out,colMeans(TS[,"C",]>0))); }
}

#Run multiple sims based on our best parm sets for different combinations of commperc, redeffifacy
nrun<-5
runlength<-175

ncperc.set<-seq((as.numeric(args[1])/100 - 0.02), (as.numeric(args[1])/100), 0.01); 
susceffect.set<-0.5

for (ii in 1:length(ncperc.set))
{
  for (jj in 1:length(susceffect.set)) {
    
    #Set % of pop who are non-compliant, reduction in NPI efficacy during a suceptible-infectious interaction when the suceptible individual is not adhering to NPIs
    ncperc=ncperc.set[ii]
    susceffect=susceffect.set[jj]

    for (run in 1:nrun) {

      #Run a sim based on a fitted parameter set
      x_region=(msim3(m3iter,parms,runlength,plotgive="TS",seed0=fit*run*1e4,InitInf=c(6:10, 71:75)));

      # #Quick plot
      # agedat.plot<-agedat[agedat$Date<="2020-07-15",]
      # dat=data.frame(ts50=1:runlength-min(which(colSums(colSums(x_region[,c("Nt1","Nt2","Nt3","Nt4","Nt5","Nt1nc","Nt2nc","Nt3nc","Nt4nc","Nt5nc"),]))>50)), Nt=colSums(colSums(x_region[,c("Nt1","Nt2","Nt3","Nt4","Nt5","Nt1nc","Nt2nc","Nt3nc","Nt4nc","Nt5nc"),])))
      # dat=dat[dat$ts50<=max(agedat.plot$ts50),]
      # plot(agedat.plot$ts50, rowSums(agedat.plot[,c( "newK.1", "newK.2", "newK.3", "newK.4", "newK.5")]), main=sprintf("ncperc= %f, susceffect = %f", ncperc, susceffect),ylim=range(c(dat$Nt-dplyr::lag(dat$Nt), rowSums(agedat.plot[,c( "newK.1", "newK.2", "newK.3", "newK.4", "newK.5")])), na.rm=T))
      # lines(dat$ts50, dat$Nt-dplyr::lag(dat$Nt), col="blue")
      
      #Create all the summary data that would normally be generated with plotgive="TS" to use for the provincial dataset
      omgs=t(apply(x_region[,"VD",], 2, function(x) meansd(x,parms[,"N"])));
      #Note: using total population sizes for weights on these age-specific values as we assume same age dist in all regions so pop in each age class is proportional to total pop
      pinfect1<-cbind(t(apply(x_region[,"probcc1",], 2, function(x) meansd(x,parms[,"N"]))), t(apply(x_region[,"probcnc1",], 2, function(x) meansd(x,parms[,"N"]))), t(apply(x_region[,"probncc1",], 2, function(x) meansd(x,parms[,"N"]))), t(apply(x_region[,"probncnc1",], 2, function(x) meansd(x,parms[,"N"]))));
      pinfect2<-cbind(t(apply(x_region[,"probcc2",], 2, function(x) meansd(x,parms[,"N"]))), t(apply(x_region[,"probcnc2",], 2, function(x) meansd(x,parms[,"N"]))), t(apply(x_region[,"probncc2",], 2, function(x) meansd(x,parms[,"N"]))), t(apply(x_region[,"probncnc2",], 2, function(x) meansd(x,parms[,"N"]))));
      pinfect3<-cbind(t(apply(x_region[,"probcc3",], 2, function(x) meansd(x,parms[,"N"]))), t(apply(x_region[,"probcnc3",], 2, function(x) meansd(x,parms[,"N"]))), t(apply(x_region[,"probncc3",], 2, function(x) meansd(x,parms[,"N"]))), t(apply(x_region[,"probncnc3",], 2, function(x) meansd(x,parms[,"N"]))));
      pinfect4<-cbind(t(apply(x_region[,"probcc4",], 2, function(x) meansd(x,parms[,"N"]))), t(apply(x_region[,"probcnc4",], 2, function(x) meansd(x,parms[,"N"]))), t(apply(x_region[,"probncc4",], 2, function(x) meansd(x,parms[,"N"]))), t(apply(x_region[,"probncnc4",], 2, function(x) meansd(x,parms[,"N"]))));
      pinfect5<-cbind(t(apply(x_region[,"probcc5",], 2, function(x) meansd(x,parms[,"N"]))), t(apply(x_region[,"probcnc5",], 2, function(x) meansd(x,parms[,"N"]))), t(apply(x_region[,"probncc5",], 2, function(x) meansd(x,parms[,"N"]))), t(apply(x_region[,"probncnc5",], 2, function(x) meansd(x,parms[,"N"]))));
      propCits=t(x_region[,"Nt",]+x_region[,"Ntnc",]); states=t(apply(x_region,c(2,3),sum));
      x_prov=cbind(states,omgs,pinfect1,pinfect2,pinfect3,pinfect4,pinfect5,propCits);
      colnames(x_prov)[147:188]<-c("meanVD", "sdVD", "meanprobcc1", "sdprobcc1", "meanprobcnc1", "sdprobcnc1", "meanprobncc1", "sdprobncc1", "meanprobncnc1", "sdprobncnc1",
                                   "meanprobcc2", "sdprobcc2", "meanprobcnc2", "sdprobcnc2", "meanprobncc2", "sdprobncc2", "meanprobncnc2", "sdprobncnc2",
                                   "meanprobcc3", "sdprobcc3", "meanprobcnc3", "sdprobcnc3", "meanprobncc3", "sdprobncc3", "meanprobncnc3", "sdprobncnc3",
                                   "meanprobcc4", "sdprobcc4", "meanprobcnc4", "sdprobcnc4", "meanprobncc4", "sdprobncc4", "meanprobncnc4", "sdprobncnc4",
                                   "meanprobcc5", "sdprobcc5", "meanprobcnc5", "sdprobcnc5", "meanprobncc5", "sdprobncc5", "meanprobncnc5", "sdprobncnc5")
      
      for (i in 1:length(x_prov[,1])) {
        
        tmp<-as.data.frame(t(x_prov[i,]))
        
        if (i==1) {
          df<-tmp
        } else {
          df<-rbind(df,tmp)
        }
      }
      
      for (i in 1:length(x_region[,1,1])) {
        
        tmp_region<-as.data.frame(t(x_region[i,,]))
        tmp_region$ts <- as.numeric(row.names(tmp_region))-1
        tmp_region$region<-i
        
        if (i==1) {
          df_region<-tmp_region
        } else {
          df_region<-rbind(df_region,tmp_region)
        }
      }
      
      #Add-ins for regional sims
      df_region$sick.c<-rowSums(df_region[,colnames(df_region) %in% c(E,Di,Da)]); df_region$know.current.c<-rowSums(df_region[,colnames(df_region) %in% c(Tk)]); df_region$recovered.c<-rowSums(df_region[,colnames(df_region) %in% c(R)]);
      df_region$sick.nc<-rowSums(df_region[,colnames(df_region) %in% c(Enc,Dinc,Danc)]); df_region$know.current.nc<-rowSums(df_region[,colnames(df_region) %in% c(Tknc)]); df_region$recovered.nc<-rowSums(df_region[,colnames(df_region) %in% c(Rnc)]);
      
      df_region$known1.c<-rowSums(df_region[, colnames(df_region) %in% c(Tk[1:4])]); df_region$known2.c<-rowSums(df_region[, colnames(df_region) %in% c(Tk[5:8])]);df_region$known3.c<-rowSums(df_region[, colnames(df_region) %in% c(Tk[9:12])]);df_region$known4.c<-rowSums(df_region[, colnames(df_region) %in% c(Tk[13:16])]);df_region$known5.c<-rowSums(df_region[, colnames(df_region) %in% c(Tk[17:20])]);
      df_region$known1.nc<-rowSums(df_region[, colnames(df_region) %in% c(Tknc[1:4])]); df_region$known2.nc<-rowSums(df_region[, colnames(df_region) %in% c(Tknc[5:8])]);df_region$known3.nc<-rowSums(df_region[, colnames(df_region) %in% c(Tknc[9:12])]);df_region$known4.nc<-rowSums(df_region[, colnames(df_region) %in% c(Tknc[13:16])]);df_region$known5.nc<-rowSums(df_region[, colnames(df_region) %in% c(Tknc[17:20])]);
      df_region$sick1.c<-rowSums(df_region[, colnames(df_region) %in% c(sick1)]); df_region$sick2.c<-rowSums(df_region[, colnames(df_region) %in% c(sick2)]); df_region$sick3.c<-rowSums(df_region[, colnames(df_region) %in% c(sick3)]); df_region$sick4.c<-rowSums(df_region[, colnames(df_region) %in% c(sick4)]); df_region$sick5.c<-rowSums(df_region[, colnames(df_region) %in% c(sick5)]);
      df_region$sick1.nc<-rowSums(df_region[, colnames(df_region) %in% c(sick1nc)]); df_region$sick2.nc<-rowSums(df_region[, colnames(df_region) %in% c(sick2nc)]); df_region$sick3.nc<-rowSums(df_region[, colnames(df_region) %in% c(sick3nc)]); df_region$sick4.nc<-rowSums(df_region[, colnames(df_region) %in% c(sick4nc)]); df_region$sick5.nc<-rowSums(df_region[, colnames(df_region) %in% c(sick5nc)]);
      df_region$fit<-unique(fitparms$fit)[fit]
      df_region$run<-run
      df_region$LL<-unique(fitparms$LL[fitparms$fit==unique(fitparms$fit)[fit]])
      
      df_region <- df_region %>% group_by(fit, run, region) %>% mutate(new.K.c=Nt - lag(Nt), new.K.nc=Ntnc - lag(Ntnc),
                                                                       new.K.1.c=Nt1 - lag(Nt1), new.K.2.c=Nt2 - lag(Nt2), new.K.3.c=Nt3 - lag(Nt3), new.K.4.c=Nt4 - lag(Nt4), new.K.5.c=Nt5 - lag(Nt5),
                                                                       new.K.1.nc=Nt1nc - lag(Nt1nc), new.K.2.nc=Nt2nc - lag(Nt2nc), new.K.3.nc=Nt3nc - lag(Nt3nc), new.K.4.nc=Nt4nc - lag(Nt4nc), new.K.5.nc=Nt5nc - lag(Nt5nc),
                                                                       new.sick.c=(sick.c + recovered.c) -  lag(sick.c + recovered.c), new.sick.nc=(sick.nc + recovered.nc) -  lag(sick.nc + recovered.nc),
                                                                       new.sick.1.c=(sick1.c + R1) - lag(sick1.c  + R1), new.sick.2.c=(sick2.c + R2) - lag(sick2.c  + R2), new.sick.3.c=(sick3.c  + R3) - lag(sick3.c  + R3), new.sick.4.c=(sick4.c + R4) - lag(sick4.c + R4), new.sick.5.c=(sick5.c + R5) - lag(sick5.c  + R5),
                                                                       new.sick.1.nc=(sick1.nc + R1nc) - lag(sick1.nc  + R1nc), new.sick.2.nc=(sick2.nc + R2nc) - lag(sick2.nc  + R2nc), new.sick.3.nc=(sick3.nc  + R3nc) - lag(sick3.nc  + R3nc), new.sick.4.nc=(sick4.nc + R4nc) - lag(sick4.nc + R4nc), new.sick.5.nc=(sick5.nc + R5nc) - lag(sick5.nc  + R5nc))
   
      #Set all first step "new" to equal total at that time since we dont have a previous ts for comparison
      df_region$new.K.c[df_region$ts==min(df_region$ts)]<-df_region$Nt[df_region$ts==min(df_region$ts)]
      df_region$new.K.1.c[df_region$ts==min(df_region$ts)]<-df_region$Nt1[df_region$ts==min(df_region$ts)]
      df_region$new.K.2.c[df_region$ts==min(df_region$ts)]<-df_region$Nt2[df_region$ts==min(df_region$ts)]
      df_region$new.K.3.c[df_region$ts==min(df_region$ts)]<-df_region$Nt3[df_region$ts==min(df_region$ts)]
      df_region$new.K.4.c[df_region$ts==min(df_region$ts)]<-df_region$Nt4[df_region$ts==min(df_region$ts)]
      df_region$new.K.5.c[df_region$ts==min(df_region$ts)]<-df_region$Nt5[df_region$ts==min(df_region$ts)]
      df_region$new.sick.c[df_region$ts==min(df_region$ts)]<-df_region$sick.c[df_region$ts==min(df_region$ts)]
      df_region$new.sick.1.c[df_region$ts==min(df_region$ts)]<-df_region$sick1.c[df_region$ts==min(df_region$ts)]
      df_region$new.sick.2.c[df_region$ts==min(df_region$ts)]<-df_region$sick2.c[df_region$ts==min(df_region$ts)]
      df_region$new.sick.3.c[df_region$ts==min(df_region$ts)]<-df_region$sick3.c[df_region$ts==min(df_region$ts)]
      df_region$new.sick.4.c[df_region$ts==min(df_region$ts)]<-df_region$sick4.c[df_region$ts==min(df_region$ts)]
      df_region$new.sick.5.c[df_region$ts==min(df_region$ts)]<-df_region$sick5.c[df_region$ts==min(df_region$ts)]
      
      df_region$new.K.nc[df_region$ts==min(df_region$ts)]<-df_region$Ntnc[df_region$ts==min(df_region$ts)]
      df_region$new.K.1.nc[df_region$ts==min(df_region$ts)]<-df_region$Nt1nc[df_region$ts==min(df_region$ts)]
      df_region$new.K.2.nc[df_region$ts==min(df_region$ts)]<-df_region$Nt2nc[df_region$ts==min(df_region$ts)]
      df_region$new.K.3.nc[df_region$ts==min(df_region$ts)]<-df_region$Nt3nc[df_region$ts==min(df_region$ts)]
      df_region$new.K.4.nc[df_region$ts==min(df_region$ts)]<-df_region$Nt4nc[df_region$ts==min(df_region$ts)]
      df_region$new.K.5.nc[df_region$ts==min(df_region$ts)]<-df_region$Nt5nc[df_region$ts==min(df_region$ts)]
      df_region$new.sick.nc[df_region$ts==min(df_region$ts)]<-df_region$sick.nc[df_region$ts==min(df_region$ts)]
      df_region$new.sick.1.nc[df_region$ts==min(df_region$ts)]<-df_region$sick1.nc[df_region$ts==min(df_region$ts)]
      df_region$new.sick.2.nc[df_region$ts==min(df_region$ts)]<-df_region$sick2.nc[df_region$ts==min(df_region$ts)]
      df_region$new.sick.3.nc[df_region$ts==min(df_region$ts)]<-df_region$sick3.nc[df_region$ts==min(df_region$ts)]
      df_region$new.sick.4.nc[df_region$ts==min(df_region$ts)]<-df_region$sick4.nc[df_region$ts==min(df_region$ts)]
      df_region$new.sick.5.nc[df_region$ts==min(df_region$ts)]<-df_region$sick5.nc[df_region$ts==min(df_region$ts)]
      
      #Add ins for provincial sims
      df$ts<-0:(length(x_prov[,1])-1)
      
      df$sick.c<-rowSums(df[,colnames(df) %in% c(E,Di,Da)]); df$know.current.c<-rowSums(df[,colnames(df) %in% c(Tk)]); df$recovered.c<-rowSums(df[,colnames(df) %in% c(R)])
      df$sick.nc<-rowSums(df[,colnames(df) %in% c(Enc,Dinc,Danc)]); df$know.current.nc<-rowSums(df[,colnames(df) %in% c(Tknc)]); df$recovered.nc<-rowSums(df[,colnames(df) %in% c(Rnc)])
      
      df$known1.c<-rowSums(df[, colnames(df) %in% c(Tk[1:4])]);df$known2.c<-rowSums(df[, colnames(df) %in% c(Tk[5:8])]);df$known3.c<-rowSums(df[, colnames(df) %in% c(Tk[9:12])]);df$known4.c<-rowSums(df[, colnames(df) %in% c(Tk[13:16])]);df$known5.c<-rowSums(df[, colnames(df) %in% c(Tk[17:20])]);
      df$known1.nc<-rowSums(df[, colnames(df) %in% c(Tknc[1:4])]);df$known2.nc<-rowSums(df[, colnames(df) %in% c(Tknc[5:8])]);df$known3.nc<-rowSums(df[, colnames(df) %in% c(Tknc[9:12])]);df$known4.nc<-rowSums(df[, colnames(df) %in% c(Tknc[13:16])]);df$known5.nc<-rowSums(df[, colnames(df) %in% c(Tknc[17:20])]);
      df$sick1.c<-rowSums(df[, colnames(df) %in% c(sick1)]); df$sick2.c<-rowSums(df[, colnames(df) %in% c(sick2)]); df$sick3.c<-rowSums(df[, colnames(df) %in% c(sick3)]); df$sick4.c<-rowSums(df[, colnames(df) %in% c(sick4)]); df$sick5.c<-rowSums(df[, colnames(df) %in% c(sick5)]);
      df$sick1.nc<-rowSums(df[, colnames(df) %in% c(sick1nc)]); df$sick2.nc<-rowSums(df[, colnames(df) %in% c(sick2nc)]); df$sick3.nc<-rowSums(df[, colnames(df) %in% c(sick3nc)]); df$sick4.nc<-rowSums(df[, colnames(df) %in% c(sick4nc)]); df$sick5.nc<-rowSums(df[, colnames(df) %in% c(sick5nc)]);
      df$fit<-unique(fitparms$fit)[fit]
      df$run<-run
      df$LL<-unique(fitparms$LL[fitparms$fit==unique(fitparms$fit)[fit]])
      
      df$new.K.c<-df$Nt - lag(df$Nt); df$new.K.nc<-df$Ntnc - lag(df$Ntnc);
      df$new.K.1.c<-df$Nt1 - lag(df$Nt1); df$new.K.2.c<-df$Nt2- lag(df$Nt2); df$new.K.3.c<-df$Nt3 - lag(df$Nt3); df$new.K.4.c<-df$Nt4 - lag(df$Nt4); df$new.K.5.c<-df$Nt5 - lag(df$Nt5);
      df$new.K.1.nc<-df$Nt1nc - lag(df$Nt1nc); df$new.K.2.nc<-df$Nt2nc - lag(df$Nt2nc); df$new.K.3.nc<-df$Nt3nc - lag(df$Nt3nc); df$new.K.4.nc<-df$Nt4nc - lag(df$Nt4nc); df$new.K.5.nc<-df$Nt5nc - lag(df$Nt5nc);
      df$new.K.c[1]<-df$Nt[1]; df$new.K.1.c[1]<-df$Nt1[1]; df$new.K.2.c[1]<-df$Nt2[1]; df$new.K.3.c[1]<-df$Nt3[1]; df$new.K.4.c[1]<-df$Nt4[1]; df$new.K.5.c[1]<-df$Nt5[1];
      df$new.K.nc[1]<-df$Ntnc[1]; df$new.K.1.nc[1]<-df$Nt1nc[1]; df$new.K.2.nc[1]<-df$Nt2nc[1]; df$new.K.3.nc[1]<-df$Nt3nc[1]; df$new.K.4.nc[1]<-df$Nt4nc[1]; df$new.K.5.nc[1]<-df$Nt5nc[1];
      
      df$new.sick.c<-(df$sick.c + df$recovered.c) -  lag(df$sick.c + df$recovered.c); df$new.sick.nc<-(df$sick.nc + df$recovered.nc) -  lag(df$sick.nc + df$recovered.nc);
      df$new.sick.1.c<-(df$sick1.c + df$R1) - lag(df$sick1.c  + df$R1); df$new.sick.2.c<-(df$sick2.c  + df$R2) - lag(df$sick2.c  + df$R2); df$new.sick.3.c<-(df$sick3.c  + df$R3) - lag(df$sick3.c  + df$R3); df$new.sick.4.c<-(df$sick4.c  + df$R4) - lag(df$sick4.c  + df$R4); df$new.sick.5.c<-(df$sick5.c  + df$R5) - lag(df$sick5.c  + df$R5);
      df$new.sick.1.nc<-(df$sick1.nc + df$R1nc) - lag(df$sick1.nc  + df$R1nc); df$new.sick.2.nc<-(df$sick2.nc + df$R2nc) - lag(df$sick2.nc + df$R2nc); df$new.sick.3.nc<-(df$sick3.nc  + df$R3nc) - lag(df$sick3.nc  + df$R3nc); df$new.sick.4.nc<-(df$sick4.nc + df$R4nc) - lag(df$sick4.nc + df$R4nc); df$new.sick.5.nc<-(df$sick5.nc + df$R5nc) - lag(df$sick5.nc + df$R5nc);
      df$new.sick.c[1]<-df$sick.c[1]; df$new.sick.1.c[1]<-df$sick1.c[1]; df$new.sick.2.c[1]<-df$sick2.c[1]; df$new.sick.3.c[1]<-df$sick3.c[1]; df$new.sick.4.c[1]<-df$sick4.c[1]; df$new.sick.5.c[1]<-df$sick5.c[1];
      df$new.sick.nc[1]<-df$sick.nc[1]; df$new.sick.1.nc[1]<-df$sick1.nc[1]; df$new.sick.2.nc[1]<-df$sick2.nc[1]; df$new.sick.3.nc[1]<-df$sick3.nc[1]; df$new.sick.4.nc[1]<-df$sick4.nc[1]; df$new.sick.5.nc[1]<-df$sick5.nc[1];
      
      #Tag rw data so we can match based on days since 50th positive case
      model50<-df$ts[min(which((df$Nt+df$Ntnc) >= 50))]
      df$ts50<-df$ts-model50
 
      #Tag rw data so we can match based on days since 50th positive case
      df_region.sum <- df_region %>%
        group_by(ts, fit, run) %>%
        dplyr::summarize(cases.K=sum(Nt+Ntnc))
      
      model50_region<-df_region$ts[min(which(df_region.sum$cases.K >= 50))]
      df_region$ts50<-df_region$ts-model50_region
      
      if (model50_region!=model50) { print("ts50 problem");}
      
      df.reduced<-df[,c(1:188,238:292)] #Knock out 49x49 matrix to reduce df size
      
      if (run==1) {
        df.prov.fit<-df.reduced
        df.region.fit<-df_region
      } else {
        df.prov.fit<-rbind(df.prov.fit, df.reduced)
        df.region.fit<-rbind(df.region.fit, df_region)
      }
      
    }
    
    df.fitparms<-as.data.frame(fitparms)
    df.fitparms$fit<-unique(fitparms$fit)[fit]
    df.fitparms$LL<-unique(fitparms$LL[fitparms$fit==unique(fitparms$fit)[fit]])
    
    df.prov.fit$ncperc<-ncperc
    df.region.fit$ncperc<-ncperc
    df.fitparms$ncperc<-ncperc
    
    df.prov.fit$susceffect<-susceffect
    df.region.fit$susceffect<-susceffect
    df.fitparms$susceffect<-susceffect
    
    if (goodfits==1 & ii==1 & jj==1){
      df.prov.comb<-df.prov.fit
      df.region.comb<-df.region.fit
      df.fitparms.comb<-df.fitparms
    } else {
      df.prov.comb<-rbind(df.prov.comb, df.prov.fit)
      df.region.comb<-rbind(df.region.comb, df.region.fit)
      df.fitparms.comb<-rbind(df.fitparms.comb, df.fitparms)
    }
    
    #Print progress
    print(length(unique(interaction(df.prov.comb$fit, df.prov.comb$ncperc, df.prov.comb$susceffect)))/(topnum*length(ncperc.set)*length(susceffect.set)))
    
    goodfits<-goodfits+1;
    
    
  }
}

} 

saveRDS(df.region.comb, file = sprintf("simsREGIONS_suscRAND_ncincr_%s_susceffect%.2f_%.2f_%.2f.rds", codeversion, susceffect, (as.numeric(args[1])/100 - 0.02), (as.numeric(args[1])/100) ));
saveRDS(df.prov.comb, file = sprintf("simsPROVINCE_suscRAND_ncincr_%s_susceffect%.2f_%.2f_%.2f.rds", codeversion, susceffect, (as.numeric(args[1])/100 - 0.02), (as.numeric(args[1])/100) ));
saveRDS(df.fitparms.comb, file = sprintf("simsFITPARMS_suscRAND_ncincr_%s_susceffect%.2f_%.2f_%.2f.rds", codeversion, susceffect, (as.numeric(args[1])/100 - 0.02), (as.numeric(args[1])/100) ));


tfin<-Sys.time()
print(tfin-tstart)