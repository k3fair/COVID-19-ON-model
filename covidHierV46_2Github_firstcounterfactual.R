#     Script runs simulations for counterfactuals for 10 March 2020 to 15 August 2020 (presence/absence of closures and/or NPI adherence)
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

library(abind); library(plyr);  library(dplyr); library(data.table);
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

setwd("C:/Users/jogad/OneDrive - University of Waterloo/S20/CovidModelling/OntarioProj/COVID-19-Hierarchy-model-master/SharcnetSims/")
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
codeversion<-"v46v2phuNUDGE"

# Set closure counterfactual type i.e. schools remain open ("schoolopen"), workplaces remain open ("workopen"), both remain open ("bothopen") or both are shut ("neitheropen", corresponds to what actually occured in the province)
# Set reopening type; with ("restricted") or without ("unrestricted") NPIs in schools/workplaces
# Set individual NPI adherence counterfactual type i.e. no individual adherence to NPIs ("vdOFF") or individual adherence to NPIs in response to case numbers ("vdON")
cftype<-"neitheropen"; 
vdtype<-"vdOFF"
reopeningtype<-"restricted";

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

goodfits<-1;
for (fit in 1:topnum) {
  
  parms=as.matrix(fitparms[fitparms$fit==unique(fitparms$fit)[fit],1:88]); # use best fit from fitting
  
  #Some parameter name differences from paper:
  #initial testing tauI_t0=cvTl, final testing tauI_tf=tauI, tauA=tauEAf*tauI
  # a1,a2,a3,a4,a5 are age-specific susceptibility modifiers (gamma_i, i=1,2,3,4,5 in the paper)
  # boost is L0 (impact of stay at home orders on NPI adherence)
  #Mfact scales how much travel happens compared to reality; Mfact=1 is just movement from commuting data
  #epsP allows individual NPI adherence efficacy to differ from that of closures (epsP=eps in paper)
  #Tg and Tl are the gobal and local closure thresholds, n0 is fraction initially infected
  #Msave is the travel matrix. Here Msave entries are numbers of commuters; msim3 converts them to proportions
  pops=colSums(DATA$Msave)
  
  ###Adjust to 2020 Q4 pop estimate (from https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1710000901)
  pop2020<-14733119
  popratio<-pop2020/sum(pops)
  
  popadj2020=round(popratio*pops)
  
  # Specification of closure/reopening events
  inCstart=50; #all events occur at a distance from the day total cases >=50 (March 10th)
  inClen_w1=79; #Duration of initial provincial workplace closure (March 25-June 11, Phase 2 commenced June 12)
  inClen_s1=178; #Duration of initial provincial school closure (Mar 14 - Sept 7, Schools reopen Sept 8)
  inClen_w2=46; #Duration of second provincial workplace closure, to first part of staggered reopening (Dec 26th 2020 - Feb 9th 2021)
  inClen_s2=51; #Duration of second provincial school closure, to first part of staggered reopening (Dec 21st 2020 - Feb 9th 2021)
  inBlen_s=70; #Duration of Summer break (June 30th - Sept 8th)
  
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
  schoolbreakgap=113; #Gap between day total cases >= 50 (Mar 10th) and day schools would have closed for Summer Holidays (July 1st)
  
  omggap=144 #Gap between day total cases >= 50 (Mar 10th) and day we switch to wave 2 omega value (Aug 1st)
  Resurge_w=50 #Duration of additional closures 
  Resurge_s=50
  NP=ncol(parms)-nrow(parms)
  
  #Dummy start values for closures s.t. we can do counterfactuals for first wave
  inCstart_w=50;
  inCstart_s=50;
  if (cftype %in% c("workopen", "bothopen")) {inCstart_w=1e8;} #Set so workplaces never close
  
  years<- 1 
  
  
  #Defining model state space. Tn, Tk, Da, and Di are all untested, tested, asymptomatics, and infecteds respectively.
  #Nt tracks cumulative # positive cases (including those recovered)
  #C tracks # days since last closure, but in output msim converts all positive C entries into eps
  #VD indicates level of NPI adherence from individuals
  #Sick indicates all individuals who are in any of E, P, A, I (i.e. all those either exposed or currently infectious)
  S=c("S1", "S2", "S3", "S4", "S5"); E=c("E1", "E2","E3", "E4", "E5"); R=c("R1", "R2", "R3", "R4","R5");
  Da=c("A1","Ak1","SA1","SAk1", "A2","Ak2","SA2","SAk2", "A3","Ak3","SA3","SAk3", "A4","Ak4","SA4","SAk4", "A5","Ak5","SA5","SAk5")
  Di=c("I1","Ik1","SI1","SIk1", "I2","Ik2","SI2","SIk2", "I3","Ik3","SI3","SIk3", "I4","Ik4","SI4","SIk4", "I5","Ik5","SI5","SIk5")
  Tn=c("A1","SA1","I1","SI1","A2","SA2","I2","SI2","A3","SA3","I3","SI3","A4","SA4","I4","SI4","A5","SA5","I5","SI5")
  Tk=c("Ak1","SAk1","Ik1","SIk1","Ak2","SAk2","Ik2","SIk2","Ak3","SAk3","Ik3","SIk3","Ak4","SAk4","Ik4","SIk4","Ak5","SAk5","Ik5","SIk5")
  All1<-c("S1", "E1", "R1", "A1","Ak1","SA1","SAk1", "I1","Ik1","SI1","SIk1"); All2<-c("S2", "E2", "R2", "A2","Ak2","SA2","SAk2", "I2","Ik2","SI2","SIk2"); All3<-c("S3", "E3", "R3", "A3","Ak3","SA3","SAk3", "I3","Ik3","SI3","SIk3"); All4<-c("S4", "E4", "R4", "A4","Ak4","SA4","SAk4", "I4","Ik4","SI4","SIk4"); All5<-c("S5", "E5", "R5", "A5","Ak5","SA5","SAk5", "I5","Ik5","SI5","SIk5");
  sick1<-c("E1", "A1","Ak1","SA1","SAk1", "I1","Ik1","SI1","SIk1"); sick2<-c("E2", "A2","Ak2","SA2","SAk2", "I2","Ik2","SI2","SIk2"); sick3<-c("E3", "A3","Ak3","SA3","SAk3", "I3","Ik3","SI3","SIk3"); sick4<-c("E4", "A4","Ak4","SA4","SAk4", "I4","Ik4","SI4","SIk4"); sick5<-c("E5", "A5","Ak5","SA5","SAk5", "I5","Ik5","SI5","SIk5");
  symp1<-c("I1","Ik1","SI1","SIk1"); symp2<-c("I2","Ik2","SI2","SIk2"); symp3<-c("I3","Ik3","SI3","SIk3"); symp4<-c("I4","Ik4","SI4","SIk4"); symp5<-c("I5","Ik5","SI5","SIk5");
  Snm=c(S,E,Da,Di,R,"Nt","Cw", "Cs", "C", "Nt1","Nt2","Nt3","Nt4","Nt5", "VD");
  
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
    S0k=S[,Tk]; S0k_1=S[,Tk[1:4]]; S0k_2=S[,Tk[5:8]]; S0k_3=S[,Tk[9:12]]; S0k_4=S[,Tk[13:16]]; S0k_5=S[,Tk[17:20]];
    S0n_1=S[,Tn[1:4]]; S0n_2=S[,Tn[5:8]]; S0n_3=S[,Tn[9:12]]; S0n_4=S[,Tn[13:16]]; S0n_5=S[,Tn[17:20]];
    S[,c(Tn,Tk)]=gpTransB(S[,c(Tn,Tk)],	as.vector(cbind(testPA1,testPA1,testI1,testI1,testPA2,testPA2,testI2,testI2,testPA3,testPA3,testI3,testI3,testPA4,testPA4,testI4,testI4,testPA5,testPA5,testI5,testI5)),seed+40); 
    
    
    #calculate pos, the vector of local active case prevalence
    S[,"Nt"]=S[,"Nt"]+rowSums(S[,Tk]-S0k); pos=rowSums(S[,Tk])/Ns; posglobal=(sum(S[,Tk])/sum(Ns)); 
    S[,"Nt1"]=S[,"Nt1"]+rowSums(S[,Tk[1:4]]-S0k_1);
    S[,"Nt2"]=S[,"Nt2"]+rowSums(S[,Tk[5:8]]-S0k_2);
    S[,"Nt3"]=S[,"Nt3"]+rowSums(S[,Tk[9:12]]-S0k_3);
    S[,"Nt4"]=S[,"Nt4"]+rowSums(S[,Tk[13:16]]-S0k_4);
    S[,"Nt5"]=S[,"Nt5"]+rowSums(S[,Tk[17:20]]-S0k_5);

    
    #Disease progression: zeta is the fraction of people who never show symptoms (modeled implicitly)
    zeta=0.2;
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
    M=Mc=parms[,-(1:NP)]; Mc=M[1:nrow(M),]*(1-closed_w)*(1-closed_w); diag(Mc)=diag(Mc)+1-colSums(Mc);
    ### Revamping McA to include age specific travel rates, first of each pair is for old/young, second for middle
    McA=abind(Mstay(Mc,parms[1,"atravel"]), Mc, Mstay(Mc,1-(1-parms[1,"atravel"])*(1-parms[1,"eta"])), Mstay(Mc,parms[1,"eta"]), Mstay(Mc,1-(1-parms[1,"atravel"])*(1-parms[1,"r"])), Mstay(Mc,parms[1,"r"]), Mstay(Mc,1-(1-parms[1,"eta"])*(1-parms[1,"r"])*(1-parms[1,"atravel"])), Mstay(Mc,1-(1-parms[1,"eta"])*(1-parms[1,"r"])), along=3);
    #Implement travel
    ### Do Ss for however many age classes you have, added atravel for age-class specific travel
    Ss=abind(reshfl2(S[,"S1"],Mstay(Mc,parms[1,"atravel"]),seed+20,FALSE),reshfl2(S[,"S2"],Mc,seed+21,FALSE),reshfl2(S[,"S3"],Mc,seed+22,FALSE),reshfl2(S[,"S4"],Mstay(Mc,parms[1,"atravel"]),seed+23,FALSE),reshfl2(S[,"S5"],Mstay(Mc,parms[1,"atravel"]),seed+24,FALSE),along=3);
    Rearr=apply(rbind(seed+(25:34), c(c(1,2,2,1,1),c(1,2,2,1,1),c(1,3,1,3,2,4,2,4,2,4,2,4,1,3,1,3,1,3,1,3),c(5,7,5,7,6,8,6,8,6,8,6,8,5,7,5,7,5,7,5,7)), S[,c(R,E,Da,Di)]), 2, function(x) reshfl2(x[-(1:2)],McA[,,x[2]],x[1]))
    
    
    #Age specific contacts, rows are age classes, columns are their contact age classes
    ageSpecifics_w<-cmat[,,1] ## age specific contacts for work
    ageSpecifics_s<-cmat[,,2] ## age specific contacts for school
    ageSpecifics_h<-cmat[,,3] ## age specific contacts for home
    ageSpecifics_o<-cmat[,,4] ## age specific contacts for other
    

    #Seasonal forcing
    scomp<-1+parms[,"B"]*cos((2*pi/365)*(timeinfo[1] + parms[,"phi"]))
    
    
    S[, "VD"]<-(1-exp(-parms[1,"omg"]*pos));
    if (timeinfo[4]>=workclosuregap2 + 20) {
      S[, "VD"] <- ifelse(timeinfo[4]<workclosuregap2+inClen_w2, (1-exp(-(parms[1,"omg"]*pos + parms[1,"boost"]))), (1-exp(-parms[1,"omg"]*pos))) #Stay at home, with staggered reopening
    }
    
    if (vdtype=="vdOFF") {S[, "VD"]<-0*pos;}
    
    #Base infection probabilities
    Infect1 = scomp*parms[,"a1"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[1,]) + (1-closed_s)%*%t(ageSpecifics_s[1,]) + (1-parms[1,"eps_h"]*S[, "VD"])%*%t(ageSpecifics_h[1,]) + (1-parms[1,"eps_o"]*S[, "VD"])%*%t(ageSpecifics_o[1,]))
    Infect2 = scomp*parms[,"a2"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[2,]) + (1-closed_s)%*%t(ageSpecifics_s[2,]) + (1-parms[1,"eps_h"]*S[, "VD"])%*%t(ageSpecifics_h[2,]) + (1-parms[1,"eps_o"]*S[, "VD"])%*%t(ageSpecifics_o[2,]))
    Infect3 = scomp*parms[,"a3"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[3,]) + (1-closed_s)%*%t(ageSpecifics_s[3,]) + (1-parms[1,"eps_h"]*S[, "VD"])%*%t(ageSpecifics_h[3,]) + (1-parms[1,"eps_o"]*S[, "VD"])%*%t(ageSpecifics_o[3,]))
    Infect4 = scomp*parms[,"a4"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[4,]) + (1-closed_s)%*%t(ageSpecifics_s[4,]) + (1-parms[1,"eps_h"]*S[, "VD"])%*%t(ageSpecifics_h[4,]) + (1-parms[1,"eps_o"]*S[, "VD"])%*%t(ageSpecifics_o[4,]))
    Infect5 = scomp*parms[,"a5"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[5,]) + (1-closed_s)%*%t(ageSpecifics_s[5,]) + (1-parms[1,"eps_h"]*S[, "VD"])%*%t(ageSpecifics_h[5,]) + (1-parms[1,"eps_o"]*S[, "VD"])%*%t(ageSpecifics_o[5,]))
    
    #Class-specific infection modifiers
    modK=1-parms[1,"eta"]; modS=(1/parms[1,"s"])-1; 
    # modA=cbind(Infect,modK*Infect,modS*Infect,modS*modK*Infect);
    
    #Introduce age specific modA 3-D arrays for not superspreader or known, not supersrpeader but known, superspreader not known, superspreader and known
    modA1=abind(Infect1,modK*Infect1,modS*Infect1,modS*modK*Infect1,along=3);
    modA2=abind(Infect2,modK*Infect2,modS*Infect2,modS*modK*Infect2,along=3);
    modA3=abind(Infect3,modK*Infect3,modS*Infect3,modS*modK*Infect3,along=3);
    modA4=abind(Infect4,modK*Infect4,modS*Infect4,modS*modK*Infect4,along=3);
    modA5=abind(Infect5,modK*Infect5,modS*Infect5,modS*modK*Infect5,along=3);
    
    #Flatten into 2-D arrays where order is based on age class to match Rearr which is ordered as (all pre/asympt by age class, all sympt by age class)
    modAbind1<-cbind(modA1[,1,],modA1[,2,],modA1[,3,],modA1[,4,],modA1[,5,]);
    modAbind2<-cbind(modA2[,1,],modA2[,2,],modA2[,3,],modA2[,4,],modA2[,5,]);
    modAbind3<-cbind(modA3[,1,],modA3[,2,],modA3[,3,],modA3[,4,],modA3[,5,]);
    modAbind4<-cbind(modA4[,1,],modA4[,2,],modA4[,3,],modA4[,4,],modA4[,5,]);
    modAbind5<-cbind(modA5[,1,],modA5[,2,],modA5[,3,],modA5[,4,],modA5[,5,]);
    
    Infects1=1 - apply((1-cbind(modAbind1, 2*modAbind1))^Rearr[,-(1:10)], 1, prod)
    Infects2=1 - apply((1-cbind(modAbind2, 2*modAbind2))^Rearr[,-(1:10)], 1, prod)
    Infects3=1 - apply((1-cbind(modAbind3, 2*modAbind3))^Rearr[,-(1:10)], 1, prod)
    Infects4=1 - apply((1-cbind(modAbind4, 2*modAbind4))^Rearr[,-(1:10)], 1, prod)
    Infects5=1 - apply((1-cbind(modAbind5, 2*modAbind5))^Rearr[,-(1:10)], 1, prod)
    
    #Implement infection and move susceptibles and newly exposeds back to home county 
    #For each Ss[,,i], entires in each column are individuals from the same region, with the rows showing if/where they travelled (so row sums give current total people in a region including visitors)
    S[,c("S1","E1")]=cbind(0,S[,"E1"]) + colSums(gpTransB(cbind(Ss[,,1],0*Ss[,,1]),round(Infects1,10),seed+35))
    S[,c("S2","E2")]=cbind(0,S[,"E2"]) + colSums(gpTransB(cbind(Ss[,,2],0*Ss[,,2]),round(Infects2,10),seed+36))
    S[,c("S3","E3")]=cbind(0,S[,"E3"]) + colSums(gpTransB(cbind(Ss[,,3],0*Ss[,,3]),round(Infects3,10),seed+37))
    S[,c("S4","E4")]=cbind(0,S[,"E4"]) + colSums(gpTransB(cbind(Ss[,,4],0*Ss[,,4]),round(Infects4,10),seed+38))
    S[,c("S5","E5")]=cbind(0,S[,"E5"]) + colSums(gpTransB(cbind(Ss[,,5],0*Ss[,,5]),round(Infects5,10),seed+39))
    
    return(S)
  }; FUN=m3iter
  
  #Implement initial closure and changes in testing probability
  closeinappl=function(parms,TS,tm=dim(TS)[3],delayInit=10){
    
    Nta=colSums(t(t(TS[,"Nt",])));
    
    omgs=cbind(parms[,"omg"],0);
    rampdays1=tm-(which(Nta>=inCstart)[1]+158-1) #Number of days past the beginning of the ramp-down
    parms[,"omg"]=pmax(omgs[,1]*exp(-parms[1,"zeta"]*rampdays1),omgs[,2]) #Ramp-down the omega value from initial to 2nd wave value
    if( (max(Nta)>=inCstart) & (tm-which(Nta>=inCstart)[1])<(158-1) ) { parms[,"omg"]=omgs[,1];} #NPI adherence with risk perception coeff omg_0 before ramp down begins
    if((max(Nta)<inCstart)) parms[,"omg"]=0; #No NPI adherence before start of pandemic (as inCstart is March 10th when total cases>=50)
    
    return(parms);
  }
  
  #Pulls info on what timestep it is
  gettime=function(TS,tm=dim(TS)[3]){
    
    tsNt=colSums(t(t(TS[,"Nt",])));
    #print(tsNt[tsNt>=inCstart]);
    if(max(tsNt)>=inCstart) {tepi<-length(tsNt[tsNt>=inCstart]) -1} else {tepi<--1} #Calc days into epidemic (since >=50 cases)
    if(max(tsNt)>=inCstart_s) {tstart_s<-length(tsNt[tsNt>=inCstart_s]) -1} else {tstart_s<--1} #Calc days past trigger date for school closures (==tepi for no counterfactuals)
    if(max(tsNt)>=inCstart_w) {tstart_w<-length(tsNt[tsNt>=inCstart_w]) -1} else {tstart_w<--1} #Calc days past trigger date for work closures (==tepi for no counterfactuals)
    
    #print(c(tepi, topen_w));
    
    timeinfo=c(tm, tepi, tstart_s, tstart_w)
    return(timeinfo);
  }

  #Function to implement simulations. InitInf sets which stages the initially infected people are in, InitInf=c(6:10) indicates all initial infections correspond to individuals who are exposed, but split between different age classes
  msim3=function(FUN,parms,Trun=365,seed0=11,plotgive=TRUE,InitInf=c(6:10)){
    #assign initial infections
    L=nrow(parms); nI=InitInf;
    
    set.seed(seed0);
    NsA=sapply(popadj2020, function(x) rmultinom(1,x,(c(3019645,3475990,3855065,2505545,592260)/13448505)));  # pop age demographics, StatCan Census 2016
    
    #Intial infections across age classes
    set.seed(seed0+1);
    Infs1=NI0_1=apply(cbind(NsA[1,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[1]),x[1],prob=x[2]/length(nI[1])));
    set.seed(seed0+2);
    Infs2=NI0_2=apply(cbind(NsA[2,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[2]),x[1],prob=x[2]/length(nI[2])));
    set.seed(seed0+3);
    Infs3=NI0_3=apply(cbind(NsA[3,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[3]),x[1],prob=x[2]/length(nI[3])));
    set.seed(seed0+4);
    Infs4=NI0_4=apply(cbind(NsA[4,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[4]),x[1],prob=x[2]/length(nI[4])));
    set.seed(seed0+5);
    Infs5=NI0_5=apply(cbind(NsA[5,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[5]),x[1],prob=x[2]/length(nI[5])));
    
    #Create object to store simulation results
    TS=array(0,dim=c(L,length(Snm),1)); 
    TS[,c(1,InitInf[1]),1]=cbind(NsA[1,]-NI0_1,Infs1); TS[,c(2,InitInf[2]),1]=cbind(NsA[2,]-NI0_2,Infs2); TS[,c(3,InitInf[3]),1]=cbind(NsA[3,]-NI0_3,Infs3); TS[,c(4,InitInf[4]),1]=cbind(NsA[4,]-NI0_4,Infs4);	TS[,c(5,InitInf[5]),1]=cbind(NsA[5,]-NI0_5,Infs5);
    TS[,57:59,1]=0; #Set Cw=Cs=C=0
    colnames(TS)=Snm; 
    #Modify travel probability as needed
    Mn=round(parms[1,"Mfact"]*parms[,-(1:NP)]*(1-diag(L))); diag(Mn)=parms[,"N"]-colSums(Mn); parms[,-(1:NP)]=Mn%*%diag(1/colSums(Mn));
    
    #Implement simulation
    set.seed(seed0+6); Seeds=rnorm(Trun,1e6,1e6);
    for(i in 2:Trun) TS=abind(TS,FUN(TS[,,i-1],closeinappl(parms,TS),gettime(TS),seed=Seeds[i]),along=3);  TS[,"C",]<-apply(TS[,c("Cw","Cs"),],c(1,3),sum);
    
    #Different levels of aggregation in model output. plotgive=3 or 3.5 give shortest output form (tracking only infections and costs) in integer format to reduce output size
    if(plotgive=="TS") return(TS); if(plotgive==TRUE) return(Resagg(TS,parms,plotgive==1));
    if(plotgive%in%c(3,3.5)){ TS[,"C",]=TS[,"C",]*matrix(parms[,"N"],nrow(parms),Trun); TSn=apply(TS,c(2,3),sum); out=matrix(as.integer(TSn),nrow(TSn),ncol(TSn)); if(plotgive==3.5) return(out[c(1:5,59),]); return(rbind(out,colMeans(TS[,"C",]>0))); }
    #In fitting also tracked mean and variance in proportion distancing (omgs) and proportion of infections in Toronto (fracTor)
    if(plotgive=="fit"){ omgs=t(apply(TS[,"VD",], 2, function(x) meansd(x,parms[,"N"]))); #print(c(omgs));
      propCits=t(TS[,"Nt",]); states=t(apply(TS,c(2,3),sum)); return(cbind(states,omgs,propCits)); }
  }

#Run multiple sims based on our best parm sets
nrun<-5
runlength<-365

for (run in 1:nrun) {
  
  #Run a sim based on a fitted parameter set
  x_region=(msim3(m3iter,parms,runlength,plotgive="TS",seed0=fit*run*1e4,InitInf=c(6:10)));
  
  #Quick plot
  agedat.plot<-agedat[agedat$Date<="2020-08-15",]
  dat=data.frame(ts50=1:runlength-min(which(colSums(colSums(x_region[,c("Nt1","Nt2","Nt3","Nt4","Nt5"),]))>50)), Nt=colSums(colSums(x_region[,c("Nt1","Nt2","Nt3","Nt4","Nt5"),])))
  dat=dat[dat$ts50<=max(agedat.plot$ts50),]
  plot(agedat.plot$ts50, rowSums(agedat.plot[,c( "newK.1", "newK.2", "newK.3", "newK.4", "newK.5")]),ylim=range(c(dat$Nt-dplyr::lag(dat$Nt), rowSums(agedat.plot[,c( "newK.1", "newK.2", "newK.3", "newK.4", "newK.5")])), na.rm=T))
  lines(dat$ts50, dat$Nt-dplyr::lag(dat$Nt), col="blue")

  #Create all the summary data that would normally be generated with plotgive="TS" to use for the provincial dataset
  omgs=t(apply(x_region[,"VD",], 2, function(x) meansd(x,parms[,"N"])));
  
  propCits=t(x_region[,"Nt",]); states=t(apply(x_region,c(2,3),sum));  
  x_prov=cbind(states,omgs,propCits);
  
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
  df_region$sick<-rowSums(df_region[,colnames(df_region) %in% c(E,Di,Da)]); df_region$know.current<-rowSums(df_region[,colnames(df_region) %in% Tk]); df_region$recovered<-rowSums(df_region[,colnames(df_region) %in% R])
  
  df_region$known1<-rowSums(df_region[, colnames(df_region) %in% Tk[1:4]]);df_region$known2<-rowSums(df_region[, colnames(df_region) %in% Tk[5:8]]);df_region$known3<-rowSums(df_region[, colnames(df_region) %in% Tk[9:12]]);df_region$known4<-rowSums(df_region[, colnames(df_region) %in% Tk[13:16]]);df_region$known5<-rowSums(df_region[, colnames(df_region) %in% Tk[17:20]]);
  df_region$sick1<-rowSums(df_region[, colnames(df_region) %in% sick1]); df_region$sick2<-rowSums(df_region[, colnames(df_region) %in% sick2]); df_region$sick3<-rowSums(df_region[, colnames(df_region) %in% sick3]); df_region$sick4<-rowSums(df_region[, colnames(df_region) %in% sick4]); df_region$sick5<-rowSums(df_region[, colnames(df_region) %in% sick5]);
  df_region$fit<-unique(fitparms$fit)[fit]
  df_region$run<-run
  df_region$LL<-unique(fitparms$LL[fitparms$fit==unique(fitparms$fit)[fit]])
  df_region$cftype<-cftype
  df_region$reopeningtype<-reopeningtype
  df_region$vdtype<-vdtype
  
  df_region <- df_region %>% group_by(fit, run, region) %>% mutate(new.K=Nt - lag(Nt),  
                                                                   new.K.1=Nt1 - lag(Nt1), new.K.2=Nt2 - lag(Nt2), new.K.3=Nt3 - lag(Nt3), new.K.4=Nt4 - lag(Nt4), new.K.5=Nt5 - lag(Nt5),
                                                                   new.sick=(sick + recovered) -  lag(sick + recovered),
                                                                   new.sick.1=(sick1 + R1) - lag(sick1  + R1), new.sick.2=(sick2  + R2) - lag(sick2  + R2), new.sick.3=(sick3  + R3) - lag(sick3  + R3), new.sick.4=(sick4  + R4) - lag(sick4  + R4), new.sick.5=(sick5  + R5) - lag(sick5  + R5))
  
  
  #Set all first step "new" to equal total at that time since we dont have a previous ts for comparison
  df_region$new.K[df_region$ts==min(df_region$ts)]<-df_region$Nt[df_region$ts==min(df_region$ts)]
  df_region$new.K.1[df_region$ts==min(df_region$ts)]<-df_region$Nt1[df_region$ts==min(df_region$ts)]
  df_region$new.K.2[df_region$ts==min(df_region$ts)]<-df_region$Nt2[df_region$ts==min(df_region$ts)]
  df_region$new.K.3[df_region$ts==min(df_region$ts)]<-df_region$Nt3[df_region$ts==min(df_region$ts)]
  df_region$new.K.4[df_region$ts==min(df_region$ts)]<-df_region$Nt4[df_region$ts==min(df_region$ts)]
  df_region$new.K.5[df_region$ts==min(df_region$ts)]<-df_region$Nt5[df_region$ts==min(df_region$ts)]
  df_region$new.sick[df_region$ts==min(df_region$ts)]<-df_region$sick[df_region$ts==min(df_region$ts)]
  df_region$new.sick.1[df_region$ts==min(df_region$ts)]<-df_region$sick1[df_region$ts==min(df_region$ts)]
  df_region$new.sick.2[df_region$ts==min(df_region$ts)]<-df_region$sick2[df_region$ts==min(df_region$ts)]
  df_region$new.sick.3[df_region$ts==min(df_region$ts)]<-df_region$sick3[df_region$ts==min(df_region$ts)]
  df_region$new.sick.4[df_region$ts==min(df_region$ts)]<-df_region$sick4[df_region$ts==min(df_region$ts)]
  df_region$new.sick.5[df_region$ts==min(df_region$ts)]<-df_region$sick5[df_region$ts==min(df_region$ts)]
  
  #Add ins for provincial sims
  df$ts<-0:(length(x_prov[,1])-1)
  
  df$sick<-rowSums(df[,colnames(df) %in% c(E,Di,Da)]); df$know.current<-rowSums(df[,colnames(df) %in% Tk]); df$recovered<-rowSums(df[,colnames(df) %in% R])
  
  df$known1<-rowSums(df[, colnames(df) %in% Tk[1:4]]);df$known2<-rowSums(df[, colnames(df) %in% Tk[5:8]]);df$known3<-rowSums(df[, colnames(df) %in% Tk[9:12]]);df$known4<-rowSums(df[, colnames(df) %in% Tk[13:16]]);df$known5<-rowSums(df[, colnames(df) %in% Tk[17:20]]);
  df$sick1<-rowSums(df[, colnames(df) %in% sick1]); df$sick2<-rowSums(df[, colnames(df) %in% sick2]); df$sick3<-rowSums(df[, colnames(df) %in% sick3]); df$sick4<-rowSums(df[, colnames(df) %in% sick4]); df$sick5<-rowSums(df[, colnames(df) %in% sick5]);
  df$fit<-unique(fitparms$fit)[fit]
  df$run<-run
  df$LL<-unique(fitparms$LL[fitparms$fit==unique(fitparms$fit)[fit]])
  df$cftype<-cftype
  df$reopeningtype<-reopeningtype
  df$vdtype<-vdtype
  df$new.K<-df$Nt - lag(df$Nt); 
  df$new.K.1<-df$Nt1 - lag(df$Nt1); df$new.K.2<-df$Nt2 - lag(df$Nt2); df$new.K.3<-df$Nt3 - lag(df$Nt3); df$new.K.4<-df$Nt4 - lag(df$Nt4); df$new.K.5<-df$Nt5 - lag(df$Nt5);
  df$new.K[1]<-df$Nt[1]; df$new.K.1[1]<-df$Nt1[1]; df$new.K.2[1]<-df$Nt2[1]; df$new.K.3[1]<-df$Nt3[1]; df$new.K.4[1]<-df$Nt4[1]; df$new.K.5[1]<-df$Nt5[1];
  df$new.sick<-(df$sick + df$recovered) -  lag(df$sick + df$recovered);
  df$new.sick.1<-(df$sick1 + df$R1) - lag(df$sick1  + df$R1); df$new.sick.2<-(df$sick2  + df$R2) - lag(df$sick2  + df$R2); df$new.sick.3<-(df$sick3  + df$R3) - lag(df$sick3  + df$R3); df$new.sick.4<-(df$sick4  + df$R4) - lag(df$sick4  + df$R4); df$new.sick.5<-(df$sick5  + df$R5) - lag(df$sick5  + df$R5);
  df$new.sick[1]<-df$sick[1]; df$new.sick.1[1]<-df$sick1[1]; df$new.sick.2[1]<-df$sick2[1]; df$new.sick.3[1]<-df$sick3[1]; df$new.sick.4[1]<-df$sick4[1]; df$new.sick.5[1]<-df$sick5[1];
  
  #Tag rw data so we can match based on days since 50th positive case
  model50<-df$ts[min(which(df$Nt >= 50))]
  df$ts50<-df$ts-model50
  
  #Tag rw data so we can match based on days since 50th positive case
  df_region.sum <- df_region %>%
    group_by(ts, fit, run) %>%
    dplyr::summarize(cases.K=sum(Nt))
  
  model50_region<-df_region$ts[min(which(df_region.sum$cases.K >= 50))]
  df_region$ts50<-df_region$ts-model50_region
  
  if (model50_region!=model50) { print("ts50 problem");}

  df.reduced<-df[,c(56:67,117:149)]
  df_region.reduced<-df_region[,56:99]
  
  if (run==1) {
    df.prov.fit<-df.reduced 
    df.region.fit<-df_region.reduced
  } else {
    df.prov.fit<-rbind(df.prov.fit, df.reduced)
    df.region.fit<-rbind(df.region.fit, df_region.reduced)
  }
  
}

df.fitparms<-as.data.frame(fitparms)
df.fitparms$fit<-unique(fitparms$fit)[fit]
df.fitparms$LL<-unique(fitparms$LL[fitparms$fit==unique(fitparms$fit)[fit]])
df.fitparms$cftype<-cftype
df.fitparms$reopeningtype<-reopeningtype
df.fitparms$vdtype<-vdtype

if (goodfits==1){
  df.prov.comb<-df.prov.fit
  df.region.comb<-df.region.fit
  df.fitparms.comb<-df.fitparms
} else {
  df.prov.comb<-rbind(df.prov.comb, df.prov.fit)
  df.region.comb<-rbind(df.region.comb, df.region.fit)
  df.fitparms.comb<-rbind(df.fitparms.comb, df.fitparms)
}

print(length(unique(df.prov.comb$fit)))

goodfits<-goodfits+1;

} 

saveRDS(df.region.comb, file = sprintf("simsREGIONS_all_COUNTERFACTUAL_%s_%s_%s_%s_V2.rds", codeversion, cftype, reopeningtype, vdtype));
saveRDS(df.prov.comb, file = sprintf("simsPROVINCE_all_COUNTERFACTUAL_%s_%s_%s_%s_V2.rds", codeversion, cftype, reopeningtype, vdtype));
saveRDS(df.fitparms.comb, file = sprintf("simsFITPARMS_all_COUNTERFACTUAL_%s_%s_%s_%s_V2.rds", codeversion,  cftype, reopeningtype, vdtype));

tfin<-Sys.time()
print(tfin-tstart);

