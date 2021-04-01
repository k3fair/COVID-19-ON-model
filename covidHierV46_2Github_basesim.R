#     Script aggregates parameter sets from fitting algorithm and runs base simulations
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

library(abind); library(plyr); library(tidyverse); library(beepr);

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

tstart<-Sys.time()

#Set name of model version for naming of output files
codeselect<-"v46v2phuNUDGE"

path = getwd() #Set to wherever output files from fitting are saved
#Find all output files from parameter fitting
file.names <- dir(path, pattern =sprintf("parmfit_y_%s_", codeselect), full.names=TRUE)
file.names.x <- dir(path, pattern =sprintf("parmfit_x_%s_", codeselect), full.names=TRUE)

#Counts number of fits that meet our criteria for being a reasonably good fit
goodfits<-1;

for (fit in 1:length(file.names))
{

  y=readRDS(file.names[fit]);
  x.dat=readRDS(file.names.x[fit]);
  
  if (y$objective<1035) #Throw out any parm combinations that don't give a reasonably good fit (reduces runtime by throwing out any parameter sets where the cost function value from the fitting is quite high)
  {
    print(y$objective)
    
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
    
    #Function to adjust for census region-level differences in transmission probability (xi_k values in paper):
    BetaMod=function(coefs,vals=regionid$phunum[order(regionid$census.region)],base=rep(1,49)) {base<-coefs[vals]; return(base);}
    
    #Note: some parameters named here are artefacts of a previous version of the model. They are retained here to avoid breaking analysis code which requires input files to have a specific dimension but do not impact simulations
    parms=cbind(N=popadj2020,Mfact=popratio*2,s=0.2,Tg_w=1,Tl_w=exp(-9.14),Tg_s=1,Tl_s=exp(-9.14),beta=BetaMod(rep(1,34))*exp(0)/mean(popadj2020),
                eps_w=0.5, eps_s=0.5,eps_h=0.5,eps_o=5,omg=exp(9.5),r=0.19,eta=0.8,alpha=0.4,sig=0.4,rho=0.67,n0=0.0001,
                atravel=0.5, a1=1,a2=1.75,a3=1,a4=5,a5=35, reopen_w=0.5, reopen_s=0.2, B=0.3, phi=80, zeta=0.025, 
                tau0=0.035, psi1=0.025,  psi2=0.025,  psi3=0.025,  psi4=0.025,  psi5=0.1, kappa=0.3, taumax=0.5, boost=0.5, DATA$Msave[-50,]);

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
    
    #emergencygap=3 #Gap between closure of schools on March 14th and state of emergency declared on March 17th
    omggap=144 #Gap between day total cases >= 50 (Mar 10th) and day we switch to wave 2 omega value (Aug 1st)
    Resurge_w=50 #Duration of additional closures 
    Resurge_s=50
    NP=ncol(parms)-nrow(parms)
    
    #Total cases corresponding to Mar 10 2020, our "start date"
    inCstart_w=50;
    inCstart_s=50;
    
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
      if(timeinfo[3]==-1 || timeinfo[2]<schoolclosuregap1) closed_s = 0*parms[,"eps_s"]; #We have not yet reached >=50 cases (prior to Mar 10) or >=50 cases reached, schools not yet closed (Mar 10 - Mar 13)
      if(timeinfo[3]>=schoolclosuregap1 && timeinfo[3]<schoolclosuregap1+inClen_s1) closed_s = parms[,"eps_s"]; #Initial school closure started (March 14-Sept 7)
      if(timeinfo[3]>=schoolclosuregap1+inClen_s1 && timeinfo[3]< schoolclosuregap2) closed_s =parms[1,"eps_s"]*parms[,"reopen_s"]; #Schools reopen with covid regs in place
      if (timeinfo[3]>=schoolclosuregap2) {
        closed_s <- ifelse(timeinfo[3]<schoolclosuregap2+inClen_s2, parms[,"eps_s"], parms[1,"eps_s"]*parms[,"reopen_s"]) #Christmas holiday to Boxing day lockdown to stay at home, with staggered reopening
      }
      if(timeinfo[4]==-1 || timeinfo[2]<workclosuregap1) closed_w = 0*parms[,"eps_w"]; #We have not yet reached >=50 cases (prior to Mar 10) or >=50 cases reached, workplaces not yet closed (Mar 10 - Mar 24)
      if(timeinfo[4]>=workclosuregap1 && timeinfo[4]<workclosuregap1+inClen_w1) {closed_w = parms[,"eps_w"];}#Initial workplace closure started (March 25-June 11th)
      if(timeinfo[4]>=workclosuregap1+inClen_w1 && timeinfo[4]< workclosuregap2) closed_w =parms[1,"eps_w"]*parms[,"reopen_w"]; #Workplaces reopen with covid regs in place
      if (timeinfo[4]>=workclosuregap2) {
        closed_w <- ifelse(timeinfo[4]<workclosuregap2+inClen_w2, parms[,"eps_w"], parms[1,"eps_w"]*parms[,"reopen_w"]) #Boxing day lockdown to stay at home, with staggered reopening
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

      #Base infection probabilities
      Infect1 = scomp*parms[,"a1"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[1,]) + (1-closed_s)%*%t(ageSpecifics_s[1,]) + (1-parms[1,"eps_h"]*S[, "VD"])%*%t(ageSpecifics_h[1,]) + (1-parms[1,"eps_o"]*S[, "VD"])%*%t(ageSpecifics_o[1,]))
      Infect2 = scomp*parms[,"a2"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[2,]) + (1-closed_s)%*%t(ageSpecifics_s[2,]) + (1-parms[1,"eps_h"]*S[, "VD"])%*%t(ageSpecifics_h[2,]) + (1-parms[1,"eps_o"]*S[, "VD"])%*%t(ageSpecifics_o[2,]))
      Infect3 = scomp*parms[,"a3"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[3,]) + (1-closed_s)%*%t(ageSpecifics_s[3,]) + (1-parms[1,"eps_h"]*S[, "VD"])%*%t(ageSpecifics_h[3,]) + (1-parms[1,"eps_o"]*S[, "VD"])%*%t(ageSpecifics_o[3,]))
      Infect4 = scomp*parms[,"a4"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[4,]) + (1-closed_s)%*%t(ageSpecifics_s[4,]) + (1-parms[1,"eps_h"]*S[, "VD"])%*%t(ageSpecifics_h[4,]) + (1-parms[1,"eps_o"]*S[, "VD"])%*%t(ageSpecifics_o[4,]))
      Infect5 = scomp*parms[,"a5"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[5,]) + (1-closed_s)%*%t(ageSpecifics_s[5,]) + (1-parms[1,"eps_h"]*S[, "VD"])%*%t(ageSpecifics_h[5,]) + (1-parms[1,"eps_o"]*S[, "VD"])%*%t(ageSpecifics_o[5,]))
      
      #Class-specific infection modifiers.
      modK=1-parms[1,"eta"]; modS=(1/parms[1,"s"])-1; 
      
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

      # Overall infection Pr is then 1-Pr(avoid infection by anyone)
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
      if( (max(Nta)>=inCstart) & (tm-which(Nta>=inCstart)[1])<(158-1) ) { parms[,"omg"]=omgs[,1];} #Voluntary distancing with risk perception coeff omg_0 before ramp down begins
      if((max(Nta)<inCstart)) parms[,"omg"]=0; #No voluntary distancing before start of pandemic (as inCstart is March 10th when total cases>=50)
      
      return(parms);
    }
    
    #Pulls info on what timestep it is
    gettime=function(TS,tm=dim(TS)[3]){
      
      tsNt=colSums(t(t(TS[,"Nt",])));
      if(max(tsNt)>=inCstart) {tepi<-length(tsNt[tsNt>=inCstart]) -1} else {tepi<--1} #Calc days into epidemic (since >=50 cases)
      if(max(tsNt)>=inCstart_s) {tstart_s<-length(tsNt[tsNt>=inCstart_s]) -1} else {tstart_s<--1} #Calc days past trigger date for school closures (==tepi for no counterfactuals)
      if(max(tsNt)>=inCstart_w) {tstart_w<-length(tsNt[tsNt>=inCstart_w]) -1} else {tstart_w<--1} #Calc days past trigger date for work closures (==tepi for no counterfactuals)
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
      if(plotgive=="fit"){ omgs=t(apply(TS[,"VD",], 2, function(x) meansd(x,parms[,"N"]))); propCits=t(TS[,"Nt",]); states=t(apply(TS,c(2,3),sum)); return(cbind(states,omgs,propCits)); }
    }
    
    ### Fitting procedure
    library(nloptr); library(mgcv);
    agelabels=c("0-19", "20-39", "40-59", "60-79", "80+")
    #To aggregate census divisions by PHU
    regionlabels=unique(regionid$phu[order(regionid$phunum)])
    bincase.new=function(x) {aggregate(x, list(PHU=regionid$phunum[order(regionid$census.region)]), FUN=sum)[,2];}
    
    LLfun=function(parmsTry,parmsFit=parmnames,extras,reps=5){
      
      parmsTry0=parmsTry;
      parmsTry[parmsFit=="zeta"]=parmsTry[parmsFit=="zeta"]/1e3
      parmsTry[parmsFit=="psi1"]=parmsTry[parmsFit=="psi1"]/1e3;
      parmsTry[parmsFit=="psi2"]=parmsTry[parmsFit=="psi2"]/1e3;
      parmsTry[parmsFit=="psi3"]=parmsTry[parmsFit=="psi3"]/1e3;
      parmsTry[parmsFit=="psi4"]=parmsTry[parmsFit=="psi4"]/1e3;
      parmsTry[parmsFit=="psi5"]=parmsTry[parmsFit=="psi5"]/1e3;
      parmsTry[parmsFit=="tau0"]=parmsTry[parmsFit=="tau0"]/1e3;
      parmsTry[parmsFit=="boost"]=parmsTry[parmsFit=="boost"]/1e3;
      #Fitting these parameters on a log scale:
      parmsTry[parmsFit=="omg"]=exp(parmsTry[parmsFit=="omg"]);
      parmsTry[parmsFit=="beta"]=exp(parmsTry[parmsFit=="beta"])/mean(popadj2020);
      # Add PHU-specific transmission modifiers
      news=parmsFit%in%names(parmsB0); Bmod=rep(1,49); if(length(news)>0) Bmod=BetaMod(parmsTry[news]);
    
      parms=cbind(N=popadj2020,Mfact=popratio*2,s=0.2,Tg_w=1,Tl_w=exp(-9.14),Tg_s=1,Tl_s=exp(-9.14),beta=Bmod*1.654565e-06,
                  eps_w=0.659, eps_s=0.659,eps_h=0.659,eps_o=0.659,omg=2.28e+05,r=0.19,eta=0.8,alpha=0.4,sig=0.4,rho=0.67,n0=0.0001,
                  atravel=0.794, a1=1.08,a2=1.05,a3=1,a4=1.89,a5=41.0, reopen_w=0.5, reopen_s=0.5, B=0.3, phi=50, zeta=0.025,
                  tau0=0.5, psi1=0.025,  psi2=0.025,  psi3=0.025,  psi4=0.025,  psi5=0.025, kappa=0.418, taumax=0.8, boost=0.5, DATA$Msave[-50,]);
      
      #Throw out bad guesses for transmission probabilities (highest probability exceeds 1)
      tcheck<-2*(1+parmsTry[parmsFit=="B"])*((1/parms[1,"s"])-1)*max(Bmod)*parmsTry[parmsFit=="beta"]*max(parmsTry[parmsFit %in% c("a1", "a2", "a3", "a4", "a5")])
      if(tcheck>1) {print(sprintf("Dud transmission probability: %f, %f, %f, %f, return %f", parmsTry[parmsFit=="B"], parmsTry0[parmsFit=="beta"], max(parmsTry[parmsFit %in% c("a1", "a2", "a3", "a4", "a5")]), max(Bmod), tcheck)); return(1e8);}
      
      parmsTry=pmax(parmsTry,0); parms[,parmsFit[!news]]=matrix(parmsTry[!news],nrow(parms),length(parmsFit[!news]),byrow=TRUE);
      parms[,"beta"]=Bmod*parms[,"beta"]; #Modify our current beta value by Bmod fcn with region specific transmision modifiers
      if(!"eps_s"%in%parmsFit) parms[,"eps_s"]=parms[,"eps_w"]; #if school distancing efficacy not specified, assume it equals work distancing efficacy
      if(!"eps_h"%in%parmsFit) parms[,"eps_h"]=parms[,"eps_w"]; #if home distancing efficacy not specified, assume it equals work distancing efficacy
      if(!"eps_o"%in%parmsFit) parms[,"eps_o"]=parms[,"eps_w"]; #if other distancing efficacy not specified, assume it equals work distancing efficacy
      if(!"reopen_s"%in%parmsFit) parms[,"reopen_s"]=parms[,"reopen_w"]; #if efficacy of distanced school reopening not specified, assume it equals efficacy of distanced work reopening
      if(is.na(reps)) return(parms); #if only want parameter matrix
      
      #Need to set 2nd dimension of sim to match cols of your msim3 output (115 for geog & age output, 66 for age output)
      sim=array(dim=c(extras["Trun"],116,0)); for(i in 1:reps) sim=abind(sim, msim3(m3iter,parms,extras["Trun"],plotgive="fit",seed0=i*1e3), along=3);
      lenR=extras["dayRfin"];
      FVN=5; #Was set to 7 before for the 7 geog groupings, now set to 5 for 5 age classes
      simp=array(dim=c(lenR+1,5+FVN,0));
      #For geog plot - set to 34 for 34 PHU
      simp_g=array(dim=c(lenR+1,34+FVN,0));
      NTS=simp_g[,1:34,0];
      
      #Counter for skips where simulation does not result in 50 total cases fast enough to have a timeseries of sufficent length for comparison to actual case counts for Ontario
      skipcount=0;
      
      for(i in 1:reps){
        #If with a years sim we don't get to 50 cases fast enough, skip this
        lcheck=sim[,"Nt",i][which(sim[,"Nt",i]>49)];
        if(length(lcheck)<(lenR+1) | max(sim[,"Nt",i])<50 ) {print(sprintf("bad run %i", i)); skipcount=skipcount+1; next;}
        
        trg=sim[,,i][which(sim[,"Nt",i]>49)[1]+(0:lenR),]; if(nrow(trg)<lenR | max(is.na(trg))==1) next;

        ###Geog aggretation based on PHU
        nts=t(apply(trg[,-(1:(length(Snm)+2))],1,bincase.new));
        NTS=abind(NTS,nts,along=3);

        simp=abind(simp, cbind(trg[,"Nt"], rowSums(trg[,c(Da,Di,R)])/trg[,"Nt"], trg[,length(Snm)+1:2], rowSums(trg[,Tk]), trg[, c("Nt1", "Nt2", "Nt3", "Nt4", "Nt5")]),  along=3)
      };
      
      if(reps-skipcount<=1) {print("Too few viable runs"); return(1e8);} #Return the "bad" value if at most 1 sim reached 50 cases quickly enough

      ###For geog aggregation based on PHU
      NTS=apply(NTS,2:3,diff);

      simp[,2,][simp[,2,]==Inf]=max(simp[,2,][simp[,2,]!=Inf]);
      if(length(dim(simp))<3) return(1e6);
      
      
      #Calculations for cost function based on time series fits
      tots.age=agedat[cumsum(rowSums(agedat[,3:7]))>50,3:7]; to.age=rowSums(tots.age); #Fitting to observed cases by ageclass, no longer drop days as we remove "grey" days in pre-processing
      tots.region=regiondat[cumsum(rowSums(regiondat[,2:35]))>50,2:35]; to.region=rowSums(tots.region); #Fitting to observed cases by ageclass, no longer drop days as we remove "grey" days in pre-processing
      mobility=omgdat$distancing #select residential mobility data, convert from percentages
      tp.age=apply(simp[1:(length(to.age)+1),1,],2,diff); tpts.age=apply(simp[1:(length(to.age)+1),6:10,],2:3,diff); #Fitting modeled new cases by ageclass
      tpts.region=NTS[1:length(to.region),,]; #Fitting modelled new cases by region
      
      tpts_DV=100*simp[1:length(mobility),3,] #Pull out mobility data for days past 50th case (keep only provincial means, not sds)
      tpts_RV=simp[1:length(to.age),2,] #Pull out testing ratio data for days past 50th case

      #Age fitting
      Tlik0.age=apply(abind(tots.age,tpts.age,along=3)[1:length(to.age),,], 1:2, function(x) (x[1]-x[-1])^2)
      for (i in 1:(reps-skipcount)) { for (j in 1:ncol(tots.age))  {Tlik0.age[i,,j]<-Tlik0.age[i,,j]/(colMeans(tots.age)[j])^2}}
      Tlik.age=sum(apply(Tlik0.age, c(2,3), mean))

      #Region fitting 
      Tlik0.region=apply(abind(tots.region,tpts.region,along=3)[1:length(to.region),,], 1:2, function(x) (x[1]-x[-1])^2)
      for (i in 1:(reps-skipcount)) { for (j in 1:ncol(tots.region))  {Tlik0.region[i,,j]<-Tlik0.region[i,,j]/(colMeans(tots.region)[j])^2}}
      Tlik.region=sum(apply(Tlik0.region, c(2,3), mean))/50
      
      #Mobility
      Mlik0=apply(abind(mobility,tpts_DV,along=2)[1:length(mobility),], 1, function(x) (x[1]-x[-1])^2)
      for (i in 1:(reps-skipcount)) {Mlik0[i,]<-Mlik0[i,]/(mean(mobility))^2}
      Mlik=sum(apply(Mlik0, 2, mean))
 
      #Calculate likelihood of testing ratio (of positive cases to total infections) and discretionary distancing level
      RV=meansd(tpts_RV[20:40,])*c(1,2) #Check testing ratio over the 20th-40th day after cases exceed 50
      sta=round(RV,2);
      Rlik=dnorm(DATA$testRat,RV[1],RV[2],log=TRUE); 
      
      LL=sum(Tlik.age,Tlik.region,4*abs(Rlik),5*(RV[2]>16),Mlik);
      stats=round(c(parmsTry0,RV[1],Tlik.age,Tlik.region,Rlik,Mlik,LL,min(c(1e8,REPORT[,ncol(REPORT)-1]),na.rm=TRUE)),2); assign("REPORT",rbind(REPORT,stats),.GlobalEnv); print(tail(REPORT,1));
      
      
      if((LL<tail(stats,1) & plotgive==TRUE) | plotgive==99){
        par(mfrow=c(2,4));
        mxplt.age=length(to.age); mxplt_mob=length(mobility); 
        matplot(tpts_DV[1:mxplt_mob,]/100,type="l",xlab="Day since 50th positive case",ylab="Individual NPI adherence",lwd=1.5,lty=1,col=1,ylim=range(c(tpts_DV/100, mobility/100))); points(mobility/100,pch=16,col=2);
        matplot(tp.age[1:mxplt.age,],type="l",xlab="Day since 50th positive case",ylab="Daily reported cases",lwd=1.5,lty=1,col=1,main=paste0(sta,collapse="_"),ylim=range(to.age)); points(to.age,pch=16,col=2);
        for(i in c(1:ncol(tots.age))){ matplot(tpts.age[1:mxplt.age,i,],type="l",lwd=2,lty=1,col=8,ylab="Daily reported cases",main=agelabels[i],ylim=c(0,max(tots.age[,i]))); points(tots.age[,i],col=2,pch=16); };
        par(mfrow=c(5,7)); mxplt.region=length(to.region);
        for(i in c(1:ncol(tots.region))){ matplot(tpts.region[1:mxplt.region,i,],type="l",lwd=2,lty=1,col=8,ylab="Daily reported cases",main=regionlabels[i],ylim=c(0,max(tots.region[,i]))); points(tots.region[,i],col=2,pch=16); };
        }
      return(pmin(LL,1e8));
    }
    
    #Fitting implementation:
    finday<-max(nrow(regiondat[cumsum(rowSums(regiondat[,2:35]))>50,2:35]), nrow(agedat[cumsum(rowSums(agedat[,3:7]))>50,]), nrow(omgdat))
    nrep=5; plotgive=FALSE;  #Set plotgive to "TRUE" to see fit of model to cases by age classes and cases by PHU
    extras=c(Trun=400,dayTfin=64,dayRfin=finday,model=3.5);
    parmsB0=c(b1=1,b2=1,b3=1,b4=1,b5=1,b6=1,b7=1,b8=1,b9=1,b10=1,b11=1,b12=1,b13=1,b14=1,b15=1,b16=1,b17=1,b18=1,b19=1,b20=1,
              b21=1,b22=1,b23=1,b24=1,b25=1,b26=1,b27=1,b28=1,b29=1,b30=1,b31=1,b32=1,b33=1,b34=1);
    pfit=c(1:60)
    fitnames=c("omg","eps_w","eps_s","eps_h","eps_o","atravel", "a1","a2","a3", "a4","a5",  "reopen_w", "reopen_s", "B", "phi", "beta", "zeta",
               "tau0", "psi1","psi2","psi3", "psi4", "psi5", "kappa","taumax", "boost",
               "b1","b2","b3","b4","b5","b6","b7", "b8", "b9", "b10", "b11", "b12", "b13", "b14", "b15", "b16", "b17","b18", "b19", "b20",
               "b21","b22","b23","b24","b25","b26","b27","b28","b29","b30","b31","b32","b33","b34");
    parmsFit=fitnames[pfit];
    #Matrix to store fitting results
    REPORT=matrix(nrow=0,ncol=length(pfit)+7); colnames(REPORT)=c(fitnames[pfit],"RV","Tlik.age","Tlik.region", "Rlik","Mlik","LL","pbLL");
    #Setting initial values and counstraints on parameters
    lowerBex=c(8,0.3,0.3,0.3,0.3,0.25,rep(0.1,5), c(0,0), c(0,1), -5, 0, 0, rep(0,5), 0, 0.3, 0, rep(0,34)); upperBex=c(12,0.9,0.9,0.9,0.9,0.75, rep(10,5), c(0.75,0.75), c(0.5,150), 10, 75, 100, rep(125,5), 0.5, 0.8, 750, rep(30,34));
    bval.guess<-c(17.43,14.27,21.24,3.18,28.88,21.68,20.12,23.89,3.61,4,14.39,27.65,12.27,19.93,9.65,4.86,5.23,25.01,28.73,2.02,1.83,14.24,22.52,3.71,16,4.25,18.03,11.53,14.76,22.5,0.8,7.98,6.77,1.89)
    parmStr0=c(9.61,0.65,0.76,0.51,0.54,0.74,0.41,0.59,0.59,1.79,7.17,0.71,0.69,0.13,93,-1.28,18.69,49.77,12.92,9.87,10.18,13.36,71.22,0.31,0.46, 500, bval.guess)
    lowerBex=lowerBex[pfit]; upperBex=upperBex[pfit];
    
parmsTry=y$solution
 
if (exists("fullparm.storage")==FALSE) {
  fullparm.storage<-c(y$solution, y$objective)
} else {
  fullparm.storage<-rbind(fullparm.storage, c(y$solution, y$objective))
}
    
#Re-run LLfun to see plot of fit, and to verify all parameter input is working correctly
plotgive="TRUE"
llinfo=LLfun(parmsTry,parmsFit=parmsFit,extras=extras,reps=nrep) 
if (abs(y$objective-llinfo)>1) {print("PROBLEM WITH FIT");} 


# Final parm set, converted from y$solution for input into visualization code
news=parmsFit%in%names(parmsB0);
fitparms=cbind(N=popadj2020,Mfact=popratio*2,s=0.2,Tg_w=1,Tl_w=exp(-9.14),Tg_s=1,Tl_s=exp(-9.14),beta=BetaMod(parmsTry[news])*exp(y$solution[16])/mean(popadj2020),
               eps_w=y$solution[2],eps_s=y$solution[3],eps_h=y$solution[4],eps_o=y$solution[5],omg=exp(y$solution[1]),r=0.19,eta=0.8,alpha=0.4,sig=0.4,rho=0.67,n0=0.0001,
               atravel=y$solution[6], a1=y$solution[7],a2=y$solution[8],a3=y$solution[9],a4=y$solution[10],a5=y$solution[11], reopen_w=y$solution[12], reopen_s=y$solution[13],
               B=y$solution[14], phi=y$solution[15], zeta=y$solution[17]/1e3,
               tau0=y$solution[18]/1e3,
               psi1=y$solution[19]/1e3,  psi2=y$solution[20]/1e3,  psi3=y$solution[21]/1e3,  psi4=y$solution[22]/1e3,  psi5=y$solution[23]/1e3,
               kappa=y$solution[24], taumax=y$solution[25], boost=y$solution[26]/1e3, DATA$Msave[-50,])

parms=fitparms; # use best fit from fitting
nrun<-5

for (run in 1:nrun) {

  #Run a sim based on a fitted parameter set
  runlength<-375
  x_region=(msim3(m3iter,parms,runlength,plotgive="TS",seed0=fit*run*1e4,InitInf=c(6:10)));

  #Quick plot showing timeseries of total cases to visually verify code is running correctly
  par(mfrow=c(1,1));
  agedat.plot<-agedat[agedat$Date<"2021-02-17",]
  dat=data.frame(ts50=1:runlength-min(which(colSums(colSums(x_region[,c("Nt1","Nt2","Nt3","Nt4","Nt5"),]))>50)), Nt=colSums(colSums(x_region[,c("Nt1","Nt2","Nt3","Nt4","Nt5"),])))
  plot(dat$ts50, dat$Nt-dplyr::lag(dat$Nt), col="blue", type="l")
  lines(agedat.plot$ts50, rowSums(agedat.plot[,c( "newK.1", "newK.2", "newK.3", "newK.4", "newK.5")]))

  
  ##### Save all output from simulation runs
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
  df_region$fit<-unlist(strsplit(strsplit(file.names[fit], "_")[[1]][5], ".rds"));
  df_region$run<-run
  df_region$LL<-y$objective

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
  df$fit<-unlist(strsplit(strsplit(file.names[fit], "_")[[1]][5], ".rds"))
  df$run<-run
  df$LL<-y$objective
  df$new.K<-df$Nt - lag(df$Nt);
  df$new.K.1<-df$Nt1 - lag(df$Nt1); df$new.K.2<-df$Nt2 - lag(df$Nt2); df$new.K.3<-df$Nt3 - lag(df$Nt3); df$new.K.4<-df$Nt4 - lag(df$Nt4); df$new.K.5<-df$Nt5 - lag(df$Nt5);
  df$new.K[1]<-df$Nt[1]; df$new.K.1[1]<-df$Nt1[1]; df$new.K.2[1]<-df$Nt2[1]; df$new.K.3[1]<-df$Nt3[1]; df$new.K.4[1]<-df$Nt4[1]; df$new.K.5[1]<-df$Nt5[1];
  df$new.sick<-(df$sick + df$recovered) -  lag(df$sick + df$recovered);
  df$new.sick.1<-(df$sick1 + df$R1) - lag(df$sick1  + df$R1); df$new.sick.2<-(df$sick2  + df$R2) - lag(df$sick2  + df$R2); df$new.sick.3<-(df$sick3  + df$R3) - lag(df$sick3  + df$R3); df$new.sick.4<-(df$sick4  + df$R4) - lag(df$sick4  + df$R4); df$new.sick.5<-(df$sick5  + df$R5) - lag(df$sick5  + df$R5);
  df$new.sick[1]<-df$sick[1]; df$new.sick.1[1]<-df$sick1[1]; df$new.sick.2[1]<-df$sick2[1]; df$new.sick.3[1]<-df$sick3[1]; df$new.sick.4[1]<-df$sick4[1]; df$new.sick.5[1]<-df$sick5[1];

  #Tag actual data so we can match based on days since 50th positive case
  model50<-df$ts[min(which(df$Nt >= 50))]
  df$ts50<-df$ts-model50

  #Tag actual data so we can match based on days since 50th positive case
  df_region.sum <- df_region %>%
    group_by(ts, fit, run) %>%
    dplyr::summarize(cases.K=sum(Nt))

  model50_region<-df_region$ts[min(which(df_region.sum$cases.K >= 50))]
  df_region$ts50<-df_region$ts-model50_region

  if (model50_region!=model50) { print("ts50 problem");}

  #Throw out columns not needed for viz to reduce file size
  df.reduced<-df[,c(56:67,117:146)]
  df_region.reduced<-df_region[,56:96]

  if (run==1) {
    df.prov.fit<-df.reduced 
    df.region.fit<-df_region.reduced
  } else {
    df.prov.fit<-rbind(df.prov.fit, df.reduced)
    df.region.fit<-rbind(df.region.fit, df_region.reduced)
  }

}

df.fitparms<-as.data.frame(fitparms)
df.fitparms$fit<-unlist(strsplit(strsplit(file.names[fit], "_")[[1]][5], ".rds"))
df.fitparms$LL<-y$objective

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

 } }


saveRDS(df.region.comb, file = sprintf("simsREGIONS_%s_V2.rds", codeselect));
saveRDS(df.prov.comb, file = sprintf("simsPROVINCE_%s_V2.rds", codeselect));
saveRDS(df.fitparms.comb, file = sprintf("simsFITPARMS_%s_V2.rds", codeselect));

tfin<-Sys.time()

print(tfin-tstart)
