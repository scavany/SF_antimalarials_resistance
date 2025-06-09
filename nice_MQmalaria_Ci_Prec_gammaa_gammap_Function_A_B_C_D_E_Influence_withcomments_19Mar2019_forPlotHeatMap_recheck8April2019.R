################ MEDICINE QUALITY MODELLING ############################
## A model for malaria treatment with varying DOSE of medicine quality ##
#########################################################################

### call deSolve 
library(deSolve)

## Remove objects
rm(list = ls())

### size of array: A=2 represent the two array (two layers: sensitive, resistance) 
A <- 2

### define simulation model times
maxt <- 35         # simulation time for model 
## (NOTE: the longer of the simulation time the higher amplitude of resistance 
## but not expand the range of the dose for resistance)
starttime <- -10   # define 10 years allow for Stability 
dt <- 1/12       # time step
times <- seq(starttime, maxt, by = dt)

###---------------------------------------------------BEGIN-----parameters <- list( ---------------------------------------------------------###
### define MODEL PARAMETERS, using list (to control the variation of the variable types and sizes)
### most parameters were set from the experiance and the review, some of them were set at zero as starting values
parameters <- list(
  N=1000,                      # population size
  mui=1/50,                    # birth rate
  muo=1/50,                    # death rate
  R0=rep(2.25,A),              # basic reproduction number (set at low level)
  ps = c(0.9,0.9),             # proportion of non-immune cases that become clinical
  pr = 0.1,                    # proportion of all immune cases that become clinical
  wait_treat = 1 / 24 / 60, # 0.5,              # average wait time before treatment after clinical symptoms (days)
  omega=1/2,                   # rate of loss of immunity = 1/(average duration of immunity)
  nuC=365/10,                  # rate of loss of clinical symptoms
  nuA=365/60,                  # rate of loss of detectibility by microscopy
  nuU=365/120,                 # rate of recovery from sub-patent infection
  t_treat = 5,                 # year treatment starts for the population
  c1max=c(0.99,0.00),          # proprtion who cure on day 1 of treatment for sensitive and resistance
  c2max=c(0.99,0.00),          # proprtion who cure on day 2 of treatment for sensitive and resistance
  c3max=c(0.99,0.00),          # proprtion who cure on day 3 of treatment for sensitive and resistance
  cpmax=c(0.99,0.60),          # proprtion who cure during partner drug treatment for sensitive and resistance
  nupmax=365/2,                # clearance by partner drug after day 3 of treatment
  rho = 365/20,                # rate of recrudescence
  nu1 = 365/1,                 # 1 / total years from day 1 to day 2
  nu2 = 365/1,                 # 1 / total years from day 2 to day 3
  nu3 = 365/1,                 # 1 / total years from day 3 to day p (here assume day p = 4)
  precmax = 0.5,#1,            # prob of recrudesence given clearance if zero drug
  precmin = 0.05,#0,           # proportion who recrudesce
  thetamax = 0.77,             # proportion who are clinical on early failure
  nuDmin=c(365/27,365/10),     # rate of elimination of perfect drug after 3 days of artmeisinin ACT for sensitive and resistance
  amp=0.7,                     # relative amplitude of seasonal forcing
  phi=0,                       # week of peak in seasonal forcing
  sensC=0.95,                  # sensitivity of diganostic test to identify clinical malaria
  sensA=0.50,                  # sensitivity of diganostic test to identify asym malaria
  sensU=0.00,                  # sensitivity of diganostic test to identify undet malaria
  kA=0.7,                      # relative infectiousness of asymptomatics compared to clinicals
  kU=0.3,                      # relative infectiousness of undetected compared to clinicals
  m1=1/10^6,                   # probability of novel mutation and fixation of a resistance mutation on day1
  m2=1/10^6,                   # probability of novel mutation and fixation of a resistance mutation on day2 
  m3=1/10^6,                   # probability of novel mutation and fixation of a resistance mutation on day3 
  mp=1/10^6,                   # probability of novel mutation and fixation of a resistance mutation on day4+
  pn=1,                        # default parasitemia for sensitive/resistance
  Ds=12,                       # perfect dose (ml/kg)   (NOTE:will get ...Warning..Internal through using ODE aroung high dose vlaues (here is >=10 ml/kg))
  ## however, it is not the problem as the resistance peak around low range of doses (there is no resistance/sensitive cases above 6 ml/kg) 
  ## WHO reference for the Artemether Lumefrantrine perfect drug = 12 ml/kg
  f=0,                         # switch for follow up cohort (if f=1 meaning the switch for the cohort is on, using for varidation with the follow up cohort)
  # (with early and late treatment failure) 
  # 1) Rogers, W.O., et al. 2) Bethell  D et al. 3) GO Gbotosho, et al.
  dose=0 ,                     # starting values of dose of artmeisinin ACT drug
  kRes2Sens=rep(0, A),         # relative sensitive parasite out-complete resistance parasite when no drug
  kSens2Res=rep(0, A),         # relative resistance parasite out-complete sensitive parasite with drug
  jSens2Res1=rep(0,A),         # relative novel mutation occurs during treatment at day1 
  jSens2Res2=rep(0,A),         # relative novel mutation occurs during treatment at day2 
  jSens2Res3=rep(0,A),         # relative novel mutation occurs during treatment at day3
  jSens2Resp=rep(0,A),         # relative novel mutation occurs during treatment at day4+ 
  jSens2ResP=rep(0,A),         # relative novel mutation occurs due to remaining drug in Pf-, drug+ stage
  jSens2ResR=rep(0,A),         # relative novel mutation occurs during immune (uninfected) stage ??????
  jSens2ResS=rep(0,A),         # relative novel mutation occurs during susceptible (uninfected) stage ??????
  lam=rep(0, A)                # ?                
)
###---------------------------------------------------END-----parameters <- list( ---------------------------------------------------------###
parameters             #with A=2
parameters$pn          #for checking parameters list above 
parameters$jSens2ResS  #for checking parameters list above 



###---------------------------------------------------BEGIN-----define MODEL INITIAL CONDITIONS--------------------------------------------###
### define MODEL INITIAL CONDITIONS
initS  <- rep(1,A)*0.9*parameters$N
initT1 <- rep(1,A)*0.0
initT2 <- rep(1,A)*0.0
initT3 <- rep(1,A)*0.0
initTp <- rep(1,A)*0.0
initTr <- rep(1,A)*0.0
initIC1 <- rep(1,A)*0.0
initIA1 <- rep(1,A)*0.05*parameters$N
initIU1 <- rep(1,A)*0.0
initP <- rep(1,A)*0.0
initR <- (rep(1,A)*0.0)+parameters$N-initS-initIA1
initCumInc <-rep(1,A)*0.0
initFail <- rep(1,A)*0.0
initpositiveDay3up <-rep(1,A)*0.0
initpositiveDay1up <-rep(1,A)*0.0
initCumIC1 <-rep(1,A)*0.0
initCliFail <-rep(1,A)*0.0
# initial conditions for the main solution vector
X <-c(initS,initT1,initT2,initT3,initTp,initTr,initIC1,initIA1,initIU1,initP,initR,initCumInc,initFail,initpositiveDay3up,initpositiveDay1up,initCumIC1,initCliFail)   

# define the indices for each variable of the main soultion vector
Sindex<-1:A
T1index<-(A+1):(2*A)
T2index<-(2*A+1):(3*A)
T3index<-(3*A+1):(4*A)
Tpindex<-(4*A+1):(5*A)
Trindex<-(5*A+1):(6*A)
IC1index<-(6*A+1):(7*A)
IA1index<-(7*A+1):(8*A)
IU1index<-(8*A+1):(9*A)
Pindex<-(9*A+1):(10*A)
Rindex<-(10*A+1):(11*A)
CumIncindex<-(11*A+1):(12*A)
Failindex<-(12*A+1):(13*A)
positiveDay3upindex<-(13*A+1):(14*A)
positiveDay1upindex<-(14*A+1):(15*A)
CumIC1index<-(15*A+1):(16*A)
CliFailindex<-(16*A+1):(17*A)
###---------------------------------------------------END-----define MODEL INITIAL CONDITIONS----------------------------------------------###



###---------------------------------------------------BEGIN-----MedQual<-function( --------------------------------------------------------###
### set up a function to solve the equations
MedQual<-function(t, X, parameters) {
  with(as.list(c(X, parameters)), {
    ## define the indices for each parameter of the main soultion vector
    S <- X[Sindex]                                # S   = Uninfected Non-immune
    T1 <- X[T1index]                              # T1  = Treated Day 1 Pf+ve
    T2 <- X[T2index]                              # T2  = Treated Day 2 Pf+ve

    T3 <- X[T3index]                              # T3  = Treated Day 3 Pf+ve
    Tp <- X[Tpindex]                              # Tp  = Treated Day >=4 Pf+ve
    Tr <- X[Trindex]                              # Tr  = Pf-ve Drug+ve ***I am quite sure Tr and P are the wrong way around (Tr is pre-recrudescent, and P is Pf-ve, drug+ve)***
    IC1 <- X[IC1index]                            # IC1 = Inifcted Clinical
    IA1 <- X[IA1index]                            # IA1 = Inifcted Asymtomatic
    IU1 <- X[IU1index]                            # IU1 = Inifcted Sub-Microscopic
    P <- X[Pindex]                                # P   = Pre Recrudescent ***Likely wrong. See Tr***
    R <- X[Rindex]                                # R   = Uninfected Immune
    CumInc <- X[CumIncindex]                      # CumInc = cumulative incidence (= treat*IC1, as only those clinical who will seek treatment and then we can seen in the system)
    Fail <- X[Failindex]                          # for checking, with follow up cohort (with early and late treatment failure), if f=1 then  dFail<- fw*Tp+sensC*fw*IC1+sensA*fw*IA1+sensU*fw*IU1 
    positiveDay3up <- X[positiveDay3upindex]      # positiveDay3up = late treatment failure (consider from T3+Tp)  
    positiveDay1up <- X[positiveDay1upindex]      # positiveDay1up = early treatment failure (consider from T1+T2+T3+Tp )
    CumIC1 <- X[CumIC1index]  
    CliFail <- X[CliFailindex]
    ## define variables that related to the parameters list
    N <- S + T1 + T2 + T3 + Tp + Tr + IC1 + IA1 + IU1 + P + R
    seas<-1+amp*cos(2*pi*(t-phi))                               # sesonality
    nu <- 1/((1/nuC)+(1/nuA)+(1/nuU))                
    I <- T1+T2+T3+Tp+Tr+IC1+IA1+IU1
    beta <- R0*(muo+nu)*seas                                      
    lam <- beta*seas*(IC1+kA*(IA1+T1+T2+T3+Tp)+kU*(IU1+Tr))/N   # force of infection 
    nuD <- nuDmin / (dose / parameters$Ds + 0.01 / 365)                # clearance rate on treatment day 4+ 
    # print(nuD)
    nup <- (dose/parameters$Ds)*nupmax            # clearance rate on treatment day 4+ (partner-drug) - SHOULDN'T THIS TAKE A SIMILAR FORM to nuD?
    # print(nup)
    theta <- 1/(1+(3/365)*nuC)                                    # ?
    treat <- (t>=t_treat)*365/(wait_treat)                      # rate of starting treatment, once passed t_treat
    tf <- 7/365                                                   # 7 days after treatment
    rf <- 52                                                    # weekly follow up
    fw <- ((t>tf)*rf)*f                                         # ?
    ## from parameters list above (which set from the experiance and the review)
    ## NOW we use real parameters values of for proportion who still positive on day 1,2,3 of treatment for sensitive from paper Lancet (as base line of the sensitive) 
    ## form the estimations there is no diffferent between Africa&Asia --> this should result from the pool study (which most of the studies run during the sensitive periods)    
    r <- c(0.96,0.92,0.76,0.5)                                    # OR form Lancet paper, take part of "Artemether Lumefrantrine Africa" that ajusted with baseline parasitamia, age <1
    ## there is NO data of the cp (proportion who still positive during partner drug treatment), then we MADE IT UP!, the OR= 0.5 for cp
    africanD1 <-sum(exp(-5.25),3.39,1.7,0)                      # list of intercpt (in the form of: b0+m1X1+m3X3+m4X4) from logistic regression ,
    africanD2 <-sum(exp(-7.35),2.28,1.43,0)                     # 
    africanD3 <-sum(exp(-8.63),1.85,4.19,0)                     # while the OR of dose ml/kg persent in the r constant
    africanD4 <-sum(exp(-9.00),1.00,5.00,0)                     # for day4+ we made it up!
    k.africa <-c(africanD1,africanD2,africanD3,africanD4)       # vector of the estimation of proportion who still positive on day 1,2,3 of treatment
    
    p.africa<-matrix(0,nrow=length(dose),ncol=4)  # define empty vector of "Probability still positive"
    c.africa<-matrix(0,nrow=length(dose),ncol=4)  # define empty vector of "Probability of Clearance"
    
    ##-----loop for function of proportion who still positive on day 1,2,3,4+ of treatment with different value of doses -----## 
    for (i in 1:4){
      ## k[i]<-(exp(-i/delta))/(1-exp(-i/delta))
      p.africa[,i]<-(k.africa[i]*r[i]^dose)/(1+(k.africa[i]*r[i]^dose))
    }
    # print(paste("p.africa:",p.africa))
    ##----- end loop for function of proportion ------------------------------------------------------------------------------##
    
    ## assigned value for "Probability of Clearance" on day 1,2,3,4+
    c.africa[,1]<-1-p.africa[,1]
    c.africa[,2]<-1-p.africa[,2]/p.africa[,1]
    c.africa[,3]<-1-p.africa[,3]/p.africa[,2]
    c.africa[,4]<-1-p.africa[,4]/p.africa[,3]
    # print(paste("c.africa:",c.africa))
    ## then assigned value for "Probability of Clearance" on day 1,2,3,4+ back to the parameters list for sensitive only (mean keep the old values of the resistance)
    c1max[1] <- c.africa[,1]<-1-p.africa[,1]
    c2max[1] <- c.africa[,2]<-1-p.africa[,2]/p.africa[,1]
    c3max[1] <- c.africa[,3]<-1-p.africa[,3]/p.africa[,2]
    cpmax[1] <- c.africa[,4]<-1-p.africa[,4]/p.africa[,4]  
    
    ## here we assume precmin follows a sigmoidal relationship with lfp
    ## kprecmin<-0.7
    ## prec0<-0.025               # the lowest possible probability of recrudescnce
    ## fp0<--2                    # a 50:50 chance of recrudescnce should happen roughly at the limit of detection which is at a parasite load of 0.01% of the original load on admission
    ## ## with durp,durp,lp0 as following
    ## durp<- 10                  # days that partner is present in the blood
    ## dura<- 3                   # days that arteminsin is present in the blood
    ## lp0<-  2                   # initial log10 parasite load in whole body
    
    ## ## here we assume that degree of resistance to artemisin is gammaa then
    gammaa <- 0.75           # degree of resitance to artmeisin (currently at about 0.3)
    gammap <- 0.75           # degree of resistance to parter drug
    ## mlamax<-3                  # sensitive rate of reduction of parasites in the body by artemisinin
    ## mlpmax<-1                  # sensitive rate of reduction of parasites in the body by the partner
    
    ## #?? there are formulae for the duration of each drug in the body (dura and durp) based on the dose
    
    ## ## the formula for precmin is then given by
    ## lfp <- lp0 - (mlamax * (1 - gammaa) + mlpmax * (1 - gammap)) * dura - mlpmax * (1-gammap) * (durp-dura)
    ## precmin <- prec0 + (1 - prec0) / (1 + exp(-kprecmin * (lfp - fp0)))
    prec <- (1 - dose / Ds) * precmax + dose / Ds * precmin      # the proportion who become pre-recrudescent
    # print(paste("prec:",prec))

    
    ## incoperate degree of resistance in the system 
    c1 <- dose / parameters$Ds * c1max
    c2 <- dose / parameters$Ds * c2max
    c3 <- dose / parameters$Ds * c3max
    cp <- dose / parameters$Ds * cpmax
    c1[2] <- (1 - gammaa) * c1[1]   # gammaa = degree of resistance for day1,2,3
    c2[2] <- (1 - gammaa) * c2[1]
    c3[2] <- (1 - gammaa) * c3[1]
    cp[2] <- (1 - gammap) * cp[1]   # gammap = degree of resistance for day4+
    ## print(paste("c1:",c1))
    ## print(paste("c2:",c2))
    ## print(paste("c3:",c3))
    ## print(paste("cp:",cp))
    
    ################ begin INFLUENCES PART #############################################################
    ###------------ there are five parts of influenceS:A,B,C,D,E ------------------------------###
    ##############################################################################################
    
    ###------------ part A  ---------------------------------------------------------------------###  
    ## cross-immunity INFLUENCES
    ## giving index of ps and recalculate ps for sensitive&resistance layers
    ## Also, am I correct in thinking that we are assuming the layers are independent?
    ps[1]<-ps[1]*(1-R[2]/N[2])+pr*(R[2]/N[2])
    ps[2]<-ps[2]*(1-R[1]/N[1])+pr*(R[1]/N[1])

    ###------------ part B + D ------------------------------------------------------------------### 
    ## sensitive Out-Cmplete Resistance WhenNoDrug  INFLUENCES
    ## 1) the present of sensitive parasites will prevent infection by resistant parasites
    lam[2] <- lam[2]*(1-(1/N[2])*((nuD[1]/(365+nuD[1]))*T1[1]+(2*nuD[1]/(365+2*nuD[1]))*T2[1]+(3*nuD[1]/(365+3*nuD[1]))*T3[1]+IC1[1]+IA1[1]+IU1[1]+Tr[1]+(nuD[1]/nuD[2]*P[1])))
    ## Add this to the manuscript
    ## 2) if a resistance is co-infected with a sensitive infection, then the sensitive infection replaces the resistance one.
    kRes2Sens[1] <- 0
    kRes2Sens[2] <- ((S[1]+R[1])/N[1])*lam[1]
    ##kRes2Sens[2] <- 0
    
    ###------------ part C + D ------------------------------------------------------------------### 
    ## resistence Out-Compete Sensitive Presence Drug  INFLUENCES
    ## 1) the present of resistance parasites will prevent infection by sensitive parasites
    lam[1] <-lam[1]*(1-(1/N[1]) * ((365/(365+nuD[1]))*T1[2] +(365/(365+2*nuD[1]))*T2[2] +(365/(365+3*nuD[1]))*T3[2] +Tp[2] +P[2]))
    ## 2) if a sensitive infection is co-infected with a resistance infection, then the resistance infection replaces the sensitive one.
    kSens2Res[1] <- ((S[2]+R[2])/N[2])*lam[2]
    kSens2Res[2] <- 0 
    
    ##------------ part E  ---------------------------------------------------------------------## 
    ## novel mutation occurs during treatment INFLUENCES
    ##  modified by mi (i=1,2,3,p)
    jSens2Res1[1] <- -365*(365/(365+nuD[1]))*m1*T1[1]/N[2]*(S[2]+R[2])
    jSens2Res1[2] <- 365*(365/(365+nuD[1]))*m1*T1[1]/N[2]*(S[2]+R[2])
    jSens2Res2[1] <- -365*(365/(365+2*nuD[1]))*m2*T2[1]/N[2]*(S[2]+R[2])
    jSens2Res2[2] <- 365*(365/(365+2*nuD[1]))*m2*T2[1]/N[2]*(S[2]+R[2])
    jSens2Res3[1] <- -365*(365/(365+3*nuD[1]))*m3*T3[1]/N[2]*(S[2]+R[2])
    jSens2Res3[2] <- 365*(365/(365+3*nuD[1]))*m3*T3[1]/N[2]*(S[2]+R[2])
    jSens2Resp[1] <- -365*mp*Tp[1]/N[2]*(S[2]+R[2])
    jSens2Resp[2] <- 365*mp*Tp[1]/N[2]*(S[2]+R[2])
    jSens2ResR[1] <- 0
    jSens2ResR[2] <- -365*(365/(365+nuD[1]))*m1*T1[1]/N[2] - 365*(365/(365+2*nuD[1]))*m2*T2[1]/N[2] - 365*(365/(365+3*nuD[1]))*m3*T3[1]/N[2] - 365*mp*Tp[1]/N[2]
    jSens2ResS[1] <- 0
    jSens2ResS[2] <- -365*(365/(365+nuD[1]))*m1*T1[1]/N[2] - 365*(365/(365+2*nuD[1]))*m2*T2[1]/N[2] - 365*(365/(365+3*nuD[1]))*m3*T3[1]/N[2] - 365*mp*Tp[1]/N[2]
    jSens2ResP[1] <- 365*(365/(365+nuD[1]))*m1*T1[1]/N[2]*(S[2]+R[2]) + 365*(365/(365+2*nuD[1]))*m2*T2[1]/N[2]*(S[2]+R[2]) + 365*(365/(365+3*nuD[1]))*m3*T3[1]/N[2]*(S[2]+R[2]) + 365*mp*Tp[1]/N[2]*(S[2]+R[2])
    jSens2ResP[2] <- 0
    ################ end INFLUENCES PART #############################################################
    
    
    ## define ------ rate of change fo the system
    dS <- mui*N +omega*R -lam*S -muo*S +jSens2ResS*S
    
    dT1 <- treat*IC1 -nu1*T1 -(nuD[1]/(365+nuD[1]))*kRes2Sens*T1 -(365/(365+nuD[1]))*kSens2Res*T1 -muo*T1 +jSens2Res1
    dT2 <- (1-c1)*nu1*T1 -nu2*T2 -(2*nuD[1]/(365+2*nuD[1]))*kRes2Sens*T2 -(365/(365+2*nuD[1]))*kSens2Res*T2 -muo*T2 +jSens2Res2
    dT3 <- (1-c2)*nu2*T2 -nu3*T3 -(3*nuD[1]/(365+3*nuD[1]))*kRes2Sens*T3 -(365/(365+3*nuD[1]))*kSens2Res*T3 -muo*T3 +jSens2Res3
    dTp <- (1-theta)*(1-c3)*nu3*T3 -(nuD+nup)*Tp  -kSens2Res*Tp -muo*Tp +jSens2Resp
    
    dTr <- prec*c1*nu1*T1 +prec*c2*nu2*T2 +prec*c3*nu3*T3 +prec*nuD*Tp -rho*Tr -kRes2Sens*Tr -muo*Tr
    dIC1 <- ps*lam*S +pr*lam*R +theta*(1-c3)*nu3*T3 +rho*Tr -treat*IC1 -nuC*IC1 -kRes2Sens*IC1 -muo*IC1
    dIA1 <- lam*(1-ps)*S +lam*(1-pr)*R +(1-prec)*nuD*Tp +nuC*IC1 -nuA*IA1 -kRes2Sens*IA1 -muo*IA1
    dIU1 <- nuA*IA1 -nuU*IU1 -kRes2Sens*IU1 -muo*IU1
    
    dP <- (1-prec)*c1*nu1*T1 +(1-prec)*c2*nu2*T2 +(1-prec)*c3*nu3*T3 +nup*Tp -nuD*P -muo*P +jSens2ResP
    dR <- nuU*IU1 +nuD*P -omega*R -lam*R +kRes2Sens*(IC1+IA1+IU1+Tr)+(((nuD[1]/(365+nuD[1]))*kRes2Sens*T1)+((2*nuD[1]/(365+2*nuD[1]))*kRes2Sens*T2)+((3*nuD[1]/(365+3*nuD[1]))*kRes2Sens*T3)) +kSens2Res*Tp + ( (365/(365+nuD[1]))*kSens2Res*T1 + (365/(365+2*nuD[1]))*kSens2Res*T2 + (365/(365+3*nuD[1]))*kSens2Res*T3) -muo*R +jSens2ResR*R
    
    dCumInc <- treat*IC1

    dFail <- fw * Tp + sensC * fw * IC1 + sensA * fw * IA1 + sensU * fw * IU1
    dpositiveDay3up <- T3 + Tp
    dpositiveDay1up <-T1 + T2 + T3 + Tp
    dCumIC1 <- ps * lam * S + pr * lam * R + theta * (1 - c3) * nu3 * T3 + rho * Tr
    dCliFail <- theta * (1 - c3) * nu3 * T3 + rho * Tr
    
    ## end  ------ rate of change fo the system

    # if (t < 10 & t > 9.5) browser()
    # return the rate of change
    list(c(
      dS, 
      dT1, dT2, dT3, dTp, dTr, 
      dIC1, dIA1, dIU1, 
      dP, dR, dCumInc, dFail, dpositiveDay3up, dpositiveDay1up, dCumIC1, dCliFail
    ))
  }
  #---END------ with(as.list(c(X, parameters)), {
  )
  #---END------ with(as.list(    
}
###---------------------------------------------------END-----MedQual<-function( ----------------------------------------------------------###



## run and plot for a range of dose 
nq <- 100
qq <- rep(0,(nq+1))
dose.vec <- rep(0,(nq+1))
cinc <- rep(0,(nq+1))
cincres <- rep(0,(nq+1))
clinc <- rep(0,(nq+1))
clincres <- rep(0,(nq+1))
cliFail <- rep(0,(nq+1))
##-------- loop over dose values --------------------------------------------------------------------------------------------------##
plotS1 <- FALSE
pdf(
  ifelse(plotS1, "figs/S1_Fig.pdf", "figs/S2_Fig.pdf"),
  width = 6.27, height = 10.2
)
par(mfrow = c(6, 4))
par(mar=c(2.1, 2.1, 3.1, 1.1), oma = c(2, 2, 0, 0))
for (i in 1:(nq+1)){
  qq[i] <- (i - 1) / nq
  dose.vec[i] <- qq[i] * parameters$Ds                                  
  # parameters["dose"]<-qq[i]
  parameters$dose <- qq[i] * parameters$Ds                                  
  # ptm <- proc.time()         # Start the clock!
  # run the model
  out <- ode(y = X, times = times, func = MedQual, parms = parameters, method="vode", maxsteps = 5000)   # run the model
  # proc.time() - ptm                                                                                     # Stop the clock
  
  ##-------- post process -----------------------------------------------------------###
  t <- out[, 1]
  S <- out[, 1 + Sindex]
  T1 <- out[, 1 + T1index]
  T2 <- out[, 1 + T2index]
  T3 <- out[, 1 + T3index]
  Tp <- out[, 1 + Tpindex]
  Tr <- out[, 1 + Trindex]
  IC1 <- out[, 1 + IC1index]
  IA1 <- out[, 1 + IA1index]
  IU1 <- out[, 1 + IU1index]
  P <- out[, 1 + Pindex]
  R <- out[, 1 + Rindex]
  CumInc <- out[, 1 + CumIncindex]
  Fail <- out[, 1 + Failindex]
  positiveDay3up <- out[, 1 + positiveDay3upindex]
  positiveDay1up <- out[, 1 + positiveDay1upindex]
  CumIC1 <- out[, 1 + CumIC1index]
  CliFail <- out[, 1 + CliFailindex]
  pop <- S + T1 + T2 + T3 + Tp + Tr + IC1 + IA1 + IU1 + P + R                     # population
  ## print(pop)                                                 # for checking 
  
  ## calculate incidence
  nt <- length(t)
  inc_month <- 0*CumInc
  inc_month[2:nt,] <- CumInc[2:nt,]-CumInc[1:(nt-1),]
  inc_month <- 1000*inc_month/pop
  totinc_month <- rowSums(inc_month)
  presinc_month <- 100*inc_month[,2]/totinc_month
  
  ### calculate clinical incidence
  ##clinc_month <- 0*IC1
  ##clinc_month[2:nt,] <- IC1[2:nt,]-IC1[1:(nt-1),]
  ##clinc_month <- 1000*clinc_month/pop
  ##totclinc_month<- rowSums(clinc_month)
  ##presclinc_month<-100*clinc_month[,2]/totclinc_month
  
  ## calculate incidence
  nt <- length(t)
  clinc_month <- 0*CumIC1
  clinc_month[2:nt,] <- CumIC1[2:nt,]-CumIC1[1:(nt-1),]
  clinc_month <- 1000*clinc_month/pop
  totclinc_month <- rowSums(clinc_month)
  presclinc_month <- 100*clinc_month[,2]/totclinc_month
  
  
  ## calculate prevalence
  prev <- 100*(T1+T2+T3+Tp+parameters$sensC*IC1+parameters$sensA*IA1)/pop
  totprev <- rowSums(prev)
  prevr <- prev[,2]-prev[,1]*prev[,2]/100
  presprev <- 100*prev[,2]/totprev

  if ((i - 1) %% 5 == 0) {
    if (plotS1) {
      plot(
        t, totprev,
        type = "l", col = "black", lwd = 2, xlim = c(0, maxt),
        main = "", xlab = "", ylab = ""
      )
      lines(t,prev[,2],col="red")
      if ((i - 1) %in% seq(0, 100, 20)) mtext("% prevalence", 2, 3)
    } else {
      plot(t,presprev,type="l",col="red",lwd=2,ylim=c(0,100),xlim=c(0,maxt),
           xlab="",ylab="")
      abline(h = 50, col = "grey") 
      if ((i - 1) %in% seq(0, 100, 20)) mtext("% resistant", 2, 3)
    }
    mtext(paste0(100 * qq[i], " % API"), 3, 0, font = 2)
    if ((i - 1) %in% seq(85, 100, 5)) mtext("time (years)", 1, 3)
  }
  
  
  
  ## plot 4 plots for each dose level: 1) incidence of total and resistance (cases per 1000 per month), 2) prevalence of total and resistance (%),  
  ## 3) incidence of resistance (%), 4) prevalence of resistance (%) 
  pdf(
    paste0("figs/temporal_dynamics_dose_", qq[i], ".pdf"),
    width = 6.2, height = 4.3, pointsize = 10
  )
  par(mfrow=c(2,2))
  par(mar=c(4.1,4.1,4.1,1.1))
  plot(t,totinc_month,type="l",col="black",lwd=2,xlim=c(0,maxt),
       main="Incidence",xlab="time (years)",ylab="cases per 1000 per month")
  lines(t,inc_month[,2],col="red")
  
  plot(t,totprev,type="l",col="black",lwd=2,xlim=c(0,maxt),
       main="Prevalence",xlab="time (years)",ylab="% prevalence")
  lines(t,prev[,2],col="red")
  # par(mar=c(4.1,4.1,1.1,1.1))
  plot(t,presinc_month,type="l",col="red",lwd=2,ylim=c(0,100),xlim=c(0,maxt),
       xlab="time (years)",ylab="% resistant")
  lines(50*pop[,1]/parameters$N,type="l",col="grey")
  
  plot(t,presprev,type="l",col="red",lwd=2,ylim=c(0,100),xlim=c(0,maxt),
       xlab="time (years)",ylab="% resistant")
  lines(50*pop[,1]/parameters$N,type="l",col="grey")
  dev.off()

  pdf(
    paste0("figs/temporal_dynamics_dose_prevonly_", qq[i], ".pdf"),
    width = 6.2, height = 4.3, pointsize = 10
  )
  par(mfrow=c(2, 1))
  par(mar=c(4.1,4.1,4.1,1.1))
  plot(t,totprev,type="l",col="black",lwd=2,xlim=c(0,maxt),
       main="Prevalence",xlab="time (years)",ylab="% prevalence")
  lines(t,prev[,2],col="red")
  plot(t,presprev,type="l",col="red",lwd=2,ylim=c(0,100),xlim=c(0,maxt),
       xlab="time (years)",ylab="% resistant")
  # lines(50*pop[,1]/parameters$N,type="l",col="grey")
  abline(h = 50, col = "grey")
  dev.off()
  print(t[which(presprev > 50)[1]])
  
  ## save key output (cumulative incidence)
  cinc[i]<-sum(CumInc[length(CumInc[,2]),])-(prod(CumInc[length(CumInc[,2]),]))/((maxt-parameters$t_treat)*1000)
  cincres[i]<-CumInc[length(CumInc[,2]),2]
  clinc[i]<-sum(CumIC1[length(CumIC1[,2]),])-(prod(CumIC1[length(CumIC1[,2]),]))/((maxt-parameters$t_treat)*1000)
  clincres[i]<-CumIC1[length(CumIC1[,2]),2]
  cliFail[i]<-sum(CliFail[length(CliFail[,2]),])-(prod(CliFail[length(CliFail[,2]),]))/((maxt-parameters$t_treat)*1000)
}
##----- end loop ---------------------------------------------------------------------------------------------------------##
dev.off()

write.csv(
  cinc,
  file = paste0(
    "tt_cinc_R0_", parameters$R0[1],
    "_AWT_", parameters$wait_treat,
    "_gammaa_", parameters$gammaa,
    "_gammap_", parameters$gammap,
    ".csv"
  )
)
write.csv(
  cincres,
  file = paste0(
    "tt_cincresR0_", parameters$R0[1],
    "_AWT_", parameters$wait_treat,
    "_gammaa_", parameters$gammaa,
    "_gammap_", parameters$gammap,
    ".csv"
  )
)
write.csv(
  clinc,
  file = paste0(
    "tt_clinc_R0_", parameters$R0[1],
    "_AWT_", parameters$wait_treat,
    "_gammaa_", parameters$gammaa,
    "_gammap_", parameters$gammap,
    ".csv"
  )
)
write.csv(
  clincres,
  file = paste0(
    "tt_clincres_R0_", parameters$R0[1],
    "_AWT_", parameters$wait_treat,
    "_gammaa_", parameters$gammaa,
    "_gammap_", parameters$gammap,
    ".csv"
  )
)
write.csv(
  cliFail,
  file = paste0(
    "tt_cliFail_R0_", parameters$R0[1],
    "_AWT_", parameters$wait_treat,
    "_gammaa_", parameters$gammaa,
    "_gammap_", parameters$gammap,
    ".csv"
  )
)


############################  END  ##########################################################################################

## NOTE TO LISA: I have not really changed the model - just the plotting


## Plot cumulative incidence vs dose
pdf(
  paste0("figs/cumulative_inc_by_dose_Tx_delay_", parameters$wait_treat, ".pdf"),
  width = 3, height = 2, pointsize = 10
)
par(mar=c(4.1,5.1,1.1,1.1))
plot(
  qq, cinc - cincres,
  type = "l",
  lwd = 2, col = "black",
  xlab = "% active ingredient",
  ylab = "Cumulative\nincidence",
  ylim = c(0, max(cinc))
)
lines(
  qq, cincres,
  lwd = 2, col = "red", lty = 2
)
legend(
  "topright",
  bty = "n",
  legend = c("Sensitive", "Resistant"),
  col = c("black", "red"),
  lty = c(1,2),
  lwd = 2
)
dev.off()
print(qq[which.max(cincres)])
print(qq[which(cincres > 100)])

parnames <- c(
  "S", "T1", "T2", "T3", "Tp", "Tr", "IC1", "IA1", "IU1", "P", "R",
  "CumInc", "Fail", "positiveDay3up", "positiveDay1up", "CumIC1", "CliFail"
)

## pdf("test.pdf", width = 4, height = 34)
## par(mfrow = c(length(parnames),2))
## for (ii in seq_along(parnames)) {
##   plot(out90[,1], out90[,ii * 2],type="l", main = parnames[ii], ylim = c(0,max(out90[,ii * 2], out100[,ii * 2])))
##   plot(out100[,1], out100[,ii * 2],type="l", main = parnames[ii], ylim = c(0,max(out90[,ii * 2], out100[,ii * 2])))
## }
## dev.off()

## par(mfrow = c(4,2))
## plot(out[, 1], out[, 1 * 2],type="l", main = parnames[1], ylim = c(0, 1000))
## plot(out2[, 1], out2[, 1 * 2],type="l", main = parnames[1], ylim = c(0, 1000))
## plot(out[, 1], out[, 6 * 2],type="l", main = parnames[6], ylim = c(0, 1000))
## plot(out2[, 1], out2[, 6 * 2],type="l", main = parnames[6], ylim = c(0, 1000))
## plot(out[, 1], out[, 11 * 2],type="l", main = parnames[11], ylim = c(0, 1000))
## plot(out2[, 1], out2[, 11 * 2],type="l", main = parnames[11], ylim = c(0, 1000))
## plot(
##   out[, 1], out[, 1 * 2] + out[, 11 * 2],
##   type = "l",
##   main = paste0(parnames[1], " + ", parnames[6], " + ", parnames[11]),
##   ylim = c(0, 1000)
## )
## plot(
##   out2[, 1], out2[, 1 * 2] + out2[, 11 * 2],
##   type = "l",
##   main = paste0(parnames[1], " + ", parnames[6], " + ", parnames[11]),
##   ylim = c(0, 1000)
## )
