## 

library(foreach)
library(doParallel)
library(mvtnorm)
library(simsurv)
library(pec)
library(DescTools)
library(zoo)
library(dplyr)
library(survival)
library(quantmod)
library(tidyr)

library(SuperLearner)
library(survSuperLearner)
library(xgboost)
# library(bartMachine)
library(KernelKnn)
library(earth)
library(fdrtool)
library(simecol)
library(MASS)
library(grf)

n <- 4000
nsim <- 200 # 1000
time.grid <- 0.01
admin.cens <- 10
taus <- c(1,2,3,4)
confounders <- c("x1","x2","x3","x4","x5","x6")
cens.sl.lib <- event.sl.lib <- c("survSL.km","survSL.coxph","survSL.loglogreg","survSL.rfsrc","survSL.expreg","survSL.weibreg")
a.sl.lib <- c("SL.mean","SL.glm","SL.nnet","SL.kernelKnn","SL.rpartPrune","SL.xgboost","SL.ranger","SL.step","SL.gam","SL.glmnet","SL.earth")
weights.sl.lib <- c("SL.glm","SL.nnet","SL.earth","SL.gam","SL.glmnet","SL.ranger","SL.xgboost")

registerDoParallel(cores=102)

## continuous
# IPCW transformation
cut_ipcw_rmtlj <- function(id, a, time, event, tau, G.a0, G.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  ind <- max(which(s<tau))
  
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  
  RMTLj.obs <- do.call(cbind, lapply(1:ns, function(u)  (s[u]-pmin(time, s[u]))*(event == cause)))
  
  term1.a0 <- RMTLj.obs/matrix(rowSums(G.a0*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  # term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a0
  # term3.a0 <- t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0*(dNct - Yt*G.dHazard.a0)/(G.a0*S.a0), 1, cumsum))-
  #   matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0*(dNct - Yt*G.dHazard.a0)/(G.a0*S.a0), 1, cumsum))
  # term4.a0 <- RMTLj.a0*t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  # term5.a0 <- t(apply(RMTLj.a0*(dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  cut.a0 <- term1.a0 # +term3.a0+term4.a0-term5.a0
  cut.a0 <- cut.a0[, ind]+(tau-s[ind])*(cut.a0[, ind+1]-cut.a0[, ind])/(s[ind+1]-s[ind])
  
  term1.a1 <- RMTLj.obs/matrix(rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  # term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a1
  # term3.a1 <- t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a1*(dNct - Yt*G.dHazard.a1)/(G.a1*S.a1), 1, cumsum))-
  #   matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a1*(dNct - Yt*G.dHazard.a1)/(G.a1*S.a1), 1, cumsum))
  # term4.a1 <- RMTLj.a1*t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  # term5.a1 <- t(apply(RMTLj.a1*(dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  cut.a1 <- term1.a1 # +term3.a1+term4.a1-term5.a1
  cut.a1 <- cut.a1[, ind]+(tau-s[ind])*(cut.a1[, ind+1]-cut.a1[, ind])/(s[ind+1]-s[ind])
  
  return(cut.a0*(1-a)+cut.a1*a)
}

# BJ transformation
cut_bj_rmtlj <- function(id, a, time, event, tau, S.a0, S.a1, Sj.a0, Sj.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  ind <- max(which(s<tau))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  S.a0.tau <- S.a0[, ind]+(tau-s[ind])*(S.a0[, ind+1]-S.a0[, ind])/(s[ind+1]-s[ind])
  S.a1.tau <- S.a1[, ind]+(tau-s[ind])*(S.a1[, ind+1]-S.a1[, ind])/(s[ind+1]-s[ind])
  S.a0.Ttilde <- rowSums(S.a0*do.call(cbind, lapply(1:ns, function(u)  ifelse(time == s[u], 1, 0))))
  S.a1.Ttilde <- rowSums(S.a1*do.call(cbind, lapply(1:ns, function(u)  ifelse(time == s[u], 1, 0))))
  S.a0.tau.min.Ttilde <- (time>tau)*S.a0.tau+(time<=tau)*S.a0.Ttilde
  S.a1.tau.min.Ttilde <- (time>tau)*S.a1.tau+(time<=tau)*S.a1.Ttilde
  
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  Fj.a0.tau <- Fj.a0[, ind]+(tau-s[ind])*(Fj.a0[, ind+1]-Fj.a0[, ind])/(s[ind+1]-s[ind])
  Fj.a1.tau <- Fj.a1[, ind]+(tau-s[ind])*(Fj.a1[, ind+1]-Fj.a1[, ind])/(s[ind+1]-s[ind])
  Fj.a0.Ttilde <- rowSums(Fj.a0*do.call(cbind, lapply(1:ns, function(u)  ifelse(time == s[u], 1, 0))))
  Fj.a1.Ttilde <- rowSums(Fj.a1*do.call(cbind, lapply(1:ns, function(u)  ifelse(time == s[u], 1, 0))))
  Fj.a0.tau.min.Ttilde <- (time>tau)*Fj.a0.tau+(time<=tau)*Fj.a0.Ttilde
  Fj.a1.tau.min.Ttilde <- (time>tau)*Fj.a1.tau+(time<=tau)*Fj.a1.Ttilde
  
  RMTLj.a0 <- t(apply(Fj.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a1 <- t(apply(Fj.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a0.tau <- RMTLj.a0[, ind]+(tau-s[ind])*(RMTLj.a0[, ind+1]-RMTLj.a0[, ind])/(s[ind+1]-s[ind])
  RMTLj.a1.tau <- RMTLj.a1[, ind]+(tau-s[ind])*(RMTLj.a1[, ind+1]-RMTLj.a1[, ind])/(s[ind+1]-s[ind])
  RMTLj.a0.Ttilde <- rowSums(RMTLj.a0*do.call(cbind, lapply(1:ns, function(u)  ifelse(time == s[u], 1, 0))))
  RMTLj.a1.Ttilde <- rowSums(RMTLj.a1*do.call(cbind, lapply(1:ns, function(u)  ifelse(time == s[u], 1, 0))))
  RMTLj.a0.tau.min.Ttilde <- (time>tau)*RMTLj.a0.tau+(time<=tau)*RMTLj.a0.Ttilde
  RMTLj.a1.tau.min.Ttilde <- (time>tau)*RMTLj.a1.tau+(time<=tau)*RMTLj.a1.Ttilde
  
  cut.a0 <- (tau-pmin(time, tau))*(event==cause)+(1-(event!=0)*(time<=tau)-(time>tau))*(RMTLj.a0.tau-RMTLj.a0.tau.min.Ttilde-Fj.a0.tau.min.Ttilde*(tau-pmin(time, tau)))/S.a0.tau.min.Ttilde
  cut.a1 <- (tau-pmin(time, tau))*(event==cause)+(1-(event!=0)*(time<=tau)-(time>tau))*(RMTLj.a1.tau-RMTLj.a1.tau.min.Ttilde-Fj.a1.tau.min.Ttilde*(tau-pmin(time, tau)))/S.a1.tau.min.Ttilde
  
  return(cut.a0*(1-a)+cut.a1*a)
}

# doubly robust censoring unbiased transformation
cut_dr_rmtlj <- function(id, a, time, event, tau, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sj.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  ind <- max(which(s<tau))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  RMTLj.a0 <- t(apply(Fj.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a1 <- t(apply(Fj.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  RMTLj.obs <- do.call(cbind, lapply(1:ns, function(u)  (s[u]-pmin(time, s[u]))*(event == cause)))
  
  term1.a0 <- RMTLj.obs/matrix(rowSums(G.a0*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  # term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a0
  term3.a0 <- t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0*(dNct - Yt*G.dHazard.a0)/(G.a0*S.a0), 1, cumsum))-
    matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0*(dNct - Yt*G.dHazard.a0)/(G.a0*S.a0), 1, cumsum))
  term4.a0 <- RMTLj.a0*t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  term5.a0 <- t(apply(RMTLj.a0*(dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  cut.a0 <- term1.a0+term3.a0+term4.a0-term5.a0 # term2.a0
  cut.a0 <- cut.a0[, ind]+(tau-s[ind])*(cut.a0[, ind+1]-cut.a0[, ind])/(s[ind+1]-s[ind])
  
  term1.a1 <- RMTLj.obs/matrix(rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  # term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a1
  term3.a1 <- t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a1*(dNct - Yt*G.dHazard.a1)/(G.a1*S.a1), 1, cumsum))-
    matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a1*(dNct - Yt*G.dHazard.a1)/(G.a1*S.a1), 1, cumsum))
  term4.a1 <- RMTLj.a1*t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  term5.a1 <- t(apply(RMTLj.a1*(dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  cut.a1 <- term1.a1+term3.a1+term4.a1-term5.a1 # term2.a1
  cut.a1 <- cut.a1[, ind]+(tau-s[ind])*(cut.a1[, ind+1]-cut.a1[, ind])/(s[ind+1]-s[ind])
  
  return(cut.a0*(1-a)+cut.a1*a)
}


# separable direct
cut_dr_sep_direct_astar1_rmtlj <- function(id, a, time, event, tau, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sjbar.a0, Sj.a1, Sjbar.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  ind <- max(which(s<tau))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  S.a <- S.a0*(1-a)+S.a1*a
  
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  G.a <- G.a0*(1-a)+G.a1*a
  
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))
  Sj.a <- Sj.a0*(1-a)+Sj.a1*a
  
  Sjbar.a0 <- t(na.locf(t(ifelse(Sjbar.a0 < 1e-3, 1e-3, Sjbar.a0))))
  Sjbar.a1 <- t(na.locf(t(ifelse(Sjbar.a1 < 1e-3, 1e-3, Sjbar.a1))))
  Sjbar.a <- Sjbar.a0*(1-a)+Sjbar.a1*a
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event != 0)))
  dNjt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  dNjbart <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event > 1)))
  
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  dMjt.a0 <- dNjt-Yt*Fj.dHazard.a0
  dMjt.a1 <- dNjt-Yt*Fj.dHazard.a1
  dMjt.a <- dMjt.a0*(1-a)+dMjt.a1*a
  
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  dMt.a0 <- dNt-Yt*S.dHazard.a0
  dMt.a1 <- dNt-Yt*S.dHazard.a1
  dMt.a <- dMt.a0*(1-a)+dMt.a1*a
  
  Fjbar.dHazard.a0 <- t(apply(cbind(0, -log(Sjbar.a0)), 1, diff))
  Fjbar.dHazard.a1 <- t(apply(cbind(0, -log(Sjbar.a1)), 1, diff))
  dMjbart.a0 <- dNjbart-Yt*Fjbar.dHazard.a0
  dMjbart.a1 <- dNjbart-Yt*Fjbar.dHazard.a1
  dMjbart.a <- dMjbart.a0*(1-a)+dMjbart.a1*a
  
  # G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  # G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  
  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  Fj.a <- Fj.a0*(1-a)+Fj.a1*a
  
  Fj.a0.a1 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (Sj.a0*Sjbar.a1))[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum)) # not Fj.dHazard.a1!
  
  RMTLj.a0 <- t(apply(Fj.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a1 <- t(apply(Fj.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a <- RMTLj.a0*(1-a)+RMTLj.a1*a
  
  RMTLj.a0.a1 <- t(apply(Fj.a0.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  
  #
  term1 <- RMTLj.a0.a1*(1-a)+RMTLj.a1*a
  
  term2 <- (matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Sjbar.a1*dMjt.a/(Sjbar.a*G.a), 1, cumsum))
            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Sjbar.a1*dMjt.a/(Sjbar.a*G.a), 1, cumsum)))
  
  term3 <- -(RMTLj.a*t(apply(Sjbar.a1*dMt.a/(Sjbar.a*S.a*G.a), 1, cumsum))-t(apply(RMTLj.a*Sjbar.a1*dMt.a/(Sjbar.a*S.a*G.a), 1, cumsum)))
  
  term4 <- (matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Sjbar.a1*Fj.a*dMt.a/(Sjbar.a*S.a*G.a), 1, cumsum))
            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Sjbar.a1*Fj.a*dMt.a/(Sjbar.a*S.a*G.a), 1, cumsum)))
  
  term5 <- (a==0)*(RMTLj.a0.a1*t(apply(dMjbart.a/(S.a*G.a), 1, cumsum))-t(apply(RMTLj.a0.a1*dMjbart.a/(S.a*G.a), 1, cumsum)))
  
  term6 <- -(a==0)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0.a1*dMjbart.a/(S.a*G.a), 1, cumsum))
                    -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0.a1*dMjbart.a/(S.a*G.a), 1, cumsum)))
  
  cut <- term1+term2+term3+term4+term5+term6
  cut <- cut[, ind]+(tau-s[ind])*(cut[, ind+1]-cut[, ind])/(s[ind+1]-s[ind])
  
  return(cut)
}


# separable indirect
cut_dr_sep_indirect_astar1_rmtlj <- function(id, a, time, event, tau, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sjbar.a0, Sj.a1, Sjbar.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  ind <- max(which(s<tau))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  S.a <- S.a0*(1-a)+S.a1*a
  
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  G.a <- G.a0*(1-a)+G.a1*a
  
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))
  Sj.a <- Sj.a0*(1-a)+Sj.a1*a
  
  Sjbar.a0 <- t(na.locf(t(ifelse(Sjbar.a0 < 1e-3, 1e-3, Sjbar.a0))))
  Sjbar.a1 <- t(na.locf(t(ifelse(Sjbar.a1 < 1e-3, 1e-3, Sjbar.a1))))
  Sjbar.a <- Sjbar.a0*(1-a)+Sjbar.a1*a
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event != 0)))
  dNjt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  dNjbart <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event > 1)))
  
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  dMjt.a0 <- dNjt-Yt*Fj.dHazard.a0
  dMjt.a1 <- dNjt-Yt*Fj.dHazard.a1
  dMjt.a <- dMjt.a0*(1-a)+dMjt.a1*a
  
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  dMt.a0 <- dNt-Yt*S.dHazard.a0
  dMt.a1 <- dNt-Yt*S.dHazard.a1
  dMt.a <- dMt.a0*(1-a)+dMt.a1*a
  
  Fjbar.dHazard.a0 <- t(apply(cbind(0, -log(Sjbar.a0)), 1, diff))
  Fjbar.dHazard.a1 <- t(apply(cbind(0, -log(Sjbar.a1)), 1, diff))
  dMjbart.a0 <- dNjbart-Yt*Fjbar.dHazard.a0
  dMjbart.a1 <- dNjbart-Yt*Fjbar.dHazard.a1
  dMjbart.a <- dMjbart.a0*(1-a)+dMjbart.a1*a
  
  # G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  # G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  
  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  Fj.a <- Fj.a0*(1-a)+Fj.a1*a
  
  Fj.a0.a1 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (Sj.a0*Sjbar.a1))[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum)) # not Fj.dHazard.a1!
  
  RMTLj.a0 <- t(apply(Fj.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a1 <- t(apply(Fj.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a <- RMTLj.a0*(1-a)+RMTLj.a1*a
  
  RMTLj.a0.a1 <- t(apply(Fj.a0.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  
  #
  term1 <- RMTLj.a0*(1-a)+RMTLj.a0.a1*a
  
  term2 <- (matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(dMjt.a/(G.a)*(1-Sjbar.a1/Sjbar.a), 1, cumsum))
            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*dMjt.a/(G.a)*(1-Sjbar.a1/Sjbar.a), 1, cumsum)))
  
  term3 <- -(RMTLj.a*t(apply(dMt.a/(S.a*G.a)*(1-Sjbar.a1/Sjbar.a), 1, cumsum))-t(apply(RMTLj.a*dMt.a/(S.a*G.a)*(1-Sjbar.a1/Sjbar.a), 1, cumsum)))
  
  term4 <- (matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a*dMt.a/(S.a*G.a)*(1-Sjbar.a1/Sjbar.a), 1, cumsum))
            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a*dMt.a/(S.a*G.a)*(1-Sjbar.a1/Sjbar.a), 1, cumsum)))
  
  term5 <- -(a==0)*(RMTLj.a0.a1*t(apply(dMjbart.a/(S.a*G.a), 1, cumsum))-t(apply(RMTLj.a0.a1*dMjbart.a/(S.a*G.a), 1, cumsum)))
  
  term6 <- (a==0)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0.a1*dMjbart.a/(S.a*G.a), 1, cumsum))
                   -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0.a1*dMjbart.a/(S.a*G.a), 1, cumsum)))
  
  cut <- term1+term2+term3+term4+term5+term6
  cut <- cut[, ind]+(tau-s[ind])*(cut[, ind+1]-cut[, ind])/(s[ind+1]-s[ind])
  
  return(cut)
}


# total
ueif_dr_rmtlj_total_mt_ht <- function(id, a, time, event, tau, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sj.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  ind <- max(which(s<tau))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  RMTLj.a0 <- t(apply(Fj.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a1 <- t(apply(Fj.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNjt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == cause)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event != 0)))
  
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply((dNjt - Yt*Fj.dHazard.a0)/(G.a0), 1, cumsum))-
                                                              t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNjt - Yt*Fj.dHazard.a0)/(G.a0), 1, cumsum)))
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a0
  term3.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a0*t(apply((dNt - Yt*S.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  term4.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMTLj.a0*(dNt - Yt*S.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  term5.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0*(dNt - Yt*S.dHazard.a0)/(S.a0*G.a0), 1, cumsum))-
                                                              t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0*(dNt - Yt*S.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  ueif.a0 <- term1.a0+term2.a0-term3.a0+term4.a0+term5.a0
  ueif.a0 <- ueif.a0[, ind]+(tau-s[ind])*(ueif.a0[, ind+1]-ueif.a0[, ind])/(s[ind+1]-s[ind])
  
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply((dNjt - Yt*Fj.dHazard.a1)/(G.a1), 1, cumsum))-
                                                              t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNjt - Yt*Fj.dHazard.a1)/(G.a1), 1, cumsum)))
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a1*t(apply((dNt - Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMTLj.a1*(dNt - Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  term5.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a1*(dNt - Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum))-
                                                              t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a1*(dNt - Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  ueif.a1 <- term1.a1+term2.a1-term3.a1+term4.a1+term5.a1
  ueif.a1 <- ueif.a1[, ind]+(tau-s[ind])*(ueif.a1[, ind+1]-ueif.a1[, ind])/(s[ind+1]-s[ind])
  
  return(ueif.a1-ueif.a0)
}


## astar=1
# separable direct
ueif_sep_direct_astar1_rmtlj <- function(id, a, ps, time, event, tau, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sjbar.a0, Sj.a1, Sjbar.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  ind <- max(which(s<tau))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))
  Sjbar.a0 <- t(na.locf(t(ifelse(Sjbar.a0 < 1e-3, 1e-3, Sjbar.a0))))
  Sjbar.a1 <- t(na.locf(t(ifelse(Sjbar.a1 < 1e-3, 1e-3, Sjbar.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  Fjbar.dHazard.a0 <- t(apply(cbind(0, -log(Sjbar.a0)), 1, diff))
  Fjbar.dHazard.a1 <- t(apply(cbind(0, -log(Sjbar.a1)), 1, diff))
  
  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  Fj.a0.a1 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (Sj.a0*Sjbar.a1))[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum)) # not Fj.dHazard.a1!
  
  RMTLj.a0 <- t(apply(Fj.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a1 <- t(apply(Fj.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a0.a1 <- t(apply(Fj.a0.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event != 0)))
  dNjt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  dNjbart <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event > 1)))
  
  # a=0
  term1.a0 <- RMTLj.a0.a1
  term2.a0 <- (matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Sjbar.a1*(dNjt-Yt*Fj.dHazard.a0)/(Sjbar.a0*G.a0), 1, cumsum))
               -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Sjbar.a1*(dNjt-Yt*Fj.dHazard.a0)/(Sjbar.a0*G.a0), 1, cumsum)))
  term3.a0 <- -(RMTLj.a0*t(apply(Sjbar.a1*(dNt-Yt*S.dHazard.a0)/(Sjbar.a0*S.a0*G.a0), 1, cumsum))
                -t(apply(RMTLj.a0*Sjbar.a1*(dNt-Yt*S.dHazard.a0)/(Sjbar.a0*S.a0*G.a0), 1, cumsum)))
  term4.a0 <- (matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Sjbar.a1*Fj.a0*(dNt-Yt*S.dHazard.a0)/(Sjbar.a0*S.a0*G.a0), 1, cumsum))
               -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Sjbar.a1*Fj.a0*(dNt-Yt*S.dHazard.a0)/(Sjbar.a0*S.a0*G.a0), 1, cumsum)))
  term5.a0 <- (RMTLj.a0.a1*t(apply((dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
               -t(apply(RMTLj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  term6.a0 <- -(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  
  ueif.a0 <- term1.a0+matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(term2.a0+term3.a0+term4.a0)+matrix(bw*(a==0)-bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(term5.a0+term6.a0)
  ueif.a0 <- ueif.a0[, ind]+(tau-s[ind])*(ueif.a0[, ind+1]-ueif.a0[, ind])/(s[ind+1]-s[ind])
  
  # a=1
  term1.a1 <- RMTLj.a1
  term2.a1 <- (matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply((dNjt-Yt*Fj.dHazard.a1)/(G.a1), 1, cumsum))
               -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNjt-Yt*Fj.dHazard.a1)/(G.a1), 1, cumsum)))
  term3.a1 <- -(RMTLj.a1*t(apply((dNt-Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
                -t(apply(RMTLj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  term4.a1 <- (matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
               -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))

  ueif.a1 <- term1.a1+matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(term2.a1+term3.a1+term4.a1)
  ueif.a1 <- ueif.a1[, ind]+(tau-s[ind])*(ueif.a1[, ind+1]-ueif.a1[, ind])/(s[ind+1]-s[ind])
  
  return(ueif.a1-ueif.a0)
}


# separable indirect
ueif_sep_indirect_astar1_rmtlj <- function(id, a, ps, time, event, tau, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sjbar.a0, Sj.a1, Sjbar.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  ind <- max(which(s<tau))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))
  Sjbar.a0 <- t(na.locf(t(ifelse(Sjbar.a0 < 1e-3, 1e-3, Sjbar.a0))))
  Sjbar.a1 <- t(na.locf(t(ifelse(Sjbar.a1 < 1e-3, 1e-3, Sjbar.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  Fjbar.dHazard.a0 <- t(apply(cbind(0, -log(Sjbar.a0)), 1, diff))
  Fjbar.dHazard.a1 <- t(apply(cbind(0, -log(Sjbar.a1)), 1, diff))
  
  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  Fj.a0.a1 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (Sj.a0*Sjbar.a1))[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum)) # not Fj.dHazard.a0!
  
  RMTLj.a0 <- t(apply(Fj.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a1 <- t(apply(Fj.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a0.a1 <- t(apply(Fj.a0.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event != 0)))
  dNjt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  dNjbart <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event > 1)))
  
  # a=0
  term1.a0 <- RMTLj.a0
  term2.a0 <- (matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply((dNjt-Yt*Fj.dHazard.a0)/(G.a0)*(1-Sjbar.a1/Sjbar.a0), 1, cumsum))
               -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNjt-Yt*Fj.dHazard.a0)/(G.a0)*(1-Sjbar.a1/Sjbar.a0), 1, cumsum)))
  term3.a0 <- -(RMTLj.a0*t(apply((dNt-Yt*S.dHazard.a0)/(S.a0*G.a0)*(1-Sjbar.a1/Sjbar.a0), 1, cumsum))
                -t(apply(RMTLj.a0*(dNt-Yt*S.dHazard.a0)/(S.a0*G.a0)*(1-Sjbar.a1/Sjbar.a0), 1, cumsum)))
  term4.a0 <- (matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0*(dNt-Yt*S.dHazard.a0)/(S.a0*G.a0)*(1-Sjbar.a1/Sjbar.a0), 1, cumsum))
               -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0*(dNt-Yt*S.dHazard.a0)/(S.a0*G.a0)*(1-Sjbar.a1/Sjbar.a0), 1, cumsum)))
  term5.a0 <- -(RMTLj.a0.a1*t(apply((dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                -t(apply(RMTLj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  term6.a0 <- (matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
               -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  
  ueif.a0 <- term1.a0+matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(term2.a0+term3.a0+term4.a0)+matrix(bw*(a==0)-bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(term5.a0+term6.a0)
  ueif.a0 <- ueif.a0[, ind]+(tau-s[ind])*(ueif.a0[, ind+1]-ueif.a0[, ind])/(s[ind+1]-s[ind])
  
  # a=1
  term1.a1 <- RMTLj.a0.a1
  term2.a1 <- (matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply((dNjt-Yt*Fj.dHazard.a1)/(G.a1)*(1-Sjbar.a1/Sjbar.a1), 1, cumsum))
               -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNjt-Yt*Fj.dHazard.a1)/(G.a1)*(1-Sjbar.a1/Sjbar.a1), 1, cumsum)))
  term3.a1 <- -(RMTLj.a1*t(apply((dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(1-Sjbar.a1/Sjbar.a1), 1, cumsum))
                -t(apply(RMTLj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(1-Sjbar.a1/Sjbar.a1), 1, cumsum)))
  term4.a1 <- (matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(1-Sjbar.a1/Sjbar.a1), 1, cumsum))
               -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(1-Sjbar.a1/Sjbar.a1), 1, cumsum)))

  ueif.a1 <- term1.a1+matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(term2.a1+term3.a1+term4.a1)
  ueif.a1 <- ueif.a1[, ind]+(tau-s[ind])*(ueif.a1[, ind+1]-ueif.a1[, ind])/(s[ind+1]-s[ind])
  
  return(ueif.a1-ueif.a0)
}

lognormal <- function(data.df, betas, sigma, n)
{return(exp(sigma*qnorm(runif(n))+as.numeric((as.matrix(data.df) %*% as.matrix(betas)))))}

loglogistic <- function(data.df, betas, gamma, n)
{return(exp(gamma*(log(runif(n))-log(1-runif(n)))+as.numeric((as.matrix(data.df) %*% as.matrix(betas)))))}

# aft <- function(data.df, betas, sigma, n)
# {return(exp(sigma*rnorm(n, mean=0, sd=1)+as.numeric((as.matrix(data.df) %*% as.matrix(betas)))))}

exponenetial <- function(data.df, lambdas, betas)
{return(-log(runif(n=nrow(data.df), min=0, max=1))/(lambdas*exp(as.matrix(data.df) %*% as.matrix(betas)))[, 1])}

dgp4_competing_nph <- function(n, admin.cens, time.grid, tau){
  corr.mat <- matrix(0.5, nrow=3, ncol=3)
  diag(corr.mat) <- 1
  x <- mvrnorm(n=n, mu=rep(0,3), Sigma=corr.mat)
  x1 <- x[,1]
  x2 <- x[,2]
  x3 <- x[,3]
  x4 <- rbinom(n,1,0.5)
  x5 <- rbinom(n,1,0.5)
  x6 <- rbinom(n,1,0.5)
  ps <- 1/(1+exp(-(0.3+0.2*x1+0.3*x2+0.3*x3-0.2*x4-0.3*x5-0.2*x6))) # 1/(1+exp(-(4.3+x1+1.5*x2+1.5*x3-x4-1.5*x5-x6))) # 1/(1+exp(-(-1+x1+1.5*x2+1.5*x3-x4-1.5*x5-x6)))
  a <- rbinom(n, 1, ps) # -1.5+0.9*x1+1.2*x2+1.2*x3+1.2*x4
  data.df <- data.frame(a, x1, x2, x3, x4, x5, x6, ow.tilt.true=ps*(1-ps))
  
  Tj1a0.coef <- c(x0=-0.1, x1=0.1, x2=-0.2, x3=0.2, x4=0.1, x5=0.8, x6=-0.2)
  Tj1a1.coef <- c(x0=-0.3, x1=0.2, x2=-0.1, x3=0.4, x4=0.2, x5=0.3, x6=0.4)
  Tj2a0.coef <- c(x0=0, x1=0.1, x2=-0.3, x3=0.1, x4=0.2, x5=0.4, x6=0.3)
  Tj2a1.coef <- c(x0=-0.2, x1=0.2, x2=0.1, x3=-0.2, x4=0.4, x5=0.3, x6=0.6)
  
  data.df$Tj1.aj0.ajbar0 <- ceiling(exponenetial(data.df=data.frame(x0=1, x1, x2, x3, x4, x5, x6), lambdas=0.12, betas=Tj1a0.coef)*100)/100
  data.df$Tj1.aj0.ajbar1 <- ceiling(exponenetial(data.df=data.frame(x0=1, x1, x2, x3, x4, x5, x6), lambdas=0.12, betas=Tj1a0.coef)*100)/100
  data.df$Tj1.aj1.ajbar0 <- ceiling(exponenetial(data.df=data.frame(x0=1, x1, x2, x3, x4, x5, x6), lambdas=0.15, betas=Tj1a1.coef)*100)/100
  data.df$Tj1.aj1.ajbar1 <- ceiling(exponenetial(data.df=data.frame(x0=1, x1, x2, x3, x4, x5, x6), lambdas=0.15, betas=Tj1a1.coef)*100)/100
  
  data.df$Tj2.aj0.ajbar0 <- ceiling(exponenetial(data.df=data.frame(x0=1, x1, x2, x3, x4, x5, x6), lambdas=0.1, betas=Tj2a0.coef)*100)/100
  data.df$Tj2.aj1.ajbar0 <- ceiling(exponenetial(data.df=data.frame(x0=1, x1, x2, x3, x4, x5, x6), lambdas=0.1, betas=Tj2a0.coef)*100)/100
  data.df$Tj2.aj0.ajbar1 <- ceiling(exponenetial(data.df=data.frame(x0=1, x1, x2, x3, x4, x5, x6), lambdas=0.08, betas=Tj2a1.coef)*100)/100
  data.df$Tj2.aj1.ajbar1 <- ceiling(exponenetial(data.df=data.frame(x0=1, x1, x2, x3, x4, x5, x6), lambdas=0.08, betas=Tj2a1.coef)*100)/100
  
  data.df$Tj1 <- (1-a)*data.df$Tj1.aj0.ajbar0 + a*data.df$Tj1.aj1.ajbar1
  data.df$Tj2 <- (1-a)*data.df$Tj2.aj0.ajbar0 + a*data.df$Tj2.aj1.ajbar1
  
  # data.df$Tj1a0 <- loglogistic(data.df=data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas=Tj1a0.coef, gamma=0.2, n=n)
  # data.df$Tj1a1 <- lognormal(data.df=data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas=Tj1a1.coef, sigma=1, n=n)
  # data.df$Tj2a0 <- lognormal(data.df=data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas=Tj2a0.coef, sigma=1, n=n)
  # data.df$Tj2a1 <- loglogistic(data.df=data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas=Tj2a1.coef, gamma=0.4, n=n)
  
  data.df$C[data.df$a==0] <- lognormal(data.df=data.frame(x0=1, x1=x1[a==0], x2=x2[a==0], x3=x3[a==0], x4=x4[a==0], x5=x5[a==0], x6=x6[a==0]), betas=c(x0=2.5, x1=0.6, x2=-0.4, x3=0.7, x4=1.5, x5=1.2, x6=1.6), sigma=0.8, n=sum(a==0))
  data.df$C[data.df$a==1] <- loglogistic(data.df=data.frame(x0=1, x1=x1[a==1], x2=x2[a==1], x3=x3[a==1], x4=x4[a==1], x5=x5[a==1], x6=x6[a==1]), betas=c(x0=2, x1=0.6, x2=0.8, x3=0.5, x4=1.2, x5=1.6, x6=1.2), gamma=0.8, n=sum(a==1))
  
  data.df$C <- pmin(data.df$C, admin.cens, na.rm=TRUE)
  data.df$end.time <- ceiling(pmin(data.df$C, data.df$Tj1, data.df$Tj2)/time.grid)*time.grid # data.df$Tj2, ceil(pmin(data.df$C,data.df$Tj1,data.df$Tj2)*365.25) # round(pmin(data.df$C,data.df$Tj1,data.df$Tj2),digits=2)
  data.df$event <- with(data.df, ifelse(pmin(Tj1,Tj2)>C,0, ifelse(Tj1<Tj2,1,2)))
  data.df <- data.df[order(data.df$end.time,-data.df$event),]
  data.df$id <- 1:n
  
  data.rep.df <- do.call("rbind", replicate(1000, data.df, simplify = FALSE))
  data.rep.df$Tj1.aj0.ajbar0 <- ceiling(exponenetial(data.df=data.frame(x0=1, data.rep.df[, c("x1","x2","x3","x4","x5","x6")]), lambdas=0.12, betas=Tj1a0.coef)*100)/100
  data.rep.df$Tj1.aj0.ajbar1 <- ceiling(exponenetial(data.df=data.frame(x0=1, data.rep.df[, c("x1","x2","x3","x4","x5","x6")]), lambdas=0.12, betas=Tj1a0.coef)*100)/100
  data.rep.df$Tj1.aj1.ajbar0 <- ceiling(exponenetial(data.df=data.frame(x0=1, data.rep.df[, c("x1","x2","x3","x4","x5","x6")]), lambdas=0.15, betas=Tj1a1.coef)*100)/100
  data.rep.df$Tj1.aj1.ajbar1 <- ceiling(exponenetial(data.df=data.frame(x0=1, data.rep.df[, c("x1","x2","x3","x4","x5","x6")]), lambdas=0.15, betas=Tj1a1.coef)*100)/100
  
  data.rep.df$Tj2.aj0.ajbar0 <- ceiling(exponenetial(data.df=data.frame(x0=1, data.rep.df[, c("x1","x2","x3","x4","x5","x6")]), lambdas=0.1, betas=Tj2a0.coef)*100)/100
  data.rep.df$Tj2.aj1.ajbar0 <- ceiling(exponenetial(data.df=data.frame(x0=1, data.rep.df[, c("x1","x2","x3","x4","x5","x6")]), lambdas=0.1, betas=Tj2a0.coef)*100)/100
  data.rep.df$Tj2.aj0.ajbar1 <- ceiling(exponenetial(data.df=data.frame(x0=1, data.rep.df[, c("x1","x2","x3","x4","x5","x6")]), lambdas=0.08, betas=Tj2a1.coef)*100)/100
  data.rep.df$Tj2.aj1.ajbar1 <- ceiling(exponenetial(data.df=data.frame(x0=1, data.rep.df[, c("x1","x2","x3","x4","x5","x6")]), lambdas=0.08, betas=Tj2a1.coef)*100)/100
  
  # data.rep.df$Tj1a0 <- loglogistic(data.df=data.frame(x0=1, data.rep.df[, c("x1","x2","x3","x4","x5","x6")]), betas=Tj1a0.coef, gamma=0.2, n=1000*n)
  # data.rep.df$Tj1a1 <- lognormal(data.df=data.frame(x0=1, data.rep.df[, c("x1","x2","x3","x4","x5","x6")]), betas=Tj1a1.coef, sigma=1, n=1000*n)
  # data.rep.df$Tj2a0 <- lognormal(data.df=data.frame(x0=1, data.rep.df[, c("x1","x2","x3","x4","x5","x6")]), betas=Tj2a0.coef, sigma=1, n=1000*n)
  # data.rep.df$Tj2a1 <- loglogistic(data.df=data.frame(x0=1, data.rep.df[, c("x1","x2","x3","x4","x5","x6")]), betas=Tj2a1.coef, gamma=0.2, n=1000*n)
  
  data.rep.df$T.aj0.ajbar0 <- pmin(data.rep.df$Tj1.aj0.ajbar0, data.rep.df$Tj2.aj0.ajbar0) # round(pmin(data.df$Tj1a0,data.df$Tj2a0),digits=3)
  data.rep.df$T.aj0.ajbar1 <- pmin(data.rep.df$Tj1.aj0.ajbar1, data.rep.df$Tj2.aj0.ajbar1)
  data.rep.df$T.aj1.ajbar0 <- pmin(data.rep.df$Tj1.aj1.ajbar0, data.rep.df$Tj2.aj1.ajbar0)
  data.rep.df$T.aj1.ajbar1 <- pmin(data.rep.df$Tj1.aj1.ajbar1, data.rep.df$Tj2.aj1.ajbar1)
  
  data.rep.df$J.aj0.ajbar0 <- with(data.rep.df, ifelse(Tj1.aj0.ajbar0<Tj2.aj0.ajbar0, 1, 2))
  data.rep.df$J.aj0.ajbar1 <- with(data.rep.df, ifelse(Tj1.aj0.ajbar1<Tj2.aj0.ajbar1, 1, 2))
  data.rep.df$J.aj1.ajbar0 <- with(data.rep.df, ifelse(Tj1.aj1.ajbar0<Tj2.aj1.ajbar0, 1, 2))
  data.rep.df$J.aj1.ajbar1 <- with(data.rep.df, ifelse(Tj1.aj1.ajbar1<Tj2.aj1.ajbar1, 1, 2))
  
  data.rep.df$Tj1.aj0.ajbar0.lost.tau <- with(data.rep.df, (tau-pmin(T.aj0.ajbar0, tau))*(J.aj0.ajbar0==1))
  data.rep.df$Tj1.aj0.ajbar1.lost.tau <- with(data.rep.df, (tau-pmin(T.aj0.ajbar1, tau))*(J.aj0.ajbar1==1))
  data.rep.df$Tj1.aj1.ajbar0.lost.tau <- with(data.rep.df, (tau-pmin(T.aj1.ajbar0, tau))*(J.aj1.ajbar0==1))
  data.rep.df$Tj1.aj1.ajbar1.lost.tau <- with(data.rep.df, (tau-pmin(T.aj1.ajbar1, tau))*(J.aj1.ajbar1==1))
  
  data.df <- data.rep.df %>%
    group_by(id) %>%
    summarise_at(vars(Tj1.aj0.ajbar0.lost.tau, Tj1.aj0.ajbar1.lost.tau, Tj1.aj1.ajbar0.lost.tau, Tj1.aj1.ajbar1.lost.tau), list(name = mean)) %>%
    rename(rmtlj.aj0.ajbar0.cond.true=Tj1.aj0.ajbar0.lost.tau_name, rmtlj.aj0.ajbar1.cond.true=Tj1.aj0.ajbar1.lost.tau_name, rmtlj.aj1.ajbar0.cond.true=Tj1.aj1.ajbar0.lost.tau_name, rmtlj.aj1.ajbar1.cond.true=Tj1.aj1.ajbar1.lost.tau_name) %>%
    merge(data.df, by="id")
  
  data.df$rmtlj.sep.direct.cate.true <- data.df$rmtlj.aj1.ajbar1.cond.true-data.df$rmtlj.aj0.ajbar1.cond.true
  data.df$rmtlj.sep.indirect.cate.true <- data.df$rmtlj.aj0.ajbar1.cond.true-data.df$rmtlj.aj0.ajbar0.cond.true
  data.df$rmtlj.total.cate.true <- data.df$rmtlj.aj1.ajbar1.cond.true-data.df$rmtlj.aj0.ajbar0.cond.true
  
  # data.df$rmtlj.a0.ite.true <- with(data.df, (tau-pmin(Tj1a0, tau))*(Tj1a0<Tj2a0))
  # data.df$rmtlj.a1.ite.true <- with(data.df, (tau-pmin(Tj1a1, tau))*(Tj1a1<Tj2a1))
  # data.df$rmtlj.diff.ite.true <- data.df$rmtlj.a1.ite.true-data.df$rmtlj.a0.ite.true
  return(data.df)
}

setwd("competing1")

start.time <- Sys.time()
foreach(i=1:nsim, .combine=c, .errorhandling="remove") %dopar% {
  sim.df <- dgp4_competing_nph(n, admin.cens, time.grid, tau=taus[4]) # 1-sum(sim.df$event)/n
  
  ## first split
  train.ind <- sample(1:n, n/2, replace=FALSE) # sample training set
  train.df <- sim.df[train.ind,]
  estimation.df <- sim.df[setdiff(1:n, train.ind),]
  
  while(max(train.df$end.time[train.df$a==0 & train.df$event==0])<taus[4] | max(train.df$end.time[train.df$a==1 & train.df$event==0])<taus[4] |
        max(train.df$end.time[train.df$a==0 & train.df$event==1])<taus[4] | max(train.df$end.time[train.df$a==1 & train.df$event==1])<taus[4] |
        max(estimation.df$end.time[estimation.df$a==0 & estimation.df$event==0])<taus[4] | max(estimation.df$end.time[estimation.df$a==1 & estimation.df$event==0])<taus[4] |
        max(estimation.df$end.time[estimation.df$a==0 & estimation.df$event==1])<taus[4] | max(estimation.df$end.time[estimation.df$a==1 & estimation.df$event==1])<taus[4] |
        max(estimation.df$end.time[estimation.df$a==0 & estimation.df$event==2])<taus[4] | max(estimation.df$end.time[estimation.df$a==1 & estimation.df$event==2])<taus[4]) {
    sim.df <- dgp4_competing_nph(2*n, admin.cens, time.grid)
    train.ind <- sample(1:n, n/2, replace=FALSE) # sample training set
    train.df <- sim.df[train.ind,]
    estimation.df <- sim.df[setdiff(1:n, train.ind),]}
  
  # ps
  train.df$naive <- 1
  train.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = estimation.df$a, X = estimation.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred) # train.df$ps <- predict(glm(a~x1+x2+x3+x4,data=train.df,family=binomial()),type="response")
  train.df$iptw <- with(train.df, 1/(ps*a+(1-ps)*(1-a)))
  train.df$ow <- with(train.df, ifelse(a, 1-ps, ps))
  train.df$ow.tilt.hat <- with(train.df, ps*(1-ps))
  
  estimation.df$naive <- 1
  estimation.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = train.df$a, X = train.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred) # estimation.df$ps <- predict(glm(a~x1+x2+x3+x4,data=estimation.df,family=binomial()),type="response")
  estimation.df$iptw <- with(estimation.df, 1/(ps*a+(1-ps)*(1-a)))
  estimation.df$ow <- with(estimation.df, ifelse(a, 1-ps, ps))
  estimation.df$ow.tilt.hat <- with(estimation.df, ps*(1-ps))
  
  train.s <- sort(unique(train.df$end.time))
  train.ns <- n_distinct(train.df$end.time)
  train.ds.mat <- matrix(diff(c(0, train.s)), nrow=n/2, ncol=n_distinct(train.df$end.time), byrow=TRUE)
  train.ind <- max(which(train.s<taus[4]))
  
  estimation.s <- sort(unique(estimation.df$end.time))
  estimation.ns <- n_distinct(estimation.df$end.time)
  estimation.ds.mat <- matrix(diff(c(0, estimation.s)), nrow=n/2, ncol=n_distinct(estimation.df$end.time), byrow=TRUE)
  estimation.ind <- max(which(estimation.s<taus[4]))
  
  ## comparsion
  # single learner
  train.S.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time, event = 1*(estimation.df$event!=0), X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                      newX = data.frame(a=0, train.df[, confounders], 0*train.df[, confounders]), new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  train.S.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time, event = 1*(estimation.df$event!=0), X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                      newX = data.frame(a=1, train.df[, confounders], 1*train.df[, confounders]), new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  
  train.Sj.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time, event = 1*(estimation.df$event==1), X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                       newX = data.frame(a=0, train.df[, confounders], 0*train.df[, confounders]), new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  train.Sj.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time, event = 1*(estimation.df$event==1), X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                       newX = data.frame(a=1, train.df[, confounders], 1*train.df[, confounders]), new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  
  train.Sjbar.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time, event = 1*(estimation.df$event>1), X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                          newX = data.frame(a=0, train.df[, confounders], 0*train.df[, confounders]), new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  train.Sjbar.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time, event = 1*(estimation.df$event>1), X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                          newX = data.frame(a=1, train.df[, confounders], 1*train.df[, confounders]), new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  
  estimation.S.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time, event = 1*(train.df$event!=0), X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                           newX = data.frame(a=0, estimation.df[, confounders], 0*estimation.df[, confounders]), new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  estimation.S.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time, event = 1*(train.df$event!=0), X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                           newX = data.frame(a=1, estimation.df[, confounders], 1*estimation.df[, confounders]), new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  
  estimation.Sj.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time, event = 1*(train.df$event==1), X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                            newX = data.frame(a=0, estimation.df[, confounders], 0*estimation.df[, confounders]), new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  estimation.Sj.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time, event = 1*(train.df$event==1), X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                            newX = data.frame(a=1, estimation.df[, confounders], 1*estimation.df[, confounders]), new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  
  estimation.Sjbar.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time, event = 1*(train.df$event>1), X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                               newX = data.frame(a=0, estimation.df[, confounders], 0*estimation.df[, confounders]), new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  estimation.Sjbar.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time, event = 1*(train.df$event>1), X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                               newX = data.frame(a=1, estimation.df[, confounders], 1*estimation.df[, confounders]), new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  
  # total
  train.Fj.a0 <- t(apply(cbind(1, train.S.a0[, 1:(train.ns-1)])*t(apply(cbind(0, -log(train.Sj.a0)), 1, diff)), 1, cumsum))
  train.Fj.a1 <- t(apply(cbind(1, train.S.a1[, 1:(train.ns-1)])*t(apply(cbind(0, -log(train.Sj.a1)), 1, diff)), 1, cumsum))
  
  estimation.Fj.a0 <- t(apply(cbind(1, estimation.S.a0[, 1:(estimation.ns-1)])*t(apply(cbind(0, -log(estimation.Sj.a0)), 1, diff)), 1, cumsum))
  estimation.Fj.a1 <- t(apply(cbind(1, estimation.S.a1[, 1:(estimation.ns-1)])*t(apply(cbind(0, -log(estimation.Sj.a1)), 1, diff)), 1, cumsum))
  
  train.rmtlj.diff <- t(apply(train.ds.mat*train.Fj.a1, 1, cumsum))-t(apply(train.ds.mat*train.Fj.a0, 1, cumsum))
  train.df$rmtlj.single.learner.cate.hat <- train.rmtlj.diff[, train.ind]+(taus[4]-train.s[train.ind])*(train.rmtlj.diff[, train.ind+1]-train.rmtlj.diff[, train.ind])/(train.s[train.ind+1]-train.s[train.ind])
  
  estimation.rmtlj.diff <- t(apply(estimation.ds.mat*estimation.Fj.a1, 1, cumsum))-t(apply(estimation.ds.mat*estimation.Fj.a0, 1, cumsum))
  estimation.df$rmtlj.single.learner.cate.hat <- estimation.rmtlj.diff[, estimation.ind]+(taus[4]-estimation.s[estimation.ind])*(estimation.rmtlj.diff[, estimation.ind+1]-estimation.rmtlj.diff[, estimation.ind])/(estimation.s[estimation.ind+1]-estimation.s[estimation.ind])
  
  # separable direct a_jbar=1
  train.Fj.a0.a1 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (train.Sj.a0*train.Sjbar.a1))[, 1:(train.ns-1)])*t(apply(cbind(0, -log(train.Sj.a0)), 1, diff)), 1, cumsum))
  train.Fj.a1.a1 <- train.Fj.a1
  
  estimation.Fj.a0.a1 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (estimation.Sj.a0*estimation.Sjbar.a1))[, 1:(estimation.ns-1)])*t(apply(cbind(0, -log(estimation.Sj.a0)), 1, diff)), 1, cumsum))
  estimation.Fj.a1.a1 <- estimation.Fj.a1
  
  train.rmtlj.sep.direct.diff <- t(apply(train.ds.mat*train.Fj.a1.a1, 1, cumsum))-t(apply(train.ds.mat*train.Fj.a0.a1, 1, cumsum))
  train.df$rmtlj.sep.direct.single.learner.cate.hat <- train.rmtlj.sep.direct.diff[, train.ind]+(taus[4]-train.s[train.ind])*(train.rmtlj.sep.direct.diff[, train.ind+1]-train.rmtlj.sep.direct.diff[, train.ind])/(train.s[train.ind+1]-train.s[train.ind])
  
  estimation.rmtlj.sep.direct.diff <- t(apply(estimation.ds.mat*estimation.Fj.a1.a1, 1, cumsum))-t(apply(estimation.ds.mat*estimation.Fj.a0.a1, 1, cumsum))
  estimation.df$rmtlj.sep.direct.single.learner.cate.hat <- estimation.rmtlj.sep.direct.diff[, estimation.ind]+(taus[4]-estimation.s[estimation.ind])*(estimation.rmtlj.sep.direct.diff[, estimation.ind+1]-estimation.rmtlj.sep.direct.diff[, estimation.ind])/(estimation.s[estimation.ind+1]-estimation.s[estimation.ind])
  
  # separable indirect a_j=0
  train.df$rmtlj.sep.indirect.single.learner.cate.hat <- train.df$rmtlj.single.learner.cate.hat-train.df$rmtlj.sep.direct.single.learner.cate.hat
  estimation.df$rmtlj.sep.indirect.single.learner.cate.hat <- estimation.df$rmtlj.single.learner.cate.hat-estimation.df$rmtlj.sep.direct.single.learner.cate.hat
  
  # two survSuperLearner
  ## train
  train.Sj.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time[estimation.df$a==0], event = 1*(estimation.df$event[estimation.df$a==0]==1), X = estimation.df[estimation.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                       newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  train.Sj.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time[estimation.df$a==1], event = 1*(estimation.df$event[estimation.df$a==1]==1), X = estimation.df[estimation.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                       newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  
  train.Sjbar.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time[estimation.df$a==0], event = 1*(estimation.df$event[estimation.df$a==0]>1), X = estimation.df[estimation.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                          newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  train.Sjbar.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time[estimation.df$a==1], event = 1*(estimation.df$event[estimation.df$a==1]>1), X = estimation.df[estimation.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                          newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  
  train.all.survival.super.a0 <- survSuperLearner(time = estimation.df$end.time[estimation.df$a==0], event = 1*(estimation.df$event[estimation.df$a==0]!=0), X = estimation.df[estimation.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                  newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
  train.G.a0 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a0$cens.SL.predict)
  train.S.a0 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a0$event.SL.predict)
  
  train.all.survival.super.a1 <- survSuperLearner(time = estimation.df$end.time[estimation.df$a==1], event = 1*(estimation.df$event[estimation.df$a==1]!=0), X = estimation.df[estimation.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                  newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
  train.G.a1 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a1$cens.SL.predict)
  train.S.a1 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a1$event.SL.predict)
  
  ## estimation
  estimation.Sj.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time[train.df$a==0], event = 1*(train.df$event[train.df$a==0]==1), X = train.df[train.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                            newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  estimation.Sj.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time[train.df$a==1], event = 1*(train.df$event[train.df$a==1]==1), X = train.df[train.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                            newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  
  estimation.Sjbar.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time[train.df$a==0], event = 1*(train.df$event[train.df$a==0]>1), X = train.df[train.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                               newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  estimation.Sjbar.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time[train.df$a==1], event = 1*(train.df$event[train.df$a==1]>1), X = train.df[train.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                               newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  estimation.all.survival.super.a0 <- survSuperLearner(time = train.df$end.time[train.df$a==0], event = 1*(train.df$event[train.df$a==0]!=0), X = train.df[train.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                       newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
  estimation.G.a0 <- Winsorize(minval = 1e-3, maxval = 1, estimation.all.survival.super.a0$cens.SL.predict)
  estimation.S.a0 <- Winsorize(minval = 1e-3, maxval = 1, estimation.all.survival.super.a0$event.SL.predict)
  
  estimation.all.survival.super.a1 <- survSuperLearner(time = train.df$end.time[train.df$a==1], event = 1*(train.df$event[train.df$a==1]!=0), X = train.df[train.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                       newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
  estimation.G.a1 <- Winsorize(minval = 1e-3, maxval = 1, estimation.all.survival.super.a1$cens.SL.predict)
  estimation.S.a1 <- Winsorize(minval = 1e-3, maxval = 1, estimation.all.survival.super.a1$event.SL.predict)
  
  train.S.a0 <- t(na.locf(t(ifelse(train.S.a0 < 1e-3, 1e-3, train.S.a0))))
  train.S.a1 <- t(na.locf(t(ifelse(train.S.a1 < 1e-3, 1e-3, train.S.a1))))
  train.G.a0 <- t(na.locf(t(ifelse(train.G.a0 < 1e-3, 1e-3, train.G.a0))))
  train.G.a1 <- t(na.locf(t(ifelse(train.G.a1 < 1e-3, 1e-3, train.G.a1))))
  train.Sj.a0 <- t(na.locf(t(ifelse(train.Sj.a0 < 1e-3, 1e-3, train.Sj.a0))))
  train.Sj.a1 <- t(na.locf(t(ifelse(train.Sj.a1 < 1e-3, 1e-3, train.Sj.a1))))
  train.Sjbar.a0 <- t(na.locf(t(ifelse(train.Sjbar.a0 < 1e-3, 1e-3, train.Sjbar.a0))))
  train.Sjbar.a1 <- t(na.locf(t(ifelse(train.Sjbar.a1 < 1e-3, 1e-3, train.Sjbar.a1))))
  
  estimation.S.a0 <- t(na.locf(t(ifelse(estimation.S.a0 < 1e-3, 1e-3, estimation.S.a0))))
  estimation.S.a1 <- t(na.locf(t(ifelse(estimation.S.a1 < 1e-3, 1e-3, estimation.S.a1))))
  estimation.G.a0 <- t(na.locf(t(ifelse(estimation.G.a0 < 1e-3, 1e-3, estimation.G.a0))))
  estimation.G.a1 <- t(na.locf(t(ifelse(estimation.G.a1 < 1e-3, 1e-3, estimation.G.a1))))
  estimation.Sj.a0 <- t(na.locf(t(ifelse(estimation.Sj.a0 < 1e-3, 1e-3, estimation.Sj.a0))))
  estimation.Sj.a1 <- t(na.locf(t(ifelse(estimation.Sj.a1 < 1e-3, 1e-3, estimation.Sj.a1))))
  estimation.Sjbar.a0 <- t(na.locf(t(ifelse(estimation.Sjbar.a0 < 1e-3, 1e-3, estimation.Sjbar.a0))))
  estimation.Sjbar.a1 <- t(na.locf(t(ifelse(estimation.Sjbar.a1 < 1e-3, 1e-3, estimation.Sjbar.a1))))
  
  # total effect
  train.Fj.a0 <- t(apply(cbind(1, train.S.a0[, 1:(train.ns-1)])*t(apply(cbind(0, -log(train.Sj.a0)), 1, diff)), 1, cumsum))
  train.Fj.a1 <- t(apply(cbind(1, train.S.a1[, 1:(train.ns-1)])*t(apply(cbind(0, -log(train.Sj.a1)), 1, diff)), 1, cumsum))
  
  estimation.Fj.a0 <- t(apply(cbind(1, estimation.S.a0[, 1:(estimation.ns-1)])*t(apply(cbind(0, -log(estimation.Sj.a0)), 1, diff)), 1, cumsum))
  estimation.Fj.a1 <- t(apply(cbind(1, estimation.S.a1[, 1:(estimation.ns-1)])*t(apply(cbind(0, -log(estimation.Sj.a1)), 1, diff)), 1, cumsum))
  
  train.rmtlj.diff <- t(apply(train.ds.mat*train.Fj.a1, 1, cumsum))-t(apply(train.ds.mat*train.Fj.a0, 1, cumsum))
  train.df$rmtlj.two.learners.cate.hat <- train.rmtlj.diff[, train.ind]+(taus[4]-train.s[train.ind])*(train.rmtlj.diff[, train.ind+1]-train.rmtlj.diff[, train.ind])/(train.s[train.ind+1]-train.s[train.ind])
  
  estimation.rmtlj.diff <- t(apply(estimation.ds.mat*estimation.Fj.a1, 1, cumsum))-t(apply(estimation.ds.mat*estimation.Fj.a0, 1, cumsum))
  estimation.df$rmtlj.two.learners.cate.hat <- estimation.rmtlj.diff[, estimation.ind]+(taus[4]-estimation.s[estimation.ind])*(estimation.rmtlj.diff[, estimation.ind+1]-estimation.rmtlj.diff[, estimation.ind])/(estimation.s[estimation.ind+1]-estimation.s[estimation.ind])
  
  # separable direct a_jbar=1
  train.Fj.a0.a1 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (train.Sj.a0*train.Sjbar.a1))[, 1:(train.ns-1)])*t(apply(cbind(0, -log(train.Sj.a0)), 1, diff)), 1, cumsum))
  train.Fj.a1.a1 <- train.Fj.a1
  
  estimation.Fj.a0.a1 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (estimation.Sj.a0*estimation.Sjbar.a1))[, 1:(estimation.ns-1)])*t(apply(cbind(0, -log(estimation.Sj.a0)), 1, diff)), 1, cumsum))
  estimation.Fj.a1.a1 <- estimation.Fj.a1
  
  train.rmtlj.sep.direct.diff <- t(apply(train.ds.mat*train.Fj.a1.a1, 1, cumsum))-t(apply(train.ds.mat*train.Fj.a0.a1, 1, cumsum))
  train.df$rmtlj.sep.direct.two.learners.cate.hat <- train.rmtlj.sep.direct.diff[, train.ind]+(taus[4]-train.s[train.ind])*(train.rmtlj.sep.direct.diff[, train.ind+1]-train.rmtlj.sep.direct.diff[, train.ind])/(train.s[train.ind+1]-train.s[train.ind])
  
  estimation.rmtlj.sep.direct.diff <- t(apply(estimation.ds.mat*estimation.Fj.a1.a1, 1, cumsum))-t(apply(estimation.ds.mat*estimation.Fj.a0.a1, 1, cumsum))
  estimation.df$rmtlj.sep.direct.two.learners.cate.hat <- estimation.rmtlj.sep.direct.diff[, estimation.ind]+(taus[4]-estimation.s[estimation.ind])*(estimation.rmtlj.sep.direct.diff[, estimation.ind+1]-estimation.rmtlj.sep.direct.diff[, estimation.ind])/(estimation.s[estimation.ind+1]-estimation.s[estimation.ind])
  
  # separable indirect a_j=0
  train.df$rmtlj.sep.indirect.two.learners.cate.hat <- train.df$rmtlj.two.learners.cate.hat-train.df$rmtlj.sep.direct.two.learners.cate.hat
  estimation.df$rmtlj.sep.indirect.two.learners.cate.hat <- estimation.df$rmtlj.two.learners.cate.hat-estimation.df$rmtlj.sep.direct.two.learners.cate.hat
  
  # # causal survival forest
  # train.df$csf.rmtlj.cate.hat <- predict(causal_survival_forest(X=data.frame(estimation.df[, confounders]), Y=estimation.df$end.time, W=estimation.df$a, D=1*(estimation.df$event!=0), target="rmtlj", horizon=taus[4]), newdata=train.df[, confounders])$predictions
  # estimation.df$csf.rmtlj.cate.hat <- predict(causal_survival_forest(X=data.frame(train.df[, confounders]), Y=train.df$end.time, W=train.df$a, D=1*(train.df$event!=0), target="rmtlj", horizon=taus[4]), newdata=estimation.df[, confounders])$predictions
  # 
  # sim.df <- rbind(train.df, estimation.df)
  # sim.df$csf.rmtlj.cate.hat <- Winsorize(minval = -taus[4], maxval = taus[4], sim.df$csf.rmtlj.cate.hat)
  
  ## transform outcome
  ## total
  # eif
  train.df$ueif.diff <- ueif_dr_rmtlj_total_mt_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, tau=taus[4], bw=train.df$iptw, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sj.a1=train.Sj.a1, freq.time=NULL, admin.cens)
  estimation.df$ueif.diff <- ueif_dr_rmtlj_total_mt_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, tau=taus[4], bw=estimation.df$iptw, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sj.a1=estimation.Sj.a1, freq.time=NULL, admin.cens)
  
  # ipcw
  train.df$ipcw.rmtlj <- cut_ipcw_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, tau=taus[4], G.a0=train.G.a0, G.a1=train.G.a1, freq.time=NULL, admin.cens)
  estimation.df$ipcw.rmtlj <- cut_ipcw_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, tau=taus[4], G.a0=estimation.G.a0, G.a1=estimation.G.a1, freq.time=NULL, admin.cens)
  
  # bj
  train.df$bj.rmtlj <- cut_bj_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, tau=taus[4], S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sj.a1=train.Sj.a1, freq.time=NULL, admin.cens)
  estimation.df$bj.rmtlj <- cut_bj_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, tau=taus[4], S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sj.a1=estimation.Sj.a1, freq.time=NULL, admin.cens)
  
  # drcut
  train.df$drcut.rmtlj <- cut_dr_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, tau=taus[4], G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sj.a1=train.Sj.a1, freq.time=NULL, admin.cens)
  estimation.df$drcut.rmtlj <- cut_dr_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, tau=taus[4], G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sj.a1=estimation.Sj.a1, freq.time=NULL, admin.cens)
  
  ## separable direct a_jbar=1
  # eif
  train.df$ueif.sep.direct.diff <- ueif_sep_direct_astar1_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, tau=taus[4], bw=train.df$iptw,
                                                                G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)
  estimation.df$ueif.sep.direct.diff <- ueif_sep_direct_astar1_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, tau=taus[4], bw=estimation.df$iptw,
                                                                     G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)
  # bj
  train.df$bj.rmtlj.sep.direct <- cut_dr_sep_direct_astar1_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, tau=taus[4], G.a0=matrix(1, nrow(train.G.a0), ncol(train.G.a0)), G.a1=matrix(1, nrow(train.G.a1), ncol(train.G.a1)), S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)
  estimation.df$bj.rmtlj.sep.direct <- cut_dr_sep_direct_astar1_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, tau=taus[4], G.a0=matrix(1, nrow(estimation.G.a0), ncol(estimation.G.a0)), G.a1=matrix(1, nrow(estimation.G.a1), ncol(estimation.G.a1)), S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)
  
  # drcut
  train.df$drcut.rmtlj.sep.direct <- cut_dr_sep_direct_astar1_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, tau=taus[4], G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)
  estimation.df$drcut.rmtlj.sep.direct <- cut_dr_sep_direct_astar1_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, tau=taus[4], G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)
  
  ## separable indirect a_j=0
  # eif
  train.df$ueif.sep.indirect.diff <- ueif_sep_indirect_astar1_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, tau=taus[4], bw=train.df$iptw,
                                                                    G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)
  estimation.df$ueif.sep.indirect.diff <- ueif_sep_indirect_astar1_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, tau=taus[4], bw=estimation.df$iptw,
                                                                         G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)
  # bj
  train.df$bj.rmtlj.sep.indirect <- cut_dr_sep_indirect_astar1_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, tau=taus[4], G.a0=matrix(1, nrow(train.G.a0), ncol(train.G.a0)), G.a1=matrix(1, nrow(train.G.a1), ncol(train.G.a1)), S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)
  estimation.df$bj.rmtlj.sep.indirect <- cut_dr_sep_indirect_astar1_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, tau=taus[4], G.a0=matrix(1, nrow(estimation.G.a0), ncol(estimation.G.a0)), G.a1=matrix(1, nrow(estimation.G.a1), ncol(estimation.G.a1)), S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)
  
  # drcut
  train.df$drcut.rmtlj.sep.indirect <- cut_dr_sep_indirect_astar1_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, tau=taus[4], G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)
  estimation.df$drcut.rmtlj.sep.indirect <- cut_dr_sep_indirect_astar1_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, tau=taus[4], G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)
  
  rm(train.G.a0, train.G.a1, train.S.a0, train.S.a1, train.Sj.a0, train.Sjbar.a0, train.Sj.a1, train.Sjbar.a1, estimation.G.a0, estimation.G.a1, estimation.S.a0, estimation.S.a1, estimation.Sj.a0, estimation.Sjbar.a0, estimation.Sj.a1, estimation.Sjbar.a1); gc()
  
  ## second split
  sim.df <- rbind(train.df, estimation.df)
  train.ind <- sample(1:n, n/2, replace=FALSE) # sample training set
  train.df <- sim.df[train.ind,]
  estimation.df <- sim.df[setdiff(1:n, train.ind),]
  
  train.df$second.split <- "train"
  estimation.df$second.split <- "estimation"
  
  # ps
  train.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = estimation.df$a, X = estimation.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred) # train.df$ps <- predict(glm(a~x1+x2+x3+x4,data=train.df,family=binomial()),type="response")
  train.df$ow <- with(train.df, ifelse(a, 1-ps, ps))
  train.df$ow.tilt.hat <- with(train.df, ps*(1-ps))
  
  estimation.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = train.df$a, X = train.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred) # estimation.df$ps <- predict(glm(a~x1+x2+x3+x4,data=estimation.df,family=binomial()),type="response")
  estimation.df$ow <- with(estimation.df, ifelse(a, 1-ps, ps))
  estimation.df$ow.tilt.hat <- with(estimation.df, ps*(1-ps))
  
  ## dml
  # total effect
  train.df$rmtlj.dml.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ueif.diff, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$rmtlj.dml.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ueif.diff, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable direct
  train.df$rmtlj.sep.direct.dml.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ueif.sep.direct.diff, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                              cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$rmtlj.sep.direct.dml.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ueif.sep.direct.diff, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable indirect
  train.df$rmtlj.sep.indirect.dml.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ueif.sep.indirect.diff, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$rmtlj.sep.indirect.dml.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ueif.sep.indirect.diff, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  ## ipcw
  # single learner
  train.df$ipcw.rmtlj.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmtlj, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, train.df[, confounders], 1*train.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmtlj, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, train.df[, confounders], 0*train.df[, confounders]))$pred
  estimation.df$ipcw.rmtlj.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmtlj, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, estimation.df[, confounders], 1*estimation.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmtlj, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, estimation.df[, confounders], 0*estimation.df[, confounders]))$pred
  # two learners
  train.df$ipcw.rmtlj.mu.a1 <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmtlj[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                 cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$ipcw.rmtlj.mu.a0 <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmtlj[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                 cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$ipcw.rmtlj.two.learners.cate.hat <- train.df$ipcw.rmtlj.mu.a1-train.df$ipcw.rmtlj.mu.a0
  
  estimation.df$ipcw.rmtlj.mu.a1 <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmtlj[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$ipcw.rmtlj.mu.a0 <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmtlj[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$ipcw.rmtlj.two.learners.cate.hat <- estimation.df$ipcw.rmtlj.mu.a1-estimation.df$ipcw.rmtlj.mu.a0
  
  ### construct outcome
  ## ipcw
  # iptw
  train.df$ipcw.rmtlj.iptw.y <- with(train.df, ipcw.rmtlj*(a/ps-(1-a)/(1-ps)))
  estimation.df$ipcw.rmtlj.iptw.y <- with(estimation.df, ipcw.rmtlj*(a/ps-(1-a)/(1-ps)))
  
  # ra
  train.df$ipcw.rmtlj.ra.y <- with(train.df, a*(ipcw.rmtlj-ipcw.rmtlj.mu.a0)+(1-a)*(ipcw.rmtlj.mu.a1-ipcw.rmtlj))
  estimation.df$ipcw.rmtlj.ra.y <- with(estimation.df, a*(ipcw.rmtlj-ipcw.rmtlj.mu.a0)+(1-a)*(ipcw.rmtlj.mu.a1-ipcw.rmtlj))
  
  # aiptw
  train.df$ipcw.rmtlj.aiptw.y <- with(train.df, ipcw.rmtlj*(a/ps-(1-a)/(1-ps))+(1-a/ps)*ipcw.rmtlj.mu.a1-(1-(1-a)/(1-ps))*ipcw.rmtlj.mu.a0)
  estimation.df$ipcw.rmtlj.aiptw.y <- with(estimation.df, ipcw.rmtlj*(a/ps-(1-a)/(1-ps))+(1-a/ps)*ipcw.rmtlj.mu.a1-(1-(1-a)/(1-ps))*ipcw.rmtlj.mu.a0)
  
  # ulearner
  train.df$ipcw.rmtlj.ulearner.y <- with(train.df, (ipcw.rmtlj-(ipcw.rmtlj.mu.a0*(1-ps)+ipcw.rmtlj.mu.a1*ps))/(a-ps))
  estimation.df$ipcw.rmtlj.ulearner.y <- with(estimation.df, (ipcw.rmtlj-(ipcw.rmtlj.mu.a0*(1-ps)+ipcw.rmtlj.mu.a1*ps))/(a-ps))
  
  ### bj
  ## single learner
  # total
  train.df$bj.rmtlj.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                 cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, train.df[, confounders], 1*train.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, train.df[, confounders], 0*train.df[, confounders]))$pred
  estimation.df$bj.rmtlj.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, estimation.df[, confounders], 1*estimation.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, estimation.df[, confounders], 0*estimation.df[, confounders]))$pred
  # separable direct
  train.df$bj.rmtlj.sep.direct.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.direct, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                            cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, train.df[, confounders], 1*train.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.direct, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, train.df[, confounders], 0*train.df[, confounders]))$pred
  estimation.df$bj.rmtlj.sep.direct.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.direct, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                                 cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, estimation.df[, confounders], 1*estimation.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.direct, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, estimation.df[, confounders], 0*estimation.df[, confounders]))$pred
  # separable indirect
  train.df$bj.rmtlj.sep.indirect.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.indirect, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                              cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, train.df[, confounders], 1*train.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.indirect, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, train.df[, confounders], 0*train.df[, confounders]))$pred
  estimation.df$bj.rmtlj.sep.indirect.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.indirect, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, estimation.df[, confounders], 1*estimation.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.indirect, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, estimation.df[, confounders], 0*estimation.df[, confounders]))$pred
  ## two learners
  # total
  train.df$bj.rmtlj.mu.a1 <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$bj.rmtlj.mu.a0 <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$bj.rmtlj.two.learners.cate.hat <- train.df$bj.rmtlj.mu.a1-train.df$bj.rmtlj.mu.a0
  
  estimation.df$bj.rmtlj.mu.a1 <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                    cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$bj.rmtlj.mu.a0 <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                    cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$bj.rmtlj.two.learners.cate.hat <- estimation.df$bj.rmtlj.mu.a1-estimation.df$bj.rmtlj.mu.a0
  
  # separable direct
  train.df$bj.rmtlj.sep.direct.mu.a1 <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.direct[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                          cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$bj.rmtlj.sep.direct.mu.a0 <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.direct[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                          cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$bj.rmtlj.sep.direct.two.learners.cate.hat <- train.df$bj.rmtlj.sep.direct.mu.a1-train.df$bj.rmtlj.sep.direct.mu.a0
  
  estimation.df$bj.rmtlj.sep.direct.mu.a1 <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.direct[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.direct.mu.a0 <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.direct[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.direct.two.learners.cate.hat <- estimation.df$bj.rmtlj.sep.direct.mu.a1-estimation.df$bj.rmtlj.sep.direct.mu.a0
  
  # separable indirect
  train.df$bj.rmtlj.sep.indirect.mu.a1 <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.indirect[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                            cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$bj.rmtlj.sep.indirect.mu.a0 <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.indirect[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                            cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$bj.rmtlj.sep.indirect.two.learners.cate.hat <- train.df$bj.rmtlj.sep.indirect.mu.a1-train.df$bj.rmtlj.sep.indirect.mu.a0
  
  estimation.df$bj.rmtlj.sep.indirect.mu.a1 <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.indirect[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                 cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.indirect.mu.a0 <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.indirect[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                 cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.indirect.two.learners.cate.hat <- estimation.df$bj.rmtlj.sep.indirect.mu.a1-estimation.df$bj.rmtlj.sep.indirect.mu.a0
  
  ## iptw
  # total
  train.df$bj.rmtlj.iptw.y <- with(train.df, bj.rmtlj*(a/ps-(1-a)/(1-ps)))
  estimation.df$bj.rmtlj.iptw.y <- with(estimation.df, bj.rmtlj*(a/ps-(1-a)/(1-ps)))
  
  # separable direct
  train.df$bj.rmtlj.sep.direct.iptw.y <- with(train.df, bj.rmtlj.sep.direct*(a/ps-(1-a)/(1-ps)))
  estimation.df$bj.rmtlj.sep.direct.iptw.y <- with(estimation.df, bj.rmtlj.sep.direct*(a/ps-(1-a)/(1-ps)))
  
  # separable indirect
  train.df$bj.rmtlj.sep.indirect.iptw.y <- with(train.df, bj.rmtlj.sep.indirect*(a/ps-(1-a)/(1-ps)))
  estimation.df$bj.rmtlj.sep.indirect.iptw.y <- with(estimation.df, bj.rmtlj.sep.indirect*(a/ps-(1-a)/(1-ps)))
  
  ## ra
  # total
  train.df$bj.rmtlj.ra.y <- with(train.df, a*(bj.rmtlj-bj.rmtlj.mu.a0)+(1-a)*(bj.rmtlj.mu.a1-bj.rmtlj))
  estimation.df$bj.rmtlj.ra.y <- with(estimation.df, a*(bj.rmtlj-bj.rmtlj.mu.a0)+(1-a)*(bj.rmtlj.mu.a1-bj.rmtlj))
  
  # separable direct
  train.df$bj.rmtlj.sep.direct.ra.y <- with(train.df, a*(bj.rmtlj.sep.direct-bj.rmtlj.sep.direct.mu.a0)+(1-a)*(bj.rmtlj.sep.direct.mu.a1-bj.rmtlj.sep.direct))
  estimation.df$bj.rmtlj.sep.direct.ra.y <- with(estimation.df, a*(bj.rmtlj.sep.direct-bj.rmtlj.sep.direct.mu.a0)+(1-a)*(bj.rmtlj.sep.direct.mu.a1-bj.rmtlj.sep.direct))
  
  # separable indirect
  train.df$bj.rmtlj.sep.indirect.ra.y <- with(train.df, a*(bj.rmtlj.sep.indirect-bj.rmtlj.sep.indirect.mu.a0)+(1-a)*(bj.rmtlj.sep.indirect.mu.a1-bj.rmtlj.sep.indirect))
  estimation.df$bj.rmtlj.sep.indirect.ra.y <- with(estimation.df, a*(bj.rmtlj.sep.indirect-bj.rmtlj.sep.indirect.mu.a0)+(1-a)*(bj.rmtlj.sep.indirect.mu.a1-bj.rmtlj.sep.indirect))
  
  ## aiptw
  # total
  train.df$bj.rmtlj.aiptw.y <- with(train.df, bj.rmtlj*(a/ps-(1-a)/(1-ps))+(1-a/ps)*bj.rmtlj.mu.a1-(1-(1-a)/(1-ps))*bj.rmtlj.mu.a0)
  estimation.df$bj.rmtlj.aiptw.y <- with(estimation.df, bj.rmtlj*(a/ps-(1-a)/(1-ps))+(1-a/ps)*bj.rmtlj.mu.a1-(1-(1-a)/(1-ps))*bj.rmtlj.mu.a0)
  
  # separable direct
  train.df$bj.rmtlj.sep.direct.aiptw.y <- with(train.df, bj.rmtlj.sep.direct*(a/ps-(1-a)/(1-ps))+(1-a/ps)*bj.rmtlj.sep.direct.mu.a1-(1-(1-a)/(1-ps))*bj.rmtlj.sep.direct.mu.a0)
  estimation.df$bj.rmtlj.sep.direct.aiptw.y <- with(estimation.df, bj.rmtlj.sep.direct*(a/ps-(1-a)/(1-ps))+(1-a/ps)*bj.rmtlj.sep.direct.mu.a1-(1-(1-a)/(1-ps))*bj.rmtlj.sep.direct.mu.a0)
  
  # separable indirect
  train.df$bj.rmtlj.sep.indirect.aiptw.y <- with(train.df, bj.rmtlj.sep.indirect*(a/ps-(1-a)/(1-ps))+(1-a/ps)*bj.rmtlj.sep.indirect.mu.a1-(1-(1-a)/(1-ps))*bj.rmtlj.sep.indirect.mu.a0)
  estimation.df$bj.rmtlj.sep.indirect.aiptw.y <- with(estimation.df, bj.rmtlj.sep.indirect*(a/ps-(1-a)/(1-ps))+(1-a/ps)*bj.rmtlj.sep.indirect.mu.a1-(1-(1-a)/(1-ps))*bj.rmtlj.sep.indirect.mu.a0)
  
  ## ulearner
  # total
  train.df$bj.rmtlj.ulearner.y <- with(train.df, (bj.rmtlj-(bj.rmtlj.mu.a0*(1-ps)+bj.rmtlj.mu.a1*ps))/(a-ps))
  estimation.df$bj.rmtlj.ulearner.y <- with(estimation.df, (bj.rmtlj-(bj.rmtlj.mu.a0*(1-ps)+bj.rmtlj.mu.a1*ps))/(a-ps))
  
  # separable direct
  train.df$bj.rmtlj.sep.direct.ulearner.y <- with(train.df, (bj.rmtlj.sep.direct-(bj.rmtlj.sep.direct.mu.a0*(1-ps)+bj.rmtlj.sep.direct.mu.a1*ps))/(a-ps))
  estimation.df$bj.rmtlj.sep.direct.ulearner.y <- with(estimation.df, (bj.rmtlj.sep.direct-(bj.rmtlj.sep.direct.mu.a0*(1-ps)+bj.rmtlj.sep.direct.mu.a1*ps))/(a-ps))
  
  # separable indirect
  train.df$bj.rmtlj.sep.indirect.ulearner.y <- with(train.df, (bj.rmtlj.sep.indirect-(bj.rmtlj.sep.indirect.mu.a0*(1-ps)+bj.rmtlj.sep.indirect.mu.a1*ps))/(a-ps))
  estimation.df$bj.rmtlj.sep.indirect.ulearner.y <- with(estimation.df, (bj.rmtlj.sep.indirect-(bj.rmtlj.sep.indirect.mu.a0*(1-ps)+bj.rmtlj.sep.indirect.mu.a1*ps))/(a-ps))
  
  ### drcut
  ## single learner
  # total
  train.df$drcut.rmtlj.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                    cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, train.df[, confounders], 1*train.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, train.df[, confounders], 0*train.df[, confounders]))$pred
  estimation.df$drcut.rmtlj.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                         cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, estimation.df[, confounders], 1*estimation.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, estimation.df[, confounders], 0*estimation.df[, confounders]))$pred
  # separable direct
  train.df$drcut.rmtlj.sep.direct.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.direct, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, train.df[, confounders], 1*train.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.direct, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, train.df[, confounders], 0*train.df[, confounders]))$pred
  estimation.df$drcut.rmtlj.sep.direct.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.direct, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                                    cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, estimation.df[, confounders], 1*estimation.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.direct, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, estimation.df[, confounders], 0*estimation.df[, confounders]))$pred
  # separable indirect
  train.df$drcut.rmtlj.sep.indirect.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.indirect, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                                 cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, train.df[, confounders], 1*train.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.indirect, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, train.df[, confounders], 0*train.df[, confounders]))$pred
  estimation.df$drcut.rmtlj.sep.indirect.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.indirect, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, estimation.df[, confounders], 1*estimation.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.indirect, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, estimation.df[, confounders], 0*estimation.df[, confounders]))$pred
  ## two learners
  # total
  train.df$drcut.rmtlj.mu.a1 <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$drcut.rmtlj.mu.a0 <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$drcut.rmtlj.two.learners.cate.hat <- train.df$drcut.rmtlj.mu.a1-train.df$drcut.rmtlj.mu.a0
  
  estimation.df$drcut.rmtlj.mu.a1 <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$drcut.rmtlj.mu.a0 <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$drcut.rmtlj.two.learners.cate.hat <- estimation.df$drcut.rmtlj.mu.a1-estimation.df$drcut.rmtlj.mu.a0
  
  # separable direct
  train.df$drcut.rmtlj.sep.direct.mu.a1 <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.direct[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$drcut.rmtlj.sep.direct.mu.a0 <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.direct[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$drcut.rmtlj.sep.direct.two.learners.cate.hat <- train.df$drcut.rmtlj.sep.direct.mu.a1-train.df$drcut.rmtlj.sep.direct.mu.a0
  
  estimation.df$drcut.rmtlj.sep.direct.mu.a1 <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.direct[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.direct.mu.a0 <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.direct[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.direct.two.learners.cate.hat <- estimation.df$drcut.rmtlj.sep.direct.mu.a1-estimation.df$drcut.rmtlj.sep.direct.mu.a0
  
  # separable indirect
  train.df$drcut.rmtlj.sep.indirect.mu.a1 <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.indirect[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$drcut.rmtlj.sep.indirect.mu.a0 <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.indirect[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$drcut.rmtlj.sep.indirect.two.learners.cate.hat <- train.df$drcut.rmtlj.sep.indirect.mu.a1-train.df$drcut.rmtlj.sep.indirect.mu.a0
  
  estimation.df$drcut.rmtlj.sep.indirect.mu.a1 <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.indirect[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                    cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.indirect.mu.a0 <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.indirect[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                    cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.indirect.two.learners.cate.hat <- estimation.df$drcut.rmtlj.sep.indirect.mu.a1-estimation.df$drcut.rmtlj.sep.indirect.mu.a0
  
  ## causal forest
  # ipcw
  train.df$ipcw.rmtlj.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(estimation.df[, confounders]), Y=estimation.df$ipcw.rmtlj, W=estimation.df$a), newdata=train.df[, confounders])$predictions
  estimation.df$ipcw.rmtlj.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(train.df[, confounders]), Y=train.df$ipcw.rmtlj, W=train.df$a), newdata=estimation.df[, confounders])$predictions
  
  ## bj
  # total
  train.df$bj.rmtlj.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(estimation.df[, confounders]), Y=estimation.df$bj.rmtlj, W=estimation.df$a), newdata=train.df[, confounders])$predictions
  estimation.df$bj.rmtlj.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(train.df[, confounders]), Y=train.df$bj.rmtlj, W=train.df$a), newdata=estimation.df[, confounders])$predictions
  
  # separable direct
  train.df$bj.rmtlj.sep.direct.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(estimation.df[, confounders]), Y=estimation.df$bj.rmtlj.sep.direct, W=estimation.df$a), newdata=train.df[, confounders])$predictions
  estimation.df$bj.rmtlj.sep.direct.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(train.df[, confounders]), Y=train.df$bj.rmtlj.sep.direct, W=train.df$a), newdata=estimation.df[, confounders])$predictions
  
  # separable indirect
  train.df$bj.rmtlj.sep.indirect.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(estimation.df[, confounders]), Y=estimation.df$bj.rmtlj.sep.indirect, W=estimation.df$a), newdata=train.df[, confounders])$predictions
  estimation.df$bj.rmtlj.sep.indirect.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(train.df[, confounders]), Y=train.df$bj.rmtlj.sep.indirect, W=train.df$a), newdata=estimation.df[, confounders])$predictions
  
  ## drcut
  # total
  train.df$drcut.rmtlj.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(estimation.df[, confounders]), Y=estimation.df$drcut.rmtlj, W=estimation.df$a), newdata=train.df[, confounders])$predictions
  estimation.df$drcut.rmtlj.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(train.df[, confounders]), Y=train.df$drcut.rmtlj, W=train.df$a), newdata=estimation.df[, confounders])$predictions
  
  # separable direct
  train.df$drcut.rmtlj.sep.direct.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(estimation.df[, confounders]), Y=estimation.df$drcut.rmtlj.sep.direct, W=estimation.df$a), newdata=train.df[, confounders])$predictions
  estimation.df$drcut.rmtlj.sep.direct.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(train.df[, confounders]), Y=train.df$drcut.rmtlj.sep.direct, W=train.df$a), newdata=estimation.df[, confounders])$predictions
  
  # separable indirect
  train.df$drcut.rmtlj.sep.indirect.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(estimation.df[, confounders]), Y=estimation.df$drcut.rmtlj.sep.indirect, W=estimation.df$a), newdata=train.df[, confounders])$predictions
  estimation.df$drcut.rmtlj.sep.indirect.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(train.df[, confounders]), Y=train.df$drcut.rmtlj.sep.indirect, W=train.df$a), newdata=estimation.df[, confounders])$predictions
  
  ## iptw
  # total
  train.df$drcut.rmtlj.iptw.y <- with(train.df, drcut.rmtlj*(a/ps-(1-a)/(1-ps)))
  estimation.df$drcut.rmtlj.iptw.y <- with(estimation.df, drcut.rmtlj*(a/ps-(1-a)/(1-ps)))
  
  # separable direct
  train.df$drcut.rmtlj.sep.direct.iptw.y <- with(train.df, drcut.rmtlj.sep.direct*(a/ps-(1-a)/(1-ps)))
  estimation.df$drcut.rmtlj.sep.direct.iptw.y <- with(estimation.df, drcut.rmtlj.sep.direct*(a/ps-(1-a)/(1-ps)))
  
  # separable indirect
  train.df$drcut.rmtlj.sep.indirect.iptw.y <- with(train.df, drcut.rmtlj.sep.indirect*(a/ps-(1-a)/(1-ps)))
  estimation.df$drcut.rmtlj.sep.indirect.iptw.y <- with(estimation.df, drcut.rmtlj.sep.indirect*(a/ps-(1-a)/(1-ps)))
  
  ## ra
  # total
  train.df$drcut.rmtlj.ra.y <- with(train.df, a*(drcut.rmtlj-drcut.rmtlj.mu.a0)+(1-a)*(drcut.rmtlj.mu.a1-drcut.rmtlj))
  estimation.df$drcut.rmtlj.ra.y <- with(estimation.df, a*(drcut.rmtlj-drcut.rmtlj.mu.a0)+(1-a)*(drcut.rmtlj.mu.a1-drcut.rmtlj))
  
  # separable direct
  train.df$drcut.rmtlj.sep.direct.ra.y <- with(train.df, a*(drcut.rmtlj.sep.direct-drcut.rmtlj.sep.direct.mu.a0)+(1-a)*(drcut.rmtlj.sep.direct.mu.a1-drcut.rmtlj.sep.direct))
  estimation.df$drcut.rmtlj.sep.direct.ra.y <- with(estimation.df, a*(drcut.rmtlj.sep.direct-drcut.rmtlj.sep.direct.mu.a0)+(1-a)*(drcut.rmtlj.sep.direct.mu.a1-drcut.rmtlj.sep.direct))
  
  # separable indirect
  train.df$drcut.rmtlj.sep.indirect.ra.y <- with(train.df, a*(drcut.rmtlj.sep.indirect-drcut.rmtlj.sep.indirect.mu.a0)+(1-a)*(drcut.rmtlj.sep.indirect.mu.a1-drcut.rmtlj.sep.indirect))
  estimation.df$drcut.rmtlj.sep.indirect.ra.y <- with(estimation.df, a*(drcut.rmtlj.sep.indirect-drcut.rmtlj.sep.indirect.mu.a0)+(1-a)*(drcut.rmtlj.sep.indirect.mu.a1-drcut.rmtlj.sep.indirect))
  
  ## aiptw
  # total
  train.df$drcut.rmtlj.aiptw.y <- with(train.df, drcut.rmtlj*(a/ps-(1-a)/(1-ps))+(1-a/ps)*drcut.rmtlj.mu.a1-(1-(1-a)/(1-ps))*drcut.rmtlj.mu.a0)
  estimation.df$drcut.rmtlj.aiptw.y <- with(estimation.df, drcut.rmtlj*(a/ps-(1-a)/(1-ps))+(1-a/ps)*drcut.rmtlj.mu.a1-(1-(1-a)/(1-ps))*drcut.rmtlj.mu.a0)
  
  # separable direct
  train.df$drcut.rmtlj.sep.direct.aiptw.y <- with(train.df, drcut.rmtlj.sep.direct*(a/ps-(1-a)/(1-ps))+(1-a/ps)*drcut.rmtlj.sep.direct.mu.a1-(1-(1-a)/(1-ps))*drcut.rmtlj.sep.direct.mu.a0)
  estimation.df$drcut.rmtlj.sep.direct.aiptw.y <- with(estimation.df, drcut.rmtlj.sep.direct*(a/ps-(1-a)/(1-ps))+(1-a/ps)*drcut.rmtlj.sep.direct.mu.a1-(1-(1-a)/(1-ps))*drcut.rmtlj.sep.direct.mu.a0)
  
  # separable indirect
  train.df$drcut.rmtlj.sep.indirect.aiptw.y <- with(train.df, drcut.rmtlj.sep.indirect*(a/ps-(1-a)/(1-ps))+(1-a/ps)*drcut.rmtlj.sep.indirect.mu.a1-(1-(1-a)/(1-ps))*drcut.rmtlj.sep.indirect.mu.a0)
  estimation.df$drcut.rmtlj.sep.indirect.aiptw.y <- with(estimation.df, drcut.rmtlj.sep.indirect*(a/ps-(1-a)/(1-ps))+(1-a/ps)*drcut.rmtlj.sep.indirect.mu.a1-(1-(1-a)/(1-ps))*drcut.rmtlj.sep.indirect.mu.a0)
  
  ## ulearner
  # total
  train.df$drcut.rmtlj.ulearner.y <- with(train.df, (drcut.rmtlj-(drcut.rmtlj.mu.a0*(1-ps)+drcut.rmtlj.mu.a1*ps))/(a-ps))
  estimation.df$drcut.rmtlj.ulearner.y <- with(estimation.df, (drcut.rmtlj-(drcut.rmtlj.mu.a0*(1-ps)+drcut.rmtlj.mu.a1*ps))/(a-ps))
  
  # separable direct
  train.df$drcut.rmtlj.sep.direct.ulearner.y <- with(train.df, (drcut.rmtlj.sep.direct-(drcut.rmtlj.sep.direct.mu.a0*(1-ps)+drcut.rmtlj.sep.direct.mu.a1*ps))/(a-ps))
  estimation.df$drcut.rmtlj.sep.direct.ulearner.y <- with(estimation.df, (drcut.rmtlj.sep.direct-(drcut.rmtlj.sep.direct.mu.a0*(1-ps)+drcut.rmtlj.sep.direct.mu.a1*ps))/(a-ps))
  
  # separable indirect
  train.df$drcut.rmtlj.sep.indirect.ulearner.y <- with(train.df, (drcut.rmtlj.sep.indirect-(drcut.rmtlj.sep.indirect.mu.a0*(1-ps)+drcut.rmtlj.sep.indirect.mu.a1*ps))/(a-ps))
  estimation.df$drcut.rmtlj.sep.indirect.ulearner.y <- with(estimation.df, (drcut.rmtlj.sep.indirect-(drcut.rmtlj.sep.indirect.mu.a0*(1-ps)+drcut.rmtlj.sep.indirect.mu.a1*ps))/(a-ps))
  
  ## third split
  sim.df <- rbind(train.df, estimation.df)
  train.ind <- sample(1:n, n/2, replace=FALSE) # sample training set
  train.df <- sim.df[train.ind,]
  estimation.df <- sim.df[setdiff(1:n, train.ind),]
  
  ## modified outcome
  # iptw
  train.df$ipcw.rmtlj.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmtlj.iptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                         cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmtlj.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmtlj.iptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                              cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # ra
  train.df$ipcw.rmtlj.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmtlj.ra.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmtlj.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmtlj.ra.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                            cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # aiptw
  train.df$ipcw.rmtlj.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmtlj.aiptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                          cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmtlj.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmtlj.aiptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # xlearner
  train.df$ipcw.rmtlj.res.a1 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$ipcw.rmtlj-estimation.df$ipcw.rmtlj.mu.a0)[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$ipcw.rmtlj.res.a0 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$ipcw.rmtlj.mu.a1-estimation.df$ipcw.rmtlj)[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmtlj.res.a1 <- predict.SuperLearner(SuperLearner(Y = (train.df$ipcw.rmtlj-train.df$ipcw.rmtlj.mu.a0)[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$ipcw.rmtlj.res.a0 <- predict.SuperLearner(SuperLearner(Y = (train.df$ipcw.rmtlj.mu.a1-train.df$ipcw.rmtlj)[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  train.df$ipcw.rmtlj.xlearner.cate.hat <- with(train.df, ps*ipcw.rmtlj.res.a0+(1-ps)*ipcw.rmtlj.res.a1)
  estimation.df$ipcw.rmtlj.xlearner.cate.hat <- with(estimation.df, ps*ipcw.rmtlj.res.a0+(1-ps)*ipcw.rmtlj.res.a1)
  
  ## modified covariate
  train.df$ipcw.rmtlj.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*ipcw.rmtlj), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmtlj.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*ipcw.rmtlj), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                            cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # efficiency augmentation
  train.df$ipcw.rmtlj.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*(ipcw.rmtlj-(ipcw.rmtlj.mu.a0*(1-ps)+ipcw.rmtlj.mu.a1*ps))), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                         cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmtlj.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*(ipcw.rmtlj-(ipcw.rmtlj.mu.a0*(1-ps)+ipcw.rmtlj.mu.a1*ps))), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                              cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # u-learner
  train.df$ipcw.rmtlj.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmtlj.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmtlj.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmtlj.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # r-learner
  train.df$ipcw.rmtlj.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmtlj.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((estimation.df$a-estimation.df$ps)^2)[,1],
                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmtlj.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmtlj.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((train.df$a-train.df$ps)^2)[,1],
                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  
  ### modified outcome
  ## iptw
  # total
  train.df$bj.rmtlj.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.iptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.iptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                            cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable direct
  train.df$bj.rmtlj.sep.direct.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.direct.iptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.direct.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.direct.iptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable indirect
  train.df$bj.rmtlj.sep.indirect.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.indirect.iptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                    cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.indirect.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.indirect.iptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                         cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  ## ra
  # total
  train.df$bj.rmtlj.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.ra.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.ra.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                          cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable direct
  train.df$bj.rmtlj.sep.direct.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.direct.ra.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.direct.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.direct.ra.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable indirect
  train.df$bj.rmtlj.sep.indirect.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.indirect.ra.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.indirect.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.indirect.ra.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  ## aiptw
  # total
  train.df$bj.rmtlj.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.aiptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.aiptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable direct
  train.df$bj.rmtlj.sep.direct.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.direct.aiptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.direct.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.direct.aiptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable indirect
  train.df$bj.rmtlj.sep.indirect.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.indirect.aiptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.indirect.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.indirect.aiptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                          cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  ## xlearner
  # total
  train.df$bj.rmtlj.res.a1 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$bj.rmtlj-estimation.df$bj.rmtlj.mu.a0)[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$bj.rmtlj.res.a0 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$bj.rmtlj.mu.a1-estimation.df$bj.rmtlj)[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.res.a1 <- predict.SuperLearner(SuperLearner(Y = (train.df$bj.rmtlj-train.df$bj.rmtlj.mu.a0)[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$bj.rmtlj.res.a0 <- predict.SuperLearner(SuperLearner(Y = (train.df$bj.rmtlj.mu.a1-train.df$bj.rmtlj)[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  train.df$bj.rmtlj.xlearner.cate.hat <- with(train.df, ps*bj.rmtlj.res.a0+(1-ps)*bj.rmtlj.res.a1)
  estimation.df$bj.rmtlj.xlearner.cate.hat <- with(estimation.df, ps*bj.rmtlj.res.a0+(1-ps)*bj.rmtlj.res.a1)
  
  # separable direct
  train.df$bj.rmtlj.sep.indirect.res.a1 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$bj.rmtlj.sep.indirect-estimation.df$bj.rmtlj.sep.indirect.mu.a0)[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$bj.rmtlj.sep.indirect.res.a0 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$bj.rmtlj.sep.indirect.mu.a1-estimation.df$bj.rmtlj.sep.indirect)[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.indirect.res.a1 <- predict.SuperLearner(SuperLearner(Y = (train.df$bj.rmtlj.sep.indirect-train.df$bj.rmtlj.sep.indirect.mu.a0)[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.indirect.res.a0 <- predict.SuperLearner(SuperLearner(Y = (train.df$bj.rmtlj.sep.indirect.mu.a1-train.df$bj.rmtlj.sep.indirect)[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  train.df$bj.rmtlj.sep.indirect.xlearner.cate.hat <- with(train.df, ps*bj.rmtlj.sep.indirect.res.a0+(1-ps)*bj.rmtlj.sep.indirect.res.a1)
  estimation.df$bj.rmtlj.sep.indirect.xlearner.cate.hat <- with(estimation.df, ps*bj.rmtlj.sep.indirect.res.a0+(1-ps)*bj.rmtlj.sep.indirect.res.a1)
  
  # separable indirect
  train.df$bj.rmtlj.sep.direct.res.a1 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$bj.rmtlj.sep.direct-estimation.df$bj.rmtlj.sep.direct.mu.a0)[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$bj.rmtlj.sep.direct.res.a0 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$bj.rmtlj.sep.direct.mu.a1-estimation.df$bj.rmtlj.sep.direct)[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.direct.res.a1 <- predict.SuperLearner(SuperLearner(Y = (train.df$bj.rmtlj.sep.direct-train.df$bj.rmtlj.sep.direct.mu.a0)[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.direct.res.a0 <- predict.SuperLearner(SuperLearner(Y = (train.df$bj.rmtlj.sep.direct.mu.a1-train.df$bj.rmtlj.sep.direct)[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  train.df$bj.rmtlj.sep.direct.xlearner.cate.hat <- with(train.df, ps*bj.rmtlj.sep.direct.res.a0+(1-ps)*bj.rmtlj.sep.direct.res.a1)
  estimation.df$bj.rmtlj.sep.direct.xlearner.cate.hat <- with(estimation.df, ps*bj.rmtlj.sep.direct.res.a0+(1-ps)*bj.rmtlj.sep.direct.res.a1)
  
  ## modified covariate
  # total
  train.df$bj.rmtlj.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*bj.rmtlj), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*bj.rmtlj), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                          cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable direct
  train.df$bj.rmtlj.sep.direct.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*bj.rmtlj.sep.direct), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.direct.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*bj.rmtlj.sep.direct), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable indirect
  train.df$bj.rmtlj.sep.indirect.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*bj.rmtlj.sep.indirect), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.indirect.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*bj.rmtlj.sep.indirect), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  ## efficiency augmentation
  # total
  train.df$bj.rmtlj.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*(bj.rmtlj-(bj.rmtlj.mu.a0*(1-ps)+bj.rmtlj.mu.a1*ps))), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*(bj.rmtlj-(bj.rmtlj.mu.a0*(1-ps)+bj.rmtlj.mu.a1*ps))), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                            cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable direct
  train.df$bj.rmtlj.sep.direct.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*(bj.rmtlj.sep.direct-(bj.rmtlj.sep.direct.mu.a0*(1-ps)+bj.rmtlj.sep.direct.mu.a1*ps))), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.direct.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*(bj.rmtlj.sep.direct-(bj.rmtlj.sep.direct.mu.a0*(1-ps)+bj.rmtlj.sep.direct.mu.a1*ps))), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable indirect
  train.df$bj.rmtlj.sep.indirect.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*(bj.rmtlj.sep.indirect-(bj.rmtlj.sep.indirect.mu.a0*(1-ps)+bj.rmtlj.sep.indirect.mu.a1*ps))), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                    cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.indirect.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*(bj.rmtlj.sep.indirect-(bj.rmtlj.sep.indirect.mu.a0*(1-ps)+bj.rmtlj.sep.indirect.mu.a1*ps))), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                         cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # u-learner
  # total
  train.df$bj.rmtlj.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable direct
  train.df$bj.rmtlj.sep.direct.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.direct.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.direct.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.direct.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable indirect
  train.df$bj.rmtlj.sep.indirect.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.indirect.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.indirect.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.indirect.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  ## r-learner
  # total
  train.df$bj.rmtlj.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((estimation.df$a-estimation.df$ps)^2)[,1],
                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((train.df$a-train.df$ps)^2)[,1],
                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable direct
  train.df$bj.rmtlj.sep.direct.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.direct.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((estimation.df$a-estimation.df$ps)^2)[,1],
                                                                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.direct.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.direct.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((train.df$a-train.df$ps)^2)[,1],
                                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable indirect
  train.df$bj.rmtlj.sep.indirect.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmtlj.sep.indirect.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((estimation.df$a-estimation.df$ps)^2)[,1],
                                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmtlj.sep.indirect.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmtlj.sep.indirect.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((train.df$a-train.df$ps)^2)[,1],
                                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  ### drcut
  ### modified outcome
  ## iptw
  # total
  train.df$drcut.rmtlj.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.iptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                          cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.iptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable direct
  train.df$drcut.rmtlj.sep.direct.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.direct.iptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.direct.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.direct.iptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                          cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable indirect
  train.df$drcut.rmtlj.sep.indirect.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.indirect.iptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.indirect.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.indirect.iptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                            cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  ## ra
  # total
  train.df$drcut.rmtlj.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.ra.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.ra.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable direct
  train.df$drcut.rmtlj.sep.direct.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.direct.ra.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.direct.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.direct.ra.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable indirect
  train.df$drcut.rmtlj.sep.indirect.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.indirect.ra.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.indirect.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.indirect.ra.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                          cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  ## aiptw
  # total
  train.df$drcut.rmtlj.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.aiptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.aiptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable direct
  train.df$drcut.rmtlj.sep.direct.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.direct.aiptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.direct.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.direct.aiptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable indirect
  train.df$drcut.rmtlj.sep.indirect.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.indirect.aiptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.indirect.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.indirect.aiptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  ## xlearner
  # total
  train.df$drcut.rmtlj.res.a1 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$drcut.rmtlj-estimation.df$drcut.rmtlj.mu.a0)[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$drcut.rmtlj.res.a0 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$drcut.rmtlj.mu.a1-estimation.df$drcut.rmtlj)[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.res.a1 <- predict.SuperLearner(SuperLearner(Y = (train.df$drcut.rmtlj-train.df$drcut.rmtlj.mu.a0)[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$drcut.rmtlj.res.a0 <- predict.SuperLearner(SuperLearner(Y = (train.df$drcut.rmtlj.mu.a1-train.df$drcut.rmtlj)[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  train.df$drcut.rmtlj.xlearner.cate.hat <- with(train.df, ps*drcut.rmtlj.res.a0+(1-ps)*drcut.rmtlj.res.a1)
  estimation.df$drcut.rmtlj.xlearner.cate.hat <- with(estimation.df, ps*drcut.rmtlj.res.a0+(1-ps)*drcut.rmtlj.res.a1)
  
  # separable direct
  train.df$drcut.rmtlj.sep.indirect.res.a1 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$drcut.rmtlj.sep.indirect-estimation.df$drcut.rmtlj.sep.indirect.mu.a0)[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$drcut.rmtlj.sep.indirect.res.a0 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$drcut.rmtlj.sep.indirect.mu.a1-estimation.df$drcut.rmtlj.sep.indirect)[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.indirect.res.a1 <- predict.SuperLearner(SuperLearner(Y = (train.df$drcut.rmtlj.sep.indirect-train.df$drcut.rmtlj.sep.indirect.mu.a0)[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.indirect.res.a0 <- predict.SuperLearner(SuperLearner(Y = (train.df$drcut.rmtlj.sep.indirect.mu.a1-train.df$drcut.rmtlj.sep.indirect)[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  train.df$drcut.rmtlj.sep.indirect.xlearner.cate.hat <- with(train.df, ps*drcut.rmtlj.sep.indirect.res.a0+(1-ps)*drcut.rmtlj.sep.indirect.res.a1)
  estimation.df$drcut.rmtlj.sep.indirect.xlearner.cate.hat <- with(estimation.df, ps*drcut.rmtlj.sep.indirect.res.a0+(1-ps)*drcut.rmtlj.sep.indirect.res.a1)
  
  # separable indirect
  train.df$drcut.rmtlj.sep.direct.res.a1 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$drcut.rmtlj.sep.direct-estimation.df$drcut.rmtlj.sep.direct.mu.a0)[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                              cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$drcut.rmtlj.sep.direct.res.a0 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$drcut.rmtlj.sep.direct.mu.a1-estimation.df$drcut.rmtlj.sep.direct)[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                              cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.direct.res.a1 <- predict.SuperLearner(SuperLearner(Y = (train.df$drcut.rmtlj.sep.direct-train.df$drcut.rmtlj.sep.direct.mu.a0)[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.direct.res.a0 <- predict.SuperLearner(SuperLearner(Y = (train.df$drcut.rmtlj.sep.direct.mu.a1-train.df$drcut.rmtlj.sep.direct)[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  train.df$drcut.rmtlj.sep.direct.xlearner.cate.hat <- with(train.df, ps*drcut.rmtlj.sep.direct.res.a0+(1-ps)*drcut.rmtlj.sep.direct.res.a1)
  estimation.df$drcut.rmtlj.sep.direct.xlearner.cate.hat <- with(estimation.df, ps*drcut.rmtlj.sep.direct.res.a0+(1-ps)*drcut.rmtlj.sep.direct.res.a1)
  
  ## modified covariate
  # total
  train.df$drcut.rmtlj.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*drcut.rmtlj), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*drcut.rmtlj), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable direct
  train.df$drcut.rmtlj.sep.direct.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*drcut.rmtlj.sep.direct), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.direct.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*drcut.rmtlj.sep.direct), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable indirect
  train.df$drcut.rmtlj.sep.indirect.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*drcut.rmtlj.sep.indirect), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.indirect.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*drcut.rmtlj.sep.indirect), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                          cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  ## efficiency augmentation
  # total
  train.df$drcut.rmtlj.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*(drcut.rmtlj-(drcut.rmtlj.mu.a0*(1-ps)+drcut.rmtlj.mu.a1*ps))), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                          cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*(drcut.rmtlj-(drcut.rmtlj.mu.a0*(1-ps)+drcut.rmtlj.mu.a1*ps))), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable direct
  train.df$drcut.rmtlj.sep.direct.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*(drcut.rmtlj.sep.direct-(drcut.rmtlj.sep.direct.mu.a0*(1-ps)+drcut.rmtlj.sep.direct.mu.a1*ps))), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.direct.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*(drcut.rmtlj.sep.direct-(drcut.rmtlj.sep.direct.mu.a0*(1-ps)+drcut.rmtlj.sep.direct.mu.a1*ps))), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                          cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable indirect
  train.df$drcut.rmtlj.sep.indirect.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*(drcut.rmtlj.sep.indirect-(drcut.rmtlj.sep.indirect.mu.a0*(1-ps)+drcut.rmtlj.sep.indirect.mu.a1*ps))), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.indirect.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*(drcut.rmtlj.sep.indirect-(drcut.rmtlj.sep.indirect.mu.a0*(1-ps)+drcut.rmtlj.sep.indirect.mu.a1*ps))), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                                            cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # u-learner
  # total
  train.df$drcut.rmtlj.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                              cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable direct
  train.df$drcut.rmtlj.sep.direct.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.direct.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                         cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.direct.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.direct.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                              cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable indirect
  train.df$drcut.rmtlj.sep.indirect.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.indirect.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.indirect.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.indirect.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  ## r-learner
  # total
  train.df$drcut.rmtlj.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((estimation.df$a-estimation.df$ps)^2)[,1],
                                                                              cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((train.df$a-train.df$ps)^2)[,1],
                                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable direct
  train.df$drcut.rmtlj.sep.direct.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.direct.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((estimation.df$a-estimation.df$ps)^2)[,1],
                                                                                         cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.direct.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.direct.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((train.df$a-train.df$ps)^2)[,1],
                                                                                              cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # separable indirect
  train.df$drcut.rmtlj.sep.indirect.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmtlj.sep.indirect.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((estimation.df$a-estimation.df$ps)^2)[,1],
                                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmtlj.sep.indirect.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmtlj.sep.indirect.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((train.df$a-train.df$ps)^2)[,1],
                                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  sim.df <- rbind(train.df, estimation.df)
  
  # ## bcf
  # # ipcw
  # sim.df$ipcw.rmtlj.bcf.cate.hat <- XBCF::getTaus(XBCF::XBCF(y=as.matrix(sim.df$ipcw.rmtlj), z=as.matrix(sim.df$a), x_con=as.matrix(sim.df[, confounders]), pihat=sim.df$ps, pcat_con=3, pcat_mod=3))
  # 
  # ## bj
  # # total
  # sim.df$bj.rmtlj.bcf.cate.hat <- XBCF::getTaus(XBCF::XBCF(y=as.matrix(sim.df$bj.rmtlj), z=as.matrix(sim.df$a), x_con=as.matrix(sim.df[, confounders]), pihat=sim.df$ps, pcat_con=3, pcat_mod=3))
  # # separable direct
  # sim.df$bj.rmtlj.sep.direct.bcf.cate.hat <- XBCF::getTaus(XBCF::XBCF(y=as.matrix(sim.df$bj.rmtlj.sep.direct), z=as.matrix(sim.df$a), x_con=as.matrix(sim.df[, confounders]), pihat=sim.df$ps, pcat_con=3, pcat_mod=3))
  # # separable indirect
  # sim.df$bj.rmtlj.sep.indirect.bcf.cate.hat <- XBCF::getTaus(XBCF::XBCF(y=as.matrix(sim.df$bj.rmtlj.sep.indirect), z=as.matrix(sim.df$a), x_con=as.matrix(sim.df[, confounders]), pihat=sim.df$ps, pcat_con=3, pcat_mod=3))
  # 
  # ## drcut
  # # total
  # sim.df$drcut.rmtlj.bcf.cate.hat <- XBCF::getTaus(XBCF::XBCF(y=as.matrix(sim.df$drcut.rmtlj), z=as.matrix(sim.df$a), x_con=as.matrix(sim.df[, confounders]), pihat=sim.df$ps, pcat_con=3, pcat_mod=3))
  # # separable direct
  # sim.df$drcut.rmtlj.sep.direct.bcf.cate.hat <- XBCF::getTaus(XBCF::XBCF(y=as.matrix(sim.df$drcut.rmtlj.sep.direct), z=as.matrix(sim.df$a), x_con=as.matrix(sim.df[, confounders]), pihat=sim.df$ps, pcat_con=3, pcat_mod=3))
  # # separable indirect
  # sim.df$drcut.rmtlj.sep.indirect.bcf.cate.hat <- XBCF::getTaus(XBCF::XBCF(y=as.matrix(sim.df$drcut.rmtlj.sep.indirect), z=as.matrix(sim.df$a), x_con=as.matrix(sim.df[, confounders]), pihat=sim.df$ps, pcat_con=3, pcat_mod=3))
  
  sim.df[, grep("cate.hat", names(sim.df), value=TRUE)] <- Winsorize(minval = -taus[4], maxval = taus[4], sim.df[, grep("cate.hat", names(sim.df), value=TRUE)])
  
  print(i)
  write.csv(sim.df, file=paste0("sim1-competing", i,".csv"), row.names = FALSE)
  NULL
}
end.time <- Sys.time()
end.time - start.time
