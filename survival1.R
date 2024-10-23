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
nsim <- 125 # 1000
time.grid <- 0.01
admin.cens <- 10
taus <- c(1,2,3,4)
confounders <- c("x1","x2","x3","x4","x5","x6")
cens.sl.lib <- event.sl.lib <- c("survSL.km","survSL.coxph","survSL.loglogreg","survSL.rfsrc","survSL.expreg","survSL.weibreg")
a.sl.lib <- c("SL.mean","SL.glm","SL.nnet","SL.kernelKnn","SL.rpartPrune","SL.xgboost","SL.ranger","SL.step","SL.gam","SL.glmnet","SL.earth")
weights.sl.lib <- c("SL.glm","SL.nnet","SL.earth","SL.gam","SL.glmnet","SL.ranger","SL.xgboost")

registerDoParallel(cores=125)

## continuous
# IPCW transformation
cut_ipcw_rmst <- function(id, a, time, event, tau, G.a0, G.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  causes <- sort(unique(event[event!=0]))
  ncauses <- length(causes)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  ind <- max(which(s<tau))
  
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  
  Delta.tau.over.G.a0 <- do.call(cbind, lapply(1:ns, function(u) ifelse(time > s[u], 1, 0)))/G.a0+matrix(event/rowSums(G.a0*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u) ifelse(time <= s[u], 1, 0)))
  Delta.tau.over.G.a1 <- do.call(cbind, lapply(1:ns, function(u) ifelse(time > s[u], 1, 0)))/G.a1+matrix(event/rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u) ifelse(time <= s[u], 1, 0)))
  # Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  # dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  
  term1.a0 <- Delta.tau.over.G.a0*pmin(matrix(time,ncol=ns,nrow=n,byrow=FALSE), matrix(s,ncol=ns,nrow=n,byrow=TRUE))
  # term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a0
  # term3.a0 <- t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNct - Yt*G.dHazard.a0)/G.a0, 1, cumsum))
  # term4.a0 <- RMST.a0*t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  # term5.a0 <- t(apply(RMST.a0*(dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  cut.a0 <- term1.a0 # +term3.a0+term4.a0-term5.a0
  cut.a0 <- cut.a0[, ind]+(tau-s[ind])*(cut.a0[, ind+1]-cut.a0[, ind])/(s[ind+1]-s[ind])
  
  term1.a1 <- Delta.tau.over.G.a1*pmin(matrix(time,ncol=ns,nrow=n,byrow=FALSE), matrix(s,ncol=ns,nrow=n,byrow=TRUE))
  # term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a1
  # term3.a1 <- t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNct - Yt*G.dHazard.a1)/G.a1, 1, cumsum))
  # term4.a1 <- RMST.a1*t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  # term5.a1 <- t(apply(RMST.a1*(dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  cut.a1 <- term1.a1 # +term3.a1+term4.a1-term5.a1
  cut.a1 <- cut.a1[, ind]+(tau-s[ind])*(cut.a1[, ind+1]-cut.a1[, ind])/(s[ind+1]-s[ind])
  
  return(cut.a0*(1-a)+cut.a1*a)
}

# BJ transformation
cut_bj_rmst <- function(id, a, time, event, tau, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  causes <- sort(unique(event[event!=0]))
  ncauses <- length(causes)
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
  
  RMST.a0 <- t(apply(S.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMST.a1 <- t(apply(S.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMST.a0.tau <- RMST.a0[, ind]+(tau-s[ind])*(RMST.a0[, ind+1]-RMST.a0[, ind])/(s[ind+1]-s[ind])
  RMST.a1.tau <- RMST.a1[, ind]+(tau-s[ind])*(RMST.a1[, ind+1]-RMST.a1[, ind])/(s[ind+1]-s[ind])
  RMST.a0.Ttilde <- rowSums(RMST.a0*do.call(cbind, lapply(1:ns, function(u)  ifelse(time == s[u], 1, 0))))
  RMST.a1.Ttilde <- rowSums(RMST.a1*do.call(cbind, lapply(1:ns, function(u)  ifelse(time == s[u], 1, 0))))
  RMST.a0.tau.min.Ttilde <- (time>tau)*RMST.a0.tau+(time<=tau)*RMST.a0.Ttilde
  RMST.a1.tau.min.Ttilde <- (time>tau)*RMST.a1.tau+(time<=tau)*RMST.a1.Ttilde
  
  cut.a0 <- pmin(time, tau)+(1-event*(time<=tau)-(time>tau))*(RMST.a0.tau-RMST.a0.tau.min.Ttilde)/S.a0.tau.min.Ttilde
  cut.a1 <- pmin(time, tau)+(1-event*(time<=tau)-(time>tau))*(RMST.a1.tau-RMST.a1.tau.min.Ttilde)/S.a1.tau.min.Ttilde
  
  return(cut.a0*(1-a)+cut.a1*a)
}

# doubly robust censoring unbiased transformation
cut_dr_rmst <- function(id, a, time, event, tau, G.a0, G.a1, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  causes <- sort(unique(event[event!=0]))
  ncauses <- length(causes)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  ind <- max(which(s<tau))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  RMST.a0 <- t(apply(S.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMST.a1 <- t(apply(S.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  
  Delta.tau.over.G.a0 <- do.call(cbind, lapply(1:ns, function(u) ifelse(time > s[u], 1, 0)))/G.a0+matrix(event/rowSums(G.a0*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u) ifelse(time <= s[u], 1, 0)))
  Delta.tau.over.G.a1 <- do.call(cbind, lapply(1:ns, function(u) ifelse(time > s[u], 1, 0)))/G.a1+matrix(event/rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u) ifelse(time <= s[u], 1, 0)))
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  
  term1.a0 <- Delta.tau.over.G.a0*pmin(matrix(time,ncol=ns,nrow=n,byrow=FALSE), matrix(s,ncol=ns,nrow=n,byrow=TRUE))
  # term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a0
  term3.a0 <- t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNct - Yt*G.dHazard.a0)/G.a0, 1, cumsum))
  term4.a0 <- RMST.a0*t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  term5.a0 <- t(apply(RMST.a0*(dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  cut.a0 <- term1.a0+term3.a0+term4.a0-term5.a0
  cut.a0 <- cut.a0[, ind]+(tau-s[ind])*(cut.a0[, ind+1]-cut.a0[, ind])/(s[ind+1]-s[ind])
  
  term1.a1 <- Delta.tau.over.G.a1*pmin(matrix(time,ncol=ns,nrow=n,byrow=FALSE), matrix(s,ncol=ns,nrow=n,byrow=TRUE))
  # term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a1
  term3.a1 <- t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNct - Yt*G.dHazard.a1)/G.a1, 1, cumsum))
  term4.a1 <- RMST.a1*t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  term5.a1 <- t(apply(RMST.a1*(dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  cut.a1 <- term1.a1+term3.a1+term4.a1-term5.a1
  cut.a1 <- cut.a1[, ind]+(tau-s[ind])*(cut.a1[, ind+1]-cut.a1[, ind])/(s[ind+1]-s[ind])
  
  return(cut.a0*(1-a)+cut.a1*a)
}

#
ueif_dr_rmst_direct_mc_ha <- function(id, a, time, event, tau, bw, tilt, G.a0, G.a1, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  causes <- sort(unique(event[event!=0]))
  ncauses <- length(causes)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  ind <- max(which(s<tau))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  RMST.a0 <- t(apply(S.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMST.a1 <- t(apply(S.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  
  Delta.tau.over.G.a0 <- do.call(cbind, lapply(1:ns, function(u) ifelse(time > s[u], 1, 0)))/G.a0+matrix(event/rowSums(G.a0*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u) ifelse(time <= s[u], 1, 0)))
  Delta.tau.over.G.a1 <- do.call(cbind, lapply(1:ns, function(u) ifelse(time > s[u], 1, 0)))/G.a1+matrix(event/rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u) ifelse(time <= s[u], 1, 0)))
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*Delta.tau.over.G.a0*pmin(matrix(time,ncol=ns,nrow=n,byrow=FALSE), matrix(s,ncol=ns,nrow=n,byrow=TRUE))
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a0
  term3.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNct - Yt*G.dHazard.a0)/G.a0, 1, cumsum))
  term4.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*RMST.a0*(t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))-1)
  term5.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMST.a0*(dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  ueif.a0 <- 1/mean(tilt)*term2.a0+1/mean(bw*(a==0))*(term1.a0+term3.a0+term4.a0-term5.a0)
  ueif.a0 <- ueif.a0[, ind]+(tau-s[ind])*(ueif.a0[, ind+1]-ueif.a0[, ind])/(s[ind+1]-s[ind])
  
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*Delta.tau.over.G.a1*pmin(matrix(time,ncol=ns,nrow=n,byrow=FALSE), matrix(s,ncol=ns,nrow=n,byrow=TRUE))
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNct - Yt*G.dHazard.a1)/G.a1, 1, cumsum))
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*RMST.a1*(t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))-1)
  term5.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMST.a1*(dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  ueif.a1 <- 1/mean(tilt)*term2.a1+1/mean(bw*(a==1))*(term1.a1+term3.a1+term4.a1-term5.a1)
  ueif.a1 <- ueif.a1[, ind]+(tau-s[ind])*(ueif.a1[, ind+1]-ueif.a1[, ind])/(s[ind+1]-s[ind])
  
  return(ueif.a1-ueif.a0) # data.frame(ueif.a0, ueif.a1)
}

lognormal <- function(data.df, betas, sigma, n)
{return(exp(sigma*qnorm(runif(n))+as.numeric((as.matrix(data.df) %*% as.matrix(betas)))))}

loglogistic <- function(data.df, betas, gamma, n)
{return(exp(gamma*(log(runif(n))-log(1-runif(n)))+as.numeric((as.matrix(data.df) %*% as.matrix(betas)))))}

# aft <- function(data.df, betas, sigma, n)
# {return(exp(sigma*rnorm(n, mean=0, sd=1)+as.numeric((as.matrix(data.df) %*% as.matrix(betas)))))}

exponenetial <- function(data.df, lambdas, betas)
{return(-log(runif(n=nrow(data.df), min=0, max=1))/(lambdas*exp(as.matrix(data.df) %*% as.matrix(betas)))[, 1])}

## proportional hazards
dgp4_survival_nph <- function(n, admin.cens, time.grid, tau){
  corr.mat <- matrix(0.5, nrow=3, ncol=3)
  diag(corr.mat) <- 1
  x <- mvrnorm(n=n, mu=rep(0,3), Sigma=corr.mat)
  x1 <- x[,1]
  x2 <- x[,2]
  x3 <- x[,3]
  x4 <- rbinom(n,1,0.5)
  x5 <- rbinom(n,1,0.5)
  x6 <- rbinom(n,1,0.5)
  ps <- 1/(1+exp(-(0.3+0.2*x1+0.3*x2+0.3*x3-0.2*x4-0.3*x5-0.2*x6))) # 1/(1+exp(-(-1+x1+1.5*x2+1.5*x3-x4-1.5*x5-x6)))
  a <- rbinom(n, 1, ps) # -1.5+0.9*x1+1.2*x2+1.2*x3+1.2*x4
  data.df <- data.frame(a, x1, x2, x3, x4, x5, x6, ow.tilt.true=ps*(1-ps))
  
  data.df$Ta0 <- ceiling(loglogistic(data.df=data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas=c(x0=0.8, x1=-0.8, x2=1, x3=0.8, x4=0.4, x5=-0.4, x6=0.8), gamma=0.2, n=n)*100)/100
  data.df$Ta1 <- ceiling(lognormal(data.df=data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas=c(x0=0.4, x1=0.6, x2=-0.8, x3=1.2, x4=0.6, x5=-0.3, x6=0.5), sigma=1, n=n)*100)/100
  data.df$Ta <- (1-a)*data.df$Ta0 + a*data.df$Ta1
  
  data.df$C[data.df$a==0] <- lognormal(data.df=data.frame(x0=1, x1=x1[a==0], x2=x2[a==0], x3=x3[a==0], x4=x4[a==0], x5=x5[a==0], x6=x6[a==0]), betas=c(x0=1.8, x1=0.6, x2=-0.8, x3=0.5, x4=0.7, x5=-0.4, x6=-0.2), sigma=0.8, n=sum(a==0))
  data.df$C[data.df$a==1] <- loglogistic(data.df=data.frame(x0=1, x1=x1[a==1], x2=x2[a==1], x3=x3[a==1], x4=x4[a==1], x5=x5[a==1], x6=x6[a==1]), betas=c(x0=2.2, x1=0.6, x2=-0.8, x3=0.5, x4=0.7, x5=0.8, x6=1.2), gamma=0.8, n=sum(a==1))
  
  data.df$C <- pmin(data.df$C, admin.cens, na.rm=TRUE)
  data.df$end.time <- ceiling(pmin(data.df$C, data.df$Ta)/time.grid)*time.grid # , data.df$Tj2, ceil(pmin(data.df$C,data.df$Tj1,data.df$Tj2)*365.25) # round(pmin(data.df$C,data.df$Tj1,data.df$Tj2),digits=2)
  data.df$event <- with(data.df,ifelse(Ta<=C, 1, 0))
  data.df <- data.df[order(data.df$end.time,-data.df$event),]
  data.df$id <- 1:n
  
  data.rep.df <- do.call("rbind", replicate(1000, data.df, simplify = FALSE))
  data.rep.df$min.Ta0.tau <- pmin(loglogistic(data.df=data.frame(x0=1, data.rep.df[, c("x1","x2","x3","x4","x5","x6")]), betas=c(x0=0.8, x1=-0.8, x2=1, x3=0.8, x4=0.4, x5=-0.4, x6=0.8), gamma=0.2, n=1000*n), tau) # gen_survtimes_exp(lambdas=0.12, covariates=as.matrix(data.frame(x0=1, data.rep.df[, confounders])), betas=c(x0=0.1, x1=0.1, x2=-0.2, x3=0.2, x4=0.1, x5=0.8, x6=-0.2))
  data.rep.df$min.Ta1.tau <- pmin(lognormal(data.df=data.frame(x0=1, data.rep.df[, c("x1","x2","x3","x4","x5","x6")]), betas=c(x0=0.4, x1=0.6, x2=-0.8, x3=1.2, x4=0.6, x5=-0.3, x6=0.5), sigma=1, n=1000*n), tau) # gen_survtimes_exp(lambdas=0.15, covariates=as.matrix(data.frame(x0=1, data.rep.df[, confounders])), betas=c(x0=0.17, x1=0.2, x2=-0.1, x3=0.4, x4=0.2, x5=0.3, x6=0.4))
  
  data.df <- data.rep.df %>%
    group_by(id) %>%
    summarise_at(vars(min.Ta0.tau, min.Ta1.tau), list(name = mean)) %>%
    rename(rmst.a0.cate.true=min.Ta0.tau_name, rmst.a1.cate.true=min.Ta1.tau_name) %>%
    merge(data.df, by="id")
  data.df$rmst.diff.cate.true <- data.df$rmst.a1.cate.true-data.df$rmst.a0.cate.true
  data.df$rmst.a0.ite.true <- pmin(data.df$Ta0, tau)
  data.df$rmst.a1.ite.true <- pmin(data.df$Ta1, tau)
  data.df$rmst.diff.ite.true <- data.df$rmst.a1.ite.true-data.df$rmst.a0.ite.true
  return(data.df)
}

setwd("survival1")

start.time <- Sys.time()
foreach(i=1:nsim, .combine=c, .errorhandling="remove") %dopar% {
  sim.df <- dgp4_survival_nph(n, admin.cens, time.grid, tau=taus[4]) # 1-sum(sim.df$event)/n
  
  ## first split
  train.ind <- sample(1:n, n/2, replace=FALSE) # sample training set
  train.df <- sim.df[train.ind,]
  estimation.df <- sim.df[setdiff(1:n, train.ind),]
  
  while(max(train.df$end.time[train.df$a==0 & train.df$event==0])<taus[4] | max(train.df$end.time[train.df$a==1 & train.df$event==0])<taus[4] |
        max(train.df$end.time[train.df$a==0 & train.df$event==1])<taus[4] | max(train.df$end.time[train.df$a==1 & train.df$event==1])<taus[4] |
        max(estimation.df$end.time[estimation.df$a==0 & estimation.df$event==0])<taus[4] | max(estimation.df$end.time[estimation.df$a==1 & estimation.df$event==0])<taus[4] |
        max(estimation.df$end.time[estimation.df$a==0 & estimation.df$event==1])<taus[4] | max(estimation.df$end.time[estimation.df$a==1 & estimation.df$event==1])<taus[4]) {
    sim.df <- dgp4_survival_nph(2*n, admin.cens, time.grid)
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
  train.ds.mat <- matrix(diff(c(0, train.s)), nrow=n/2, ncol=n_distinct(train.df$end.time), byrow=TRUE)
  train.ind <- max(which(train.s<taus[4]))
  
  estimation.s <- sort(unique(estimation.df$end.time))
  estimation.ds.mat <- matrix(diff(c(0, estimation.s)), nrow=n/2, ncol=n_distinct(estimation.df$end.time), byrow=TRUE)
  estimation.ind <- max(which(estimation.s<taus[4]))
  
  ## comparsion
  # single learner
  train.S.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time, event = 1*(estimation.df$event!=0), X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                 newX = data.frame(a=0, train.df[, confounders], 0*train.df[, confounders]), new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  train.S.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time, event = 1*(estimation.df$event!=0), X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                 newX = data.frame(a=1, train.df[, confounders], 1*train.df[, confounders]), new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  
  estimation.S.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time, event = 1*(train.df$event!=0), X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                      newX = data.frame(a=0, estimation.df[, confounders], 0*estimation.df[, confounders]), new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  estimation.S.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time, event = 1*(train.df$event!=0), X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                      newX = data.frame(a=1, estimation.df[, confounders], 1*estimation.df[, confounders]), new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
  
  train.rmst.diff <- t(apply(train.ds.mat*train.S.a1, 1, cumsum))-t(apply(train.ds.mat*train.S.a0, 1, cumsum))
  train.df$rmst.single.learner.cate.hat <- train.rmst.diff[, train.ind]+(taus[4]-train.s[train.ind])*(train.rmst.diff[, train.ind+1]-train.rmst.diff[, train.ind])/(train.s[train.ind+1]-train.s[train.ind])
  
  estimation.rmst.diff <- t(apply(estimation.ds.mat*estimation.S.a1, 1, cumsum))-t(apply(estimation.ds.mat*estimation.S.a0, 1, cumsum))
  estimation.df$rmst.single.learner.cate.hat <- estimation.rmst.diff[, estimation.ind]+(taus[4]-estimation.s[estimation.ind])*(estimation.rmst.diff[, estimation.ind+1]-estimation.rmst.diff[, estimation.ind])/(estimation.s[estimation.ind+1]-estimation.s[estimation.ind])
  
  # two survSuperLearner
  ## train
  train.all.survival.super.a0 <- survSuperLearner(time = estimation.df$end.time[estimation.df$a==0], event = 1*(estimation.df$event[estimation.df$a==0]!=0), X = estimation.df[estimation.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                  newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
  train.G.a0 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a0$cens.SL.predict)
  train.S.a0 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a0$event.SL.predict)
  
  train.all.survival.super.a1 <- survSuperLearner(time = estimation.df$end.time[estimation.df$a==1], event = 1*(estimation.df$event[estimation.df$a==1]!=0), X = estimation.df[estimation.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                  newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
  train.G.a1 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a1$cens.SL.predict)
  train.S.a1 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a1$event.SL.predict)
  
  ## estimation
  estimation.all.survival.super.a0 <- survSuperLearner(time = train.df$end.time[train.df$a==0], event = 1*(train.df$event[train.df$a==0]!=0), X = train.df[train.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                       newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
  estimation.G.a0 <- Winsorize(minval = 1e-3, maxval = 1, estimation.all.survival.super.a0$cens.SL.predict)
  estimation.S.a0 <- Winsorize(minval = 1e-3, maxval = 1, estimation.all.survival.super.a0$event.SL.predict)
  
  estimation.all.survival.super.a1 <- survSuperLearner(time = train.df$end.time[train.df$a==1], event = 1*(train.df$event[train.df$a==1]!=0), X = train.df[train.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                       newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
  estimation.G.a1 <- Winsorize(minval = 1e-3, maxval = 1, estimation.all.survival.super.a1$cens.SL.predict)
  estimation.S.a1 <- Winsorize(minval = 1e-3, maxval = 1, estimation.all.survival.super.a1$event.SL.predict)
  
  train.rmst.diff <- t(apply(train.ds.mat*train.S.a1, 1, cumsum))-t(apply(train.ds.mat*train.S.a0, 1, cumsum))
  train.df$rmst.two.learners.cate.hat <- train.rmst.diff[, train.ind]+(taus[4]-train.s[train.ind])*(train.rmst.diff[, train.ind+1]-train.rmst.diff[, train.ind])/(train.s[train.ind+1]-train.s[train.ind])
  
  estimation.rmst.diff <- t(apply(estimation.ds.mat*estimation.S.a1, 1, cumsum))-t(apply(estimation.ds.mat*estimation.S.a0, 1, cumsum))
  estimation.df$rmst.two.learners.cate.hat <- estimation.rmst.diff[, estimation.ind]+(taus[4]-estimation.s[estimation.ind])*(estimation.rmst.diff[, estimation.ind+1]-estimation.rmst.diff[, estimation.ind])/(estimation.s[estimation.ind+1]-estimation.s[estimation.ind])
  
  # causal survival forest
  train.df$rmst.csf.cate.hat <- predict(causal_survival_forest(X=data.frame(estimation.df[, confounders]), Y=estimation.df$end.time, W=estimation.df$a, D=1*(estimation.df$event!=0), target="RMST", horizon=taus[4]), newdata=train.df[, confounders])$predictions
  estimation.df$rmst.csf.cate.hat <- predict(causal_survival_forest(X=data.frame(train.df[, confounders]), Y=train.df$end.time, W=train.df$a, D=1*(train.df$event!=0), target="RMST", horizon=taus[4]), newdata=estimation.df[, confounders])$predictions
  
  ## transform outcome
  # eif
  train.df$ueif.diff <- ueif_dr_rmst_direct_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, tau=taus[4], bw=train.df$iptw, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, freq.time=NULL, admin.cens)
  estimation.df$ueif.diff <- ueif_dr_rmst_direct_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, tau=taus[4], bw=estimation.df$iptw, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, freq.time=NULL, admin.cens)
  
  # ipcw
  train.df$ipcw.rmst <- cut_ipcw_rmst(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, tau=taus[4], G.a0=train.G.a0, G.a1=train.G.a1, freq.time=NULL, admin.cens)
  estimation.df$ipcw.rmst <- cut_ipcw_rmst(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, tau=taus[4], G.a0=estimation.G.a0, G.a1=estimation.G.a1, freq.time=NULL, admin.cens)
  
  # bj
  train.df$bj.rmst <- cut_bj_rmst(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, tau=taus[4], S.a0=train.S.a0, S.a1=train.S.a1, freq.time=NULL, admin.cens)
  estimation.df$bj.rmst <- cut_bj_rmst(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, tau=taus[4], S.a0=estimation.S.a0, S.a1=estimation.S.a1, freq.time=NULL, admin.cens)
  
  # drcut
  train.df$drcut.rmst <- cut_dr_rmst(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, tau=taus[4], G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, freq.time=NULL, admin.cens)
  estimation.df$drcut.rmst <- cut_dr_rmst(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, tau=taus[4], G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, freq.time=NULL, admin.cens)
  
  rm(train.G.a0, train.G.a1, train.S.a0, train.S.a1, estimation.G.a0, estimation.G.a1, estimation.S.a0, estimation.S.a1); gc()
  
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
  
  # dml
  train.df$rmst.dml.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ueif.diff, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$rmst.dml.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ueif.diff, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  ## ipcw
  # single learner
  train.df$ipcw.rmst.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmst, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, train.df[, confounders], 1*train.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmst, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, train.df[, confounders], 0*train.df[, confounders]))$pred
  estimation.df$ipcw.rmst.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmst, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, estimation.df[, confounders], 1*estimation.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmst, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, estimation.df[, confounders], 0*estimation.df[, confounders]))$pred
  # two learners
  train.df$ipcw.rmst.mu.a1 <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmst[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$ipcw.rmst.mu.a0 <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmst[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$ipcw.rmst.two.learners.cate.hat <- train.df$ipcw.rmst.mu.a1-train.df$ipcw.rmst.mu.a0
  
  estimation.df$ipcw.rmst.mu.a1 <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmst[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$ipcw.rmst.mu.a0 <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmst[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$ipcw.rmst.two.learners.cate.hat <- estimation.df$ipcw.rmst.mu.a1-estimation.df$ipcw.rmst.mu.a0
  
  ## bj
  # single learner
  train.df$bj.rmst.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmst, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, train.df[, confounders], 1*train.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmst, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, train.df[, confounders], 0*train.df[, confounders]))$pred
  estimation.df$bj.rmst.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmst, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                     cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, estimation.df[, confounders], 1*estimation.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = train.df$bj.rmst, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, estimation.df[, confounders], 0*estimation.df[, confounders]))$pred
  # two learners
  train.df$bj.rmst.mu.a1 <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmst[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                              cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$bj.rmst.mu.a0 <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmst[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                              cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$bj.rmst.two.learners.cate.hat <- train.df$bj.rmst.mu.a1-train.df$bj.rmst.mu.a0
  
  estimation.df$bj.rmst.mu.a1 <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmst[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$bj.rmst.mu.a0 <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmst[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$bj.rmst.two.learners.cate.hat <- estimation.df$bj.rmst.mu.a1-estimation.df$bj.rmst.mu.a0
  
  ## drcut
  # single learner
  train.df$drcut.rmst.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmst, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                   cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, train.df[, confounders], 1*train.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmst, X = data.frame(estimation.df[, c("a", confounders)], estimation.df$a*estimation.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, train.df[, confounders], 0*train.df[, confounders]))$pred
  estimation.df$drcut.rmst.single.learner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmst, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=1, estimation.df[, confounders], 1*estimation.df[, confounders]))$pred-
    predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmst, X = data.frame(train.df[, c("a", confounders)], train.df$a*train.df[, confounders]), family = gaussian(), SL.library = a.sl.lib,
                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=data.frame(a=0, estimation.df[, confounders], 0*estimation.df[, confounders]))$pred
  # two learners
  train.df$drcut.rmst.mu.a1 <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmst[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                 cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$drcut.rmst.mu.a0 <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmst[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                 cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$drcut.rmst.two.learners.cate.hat <- train.df$drcut.rmst.mu.a1-train.df$drcut.rmst.mu.a0
  
  estimation.df$drcut.rmst.mu.a1 <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmst[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$drcut.rmst.mu.a0 <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmst[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$drcut.rmst.two.learners.cate.hat <- estimation.df$drcut.rmst.mu.a1-estimation.df$drcut.rmst.mu.a0
  
  # causal forest
  train.df$ipcw.rmst.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(estimation.df[, confounders]), Y=estimation.df$ipcw.rmst, W=estimation.df$a), newdata=train.df[, confounders])$predictions
  estimation.df$ipcw.rmst.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(train.df[, confounders]), Y=train.df$ipcw.rmst, W=train.df$a), newdata=estimation.df[, confounders])$predictions
  
  train.df$bj.rmst.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(estimation.df[, confounders]), Y=estimation.df$bj.rmst, W=estimation.df$a), newdata=train.df[, confounders])$predictions
  estimation.df$bj.rmst.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(train.df[, confounders]), Y=train.df$bj.rmst, W=train.df$a), newdata=estimation.df[, confounders])$predictions
  
  train.df$drcut.rmst.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(estimation.df[, confounders]), Y=estimation.df$drcut.rmst, W=estimation.df$a), newdata=train.df[, confounders])$predictions
  estimation.df$drcut.rmst.causal.forest.cate.hat <- predict(causal_forest(X=data.frame(train.df[, confounders]), Y=train.df$drcut.rmst, W=train.df$a), newdata=estimation.df[, confounders])$predictions
  
  ## ipcw construct outcome
  # iptw
  train.df$ipcw.rmst.iptw.y <- with(train.df, ipcw.rmst*(a/ps-(1-a)/(1-ps)))
  estimation.df$ipcw.rmst.iptw.y <- with(estimation.df, ipcw.rmst*(a/ps-(1-a)/(1-ps)))
  
  # ra
  train.df$ipcw.rmst.ra.y <- with(train.df, a*(ipcw.rmst-ipcw.rmst.mu.a0)+(1-a)*(ipcw.rmst.mu.a1-ipcw.rmst))
  estimation.df$ipcw.rmst.ra.y <- with(estimation.df, a*(ipcw.rmst-ipcw.rmst.mu.a0)+(1-a)*(ipcw.rmst.mu.a1-ipcw.rmst))
  
  # aiptw
  train.df$ipcw.rmst.aiptw.y <- with(train.df, ipcw.rmst*(a/ps-(1-a)/(1-ps))+(1-a/ps)*ipcw.rmst.mu.a1-(1-(1-a)/(1-ps))*ipcw.rmst.mu.a0)
  estimation.df$ipcw.rmst.aiptw.y <- with(estimation.df, ipcw.rmst*(a/ps-(1-a)/(1-ps))+(1-a/ps)*ipcw.rmst.mu.a1-(1-(1-a)/(1-ps))*ipcw.rmst.mu.a0)
  
  # ulearner
  train.df$ipcw.rmst.ulearner.y <- with(train.df, (ipcw.rmst-(ipcw.rmst.mu.a0*(1-ps)+ipcw.rmst.mu.a1*ps))/(a-ps))
  estimation.df$ipcw.rmst.ulearner.y <- with(estimation.df, (ipcw.rmst-(ipcw.rmst.mu.a0*(1-ps)+ipcw.rmst.mu.a1*ps))/(a-ps))
  
  ## bj construct outcome
  # iptw
  train.df$bj.rmst.iptw.y <- with(train.df, bj.rmst*(a/ps-(1-a)/(1-ps)))
  estimation.df$bj.rmst.iptw.y <- with(estimation.df, bj.rmst*(a/ps-(1-a)/(1-ps)))
  
  # ra
  train.df$bj.rmst.ra.y <- with(train.df, a*(bj.rmst-bj.rmst.mu.a0)+(1-a)*(bj.rmst.mu.a1-bj.rmst))
  estimation.df$bj.rmst.ra.y <- with(estimation.df, a*(bj.rmst-bj.rmst.mu.a0)+(1-a)*(bj.rmst.mu.a1-bj.rmst))
  
  # aiptw
  train.df$bj.rmst.aiptw.y <- with(train.df, bj.rmst*(a/ps-(1-a)/(1-ps))+(1-a/ps)*bj.rmst.mu.a1-(1-(1-a)/(1-ps))*bj.rmst.mu.a0)
  estimation.df$bj.rmst.aiptw.y <- with(estimation.df, bj.rmst*(a/ps-(1-a)/(1-ps))+(1-a/ps)*bj.rmst.mu.a1-(1-(1-a)/(1-ps))*bj.rmst.mu.a0)
  
  # ulearner
  train.df$bj.rmst.ulearner.y <- with(train.df, (bj.rmst-(bj.rmst.mu.a0*(1-ps)+bj.rmst.mu.a1*ps))/(a-ps))
  estimation.df$bj.rmst.ulearner.y <- with(estimation.df, (bj.rmst-(bj.rmst.mu.a0*(1-ps)+bj.rmst.mu.a1*ps))/(a-ps))
  
  ## drcut construct outcome
  # iptw
  train.df$drcut.rmst.iptw.y <- with(train.df, drcut.rmst*(a/ps-(1-a)/(1-ps)))
  estimation.df$drcut.rmst.iptw.y <- with(estimation.df, drcut.rmst*(a/ps-(1-a)/(1-ps)))
  
  # ra
  train.df$drcut.rmst.ra.y <- with(train.df, a*(drcut.rmst-drcut.rmst.mu.a0)+(1-a)*(drcut.rmst.mu.a1-drcut.rmst))
  estimation.df$drcut.rmst.ra.y <- with(estimation.df, a*(drcut.rmst-drcut.rmst.mu.a0)+(1-a)*(drcut.rmst.mu.a1-drcut.rmst))
  
  # aiptw
  train.df$drcut.rmst.aiptw.y <- with(train.df, drcut.rmst*(a/ps-(1-a)/(1-ps))+(1-a/ps)*drcut.rmst.mu.a1-(1-(1-a)/(1-ps))*drcut.rmst.mu.a0)
  estimation.df$drcut.rmst.aiptw.y <- with(estimation.df, drcut.rmst*(a/ps-(1-a)/(1-ps))+(1-a/ps)*drcut.rmst.mu.a1-(1-(1-a)/(1-ps))*drcut.rmst.mu.a0)
  
  # ulearner
  train.df$drcut.rmst.ulearner.y <- with(train.df, (drcut.rmst-(drcut.rmst.mu.a0*(1-ps)+drcut.rmst.mu.a1*ps))/(a-ps))
  estimation.df$drcut.rmst.ulearner.y <- with(estimation.df, (drcut.rmst-(drcut.rmst.mu.a0*(1-ps)+drcut.rmst.mu.a1*ps))/(a-ps))
  
  ## third split
  sim.df <- rbind(train.df, estimation.df)
  train.ind <- sample(1:n, n/2, replace=FALSE) # sample training set
  train.df <- sim.df[train.ind,]
  estimation.df <- sim.df[setdiff(1:n, train.ind),]
  
  ## ipcw modified outcome
  # iptw
  train.df$ipcw.rmst.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmst.iptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmst.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmst.iptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # ra
  train.df$ipcw.rmst.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmst.ra.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmst.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmst.ra.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # aiptw
  train.df$ipcw.rmst.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmst.aiptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                         cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmst.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmst.aiptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                              cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # xlearner
  train.df$ipcw.rmst.res.a1 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$ipcw.rmst-estimation.df$ipcw.rmst.mu.a0)[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                 cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$ipcw.rmst.res.a0 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$ipcw.rmst.mu.a1-estimation.df$ipcw.rmst)[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                 cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmst.res.a1 <- predict.SuperLearner(SuperLearner(Y = (train.df$ipcw.rmst-train.df$ipcw.rmst.mu.a0)[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$ipcw.rmst.res.a0 <- predict.SuperLearner(SuperLearner(Y = (train.df$ipcw.rmst.mu.a1-train.df$ipcw.rmst)[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  train.df$ipcw.rmst.xlearner.cate.hat <- with(train.df, ps*ipcw.rmst.res.a0+(1-ps)*ipcw.rmst.res.a1)
  estimation.df$ipcw.rmst.xlearner.cate.hat <- with(estimation.df, ps*ipcw.rmst.res.a0+(1-ps)*ipcw.rmst.res.a1)
  
  ## modified covariate
  train.df$ipcw.rmst.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*ipcw.rmst), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmst.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*ipcw.rmst), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # efficiency augmentation
  train.df$ipcw.rmst.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*(ipcw.rmst-(ipcw.rmst.mu.a0*(1-ps)+ipcw.rmst.mu.a1*ps))), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                        cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmst.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*(ipcw.rmst-(ipcw.rmst.mu.a0*(1-ps)+ipcw.rmst.mu.a1*ps))), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # u-learner
  train.df$ipcw.rmst.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmst.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                            cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmst.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmst.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                 cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # r-learner
  train.df$ipcw.rmst.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$ipcw.rmst.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((estimation.df$a-estimation.df$ps)^2)[,1],
                                                                            cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$ipcw.rmst.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$ipcw.rmst.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((train.df$a-train.df$ps)^2)[,1],
                                                                                 cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # # bcf
  # sim.df$ipcw.rmst.bcf.cate.hat <- XBCF::getTaus(XBCF::XBCF(y=as.matrix(sim.df$ipcw.rmst), z=as.matrix(sim.df$a), x_con=as.matrix(sim.df[, confounders]), pihat=sim.df$ps, pcat_con=3, pcat_mod=3))
  # while(is.null(sim.df$ipcw.rmst.bcf.cate.hat))
  # {sim.df$ipcw.rmst.bcf.cate.hat <- XBCF::getTaus(XBCF::XBCF(y=as.matrix(sim.df$ipcw.rmst), z=as.matrix(sim.df$a), x_con=as.matrix(sim.df[, confounders]), pihat=sim.df$ps, pcat_con=3, pcat_mod=3))}
  # sim.df$ipcw.rmst.bcf.cate.hat <- Winsorize(minval = -taus[4], maxval = taus[4], sim.df$ipcw.rmst.bcf.cate.hat)
  
  ## bj modified outcome
  # iptw
  train.df$bj.rmst.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmst.iptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmst.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmst.iptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # ra
  train.df$bj.rmst.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmst.ra.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                    cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmst.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmst.ra.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                         cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # aiptw
  train.df$bj.rmst.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmst.aiptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmst.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmst.aiptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                            cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # xlearner
  train.df$bj.rmst.res.a1 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$bj.rmst-estimation.df$bj.rmst.mu.a0)[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$bj.rmst.res.a0 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$bj.rmst.mu.a1-estimation.df$bj.rmst)[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmst.res.a1 <- predict.SuperLearner(SuperLearner(Y = (train.df$bj.rmst-train.df$bj.rmst.mu.a0)[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                    cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$bj.rmst.res.a0 <- predict.SuperLearner(SuperLearner(Y = (train.df$bj.rmst.mu.a1-train.df$bj.rmst)[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                    cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  train.df$bj.rmst.xlearner.cate.hat <- with(train.df, ps*bj.rmst.res.a0+(1-ps)*bj.rmst.res.a1)
  estimation.df$bj.rmst.xlearner.cate.hat <- with(estimation.df, ps*bj.rmst.res.a0+(1-ps)*bj.rmst.res.a1)
  
  ## modified covariate
  train.df$bj.rmst.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*bj.rmst), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                    cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmst.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*bj.rmst), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                         cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # efficiency augmentation
  train.df$bj.rmst.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*(bj.rmst-(bj.rmst.mu.a0*(1-ps)+bj.rmst.mu.a1*ps))), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                      cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmst.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*(bj.rmst-(bj.rmst.mu.a0*(1-ps)+bj.rmst.mu.a1*ps))), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # u-learner
  train.df$bj.rmst.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmst.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                          cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmst.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmst.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # r-learner
  train.df$bj.rmst.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$bj.rmst.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((estimation.df$a-estimation.df$ps)^2)[,1],
                                                                          cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$bj.rmst.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$bj.rmst.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((train.df$a-train.df$ps)^2)[,1],
                                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # # bcf
  # sim.df$bj.rmst.bcf.cate.hat <- XBCF::getTaus(XBCF::XBCF(y=as.matrix(sim.df$bj.rmst), z=as.matrix(sim.df$a), x_con=as.matrix(sim.df[, confounders]), pihat=sim.df$ps, pcat_con=3, pcat_mod=3))
  # while(is.null(sim.df$bj.rmst.bcf.cate.hat))
  # {sim.df$bj.rmst.bcf.cate.hat <- XBCF::getTaus(XBCF::XBCF(y=as.matrix(sim.df$bj.rmst), z=as.matrix(sim.df$a), x_con=as.matrix(sim.df[, confounders]), pihat=sim.df$ps, pcat_con=3, pcat_mod=3))}
  # sim.df$bj.rmst.bcf.cate.hat <- Winsorize(minval = -taus[4], maxval = taus[4], sim.df$bj.rmst.bcf.cate.hat)
  
  ## drcut modified outcome
  # iptw
  train.df$drcut.rmst.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmst.iptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                         cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmst.iptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmst.iptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                              cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # ra
  train.df$drcut.rmst.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmst.ra.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmst.ra.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmst.ra.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                            cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # aiptw
  train.df$drcut.rmst.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmst.aiptw.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                          cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmst.aiptw.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmst.aiptw.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # xlearner
  train.df$drcut.rmst.res.a1 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$drcut.rmst-estimation.df$drcut.rmst.mu.a0)[estimation.df$a==1], X = estimation.df[estimation.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  train.df$drcut.rmst.res.a0 <- predict.SuperLearner(SuperLearner(Y = (estimation.df$drcut.rmst.mu.a1-estimation.df$drcut.rmst)[estimation.df$a==0], X = estimation.df[estimation.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmst.res.a1 <- predict.SuperLearner(SuperLearner(Y = (train.df$drcut.rmst-train.df$drcut.rmst.mu.a0)[train.df$a==1], X = train.df[train.df$a==1, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  estimation.df$drcut.rmst.res.a0 <- predict.SuperLearner(SuperLearner(Y = (train.df$drcut.rmst.mu.a1-train.df$drcut.rmst)[train.df$a==0], X = train.df[train.df$a==0, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  train.df$drcut.rmst.xlearner.cate.hat <- with(train.df, ps*drcut.rmst.res.a0+(1-ps)*drcut.rmst.res.a1)
  estimation.df$drcut.rmst.xlearner.cate.hat <- with(estimation.df, ps*drcut.rmst.res.a0+(1-ps)*drcut.rmst.res.a1)
  
  ## modified covariate
  train.df$drcut.rmst.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*drcut.rmst), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                       cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmst.mc.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*drcut.rmst), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                            cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # efficiency augmentation
  train.df$drcut.rmst.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(estimation.df, 2*(2*a-1)*(drcut.rmst-(drcut.rmst.mu.a0*(1-ps)+drcut.rmst.mu.a1*ps))), X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(estimation.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                         cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmst.mcea.cate.hat <- predict.SuperLearner(SuperLearner(Y = with(train.df, 2*(2*a-1)*(drcut.rmst-(drcut.rmst.mu.a0*(1-ps)+drcut.rmst.mu.a1*ps))), X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = with(train.df, (2*a-1)*(a-ps)/(4*ps*(1-ps)))[,1],
                                                                              cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # u-learner
  train.df$drcut.rmst.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmst.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmst.ulearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmst.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = a.sl.lib,
                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  # r-learner
  train.df$drcut.rmst.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = estimation.df$drcut.rmst.ulearner.y, X = estimation.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((estimation.df$a-estimation.df$ps)^2)[,1],
                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred
  estimation.df$drcut.rmst.rlearner.cate.hat <- predict.SuperLearner(SuperLearner(Y = train.df$drcut.rmst.ulearner.y, X = train.df[, confounders], family = gaussian(), SL.library = weights.sl.lib, obsWeights = ((train.df$a-train.df$ps)^2)[,1],
                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred
  sim.df <- rbind(train.df, estimation.df)
  # # bcf
  # sim.df$drcut.rmst.bcf.cate.hat <- XBCF::getTaus(XBCF::XBCF(y=as.matrix(sim.df$drcut.rmst), z=as.matrix(sim.df$a), x_con=as.matrix(sim.df[, confounders]), pihat=sim.df$ps, pcat_con=3, pcat_mod=3))
  # while(is.null(sim.df$drcut.rmst.bcf.cate.hat))
  # {sim.df$drcut.rmst.bcf.cate.hat <- XBCF::getTaus(XBCF::XBCF(y=as.matrix(sim.df$drcut.rmst), z=as.matrix(sim.df$a), x_con=as.matrix(sim.df[, confounders]), pihat=sim.df$ps, pcat_con=3, pcat_mod=3))}
  # sim.df$drcut.rmst.bcf.cate.hat <- Winsorize(minval = -taus[4], maxval = taus[4], sim.df$drcut.rmst.bcf.cate.hat)
  
  sim.df[, grep("cate.hat", names(sim.df), value=TRUE)] <- Winsorize(minval = -taus[4], maxval = taus[4], sim.df[, grep("cate.hat", names(sim.df), value=TRUE)])
  
  print(i)
  write.csv(sim.df, file=paste0("sim1-survival", i,".csv"), row.names = FALSE)
  NULL
}
end.time <- Sys.time()
end.time - start.time
