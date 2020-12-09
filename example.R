rm(list=ls())
library(tidyverse)
source("helpers.R")
library(parallel)
library(doSNOW)
library(foreach)

ncl = detectCores(logical = FALSE)
cl = makeCluster(ncl)
registerDoSNOW(cl)

pi1=0.01; beta1=1; beta2=0.8; beta3=0.2; link="probit"
beta = c(beta1,beta2,beta3)
beta0 = beta0_find(pi1,beta,link)
beta_true = c(beta0,beta)

## One-marker
gamma1_hat = ParaEst_OneMarker(beta_true,link)
gamma = c(gamma1_hat,0,0); gamma.all = gamma
x.all = crossing(x1=c(-5,5),x2=c(-5,5)) %>% as.data.frame()
y.all = gamma[1] + gamma[2]*x.all$x1 + gamma[3]*x.all$x2
y.min = min(y.all)
y.max = max(y.all)
ap1 = ap_fun(gamma,pi1,beta_true,y.min,y.max,link)
auc1 = auc_fun(gamma,pi1,beta_true,y.min,y.max,link)
bs1 = bs_fun(gamma,beta_true,link)
sbs1 = 1 - bs1/(pi1*(1-pi1))

gamma2_hat = ParaEst_TwoMarker(beta_true,link)
gamma = c(gamma2_hat,0); gamma.all = rbind(gamma.all,gamma)
x.all = crossing(x1=c(-5,5),x2=c(-5,5)) %>% as.data.frame()
y.all = gamma[1] + gamma[2]*x.all$x1 + gamma[3]*x.all$x2
y.min = min(y.all)
y.max = max(y.all)
ap2 = ap_fun(gamma,pi1,beta_true,y.min,y.max,link, given="x1")
auc2 = auc_fun(gamma,pi1,beta_true,y.min,y.max,link, given="x1")
## Remark: if gamma[2] is too small, suggesting using given="x2"
bs2 = bs_fun(gamma,beta_true,link)
sbs2 = 1 - bs2/(pi1*(1-pi1))