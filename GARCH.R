setwd("~/Documents/2016Spring/FIN580 Risk Mgt/HW8")

library("fGarch")
library("stats4")
#data set
data_efa<-read.csv(file = "EFA.csv",header = TRUE, sep = ",")
data_iwm<-read.csv(file = "IWM.csv",header = TRUE, sep = ",")
#price level
efa<-data_efa$Adj.Close
iwm<-data_iwm$Adj.Close
#simple return
efa_ret<-efa[-1]/efa[-length(efa)]-1
iwm_ret<-iwm[-1]/iwm[-length(iwm)]-1

#GARCH
#efa_gch_fit<-garchFit(formula = ~garch(1,1),data = efa_ret,cond.dist = c("norm"))
#coef(efa_gch_fit)
#efa_sigma_1<-volatility(efa_gch_fit)[1]

#initial guess
efa_s1<-sd(efa_ret)
efa_s<-sd(efa_ret)
iwm_s1<-sd(iwm_ret)
iwm_s<-sd(iwm_ret)

#NGARCH parameter estimation
sigma<-rep(0,1010)
Nll<-function(a,b,s1,s,theta){
  
  if (a+b<1 && s1>0 && s){
    sigma[1] <- s1
    sum = -0.5*(log(sigma[1]^2)+R[1]^2/sigma[1]^2)
    for (t in 1:1009){
      
      sigma[t+1] <- sqrt((1-a*(1+theta^2)-b)*s^2+a*(R[t]-theta*sigma[t])^2+b*sigma[t]^2)
      sum = sum + -0.5*(log(sigma[t]^2)+R[t]^2/sigma[t]^2)
    }
    return (-sum)
  }
  else NA
}

#EFA NGARCH fit
R<-efa_ret
fit_efa_1<-mle(Nll,start = list(a = 0.1, b = 0.9 ,s1 = efa_s1 ,s= efa_s,theta = 0),method = "BFGS")
fit_efa_1
a = 0.04559882; b=0.86976272;  s1=0.01506999; s=0.01576239;  theta=1.30313525
efa_ngarch_list<-c(a,b,s1,s,theta)
efa_sigma<-rep(0,1010)
efa_sigma[1]<-s1
for(i in 2:1010)
{
  efa_sigma[i]<-sqrt((1-a*(1+theta^2)-b)*s^2+a*(efa_ret[i-1]-theta*efa_sigma[i-1])^2+b*efa_sigma[i-1]^2)
}

#IWM NGARCH fit
R<-iwm_ret
fit_iwm_1<-mle(Nll,start = list(a = 0.05, b = 0.9 ,s1 = iwm_s1 ,s= iwm_s,theta = 0),method = "BFGS")
fit_iwm_1
a=0.02965582; b=0.89680813; s1=0.01445754; s=0.01373495; theta=1.50590326 
iwm_ngarch_list<-c(a,b,s1,s,theta)
iwm_sigma<-rep(0,1010)
iwm_sigma[1]<-s1
for(i in 2:1010)
{
  iwm_sigma[i]<-sqrt((1-a*(1+theta^2)-b)*s^2+a*(iwm_ret[i-1]-theta*iwm_sigma[i-1])^2+b*iwm_sigma[i-1]^2)
}

#DCC

#normalized return
efa_z_t<-efa_ret/efa_sigma[1:1009]
iwm_z_t<-iwm_ret/iwm_sigma[1:1009]
pair_z_t<-efa_z_t*iwm_z_t

rho_12<-rep(0,1010)
q11<-rep(0,1010)
q22<-rep(0,1010)
q12<-rep(0,1010)

dcc<-function(alpha,beta)
{
  q11[1]<-1
  q22[1]<-1
  q12[1]<-mean(pair_z_t)
  rho_12[1]<-q12[1]
  sum=-0.5*(log(1-rho_12[1]^2)+(efa_z_t[1]^2+iwm_z_t[1]^2-2*rho_12[1]*pair_z_t[1])/(1-rho_12[1]^2))
  for(t in 1:1009){
  q11[t+1]<-1+alpha*(efa_z_t[t]^2-1)+beta*(q11[t]-1)
  q22[t+1]<-1+alpha*(iwm_z_t[t]^2-1)+beta*(q22[t]-1)
  q12[t+1]<-rho_12[1]+alpha*(pair_z_t[t]-rho_12[1])+beta*(q12[t]-rho_12[1])
  rho_12[t+1]<-q12[t+1]/sqrt(q11[t+1]*q22[t+1])
  sum<-sum-0.5*(log(1-rho_12[t+1]^2)+(efa_z_t[t+1]^2+iwm_z_t[t+1]^2-2*rho_12[t+1]*pair_z_t[t+1])/(1-rho_12[t+1]^2))
  }
  return (-sum)
}

fit_dcc_1<-mle(dcc,start = list(alpha = 0.05, beta = 0.9),method = "BFGS")
fit_dcc_1
alpha=0.02310267; beta=0.97769212 
dcc_list<-c(alpha,beta)
q11[1]<-1
q22[1]<-1
q12[1]<-mean(pair_z_t)
rho_12[1]<-q12[1]
for(t in 1:1009)
{
  q11[t+1]<-1+alpha*(efa_z_t[t]^2-1)+beta*(q11[t]-1)
  q22[t+1]<-1+alpha*(iwm_z_t[t]^2-1)+beta*(q22[t]-1)
  q12[t+1]<-rho_12_initial+alpha*(pair_z_t[t]-rho_12_initial)+beta*(q12[t]-rho_12_initial)
  rho_12[t+1]<-q12[t+1]/sqrt(q11[t+1]*q22[t+1])
}
plot(rho_12,type="l")

#forecast with MC
risk_free<-1.2/25200
portfolio_value<-890
for(x in 2:10000)
{
#fcst FIRST sigma
fcst_efa_sigma<-rep(0,252)
fcst_iwm_sigma<-rep(0,252)
fcst_efa_sigma[1]<-efa_sigma[1010]
fcst_iwm_sigma[1]<-iwm_sigma[1010]
#fcst unconditional z
fcst_efa_z_t<-rnorm(252,0,1)
fcst_iwm_z_t<-rnorm(252,0,1)
#fcst FIRST q, rho
fcst_q11<-rep(0,252)
fcst_q22<-rep(0,252)
fcst_q12<-rep(0,252)
fcst_rho_12<-rep(0,252)
fcst_q11[1]<-q11[1010]
fcst_q22[1]<-q22[1010]
fcst_q12[1]<-q12[1010]
fcst_rho_12[1]<-rho_12[1010]
#fcst conditional z, using rho and unconditiaonal z
fcst_efa_z_cdtl<-rep(0,252)
fcst_iwm_z_cdtl<-rep(0,252)
fcst_pair_z_cdtl<-rep(0,252)
fcst_efa_z_cdtl<-fcst_efa_z_t
fcst_iwm_z_cdtl[1]<-fcst_rho_12[1]*fcst_efa_z_t[1]+sqrt(1-fcst_rho_12[1]^2)*fcst_iwm_z_t[1]
fcst_pair_z_cdtl[1]<-fcst_efa_z_cdtl[1]*fcst_iwm_z_cdtl[1]
#fcst FIRST return, using sigma and conditional z
fcst_efa_return<-rep(0,252)
fcst_iwm_return<-rep(0,252)
fcst_efa_return[1]<-fcst_efa_sigma[1]*fcst_efa_z_cdtl[1]-0.5*fcst_efa_sigma[1]^2
fcst_iwm_return[1]<-fcst_iwm_sigma[1]*fcst_iwm_z_cdtl[1]-0.5*fcst_iwm_sigma[1]^2

#fcst 2:252
for(t in 1:251)
{
  #updata sigma
  fcst_efa_sigma[t+1]<-sqrt((1-efa_ngarch_list[1]*(1+efa_ngarch_list[5]^2)-efa_ngarch_list[2])*efa_ngarch_list[4]^2+efa_ngarch_list[1]*(fcst_efa_return[t]-efa_ngarch_list[5]*fcst_efa_sigma[t])^2+efa_ngarch_list[2]*fcst_efa_sigma[t]^2)
  fcst_iwm_sigma[t+1]<-sqrt((1-iwm_ngarch_list[1]*(1+iwm_ngarch_list[5]^2)-iwm_ngarch_list[2])*iwm_ngarch_list[4]^2+iwm_ngarch_list[1]*(fcst_iwm_return[t]-iwm_ngarch_list[5]*fcst_iwm_sigma[t])^2+iwm_ngarch_list[2]*fcst_iwm_sigma[t]^2)
  #updata q, rho
  fcst_q11[t+1]<-1+alpha*(fcst_efa_z_cdtl[t]^2-1)+beta*(fcst_q11[t]-1)
  fcst_q22[t+1]<-1+alpha*(fcst_iwm_z_cdtl[t]^2-1)+beta*(fcst_q22[t]-1)
  fcst_q12[t+1]<-rho_12[1]+alpha*(fcst_pair_z_cdtl[t]-rho_12[1])+beta*(fcst_q12[t]-rho_12[1])
  fcst_rho_12[t+1]<-fcst_q12[t+1]/sqrt(fcst_q11[t+1]*fcst_q22[t+1])
  #updata conditional z (only for IWM and pairwise)
  fcst_iwm_z_cdtl[t+1]<-fcst_rho_12[t+1]*fcst_efa_z_t[t+1]+sqrt(1-fcst_rho_12[t+1]^2)*fcst_iwm_z_t[t+1]
  fcst_pair_z_cdtl[t+1]<-fcst_efa_z_cdtl[t+1]*fcst_iwm_z_cdtl[t+1]
  #fcst return
  fcst_efa_return[t+1]<-fcst_efa_sigma[t+1]*fcst_efa_z_cdtl[t+1]+risk_free
  fcst_iwm_return[t+1]<-fcst_iwm_sigma[t+1]*fcst_iwm_z_cdtl[t+1]+risk_free
}
######  plot(fcst_efa_return)
#portfolio price
efa_maturity_price<-efa[1010]
iwm_maturity_price<-iwm[1010]
for(t in 1:252)
{
  efa_maturity_price<-efa_maturity_price*(1+fcst_efa_return[t])
  iwm_maturity_price<-iwm_maturity_price*(1+fcst_iwm_return[t])
}

if(iwm_maturity_price>=94.673 && efa_maturity_price>=55.752)
{
  portfolio_value<-(x-1)/x*portfolio_value+1/x*1000
}
else
{
  loss<-1000*1.1765*(min(efa_maturity_price/55.752*0.85-1,iwm_maturity_price/94.673*0.85-1)+0.15)
  portfolio_value<-(x-1)/x*portfolio_value+1/x*(1000+loss)
}
}

(portfolio_value+1000*0.0502)/(1+risk_free)^255
#coupon_date<-seq(from=21, to=252, by= 21)
#coupons<-10*5.02/(1+risk_free)^coupon_date
