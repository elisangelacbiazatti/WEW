require(maxLik)
require(survival)


#TTT plot#
require(AdequacyModel)
TTT(tempo, col="red", lwd=2.5, grid=TRUE, lty=2)

x0=c(rep(1,29))
n=length(tempo)
x1=grupo
y=log(tempo)


################
#log Weibull Extended Weibull
################
log.veroLWEW <- function(par){
  beta0 <- par[1]
  beta1 <- par[2]
  sigma <- par[3]
  b <- par[4]
  alpha <- par[5]
  
  mu <- beta0*x0 + beta1*x1 
  w <- (y-mu)/sigma
  
  df=(b/sigma)*exp(w)*(1+alpha*exp(w))^(1/alpha-1)*((1+alpha*exp(w))^(1/alpha)-1)^(b-1)*
    exp(-((1+alpha*exp(w))^(1/alpha)-1)^b)
  
  st<-exp(-((1+alpha*exp(w))^(1/alpha)-1)^b)
  
  lv <- censura*log(df)+(1-censura)*log(st)
  logv <- sum(lv)
  if ((sigma>0)&&(b>0)&&(alpha>0))
    return(logv)
  else return (-Inf)
}

set.seed(1729) 
V1<-maxLik(log.veroLWEW, start=c(6.519875, 1 , 4,  2,  5),method="BFGS")

summary(V1)
AIC(V1)

logveroLWEW=logLik(V1)
logveroLWEW


npwew<-5
AICw<-(-2*logveroLWEW)+(2*npwew)  
AICw

AICcw<-AICw + ((2*npwew^2+2*npwew)/(n-npwew-1))
AICcw

BICw<-(-2*logveroLWEW) + npwew*log(n)
BICw


################
#log Extended Weibull
################
log.veroLEW <- function(par){
  beta0 <- par[1]
  beta1 <- par[2]
  sigma <- par[3]
  b <- 1
  alpha <- par[4]
  
  mu <- beta0*x0 + beta1*x1 
  w <- (y-mu)/sigma
  
  df=(b/sigma)*exp(w)*(1+alpha*exp(w))^(1/alpha-1)*((1+alpha*exp(w))^(1/alpha)-1)^(b-1)*
    exp(-((1+alpha*exp(w))^(1/alpha)-1)^b)
  
  st<-exp(-((1+alpha*exp(w))^(1/alpha)-1)^b)
  
  lv <- censura*log(df)+(1-censura)*log(st)
  logv <- sum(lv)
  if ((sigma>0)&&(alpha>0))
    return(logv)
  else return (-Inf)
}
set.seed(1729)
v1EW<-maxLik(log.veroLEW, start=c(.3 ,5, 1, 5),method="BFGS")
summary(v1EW)
AIC(v1EW)

logveroLEW=logLik(v1EW)
logveroLEW

npew<-4
AICEw<-(-2*logveroLEW)+(2*npew) 
AICEw

AICcEw<-AICEw + ((2*npew^2+2*npew)/(n-npew-1))
AICcEw

BICEw<-(-2*logveroLEW) + npew*log(n)
BICEw


#TRV
logveroLEW     
logveroLWEW   
TRV<-2*(logveroLWEW-logveroLEW)
TRV
1-pchisq(TRV,1)


################
#log Weibull Weibull
################
log.veroLWW <- function(par){
  beta0 <- par[1]
  beta1 <- par[2]
  sigma <- par[3]
  b <- par[4]
  
  
  mu <- beta0*x0 + beta1*x1 
  w <- (y-mu)/sigma
  
  df=(b/sigma)*exp(w+exp(w)-(exp(exp(w))-1)^b)*(exp(exp(w))-1)^(b-1)
  
  
  st<-exp(-(exp(exp(w))-1)^b)
  
  lv <- censura*log(df)+(1-censura)*log(st)
  logv <- sum(lv)
  if ((sigma>0)&&(b>0))
    return(logv)
  else return (-Inf)
}
set.seed(1729)
v12<-maxLik(log.veroLWW, start=c(2.519875, 1.7 , 2 , 4),method="BFGS")

summary(v12)
AIC(v12)
logveroLWW=logLik(v12)

np12<-4
AIC12<-(-2*logveroLWW)+(2*np12)
AIC12

AICc12<-AIC12 + ((2*np12^2+2*np12)/(n-np12-1))
AICc12

BIC12<-(-2*logveroLWW) + np12*log(n)
BIC12




TRV1<-2*(logveroLWEW-logveroLWW)
TRV1

1-pchisq(TRV1,1)

################
#log Gamma Extended Weibull
################
log.veroLGEW <- function(par){
  beta0 <- par[1]
  beta1 <- par[2]
  sigma <- par[3]
  a <- par[4]
  alpha <- par[5]
  
  mu <- beta0*x0 + beta1*x1 
  w <- (y-mu)/sigma
  
  df=((exp(w)*(1+alpha*exp(w))^(-1/alpha-1))/(sigma*alpha^(a-1)*gamma(a)))*(log(1+alpha*exp(w)))^(a-1)  
  
  
  st<-1-pgamma((1/alpha)*log(1+alpha*exp(w)),a)
  
  lv <- censura*log(df)+(1-censura)*log(st)
  logv <- sum(lv)
  if ((sigma>0)&&(a>0)&&(alpha>0))
    return(logv)
  else return (-Inf)
}
set.seed(1729)
V2<-maxLik(log.veroLGEW, start=c(2.5, 1.849902 , 2 , 0.9, 4),method="BFGS")
summary(V2)


logveroLGEW=logLik(V2)
logveroLGEW

npgew<-5
AICg<-(-2*logveroLGEW)+(2*npgew)
AICg

AICcg<-AICg + ((2*npgew^2+2*npgew)/(n-npgew-1))
AICcg

BICg<-(-2*logveroLGEW) + npgew*log(n)
BICg


###############
#log beta Weibull 
###############
log.verob <- function(par){
  beta0 <- par[1]
  beta1 <- par[2]
  sigma <- par[3]
  a<-par[4]
  b<-par[5]  
  mu <- beta0*x0 + beta1*x1 
  w <- (y-mu)/sigma
  df <- (1/(beta(a,b)*sigma))*exp(w-b*exp(w))*(1-exp(-exp(w)))^(a-1)
  G <- (1-exp(-exp(w)))
  logv <- sum(log(df))
  if ((sigma>0)&&(a>0)&&(b>0))
    return(logv)
  else return (-Inf)
}
set.seed(1729)
V11b<-maxLik(log.verob, start=c(5,  0.056187469,  0.830053984, 0.536347759, 0.5),method="BFGS")
summary(V11b)

logVerowb<-logLik(V11b)
logVerowb

npB<-5
AICB<-(-2*logVerowb)+(2*npB)  
AICB

AICcB<-AICB + ((2*npB^2+2*npB)/(n-npB-1))
AICcB

BICB=(-2*logVerowb) + npB*log(n)
BICB



################
#log Kw Weibull
################
log.verokw <- function(par){
  beta0 <- par[1]
  beta1 <- par[2]
  sigma <- par[3]
  a <- par[4]
  b <- par[5]  
  mu <- beta0*x0 + beta1*x1 
  w <- (y-mu)/sigma
  df=((a*b)/sigma)*exp(w-exp(w))*(1-exp(-exp(w)))^(a-1)*(1-(1-exp(-exp(w)))^a)^(b-1)
  logv <- sum(log(df))
  if ((sigma>0)&&(a>0)&&(b>0))
    return(logv)
  else return (-Inf)
}
set.seed(1729)
vkw<-maxLik(log.verokw, start=c(2,  0.056187469,  1.5, 0.536347759, 0.5),method="BFGS")
summary(vkw)

logVerokw<-logLik(vkw)
logVerokw


npkw<-5
AICkw<-(-2*logVerokw)+(2*npkw)  
AICkw

AICckw<-AICkw + ((2*npkw^2+2*npkw)/(n-npkw-1))
AICckw

BICkw=(-2*logVerokw) + npkw*log(n)
BICkw


