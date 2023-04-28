#App 1

library(AdequacyModel)
library(MASS)
require(GenSA)
require(survival)


fit.sa2<- function(data,density) {
  minusllike<-function(x) -sum(log(density(data,x[1],x[2])))   
  lower <- c(0,0) 
  upper <- c(10,10)
  out <- GenSA(lower = lower, upper = upper, fn = minusllike, control=list(verbose=TRUE,maxit=7000,max.time=3))
  return(out[c("value","par","counts")])
}


fit.sa3<- function(data,density) {
  minusllike<-function(x) -sum(log(density(data,x[1],x[2],x[3]))) 
  lower <- c(0,0,0) 
  upper <- c(100,100,100)
  out <- GenSA(lower = lower, upper = upper, fn = minusllike, control=list(verbose=TRUE,maxit=10000,max.time=3))
  return(out[c("value","par","counts")])
}


fit.sa4<- function(data,density) {
  minusllike<-function(x) -sum(log(density(data,x[1],x[2],x[3],x[4]))) 
  lower <- c(0,0,0,0) 
  upper <- c(100,100,100,100)
  out <- GenSA(lower = lower, upper = upper, fn = minusllike, control=list(verbose=TRUE,max.time=2))
  return(out[c("value","par","counts")])
}

#######################

#WEW

######################
pdf.WEW = function(x,b,alpha,beta,lambda){
  a<-1
  g <- lambda*beta*x^(beta-1)*(1+alpha*lambda*x^beta)^(-1/alpha-1)
  G <- 1-(1+alpha*lambda*x^beta)^(-1/alpha)
  
  fd=a*b*g*(G^(b-1)/(1-G)^(b+1))*exp(-a*(G/(1-G))^(b))
  fd
}
set.seed(1729)
fit.sa4(dados,pdf.WEW)

pdf_wew = function(par, x){
  a<-1
  b = par[1]
  alpha = par[2]
  beta = par[3]
  lambda = par[4]
  g <- lambda*beta*x^(beta-1)*(1+alpha*lambda*x^beta)^(-1/alpha-1)
  G <- 1-(1+alpha*lambda*x^beta)^(-1/alpha)
  
  fd=a*b*g*(G^(b-1)/(1-G)^(b+1))*exp(-a*(G/(1-G))^(b))
  fd
}
cdf_wew = function(par, x){
  a<-1
  b = par[1]
  alpha = par[2]
  beta = par[3]
  lambda = par[4]
  g <- lambda*beta*x^(beta-1)*(1+alpha*lambda*x^beta)^(-1/alpha-1)
  G <- 1-(1+alpha*lambda*x^beta)^(-1/alpha)
  fa=1-exp(-a*(G/(1-G))^(b))
}
set.seed(1729)
resultswew  = goodness.fit(pdf = pdf_wew , cdf = cdf_wew ,
                           starts = c(1, 95, 35, 11),
                           data = dados, method = "SANN", domain = c(0, Inf),
                           mle = NULL);resultswew

chutes<-c(1, 95, 35, 11)
fit.gtnhWEW<- fitdistr(dados,pdf.WEW,start=list(b=chutes[1],alpha=chutes[2],beta=chutes[3], lambda=chutes[4]),control=list(ndeps=c(1e-6,1e-12,1e-8,1e-8),maxit=10000))
fit.gtnhWEW

#######################

#EXTENDED WEIBULL

###################### 

pdf.EW = function(x,alpha,beta,lambda){
  a<-1
  b<-1
  g <- lambda*beta*x^(beta-1)*(1+alpha*lambda*x^beta)^(-1/alpha-1)
  G <- 1-(1+alpha*lambda*x^beta)^(-1/alpha)
  
  fd=a*b*g*(G^(b-1)/(1-G)^(b+1))*exp(-a*(G/(1-G))^(b))
  fd
}
set.seed(1729)
fit.sa3(dados,pdf.EW)

pdf_ew = function(par, x){
  a<-1
  b = 1
  alpha = par[1]
  beta = par[2]
  lambda = par[3]
  g <- lambda*beta*x^(beta-1)*(1+alpha*lambda*x^beta)^(-1/alpha-1)
  G <- 1-(1+alpha*lambda*x^beta)^(-1/alpha)
  
  fd=a*b*g*(G^(b-1)/(1-G)^(b+1))*exp(-a*(G/(1-G))^(b))
  fd
}
cdf_ew = function(par, x){
  a<-1
  b = 1
  alpha = par[1]
  beta = par[2]
  lambda = par[3]
  g <- lambda*beta*x^(beta-1)*(1+alpha*lambda*x^beta)^(-1/alpha-1)
  G <- 1-(1+alpha*lambda*x^beta)^(-1/alpha)
  fa=1-exp(-a*(G/(1-G))^(b))
}
set.seed(1729)
resultsew  = goodness.fit(pdf = pdf_ew , cdf = cdf_ew ,
                          starts = c(100, 40.59750358 ,  0.04570209
                          ),
                          data = dados, method = "SANN", domain = c(0, Inf),
                          mle = NULL);resultsew


chutes<-c(100, 40.59750358 ,  0.04570209)
fit.gtnhEW<- fitdistr(dados,pdf.EW,start=list(alpha=chutes[1],beta=chutes[2], lambda=chutes[3]),control=list(ndeps=c(1e-6,1e-8,1e-12),maxit=10000))
fit.gtnhEW


#######################

#Weibull Weibull

###################### 

pdf.WW = function(x,b,beta,lambda){
  a<-1
  g <- lambda*beta*x^(beta-1)*exp(-lambda*x^beta)
  G<- 1-exp(-lambda*x^beta)
  
  fd=a*b*g*(G^(b-1)/(1-G)^(b+1))*exp(-a*(G/(1-G))^(b))
  fd
}
set.seed(1729)
fit.sa3(dados,pdf.WW)

pdf_ww = function(par, x){
  a<-1
  b = par[1]
  beta = par[2]
  lambda = par[3]
  
  g <- lambda*beta*x^(beta-1)*exp(-lambda*x^beta)
  G<- 1-exp(-lambda*x^beta)
  
  fd=a*b*g*(G^(b-1)/(1-G)^(b+1))*exp(-a*(G/(1-G))^(b))
  fd
}
cdf_ww = function(par, x){
  a<-1
  b = par[1]
  beta = par[2]
  lambda = par[3]
 
  g <- lambda*beta*x^(beta-1)*exp(-lambda*x^beta)
  G<- 1-exp(-lambda*x^beta)
  fa=1-exp(-a*(G/(1-G))^(b))
}
set.seed(1729)
resultsww  = goodness.fit(pdf = pdf_ww , cdf = cdf_ww ,
                          starts = c(100,  0.00838344,   0.68114562
                          ),
                          data = dados, method = "SANN", domain = c(0, 1000),
                          mle = NULL);resultsww

chutes<-c(100,  0.00838344,   0.68114562)
fit.gtnhWW<- fitdistr(dados,pdf.WW,start=list(b=chutes[1],beta=chutes[2], lambda=chutes[3]),control=list(ndeps=c(1e-8,1e-10,1e-12),maxit=10000))
fit.gtnhWW

####################### 

#W exponencial

######################

pdf.WE = function(x,b,lambda){
  a<-1
  beta=1
 
  g <- lambda*beta*x^(beta-1)*exp(-lambda*x^beta)
  G<- 1-exp(-lambda*x^beta)
  
  fd=a*b*g*(G^(b-1)/(1-G)^(b+1))*exp(-a*(G/(1-G))^(b))
  fd
}
set.seed(1729)
fit.sa2(dados,pdf.WE)

pdf_we = function(par, x){
  a<-1
 
  beta=1
  b = par[1]
  lambda = par[2]
 
  g <- lambda*beta*x^(beta-1)*exp(-lambda*x^beta)
  G<- 1-exp(-lambda*x^beta)
  
  fd=a*b*g*(G^(b-1)/(1-G)^(b+1))*exp(-a*(G/(1-G))^(b))
  fd
}
cdf_we = function(par, x){
  a<-1
  b = par[1]
  
  beta=1
  lambda = par[2]
  
  g <- lambda*beta*x^(beta-1)*exp(-lambda*x^beta)
  G<- 1-exp(-lambda*x^beta)
  
  fa=1-exp(-a*(G/(1-G))^(b))
}
set.seed(1729)
resultswe  = goodness.fit(pdf = pdf_we , cdf = cdf_we ,
                          starts = c(0.78266544, 0.06626933
                          ),
                          data = dados, method = "SANN", domain = c(0, Inf),
                          mle = NULL);resultswe



#######################

#GEW

######################

pdf.gew = function(x,a,alpha,beta,lambda){
  
  g <- lambda*beta*x^(beta-1)*(1+alpha*lambda*x^beta)^(-1/alpha-1)
  G <- 1-(1+alpha*lambda*x^beta)^(-1/alpha)
  
  fd=g/gamma(a)*(-log(1-G))^(a-1)
  fd
}

set.seed(1729)
fit.sa4(dados,pdf.gew)

pdf_gew = function(par, x){
  
  a= par[1]
  alpha = par[2]
  beta = par[3]
  lambda = par[4]
  g <- lambda*beta*x^(beta-1)*(1+alpha*lambda*x^beta)^(-1/alpha-1)
  G <- 1-(1+alpha*lambda*x^beta)^(-1/alpha)
  fd=g/gamma(a)*(-log(1-G))^(a-1)
  fd
}


cdf_gew = function(par, x){
  
  a= par[1]
  alpha = par[2]
  beta = par[3]
  lambda = par[4]
  g <- lambda*beta*x^(beta-1)*(1+alpha*lambda*x^beta)^(-1/alpha-1)
  G <- 1-(1+alpha*lambda*x^beta)^(-1/alpha)
  fa=pgamma(-log(1-G),a)
  
}

set.seed(1729)
resultsgew  = goodness.fit(pdf = pdf_gew , cdf = cdf_gew ,
                           starts = c(5,  0.8,  2, 100),
                           data = dados, method = "SANN", domain = c(0, Inf),
                           mle = NULL);resultsgew


chutes<-c(5,  0.8,  2, 100)
fit.gtnhGEW<- fitdistr(dados,pdf.gew,start=list(a=chutes[1],alpha=chutes[2],beta=chutes[3], lambda=chutes[4]),control=list(ndeps=c(1e-8,1e-6,1e-6,1e-12),maxit=10000))
fit.gtnhGEW


########################################################
#Kw-Wei 
#######################################################
kwei.pdf=function(x,a,b,k,beta){
  g <- (k/beta)*(x/beta)^(k-1)*exp(-(x/beta)^k)
  G <- 1-exp(-(x/beta)^k) 
  fd<-a*b*g*G^(a-1)*(1-G^a)^(b-1)
  fd
}
set.seed(1729)
fit.sa4(dados,kwei.pdf)

pdf_KwWei = function(par, x){
  a = par[1]
  b = par[2]
  k = par[3]
  beta = par[4]
  g <- (k/beta)*(x/beta)^(k-1)*exp(-(x/beta)^k)
  G <- 1-exp(-(x/beta)^k) 
  a*b*g*G^(a-1)*(1-G^a)^(b-1)
}

cdf_KwWei = function(par, x){
  a = par[1]
  b = par[2]
  k = par[3]
  beta = par[4]
  g <- (k/beta)*(x/beta)^(k-1)*exp(-(x/beta)^k)
  G <- 1-exp(-(x/beta)^k)
  fa<-1-(1-G^a)^b
}
set.seed(1729)
resultsKwWei = goodness.fit(pdf = pdf_KwWei, cdf = cdf_KwWei,
                            starts = c(17.691334, 100,   0.306549, 100), 
                            data = dados, method = "SANN", domain = c(0, Inf),
                            mle = NULL);resultsKwWei

chutes<-c(17.691334, 100,   0.306549, 100)
fit.gtnhkw<- fitdistr(dados,kwei.pdf,start=list(a=chutes[1],b=chutes[2],k=chutes[3], beta=chutes[4]),control=list(ndeps=c(1e-6,1e-10,1e-6,1e-12),maxit=10000))
fit.gtnhkw



#################
# Beta Weibull
#################
bw.pdf=function(x,a,b,k,beta){
  g <- (k/beta)*(x/beta)^(k-1)*exp(-(x/beta)^k)
  G <- 1-exp(-(x/beta)^k)
  fd=(1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1) 
  fd
}

set.seed(1729)
fit.sa4(dados,bw.pdf)

cdfbw=function(par,x)
{
  a=par[1]
  b=par[2]
  k = par[3]
  beta = par[4]
  g <- (k/beta)*(x/beta)^(k-1)*exp(-(x/beta)^k)
  G <- 1-exp(-(x/beta)^k)
  pbeta(G, a, b)
}
pdfbw=function(par,x)
{
  a=par[1]
  b=par[2]
  k = par[3]
  beta = par[4]
  g <- (k/beta)*(x/beta)^(k-1)*exp(-(x/beta)^k)
  G <- 1-exp(-(x/beta)^k)
  (1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1)
}
set.seed(1729)
resultsbw = goodness.fit(pdf = pdfbw, cdf = cdfbw,
                         starts = c(.01,  200,  0.01, 200),
                         data = dados, method = "SANN", domain = c(0, Inf),
                         mle = NULL);resultsbw

chutes<-c(.01,  200,  0.01, 200)
fit.gtnhbw<- fitdistr(dados,bw.pdf,start=list(a=chutes[1],b=chutes[2],k=chutes[3], beta=chutes[4]),control=list(ndeps=c(1e-6,1e-12,1e-6,1e-12),maxit=10000))
fit.gtnhbw

######################################################

#BGE OU BEE

######################################################



bee.pdf=function(x,a,b,lambda,alpha){
  g <- alpha*lambda*exp(-(lambda*x))*(1-exp(-(lambda*x)))^(alpha-1)
  G <- (1-exp(-(lambda*x)))^(alpha) 
  fd=(1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1) 
  fd
}
set.seed(1729)
fit.sa4(dados,bee.pdf)

pdf.bee=function(par,x){
  a=par[1]
  b=par[2]
  lambda = par[3]
  alpha = par[4]
  g <- alpha*lambda*exp(-(lambda*x))*(1-exp(-(lambda*x)))^(alpha-1)
  G <- (1-exp(-(lambda*x)))^(alpha) 
  fd=(1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1) 
  fd
}
cdf.bee=function(par,x){
  a=par[1]
  b=par[2]
  lambda = par[3]
  alpha = par[4]
  g <- alpha*lambda*exp(-(lambda*x))*(1-exp(-(lambda*x)))^(alpha-1)
  G <- (1-exp(-(lambda*x)))^(alpha) 
  pbeta(G, a, b)
}
set.seed(1729)
resultsbee = goodness.fit(pdf = pdf.bee, cdf = cdf.bee,
                          starts = c( 12.2066526,   2.2120488,   0.06489262, 10),
                          data = dados, method = "SANN", domain = c(0, Inf),
                          mle = NULL);resultsbee

chutes<-c(12.2066526,   2.2120488,   0.06489262, 10)
fit.gtnhbee<- fitdistr(dados,bee.pdf,start=list(a=chutes[1],b=chutes[2],lambda=chutes[3], alpha=chutes[4]),control=list(ndeps=c(1e-6,1e-6,1e-6,1e-8),maxit=10000))
fit.gtnhbee

#######################################################
#Kw-Gamma
#######################################################
kga.pdf=function(x,a,b,alpha,beta){
  g=dgamma(x,alpha,beta)
  G=pgamma(x,alpha,beta)
  fd<-a*b*g*G^(a-1)*(1-G^a)^(b-1)
  fd
}
set.seed(1729)
fit.sa4(dados,kga.pdf)

pdf_Kwga = function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  beta = par[4]
  g=dgamma(x,alpha,beta)
  G=pgamma(x,alpha,beta)
  a*b*g*G^(a-1)*(1-G^a)^(b-1)
  
}

cdf_Kwga = function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  beta = par[4]
  g=dgamma(x,alpha,beta)
  G=pgamma(x,alpha,beta)
  fa<-1-(1-G^a)^b
}
set.seed(1729)
resultsKwga = goodness.fit(pdf = pdf_Kwga, cdf = cdf_Kwga,
                           starts = c( 4.0067102, 0.2894311, 0.7111096, 0.5392743), 
                           data = dados, method = "SANN", domain = c(0, Inf),
                           mle = NULL);resultsKwga

chutes<-c( 4.0067102, 0.2894311, 0.7111096, 0.5392743)
fit.gtnhkga<- fitdistr(dados,kga.pdf,start=list(a=chutes[1],b=chutes[2],alpha=chutes[3], beta=chutes[4]),control=list(ndeps=c(1e-6,1e-6,1e-6,1e-6),maxit=10000))
fit.gtnhkga

#################
# Beta Gamma
#################
bga.pdf=function(x,a,b,alpha,beta){
  g=dgamma(x,alpha,beta)
  G=pgamma(x,alpha,beta)
  fd=(1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1) 
  fd
}
set.seed(1729)
fit.sa4(dados,bga.pdf)

cdfbga=function(par,x)
{
  a=par[1]
  b=par[2]
  alpha=par[3]
  beta=par[4]
  g=dgamma(x,alpha,beta)
  G=pgamma(x,alpha,beta)
  pbeta(G, a, b)
}
pdfbga=function(par,x)
{
  a=par[1]
  b=par[2]
  alpha=par[3]
  beta=par[4]
  g=dgamma(x,alpha,beta)
  G=pgamma(x,alpha,beta)
  (1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1)
}
set.seed(1729)
resultsbga = goodness.fit(pdf = pdfbga, cdf = cdfbga,
                          starts = c(0.03046034,  0.19825689, 40.77597868,  2.31795404),
                          data = dados, method = "SANN", domain = c(0, Inf),
                          mle = NULL);resultsbga

############
chutes<-c(0.03046034,  0.19825689, 40.77597868,  2.31795404)
fit.gtnhbga<- fitdistr(dados,bga.pdf,start=list(a=chutes[1],b=chutes[2],alpha=chutes[3], beta=chutes[4]),control=list(ndeps=c(1e-12,1e-12,1e-12,1e-10),maxit=10000))
fit.gtnhbga


#plots

truehist(dados,
         ylim=c(0,0.15),
         col = "white",ylab="f(x)",xlab = "x",nbins = 11)
curve(pdf_wew(fit.gtnhWEW$estimate,x),add=TRUE, lwd = 3, col="red",lty=1)
curve(pdf_gew(fit.gtnhGEW$estimate,x),add=TRUE, lwd = 3, col="yellow",lty=2)
curve(pdf_Kwga(fit.gtnhkga$estimate,x),add=TRUE, lwd = 3, col="blue",lty=3)
curve(pdf.bee(resultsbee$mle,x),add=TRUE, lwd = 3, col="green",lty=4)

legend(25,.15, legend = c( "WEW", "GEW", "KwGa", "BEE"),
       col = c("red","yellow","blue","green"), lwd = 3 ,lty = c(1,2,3,4) , bty ="n")


km<- survfit(Surv(dados) ~ 1) 
plot(km$time, 1-km$surv, xlab = "x", ylab="F(x)",lwd=2,type = "l")
curve(cdf_wew(fit.gtnhWEW$estimate,x),add=TRUE, lwd = 3, col="red",lty=1)
curve(cdf_gew(fit.gtnhGEW$estimate,x),add=TRUE, lwd = 3, col="yellow",lty=2)
curve(cdf_Kwga(fit.gtnhkga$estimate,x),add=TRUE, lwd = 3, col="blue",lty=3)
curve(cdf.bee(resultsbee$mle,x),add=TRUE, lwd = 3, col="green",lty=4)

legend(25,.6, legend = c( "WEW", "GEW", "KwGa", "BEE"),
       col = c("red","yellow","blue","green"), lwd = 3 ,lty = c(1,2,3,4) , bty ="n")




