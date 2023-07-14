library(glmnet);library(drgee)
source('gen_binary.R')
n=500;p=15;d=p+1
kn=log(log(n))*log(2*p)
expit<-function(x) exp(x)/(1+exp(x))

ds=gen_binary(n,p)
x <- ds[[1]];y=ds[[2]];A=ds[[3]]
X=cbind(1,x)
D=exp(-x[,2]^2-x[,1]^2+x[,3]-x[,4])

# treatment correct, using A-learning

## fitted by drgee package
fit = drgee(oformula=y~x,
            eformula=A~x+D,
            iaformula=~x, olink="logit", elink="logit",
            estimation.method="dr")

## coefficients of baseline model
beta=fit$coefficients.all[(d*5-p+1):(d*5+1)]

## blip parameter
psi0=fit$coefficients


##  ## coefficients of treatment model

alpha=fit$coefficients.all[(d*3-p):(d*3+1)]
an=expit(cbind(X,D)%*%alpha)
yn=expit(X%*%beta)
pistar=1+(1-an)*yn/(an*expit(X%*%(psi0+beta)))

## constructing the weights
w=abs(A-1/pistar)
w=c(w/sum(w)*n)

## penalty factor, also known as the adaptive weights
pfac=1/abs(c(beta[-1],psi0))
pfac[d]=0
pfac=pfac^3

## run glmnet with weights
m=glmnet(cbind(x,A,A*x), y, family = "binomial",weights = w,
         penalty.factor=pfac)

dev=(1-m$dev.ratio)*m$nulldev
s=m$df;gic=1/n *(dev+kn*s)
id=which.min(gic)
psi0=coef(m)[,id][(d+1):(2*d)]

### obtain our final estimator
psi0
