gen_binary=function(n,p){
  X=matrix(rnorm(n*p,mean = .5,sd=1),nrow = n)

  x1=X[,1];x2=X[,2]
  free=x1-.2*x2-exp(-x2^2-x1^2+X[,3]-X[,4])
  expit=function(x) exp(x)/(1+exp(x))

  q1=expit(free+1-2*x1)
  q0=expit(free)
  
  pro=q0*(q1+1-q0)/(1-(q1-q0)^2)
  
  a<-rbinom(n,1,prob = pro)
  mu=free+1*a-2*a*x1
  y<-rbinom(n,1,expit(mu))
  return(list(X,y,a))
}
