norm.reg.mcmc <- function(y,X,betastrt,betamean,betavar,s2mean,s2sd,ngibbs,no.print=FALSE){

#
#  (20070202)  Mevin Hooten  [last modified:  20110915]
#  Multiple Linear Bayesian Regression  (iid errors) w/ DIC calculations.
#
#  Example Use:
#
#  simdata.df=read.table("simdata.txt")
#  pairs(simdata.df)
#  y=as.vector(simdata.df[,1])
#  X=as.matrix(simdata.df[,2:4])
#  tmp.out=norm.reg.mcmc(y,X,solve(t(X)%*%X)%*%t(X)%*%y,c(0,0,0),100,1,100,1000)
#  X1=X[,-3]
#  tmp.out=norm.reg.mcmc(y,X1,solve(t(X1)%*%X1)%*%t(X1)%*%y,c(0,0),100,1,100,1000)
#  X2=X[,-2]
#  tmp.out=norm.reg.mcmc(y,X2,solve(t(X2)%*%X2)%*%t(X2)%*%y,c(0,0),100,1,100,1000)
#
#

###
### Subroutines 
###

invgammastrt <- function(igmn,igvar){
  q <- 2+(igmn^2)/igvar
  r <- 1/(igmn*(q-1))
  list(r=r,q=q)
}

###
### Hyperpriors
###

n=dim(X)[1]
p=dim(X)[2]
n.burn=round(.1*ngibbs)
r=invgammastrt(s2mean,(s2sd^2))$r
q=invgammastrt(s2mean,(s2sd^2))$q
beta0=matrix(betastrt,p,1)
Sig0=betavar*diag(p)

betasave=matrix(0,p,ngibbs)
s2save=rep(0,ngibbs)
Davgsave=rep(0,ngibbs)
y.pred.mean=rep(0,n)

###
### Starting Values
###

beta=matrix(betastrt,p,1)

###
### Gibbs Loop
###

for(k in 2:ngibbs){
  if(k%%100==0) cat(k," ");flush.console()

  ###
  ### Sample s2
  ###

  tmpr=(1/r+.5*t(y-X%*%beta)%*%(y-X%*%beta))^(-1)
  tmpq=n/2+q
  s2=1/rgamma(1,tmpq,,tmpr) 

  ###
  ### Sample beta
  ###

  tmpvar=solve(t(X)%*%X/s2 + solve(Sig0))
  tmpmean=tmpvar%*%(t(X)%*%y/s2 + solve(Sig0)%*%betamean)
  beta=tmpmean+t(chol(tmpvar))%*%matrix(rnorm(p),p,1)

  ###
  ### DIC Calculations 
  ###

  Davgsave[k]=-2*sum(dnorm(y,X%*%beta,sqrt(s2),log=TRUE))

  ###
  ### Posterior Predictive Calculations 
  ###

  if(k > n.burn){
    y.pred=rnorm(n,X%*%beta,sqrt(s2))
    y.pred.mean=y.pred.mean+y.pred/(ngibbs-n.burn)
  }

  ###
  ### Save Samples
  ###
  
  betasave[,k]=beta
  s2save[k]=s2

}
cat("\n");flush.console()

###
###  Calculate DIC
###

if(dim(X)[2]==1){
  postbetamn=mean(betasave[,-(1:n.burn)])
}
if(dim(X)[2]>1){
  postbetamn=apply(betasave[,-(1:n.burn)],1,mean)
}
posts2mn=mean(s2save[-(1:n.burn)])
cat("Posterior Mean for Beta:","\n")
print(postbetamn)
cat("Posterior Mean for s2:","\n")
print(posts2mn)
Dhat=-2*(sum(dnorm(y,X%*%postbetamn,sqrt(posts2mn),log=TRUE)))
Davg=mean(Davgsave[-(1:n.burn)])
pD=Davg-Dhat
DIC=2*Davg-Dhat

cat("Dhat:",Dhat,"Davg:",Davg,"pD:",pD,"DIC:",DIC,"\n")

###
###  Calculate MSE 
###

MSE=mean((y-y.pred.mean)^2)
cat("RMSE:",sqrt(MSE),"\n")

###
###  Write Output
###

list(betasave=betasave,s2save=s2save,y=y,X=X,ngibbs=ngibbs,n=n,r=r,q=q,p=p,Dhat=Dhat,Davg=Davg,pD=pD,DIC=DIC,MSE=MSE)

}
