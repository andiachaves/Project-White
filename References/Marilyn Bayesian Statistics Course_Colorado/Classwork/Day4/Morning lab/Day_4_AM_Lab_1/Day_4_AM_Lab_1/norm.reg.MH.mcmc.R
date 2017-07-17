norm.reg.MH.mcmc <- function(y,X,betastrt,betamean,betavar,s2mean,s2sd,s2.tune,ngibbs,no.print=FALSE){

#
#  (20070202)  Mevin Hooten  [last modified:  20130513]
#  Multiple Linear Bayesian Regression  (iid errors) w/ DIC calculations.
#
#  Example Use:
#
#  simdata.df=read.table("simdata.txt")
#  pairs(simdata.df)
#  y=as.vector(simdata.df[,1])
#  X=as.matrix(simdata.df[,2:4])
#  ex15.out=norm.reg.MH.mcmc(y,X,solve(t(X)%*%X)%*%t(X)%*%y,c(0,0,0),100,0,1,.1,1000)
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
### Setup Variables and Priors
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
### Start Values
###

beta=matrix(betastrt,p,1)
s2=1

###
### Gibbs Loop
###

for(k in 2:ngibbs){
  if(k%%100==0) cat(k," ");flush.console()

  ###
  ### Sample s2
  ###

  s2.star=exp(rnorm(1,log(s2),s2.tune))  
  mh.1=sum(dnorm(y,X%*%beta,sqrt(s2.star),log=TRUE))+dnorm(log(s2.star)/2,s2mean,s2sd,log=TRUE)
  mh.2=sum(dnorm(y,X%*%beta,sqrt(s2),log=TRUE))+dnorm(log(s2)/2,s2mean,s2sd,log=TRUE)
  mh=exp(mh.1-mh.2)
 
  if(mh > runif(1)){
    s2=s2.star
  }

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

postbetamn=apply(betasave[,-(1:n.burn)],1,mean)
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
