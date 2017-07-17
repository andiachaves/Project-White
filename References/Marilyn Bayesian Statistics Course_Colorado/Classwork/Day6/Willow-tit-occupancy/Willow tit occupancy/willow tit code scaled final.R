library(rjags)
rm(list=ls())

##Summarize original data to deal with missing values
setwd("/Users/Tom/Documents/NSF Statistics Workshop/2013 Course/Hobbs course materials/Occupancy lab/Willow tit occupancy/")

obs=read.csv("wtmatrix.csv",header=TRUE,sep=",",na.strings=c("NA"))
obs=obs[,c("y.1","y.2","y.3","elev","forest")]
y<-as.matrix(obs[,c("y.1","y.2","y.3")])
n<-apply(!is.na(y),1,sum)  #nice way to sum non missing, n is number of visits
M<-nrow(y)
y<-apply(y,1,sum,na.rm=TRUE) #y is number of times birds were observed

obs=as.data.frame(cbind(n,y,obs$elev,obs$forest))
names(obs)=c("number_visits", "number_detections", "elev", "forest")
write.csv(obs,file="Swiss BB data.csv")


######Student work starts here
obs=read.csv(file="Swiss BB data.csv")
#create vectors of standardized covariates
elev<-as.vector(scale(obs[,"elev"],center=TRUE))
forest<-as.vector(scale(obs[,"forest"],center=TRUE))
elev2<-elev*elev

sink("model.txt")
cat("
      model { 
      	
      	#taus for dnorm can be larger than ususal because predictors are standarized 
      	p ~ dunif(0,1)                
        b0 ~ dnorm(0,.0001) 
        b1 ~ dnorm(0,.0001) 
        b2 ~ dnorm(0,.0001) 
        b3 ~ dnorm(0,.0001) 
        logit(psi0)<-b0
        for(i in 1:length(y)){ 
           z[i] ~ dbin(psi[i],1)      
           logit(psi[i]) <- b0 + b1*elev[i] + b2*elev2[i] + b3*forest[i] 
   
           tmp[i]<-z[i]*p          
           y[i] ~ dbin(tmp[i],n[i])  
         } 
        
       #estimate optimum elevation at mean forest cover for scaled scaled values
       elev.xc<- -(1/2)*b1/b2
       #back transform to get elevation in meters
       elev.opt <- elev.xc *sd.elev + mu.elev
       for(j in 1:length(elev.x)){
       	logit(psi.elev[j]) <- b0 + b1*(elev.x[j]-mu.elev)/sd.elev + b2*((elev.x[j]-mu.elev)/sd.elev)^2
		}
     
      }#end of model

      

",fill=TRUE)
sink()



#make some elevations for predicting occupancy
elev.x=seq(500,2500,100)
data <- list ( y=obs$number_detections,n=obs$number_visits,forest=forest,elev=elev,elev2=elev2, mu.elev=mean(obs$elev), sd.elev=sd(obs$elev),elev.x=elev.x)


inits= list(
list( z=rep(1,M),b0=rnorm(1),p=runif(1),b1=rnorm(1),b2=rnorm(1),b3=rnorm(1)),
list ( z=rep(1,M),b0=rnorm(1),p=runif(1),b1=rnorm(1),b2=rnorm(1),b3=rnorm(1)) 
)


parameters <- c("p","b0","b1","b2","b3")
load.module("dic")
jm=jags.model("model.txt",data=data, inits=inits, n.chains=length(inits), n.adapt=3000)

update(jm,n.iter=5000)
z=coda.samples(jm,variable.names=parameters,n.iter=10000)
#check convergence
plot(z)
gelman.diag(z)
heidel.diag(z)
zj=jags.samples(jm,variable.names=c("psi.elev","elev.opt"),n.iter=10000)
#zdic=dic.samples(jm,n.iter=1000)

#output derived quantities

pdf(file="Bird plots.pdf", width=6, height=4)
par(mfrow=c(1,2))
psi = summary(zj$psi.elev,quantile,c(.025,.5,.975))$stat
plot(elev.x,psi[2,],typ="l",ylim=c(0,1),xlab="Elevation (m)", ylab="Probability of occupancy")
lines(elev.x,psi[1,],typ="l", lty="dashed")
lines(elev.x,psi[3,],typ="l", lty="dashed")
abline(v=1790)
text(1300,.9, "Optimum elevation", cex=.5)


#plot(density(zj$elev.opt), main="Optimum elevation at mean forest cover" ,xlab="Elevation", xlim=c(1500,2500))
#abline(v=summary(zj$elev.opt,quantile,c(.025))$stat)
hist(zj$elev.opt, breaks = 300,main="",xlab="Elevation (m)", xlim=c(1500,2500), freq=FALSE)
abline(v=summary(zj$elev.opt,quantile,c(.975))$stat, lty="dashed")
abline(v=summary(zj$elev.opt,quantile,c(.025))$stat, lty="dashed")
dev.off()