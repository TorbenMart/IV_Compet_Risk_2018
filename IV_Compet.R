library(Rcpp)


sourceCpp("est_compet.cpp") # Rcpp function (to speed things up) that does the estimation. Make sure R can access the file.

simny.Xcont.Gbin=function(n,gamG=1){ # Function to simulate data; Sim setup 1 in paper.
G=rbinom(n,1,0.5)
Sigma=matrix(c(0.25,0.15,0.15,0.25),2,2);mu=c(0,0) # Evt sæt 0.15 til -0.16667, så vil Aalen est. X's' effekt til 0.
Sigma=matrix(c(0.25,-0.16667,-0.16667,0.25),2,2);mu=c(0,0) 
library(MASS)
tmp=mvrnorm(n , mu, Sigma)
X=c(tmp[,1])+0.5+gamG*G
U=c(tmp[,2])+1.5
# gamU=0.15/0.25;-0.16667/0.25
# summary(lm(X~U+G))
ud=list()	
ud$X=X	
ud$G=G
ud$U=U	
ud$theta.true=0.3
return(ud)
}



estf=function(nsim){ # Function that does the estimation outlined in paper.

#################    Simulating data   ##########

n=nsim


del=0.01;del.indik=1

simvar=simny.Xcont.Gbin(n,0.32) #Sim setup 1 in paper, here correlation of 0.3



G=simvar$G
X=simvar$X
U=simvar$U
summary(X);summary(U)
theta.true=simvar$theta.true


res.anova=lm(X~G)
F.stat=summary(res.anova)$fstatistic[1];#print(c("F.stat",F.stat))
#print(F.stat)
beta.hat=1

alpha0=0.1; alphaG=0; alphaX=0; alphaU=0.08
beta0=0.07; betaG=0; betaX=0.1; betaU=0.1

#The two cause specific hazards
lamb1=(alpha0+alphaG*G+alphaX*X+alphaU*U) 
lamb2=(beta0+betaG*G+betaX*X+betaU*U) 
min.intens1=min(c(lamb1))
min.intens2=min(c(lamb2))
doed=rexp(n)/(lamb1+lamb2) 
summary(doed)
indik.neg=as.numeric(summary(doed)[1]<0)

eps=rbinom(n,1,lamb2/(lamb1+lamb2))+1 # Which of the causes are we seeing.

# Inducing censoring
#
cen1<-runif(n,0,3.5)
cen.tmp=rbinom(n,1,0.2)
cen1=cen1*cen.tmp
cen1[cen1==0]=10
cen.t=3.5
cen2<-rep(cen.t,n)
cen<-apply(cbind(cen1,cen2),1,min);	
time=apply(cbind(doed,cen),1,min)
status=as.numeric(cen>=doed)
stime=sort(time[status==1])
k=length(stime)
stime1=sort(time)  # All ordered time-points, also censored.
stime1.tmp=sort(time)  # All ordered time-points, also censored.

#print(cor.test(G,X,method="spearman"))
#print(c('Percent censorings:',1-sum(status/n)))


status[status>0]=eps[status>0]
status1=status*0
status1[status==1]=1
status2=status*0
status2[status==2]=1
stime=sort(time[status>0])
stime1=sort(time[status==1]);k1=length(stime1)
stime2=sort(time[status==2]);k2=length(stime2)
res=B_est_compet(time, status1,  status2, stime, G, X)

### Calculating variances
dB1=dB2=numeric(k)
dB1[1]=res$B1[1];dB2[1]=res$B2[1]
dB1[2:k]=res$B1[2:k]-res$B1[1:(k-1)];dB2[2:k]=res$B2[2:k]-res$B2[1:(k-1)]
dB.dot=B.dot=matrix(0,k,2)


G.c=G-mean(G)
I=matrix(c(1,0,0,1),2,2);b=matrix(c(1,1),2,1)
F_0t=list()
d.epsB.tmp=epsB.tmp=list()
d.epsB=epsB=list()
epsB.thet=list()
Sigma=list()
se=matrix(0,k,2)
V=V.tmp=dV=matrix(0,k,2)
status.t=as.numeric(time>=stime[1])
dN.tmp=as.numeric(time==stime[1])
dN1=matrix(status1*dN.tmp,n,1)
dN2=matrix(status2*dN.tmp,n,1)
H=status.t*G.c/sum(status.t*G.c*X)
H.deriv=status.t*G.c*(X*sum(status.t*G.c*X)-sum(status.t*G.c*X*X))/sum(status.t*G.c*X)^2
F_0t[[1]]=I+b%*%matrix(H.deriv,1,n)%*%cbind(dN1,dN2)
V.tmp[1,]=dV[1,]=matrix(H,1,n)%*%(cbind(dN1,dN2)-cbind(dB1[1]*X,dB2[1]*X) )%*%solve(F_0t[[1]])
mat.tmp=(cbind(dN1,dN2)-cbind(dB1[1]*X,dB2[1]*X) )%*%solve(F_0t[[1]])
d.epsB.tmp[[1]]=epsB.tmp[[1]]=epsB[[1]]=mat.tmp*c(H)
T1=sum(status.t*G.c*dN1);T1.dot=-sum(status.t*dN1)
N1=sum(status.t*G.c*X);N1.dot=-sum(status.t*X)
T2=sum(status.t*G.c*dN2);T2.dot=-sum(status.t*c(dN2))
N2=sum(status.t*G.c*X);N2.dot=-sum(status.t*X)
B.dot[1,]=dB.dot[1,]=c( (T1.dot*N1-T1*N1.dot)/N1^2,  (T2.dot*N2-T2*N2.dot)/N2^2  )
epsB.thet[[1]]=epsB[[1]]+cbind(G.c*B.dot[1,1],G.c*B.dot[1,2])/n
Sigma[[1]]=t(epsB.thet[[1]])%*%epsB.thet[[1]]
se[1,]=sqrt(diag(Sigma[[1]]))

for (j in 2:k){
status.t=as.numeric(time>=stime[j])
dN.tmp=as.numeric(time==stime[j])
dN1=matrix(status1*dN.tmp,n,1)
dN2=matrix(status2*dN.tmp,n,1)
H=status.t*G.c*exp(X*(res$B1[j-1]+res$B2[j-1]))/sum(status.t*G.c*X*exp(X*(res$B1[j-1]+res$B2[j-1])))
H.deriv=status.t*G.c*exp(X*(res$B1[j-1]+res$B2[j-1]))*(X*sum(status.t*G.c*X*exp(X*(res$B1[j-1]+res$B2[j-1])))-
sum(status.t*G.c*exp(X*(res$B1[j-1]+res$B2[j-1]))*X*X))/sum(status.t*G.c*X*exp(X*(res$B1[j-1]+res$B2[j-1])))^2
F_0t[[j]]=F_0t[[j-1]]%*%(I+b%*%matrix(H.deriv,1,n)%*%cbind(dN1,dN2))
dV[j,]=matrix(H,1,n)%*%(cbind(dN1,dN2)-cbind(dB1[j]*X,dB2[j]*X) )%*%F_0t[[j-1]]
dmat.tmp=(cbind(dN1,dN2)-cbind(dB1[j]*X,dB2[j]*X) )%*%solve(F_0t[[j]])
d.epsB.tmp[[j]]=dmat.tmp*H
epsB.tmp[[j]]=epsB.tmp[[j-1]]+d.epsB.tmp[[j]]
epsB[[j]]=epsB.tmp[[j]]%*%F_0t[[j]]
V.tmp[j,]=V.tmp[j-1,]+dV[j,]
V[j,]=c(matrix(V.tmp,1,2)%*%solve(F_0t[[j]]))
tmp1=exp(X*(res$B1[j-1]+res$B2[j-1]))
T1=sum(status.t*G.c*tmp1*dN1);
T1.dot=sum(status.t*tmp1*(G.c*X*(B.dot[j-1,1]+B.dot[j-1,2]) - 1)*dN1)
N1=sum(status.t*G.c*tmp1*X);
N1.dot=sum(status.t*X*tmp1*(G.c*X*(B.dot[j-1,1]+B.dot[j-1,2]) - 1) )
T2=sum(status.t*G.c*tmp1*dN2);
T2.dot=sum(status.t*tmp1*(G.c*X*(B.dot[j-1,1]+B.dot[j-1,2]) - 1)*dN2)
N2=sum(status.t*G.c*tmp1*X);
N2.dot=sum(status.t*X*tmp1*(G.c*X*(B.dot[j-1,1]+B.dot[j-1,2]) - 1) )
dB.dot[j,]=c( (T1.dot*N1-T1*N1.dot)/N1^2,  (T2.dot*N2-T2*N2.dot)/N2^2  )
B.dot[j,]=B.dot[j-1,]+dB.dot[j,]
epsB.thet[[j]]=epsB[[j]]+cbind(G.c*B.dot[j,1],G.c*B.dot[j,2])/n
Sigma[[j]]=t(epsB.thet[[j]])%*%epsB.thet[[j]]
se[j,]=sqrt(diag(Sigma[[j]]))

}

### End Calculating variances


ud=list()
ud$true.coef=c(alphaX,betaX)
ud$stime=stime
ud$B1=res$B1
ud$B2=res$B2
ud$se.B1=se[,1]
ud$se.B2=se[,2]

return(ud)

}

res=estf(3200)
names(res)

# res contains the following:
# true.coef: Data simulated from constant effects model; contains the true values.
# stime: ordered event times (both causes)
# B1: Estimated B1(t), see paper, display (2.5); same length as stime.
# B2: Estimated B2(t), see paper, display (2.5); same length as stime.
# se.B1: Estimated s.e. of \hat B1 at all event times; same length as stime.   
# se.B2: Estimated s.e. of \hat B2 at all event times; same length as stime.


#We now plot the estimated cumulative coefficients along with 95% point wise confidence bands. Both causes.

min1=min(res$B1-1.96*res$se.B1)
max1=max(res$B1+1.96*res$se.B1)
min2=min(res$B2-1.96*res$se.B2)
max2=max(res$B2+1.96*res$se.B2)

par(mfrow=c(1,2))
plot(res$stime,res$B1,type="s",xlab="Time",ylab="Estimated cumulative regression coefficient",main="Cause 1",ylim=c(min1,max1));
lines(res$stime,res$true.coef[1]*res$stime)
lines(res$stime,res$B1-1.96*res$se.B1,type="s",lty=2)
lines(res$stime,res$B1+1.96*res$se.B1,type="s",lty=2)
plot(res$stime,res$B2,type="s",xlab="Time",ylab="Estimated cumulative regression coefficient",main="Cause 2",ylim=c(min2,max2))
lines(res$stime,res$true.coef[2]*res$stime)
lines(res$stime,res$B2-1.96*res$se.B2,type="s",lty=2)
lines(res$stime,res$B2+1.96*res$se.B2,type="s",lty=2)
