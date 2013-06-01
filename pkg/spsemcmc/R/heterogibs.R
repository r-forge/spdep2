#Heteroscatastic Gibbs sampler
#set no of obs and variables
heterogibs <-function(){
n<-100
k<-3
#generate the data set
x <-matrix(rnorm(n*k),n,k)
b <-matrix(rep(1,k),k,1)
tt <-matrix(rep(1,n),n,1)
tt[51:100,1] <-1:50
y<-x%*%b+rnorm(n)*sqrt(tt)
#set no of draws, omits, create output variable
ndraw <-1100
nomit <-100
bsave <-matrix(rep(0,ndraw*k),ndraw,k)
ssave <-matrix(rep(0,ndraw*1),ndraw,1)
vsave <-matrix(rep(0,ndraw*n),ndraw,n)
#prior for b mean and variances
c <-matrix(c(1,1,1),3,1)
R <-diag(rep(1,k))
T <-diag(rep(1,k))
Q <-T
q <-Q%*%c
b0 <-qr.solve(x,y)
#setting up initial values
sige <-t((y-x%*%b0))%*%(y-x%*%b0)/(n-k)
V <-matrix(rep(1,n),n,1)
In <-matrix(rep(1,n),n,1)
rval <-4
qpq <-t(Q)%*%Q
qpv <-t(Q)%*%q
for (i in 1:ndraw) {
ys <-y*sqrt(V)
xs <-x*matrix(rep(sqrt(V),3),100,3)
xpxi <-solve(t(xs)%*%xs+as.numeric(sige)*qpq)
b <-xpxi%*%(t(xs)%*%ys+as.numeric(sige)*qpv)
b <-mvrnorm(1,b,as.numeric(sige)*xpxi)
bsave[i,] <-t(b)
e <-ys - xs%*%b
ssr <-t(e)%*%e
chi <-rchisq(1,df=n)
sige <-ssr/chi
ssave[i,] <-sige
chiv <-rchisq(n,rval+1)
vi <-(e*e/as.numeric(sige))+In*rval/chiv
V <-In/vi
vsave[i,] <-t(vi)
}
##computation of means and sd and other statistic
bhat <-mean(bsave[(nomit+1):ndraw,])
bstd <-sd(bsave[(nomit+1):ndraw,])
tstat <-bhat/bstd
smean <-mean(ssave[(nomit+1):ndraw,])
vmean <-mean(vsave[(nomit+1):ndraw,])
tout <-qt(t(tstat),n)
##output
output <-cbind(bhat,tstat,tout,smean,vmean)
output
}




