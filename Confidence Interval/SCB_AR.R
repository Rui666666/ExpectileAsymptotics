library(expectreg)
library(VGAM)


alpha = seq(from = 0.1, to =0.9, by=0.1)

simulate_ar1 = function(n, phi,family="t",sd=1, df=5) {
  Y = numeric(n)
  if(family=="t")
    for (t in 2:n) 
      Y[t] = phi*Y[t-1] + rt(1,df=df)
  else
    for(t in 2:n)
      Y[t] = phi*Y[t-1] + rnorm(1,sd=sd)
  return(Y)
}

expectile_truth = function(phi=0,family="t",df=5,sd=1,alpha=seq(0.1,0.9,0.1)){
  n = 1000
  iter = 1000
  expectile_hat = matrix(0,ncol=9,nrow=iter)
  for(i in 1:iter){
    Data = simulate_ar1(n = n,family=family,df=df,sd=sd,phi=phi)
    expectile_hat[i,] = expectile(Data, probs=alpha)
  }
  return(apply(expectile_hat, 2, mean))
}

expectile_SCB = function(Data,b,alpha){
  theta = expectile(Data,probs = alpha)
  n = length(Data)
  theta_hat = matrix(0, nrow= n-b+1,ncol=length(alpha))
  for(i in 1:(n-b+1)){
    theta_hat[i,] = expectile(Data[i:(i+b-1)],probs=alpha)
  }
  pivot = sqrt(b)*abs(theta_hat - rep(1,n-b+1)%*%t(theta))
  max_pivot = apply(pivot,1,max)
  return(quantile(max_pivot,0.95))
}

Test_SCB = function(Data, b_scale=1/3, b = NA,
                    alpha = seq(0.1,0.9,by=0.1), truth){
  n = length(Data)
  expectile_hat = expectile(Data, probs=alpha)
  if(is.na(b)) b = ceiling(n^b_scale)
  q = expectile_SCB(Data,b,alpha)
  return(all(sqrt(n)*abs(expectile_hat-truth)<q))
}

coverage = function(n=500, phi=0, family="t",sd=1,df=5,b=22,
                    b_scale=1/3,iter=10000,alpha=seq(0.1,0.9,0.1)){
  start_time = Sys.time()
  n = 500
  phi = 0.3
  cover = numeric(iter)
  coverage_rate_SCB = 0
  truth = expectile_truth(phi=phi,family=family,sd=sd,df=df,alpha=alpha)
  for(i in 1:iter){
    Data = simulate_ar1(n=n,family=family,sd=sd,df=df,phi=phi)
    cover[i] = Test_SCB(Data,b_scale=b_scale,b=b,truth=truth,alpha=alpha)
  }
  coverage_rate_SCB = mean(cover)
  end_time = Sys.time()
  print(end_time-start_time)
  return(coverage_rate_SCB)
}
coverage(iter=100)

# N(0,1) phi = 0,-.1,-.5,-.8. CR .9426,.9464,.9358,.8276
# t5 phi = 0,.1,.2,.3,.4,.5. CR .9383,.9345, 0.9327,0.9296, 0.9246, 0.9146
# t5 phi = -.1,-.2,-.3,-.4,-.5. CR  0.9384,0.9403,0.9419,0.9348,0.9273

# plot
alpha=seq(0.1,0.9,0.01)
qs=expectile_SCB(Data,b=22,alpha=alpha)
qp=c()
for(j in 1:81){
  mu = expectile(Data,probs=alpha[j])
  n = length(Data)
  hat_mu = numeric(n-b+1)
  for(i in 1:(n-b+1)){
    hat_mu[i] = expectile(Data[i:(i+b-1)],probs=alpha[j])
  }
  qp[j] = quantile(abs(sqrt(b)*(hat_mu-mu)),0.95)
}

exp_hat = expectile(Data,probs=alpha)

par(mfrow=c(1,1),mar=c(2.5,2,1,1))
plot(x= alpha,y = exp_hat, bty="l",xlab="",ylab="",type="l",ylim=c(-1.5,1))
lines(x=alpha, y=exp_hat-qs/sqrt(n),col="blue")
lines(x=alpha, y = exp_hat -qp/sqrt(n),col="orange")
lines(x=alpha, y=exp_hat+qs/sqrt(n),col="blue")
lines(x=alpha, y = exp_hat +qp/sqrt(n),col="orange")



