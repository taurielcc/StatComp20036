## -----------------------------------------------------------------------------
Y= c(10.98, 11.13, 12.51, 8.40, 9.27, 8.73, 6.36, 8.50,
7.82, 9.14, 8.24, 12.19, 11.88,9.57, 10.94, 9.58, 10.09,
8.11, 6.83, 8.88, 7.68, 8.47, 8.86, 10.36, 11.08)
X =c(35.3, 29.7, 30.8, 58.8, 61.4, 71.3, 74.4, 76.7, 70.7,
57.5, 46.4, 28.9, 28.1, 39.1, 46.8, 48.5, 59.3, 70.0,
70.0, 74.5, 72.1, 58.1, 44.6, 33.4, 28.6)
fit <- lm(Y~X)
coef <- round(fit$coef,3)
plot(X,Y)
abline(a=coef[1],b=coef[2],col='red',lwd=2)
text(55,10,paste0('Y=',coef[1],coef[2],'X'),cex=1.5)

## ----results = 'asis'---------------------------------------------------------
library(xtable)
print(xtable(head(mtcars)),type="html")

## -----------------------------------------------------------------------------
set.seed(303)
my.pareto <- function(a,b,n){
  u=runif(n)
  x=b*(1-u)^(-1/a)
  return(x)
}

x<-my.pareto(2,2,1000)#simulate a random sample from the Pareto(2, 2) distribution. 

#Graph the density histogram of the sample with the Pareto(2, 2) density superimposed for comparison.
hist(x, main = "Sample DensityHistogram vs Theoretical Density Curve for Pareto(2,2)", col = "green", prob =TRUE)
y<-seq(min(x),max(x),0.01)
lines(y,8/y^3,col="red")
  

## -----------------------------------------------------------------------------
set.seed(309)
my.epanechnikov <- function(n){
  x <- c()
  for(i in 1:n){
    u1 <- runif(1,-1,1)
    u2 <- runif(1,-1,1)
    u3 <- runif(1,-1,1)
    if(abs(u3) >= abs(u2) && abs(u3) >= abs(u1)) x[i] <- u2
    else x[i] <- u3
  }

  hist(x, main = "Sample Density Histogram vs Theoretical Density Curve", col = "green", prob =TRUE, xlim = c(-1, 1),ylim = c(0,1))
  
  domain <- seq(-1, 1, 0.01)
  f <- function(x){
    3/4 * (1 - x^2)
  }
  lines(domain, f(domain), col = "red") 
}

my.epanechnikov(1000)



## -----------------------------------------------------------------------------
empirical_sample<-function(r,beta,n){
  u=runif(n)
  y=beta*((1-u)^(-1/r)-1)
  return(y)
}
x<-empirical_sample(r=4,beta=2,n=1000)
hist(x, main = "Sample Density Histogram vs Theoretical Density Curve", col = "green", prob =TRUE)
y<-seq(min(x),max(x),0.01)
lines(y,64/(2+y)^5,col="red")


## -----------------------------------------------------------------------------
set.seed(501)
n <- 1e4
x <- runif(n, min=0, max=pi/3)
theta.hat <- mean(sin(x)) * pi/3
theta.true<- -cos(pi/3) - (-cos(0))
print(c(theta.hat,theta.true))


## -----------------------------------------------------------------------------
set.seed(507)
m = 10000
u1 = runif( as.integer(m/2) )
u.anti = 1 - u1
g.anti = ( exp(u1) + exp(u.anti) )/2
Theta.anti = mean( g.anti )
var.anti = var( g.anti )/(m/2)
u = runif(m)
Theta.MC = mean( exp(u) )
var.MC = var( exp(u) )/m
perc.reduc = 100*(var.MC - var.anti)/var.MC

#Theoretical value
Theo.MC<-((exp(1)^2-1)/2 - (exp(1)-1)^2)/m
Theo.anti<-((exp(1)^2-1)/2-2*(exp(1)-1)^2+exp(1))/m
Theo.reduc<-100*(Theo.MC - Theo.anti)/Theo.MC


## -----------------------------------------------------------------------------
cat("simple Monte Carlo method: ",Theta.MC, "\n", "antithetic variate approach: ",Theta.anti)

## -----------------------------------------------------------------------------
cat(" simple Monte Carlo method: ",var.MC, "\n", "antithetic variate approach: ", var.anti,"\n","the percent reduction in variance", perc.reduc,"%")

## -----------------------------------------------------------------------------
cat(" Theoretical value by simple MC: ",Theo.MC, "\n", "Theoretical value by antithetic variate: ",Theo.anti,"\n","the percent reduction in variance", Theo.reduc,"%")

## -----------------------------------------------------------------------------
m<-1e6
t<-rgamma(m,shape = 3/2,scale = 2)
theta.hat<-0.5*mean(t>1)
se<-sd(0.5*(t>1))
cat("theta.hat: ",theta.hat,"\n","se of g/f:",se)

## -----------------------------------------------------------------------------
m<-1e6
x<-rnorm(m)
theta.hat<-mean(x^2*(x>1))
se<-sd(x^2*(x>1))
cat("theta.hat: ",theta.hat,"\n","se of g/f:",se)

## -----------------------------------------------------------------------------
M=1e4 #number of replicates
k<-5 #number of stratum
r=M/k #replicates per stratum
N=50 #number of times to repeat the estimation
T2<-numeric(k)
est<-matrix(0,N,2)

g<-function(x){
  exp(-x-log(1+x^2))*(x>0)*(x<1)
}
f<-function(x){
  exp(-x)/(1-exp(-1))
}

for(i in 1:N){
  est[i,1]<-mean(g(runif(M)))
  for(j in 1:k){
    u<-runif(r,(j-1)/k,j/k)
    x=-log(1-u*(1-exp(-1)))
    T2[j]<-mean(g(x)/(k*f(x)))
      }
  est[i,2]<-sum(T2)
}

apply(est,2,mean)

apply(est,2,sd)


## -----------------------------------------------------------------------------
f_j<-function(x,j){
  exp(-x)/(exp(-(j-1)/5)-exp(-j/5))
}

theta.hat<-numeric(k)
se_j<-numeric(k)
for(j in 1:k){
  u<-runif(r,(j-1)/k,j/k)
  x=-log(exp(-(j-1)/5)-u*(exp(-(j-1)/5)-exp(-j/5)))
  theta.hat[j]<-mean(g(x)/f_j(x,j))
  se_j[j]<-sd(g(x)/f_j(x,j))
}

cat("theta.hat: ",sum(theta.hat),"\n","se:",sum(se_j))


## -----------------------------------------------------------------------------
m <- 1000
n <- 2000
mu <- 1
sigma <- 0.5
alpha <- 0.05
res <- replicate(m, expr={
  x <- rlnorm(n, meanlog=mu, sdlog=sigma)
  log.mean <- mean(log(x))
  log.se <- sd(log(x)) / sqrt(n)
  t.val <- qt(1 - alpha / 2, n - 1)
  
  CI <- c(log.mean, log.mean - (t.val * log.se), log.mean + (t.val * log.se))
  CI
})
cat("an empirical estimate of 95% confidence interval for the parameter mu:" ,"\n", "(", mean(res[2,]), mean(res[3,]),")")


## -----------------------------------------------------------------------------
UCL <- replicate(m, expr = {
  x <- rlnorm(n, meanlog=mu, sdlog=sigma)
  log.mean <- mean(log(x))
  log.se <- sd(log(x)) / sqrt(n)
  (log.mean -mu)/log.se
} )
CL<-mean(UCL > qt(alpha , n - 1))
CL

## -----------------------------------------------------------------------------
set.seed(123)
m<-1000
c<-numeric(m)
n<-20
for (i in 1:m) {
  x<-rchisq(n,2)
  alpha<-0.05
  u<-mean(x)+sqrt(var(x)/n)*qt(1-alpha/2,n-1)
  l<-mean(x)+sqrt(var(x)/n)*qt(alpha/2,n-1)
  if(2<=u&&2>=l) c[i]=1
}
cat("coverage probability(t-interval)=",sum(c)/m)
#The example using Chi-square interval
alpha <- 0.05
UCL <- replicate(m, expr = {
  x <- rchisq(n, df = 2)
  (n-1) * var(x) / qchisq(alpha, df = n-1)
} )
cat("\ncoverage probability(Chi-square interval)=",mean(UCL > 4))

## -----------------------------------------------------------------------------
level<-0.1 # the signifiance level
n <- 30   # sample size
m <- 2500  # reputation times

# statistics function
skewness <- function(x){
  # input data x
  # return statistics b1
  x_bar <- mean(x)
  numerator <- mean((x-x_bar)^3)
  denominator <- (mean((x-x_bar)^2))^1.5
  return(numerator / denominator)
}

# critical value for the skewness test at level=0.1
cv <- qnorm(1-level/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))

# replicate the procedure with different parameter alpha
pwr <- function(a){
  stat <- replicate(m, expr={
                    X <- rbeta(n, a, a)
                    skewness(X) })
  length(which(abs(stat) >= cv)) / m
}

# calculate power
alpha <- seq(0.5, 40, 0.5)
powers <- unlist(base::lapply(X=alpha, FUN=pwr))

# plot power vs alpha
plot(alpha, powers, type = "b", xlab = bquote(alpha), ylab='Power', ylim = c(0,1))
# plot test level line, horizon
abline(h = .1, lty = 3)
# add standard errors
se <- sqrt(powers * (1-powers) / m) 
lines(alpha, powers+se, lty = 3,col="blue")
lines(alpha, powers-se, lty = 3,col="blue")

## -----------------------------------------------------------------------------
df <- seq(1, 40) #degree of freedom for t-test
# replicate the procedure with different parameter alpha
pwr <- function(a){
  stat <- replicate(m, expr={
                    X <- rt(n, a)
                    skewness(X)})
  length(which(abs(stat) >= cv)) / m
}
powers <- unlist(base::lapply(X=df, FUN=pwr))
# plot power vs df
plot(df, powers, type = "b", xlab = 'Degree of Freedom', ylab='Power', ylim = c(0,1))
# plot test level line, horizon
abline(h = .1, lty = 3)
# add standard errors
se <- sqrt(powers * (1-powers) / m) 
lines(df, powers+se, lty = 3)
lines(df, powers-se, lty = 3)

## -----------------------------------------------------------------------------
alpha = 0.055
mu1 <- mu2 <- 0
sigma1 <- 1
sigma2 <- 1.5
sample_num <- c(15, 100, 1000)
m <- 1000
tests_F <- tests_CF <- numeric(length(sample_num))
testF <- function(x, y) {
  f_test <- var.test(x, y, alternative = "two.sided", conf.level = 1-alpha)
  return(f_test$p.value < alpha)
}
tests5 = function(x, y) {
  X<- x-mean(x)
  Y<- y-mean(y)
  out_x <- sum(X > max(Y)) + sum(X < min(Y))
  out_y <- sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(out_x, out_y)) > 5))
}

for (i in 1 : 3) {
  n <- sample_num[i]
  power_F <- mean(replicate(m, expr = {
    x <- rnorm(n, mu1, sigma1)
    y <- rnorm(n, mu2, sigma2)
    testF(x, y)
  }))
  power_CF <- mean(replicate(m, expr = {
    x <- rnorm(n, mu1, sigma1)
    y <- rnorm(n, mu2, sigma2)
    tests5(x, y)
  }))
  
  tests_F[i] <- power_F
  tests_CF[i] <- power_CF
}

#print the results
m<-matrix(c(tests_F, tests_CF),byrow=TRUE,nrow=2,
          dimnames = list(c("F test","Count five test"),c("size=15","size=100","size=1000")))
knitr::kable(m,align = "c",caption = "the power of the Count Five test vs F test")


## -----------------------------------------------------------------------------
mardia_multskew<-function(x,s=var(x)){
  x<-as.matrix(x)
  n<-nrow(x)
  #d<-ncol(x)

  #creat 1 matrix
  one<-matrix(1,n)
  
  #creat idendity matrix n*n
  id<-matrix(0,n,n)
  diag(id)<-1
  
  #creat Q matrix (I-1/n * 1[n] * 1[n]' )
  Q <- id - 1/n * one %*% t(one)

  # create g matrix 
  g <- Q %*% x %*% solve(s) %*% t(x) %*% Q
  
  b1d <- 1/(n^2) * sum(g^3)# b1d is used for skew measure
  k1<-n*b1d/6
  return(k1)
}


## -----------------------------------------------------------------------------
#repeat Examples 6.8
alpha<-0.05
n<-c(10,20,30,50,100) #sample size
d<-2
cv<-qchisq(1-alpha,df=d*(d+1)*(d+2)/6)#crit. values for chi-squared distribution

p.reject<-numeric(length(n))
m=1000

for (i in 1:length(n)) {
  skew<-numeric(m)
  for (j in 1:m) {
    data<-data.frame(x=rnorm(n[i]),y=rnorm(n[i]))
    #test decision is 1 (reject) or 0
    skew[j]<-as.integer(mardia_multskew(data)>=cv)
  }
  p.reject[i]<-mean(skew) #proportion rejected
}

p.reject


## -----------------------------------------------------------------------------
alpha<-0.1
n<-30
m=2500
eps<-c(seq(0,0.15,0.01),seq(0.15,1,0.05))
N<-length(eps)
pwr<-numeric(N)
cv<-qchisq(1-alpha,df=d*(d+1)*(d+2)/6)#crit. values for chi-squared distribution

for(j in 1:N){
  e<-eps[j]
  skew<-numeric(m)
  for (i in 1:m) {
    sigma<-sample(c(1,10),replace = TRUE,size = n,prob=c(1-e,e))
    data<-data.frame(rnorm(n,0,sigma),rnorm(n,0,sigma))
    skew[i]<-as.integer(mardia_multskew(data)>=cv)
  }
  pwr[j]<-mean(skew)
}


#plot power vs epsilon
plot(eps, pwr, type = "b",
xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(eps, pwr+se, lty = 3)
lines(eps, pwr-se, lty = 3)


## -----------------------------------------------------------------------------
library(bootstrap); attach(law)
n = length( LSAT )
R.hat = cor( LSAT,GPA )
R.jack = numeric( n )
for (i in 1:n) { R.jack[i] = cor( LSAT[-i],GPA[-i] ) }
bias.jack = (n-1)*( mean(R.jack) - R.hat )
R.bar = mean(R.jack)
se.jack = sqrt( (n-1)*mean( (R.jack-R.bar)^2 ) ) 
detach(law)
cat("jackknife estimate for bias: ",bias.jack,"\n","jackknife estimate for se: ",se.jack)

## -----------------------------------------------------------------------------
library(boot)
attach(aircondit)
x<-hours
detach(aircondit)
time.hat<-mean(x)  #MLE of 1/lambda
B=2000
set.seed(705)
time.boot<-function(x,i){
  mean(x[i])
}
boot.obj<-boot(data=x,statistic = time.boot,R=B)
print(boot.ci(boot.obj,type=c("norm","basic","perc","bca")))

## -----------------------------------------------------------------------------
data(scor,package = "bootstrap")
lamda<-eigen(cov(scor))$values
theta.hat<-lamda[1]/sum(lamda)

#Jackknife
n=nrow(scor)
theta.jack<-numeric(n)
for (i in 1:n) {
  lambda.jack<-eigen(cov(scor[-i,]))$values
  theta.jack[i]<-lambda.jack[1]/sum(lambda.jack)
}
bias.jack<-(n-1)*(mean(theta.jack)-theta.hat)
se.jack<-sqrt((n-1)*mean((theta.jack-mean(theta.jack))^2))

#print the results
cat("theta.hat: ",theta.hat,"\n","bias.jeck: ",bias.jack,"\n","se.jack: ",se.jack)

## -----------------------------------------------------------------------------
library(DAAG); 
attach(ironslag)
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n*(n-1)) # 'leave two out' has n(n-1) combinations
for (i in 1:n){
  for (j in i:n){
    if (i != j){
      y=magnetic[c(-i,-j)]
      x=chemical[c(-i,-j)]
      
      J1 <- lm(y ~ x)
      yhat11 <- J1$coef[1] + J1$coef[2] * chemical[i]
      yhat12 <- J1$coef[1] + J1$coef[2] * chemical[j]
      e1[(i-1)*n+j] <- sqrt((magnetic[i] - yhat11)^2+(magnetic[j] - yhat12)^2)
      
      J2 <- lm(y ~ x + I(x^2))
      yhat21 <- J2$coef[1] + J2$coef[2] * chemical[i] +
      J2$coef[3] * chemical[i]^2
      yhat22 <- J2$coef[1] + J2$coef[2] * chemical[j] +
      J2$coef[3] * chemical[j]^2
      e2[(i-1)*n+j] <- sqrt((magnetic[i] - yhat21)^2+(magnetic[j] - yhat22)^2)
      
      J3 <- lm(log(y) ~ x)
      logyhat31 <- J3$coef[1] + J3$coef[2] * chemical[i]
      logyhat32 <- J3$coef[1] + J3$coef[2] * chemical[j]
      yhat31 <- exp(logyhat31)
      yhat32 <- exp(logyhat32)
      e3[(i-1)*n+j] <- sqrt((magnetic[i] - yhat31)^2+(magnetic[j] - yhat32)^2)
      
      J4 <- lm(log(y) ~ log(x))
      logyhat41 <- J4$coef[1] + J4$coef[2] * log(chemical[i])
      logyhat42 <- J4$coef[1] + J4$coef[2] * log(chemical[j])
      yhat41 <- exp(logyhat41)
      yhat42 <- exp(logyhat42)
      e4[(i-1)*n+j] <- sqrt((magnetic[i] - yhat41)^2+(magnetic[j] - yhat42)^2)
    }
  }
}
detach(ironslag)
# estimates for prediction error
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}
set.seed(123)
# Count Five test permutation
count5test_permutation <- function(z) {
  n <- length(z)
  x <- z[1:(n/2)]
  y <- z[-(1:(n/2))]
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y)) 
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0) 
  return(as.integer(max(c(outx, outy)) > 5))
}
permutation <- function(z,R) {
  n <- length(z)
  out <- numeric(R)
  for (r in 1: R){
    p <- sample(1:n ,n ,replace = FALSE)
    out[r] <- count5test_permutation(z[p])
  }
  sum(out)/R
}              
n1 <- 20
n2 <- 50
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m <- 1e3
alphahat1 <- mean(replicate(m, expr={
  x <- rnorm(n1, mu1, sigma1)
  y <- rnorm(n2, mu2, sigma2)
  x <- x - mean(x) #centered by sample mean
  y <- y - mean(y)
  count5test(x, y)
}))
alphahat2 <- mean(replicate(m, expr={
  x <- rnorm(n1, mu1, sigma1)
  y <- rnorm(n2, mu2, sigma2)
  x <- x - mean(x) #centered by sample mean 
  y <- y - mean(y)
  z <- c(x,y)
  permutation(z,1000) 
})<0.05)
round(c(count5test=alphahat1,count5test_permutation=alphahat2),4)

## -----------------------------------------------------------------------------
library(RANN)
library(boot)
library(energy)
library(Ball)
library(survival)

m <- 30#times to loop
k<-3
p<-2#ncol
# mu <- 0.5
n1 <- n2 <- 50#nrow
R<-999#the number of replications in the bd.test function
n <- n1+n2
N = c(n1,n2)

Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]
  n2 <- sizes[2]
  n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0)
  z <- z[ix, ]
  NN <- nn2(data=z, k=k+1)
  block1 <- NN$nn.idx[1:n1,-1] 
  block2 <- NN$nn.idx[(n1+1):n,-1] 
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5) 
  (i1 + i2) / (k * n)
}


eqdist.nn <- function(z,sizes,k){#NN
  boot.obj <- boot(data=z,statistic=Tn,R=R,
  sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

p.values <- matrix(NA,m,3)#to store p values

for(i in 1:m) {
  x <- matrix(rnorm(n1 * p, sd = 1), ncol = p)#unequal variances
  y <- matrix(rnorm(n2 * p, sd = 1.4), ncol = p)
  z <- rbind(x, y)
  p.values[i, 1] <- eqdist.nn(z, N, k)$p.value#NN
  p.values[i, 2] <- eqdist.etest(z, sizes = N, R = R)$p.value#in the energy package
  p.values[i, 3] <- bd.test(x = x, y = y, R = 999, seed = i)$p.value#"Ball Divergence" in the ball package
}

alpha <- 0.1#confidence level
power <- apply(p.values<alpha,2,mean)#compute the number of p.values which is less than 0.1 in each column
pow <- data.frame(methods = c('NN','energy','Ball'),power)
print(pow)

## -----------------------------------------------------------------------------
for(i in 1:m) {
  x <- matrix(rnorm(n1 * p, mean = 0.4, sd = 1), ncol = p)#unequal variances and unequal expectations
  y <- matrix(rnorm(n2 * p, mean = 0, sd = 1.4), ncol = p)
  z <- rbind(x, y)
  p.values[i, 1] <- eqdist.nn(z, N, k)$p.value
  p.values[i, 2] <- eqdist.etest(z, sizes = N, R = R)$p.value
  p.values[i, 3] <- bd.test(x = x,  y = y,  R = 999,  seed = i)$p.value
}
alpha <- 0.1
power <- apply(p.values<alpha,2,mean)
pow <- data.frame(methods = c('NN','energy','Ball'),power)
print(pow)

## -----------------------------------------------------------------------------
for(i in 1:m) {
  x <- matrix(rt(n1 * p,df = 1), ncol = p)# t distribution
  y <- matrix(rnorm(n2 * p,sd = sample(c(1,1.3),size = n2*p, prob = c(0.5,0.5),replace = T)), ncol = p)#bimodel distribution
  z <- rbind(x, y)
  p.values[i, 1] <- eqdist.nn(z, N, k)$p.value
  p.values[i, 2] <- eqdist.etest(z, sizes = N, R = R)$p.value
  p.values[i, 3] <- bd.test(x = x, y = y, R = 999, seed = i)$p.value
}
alpha <- 0.01
power <- apply(p.values<alpha,2,mean)
pow <- data.frame(methods = c('NN','energy','Ball'),power)
print(pow)

## -----------------------------------------------------------------------------
n1 <- 50
n2 <- 5
n <- n1+n2
N = c(n1,n2)
for(i in 1:m) {
  x <- matrix(rnorm(n1*p,mean = 1), ncol = p)#100 samples
  y <- matrix(rnorm(n2*p,mean = 2), ncol = 2)#10 samples
  z <- rbind(x, y)
  p.values[i, 1] <- eqdist.nn(z, N, k)$p.value
  p.values[i, 2] <- eqdist.etest(z, sizes = N, R = R)$p.value
  p.values[i, 3] <- bd.test(x = x, y = y, R = 999, seed = i)$p.value
}
alpha <- 0.1
power <- apply(p.values<alpha,2,mean)
pow <- data.frame(methods = c('NN','energy','Ball'),power)
print(pow)

## -----------------------------------------------------------------------------
f_laplace<-function(x){
  0.5*exp(-abs(x))
}

rw.Metropolis<-function(sig,x0,m){
  x<-numeric(m)
  x[1]<-x0
  u<-runif(m)
  k<-0
  for(i in 2:m){
    y<-rnorm(1,x[i-1],sig)
      if(u[i]<=(f_laplace(y)/f_laplace(x[i-1])))
        x[i]<-y else{
          x[i]<-x[i-1]
          k<-k+1
        }
  }
  return(list(x=x,k=k))
}

## -----------------------------------------------------------------------------
m=2000
sig<-c(0.05,0.5,2,16)
x0<-25
rw1<-rw.Metropolis(sig[1],x0,m)
rw2<-rw.Metropolis(sig[2],x0,m)
rw3<-rw.Metropolis(sig[3],x0,m)
rw4<-rw.Metropolis(sig[4],x0,m)
accept<-1-c(rw1$k/m,rw2$k/m,rw3$k/m,rw4$k/m)  #acceptance rates
results<-rbind(sig,accept)
colnames(results)<-c("chain1","chain2","chain3","chain4")
print(results)

## -----------------------------------------------------------------------------
rw<-cbind(rw1$x,rw2$x,rw3$x,rw4$x)
for(i in 1:4){
  plot(rw[,i],type="l",xlab=bquote(sigma==.(round(sig[i],3))),ylab="X",ylim=range(rw[,i]))
}

## -----------------------------------------------------------------------------
Gelman_Rubin<-function(psi){
  psi<-as.matrix(psi) #psi[i,j] is the statistic psi(X[i,1:j]) for the chain in i-th row of X
  n<-ncol(psi)
  k<-nrow(psi)
  psi.mean<-rowMeans(psi)
  B<-n*var(psi.mean) #between variance est
  psi.w<-apply(psi,1,"var") #within variance
  W<-mean(psi.w) #within est
  v.hat<-W*(n-1)/n+(B/n) #upper variance est
  r.hat<-v.hat/W #G-R statistic
  return(r.hat)
}

sigma <- 2 #parameter of proposal distribution
n=15000 #length of chains
k <- 4 #number of chains to generate
b <- 1000 #burn-in length
x0<-c(-15,-5,15,25)


my_chains<-function(sig,x0,m){
  x<-numeric(m)
  x[1]<-x0
  u<-runif(m)
  for(i in 2:m){
    y<-rnorm(1,x[i-1],sig)
      if(u[i]<=(f_laplace(y)/f_laplace(x[i-1])))
        x[i]<-y else{
          x[i]<-x[i-1]
        }
  }
  return(x)
}

#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
X[i, ] <- my_chains(sigma,x0[i],n)

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman_Rubin(psi))

#plot psi for the four chains
for (i in 1:k)
plot(psi[i, (b+1):n], type="l",xlab=i, ylab=bquote(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
  rhat[j] <- Gelman_Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.1, lty=2)



## -----------------------------------------------------------------------------
k<-c(4:25,100,500,1000)

object<-function(a,df){
  a2<-a^2
  arg<-sqrt(a2*df/(df+1-a2))
  sk<-pt(q=arg,df=df,lower=F)
  
  arg=sqrt(a2*(df-1)/(df-a2))
  skm1<-pt(q=arg,df=df-1,lower=F)
  
  return(sk-skm1)
}

for(i in 1:length(k)){
  print(c(as.integer(k[i]),uniroot(object,lower=1,upper=2,df=k[i])$root))
  
}


## -----------------------------------------------------------------------------
blood<-function(p,q,na,nb,nab,noo,N){
  n<-sum(na,nb,nab,noo)
  pt<-qt<-rt<-l_obs<-numeric(0)
  r<-1-p-q
  l_obs[1]<-na*log(p^2+2*p*r)+nb*log(q^2+2*q*r)+noo*log(r^2)+nab*log(2*p*q)
  
  pt[1]<-p
  qt[1]<-q
  rt[1]<-1-p-q
  for(i in 1:N){
    #E-step
    p.old<-pt[i]
    q.old<-qt[i]
    r.old<-rt[i]
    
    naa_t<-na*p.old^2/(p.old^2+2*p.old*r.old)
    nao_t<-na*2*p.old*r.old/(p.old^2+2*p.old*r.old)
    nbb_t<-nb*q.old^2/(q.old^2+2*q.old*r.old)
    nbo_t<-nb*2*q.old*r.old/(q.old^2+2*q.old*r.old)
    
    #M-step: update p,q,r
    pt[i+1]<-(2*naa_t+nao_t+nab)/(2*n)
    qt[i+1]<-(2*nbb_t+nbo_t+nab)/(2*n)
    rt[i+1]<-(2*noo+nao_t+nbo_t)/(2*n)
  
    # calculate the corresponding log-maximum likelihood values (for observed data)
    l_obs[i+1]<-na*log(pt[i+1]^2+2*pt[i+1]*rt[i+1])+nb*log(qt[i+1]^2+2*qt[i+1]*rt[i+1])+noo*log(rt[i+1]^2)+nab*log(2*pt[i+1]*qt[i+1])
    
    
    if(abs(pt[i+1]-p.old)/p.old < tol) break
  }
  return(list(pt=pt,qt=qt,rt=rt,l_obs=l_obs))
}


tol <- .Machine$double.eps^0.5
na<-444;nb<-132;noo<-361;nab<-63
N <- 1000 #max. number of iterations
p<-1/3;q<-1/3  #initial est. for p,q
my_EM_results<-blood(p,q,na,nb,nab,noo,N)
p_t<-my_EM_results$pt;q_t<-my_EM_results$qt;
r_t<-my_EM_results$rt;l_obs<-my_EM_results$l_obs
print(cbind(p_t,q_t,r_t,l_obs))  
plot(1:length(l_obs), l_obs, "b", xlab="Iter. number", main="Log-likelihood for observed data")


## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

## -----------------------------------------------------------------------------
set.seed(1)
attach(mtcars)

## -----------------------------------------------------------------------------
# loop
for (x in formulas){
   print(lm(x))
}

## -----------------------------------------------------------------------------
# lapply
lapply(formulas, lm)
detach()

## -----------------------------------------------------------------------------
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)

## -----------------------------------------------------------------------------
sapply(trials,function(x){
  return(round(x$p.value,4))
})

## -----------------------------------------------------------------------------
##get rid of the anonymous function
round(sapply(trials, "[[", "p.value"),4)

## -----------------------------------------------------------------------------
library(parallel)
mcvMap <- function(f, FUN.VALUE , ...) {
    out <- mcMap(f, ...)
    vapply(out, identity, FUN.VALUE)
}

## -----------------------------------------------------------------------------
data <- matrix(rnorm(20, 0, 10), nrow = 4)
x <- as.data.frame(data)
answer1 <- Map("/",x,vapply(x,mean,c(1)))
answer2 <- lapply(x,function(data){data/(mean(data))})
print(data.frame(unlist(answer1),unlist(answer2)))

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)

set.seed(123)
#define cpp function
cppFunction('List rwMC(double sigma, double x0, int N) {
     NumericVector x(N);
     x[0]=x0;
     auto u = runif(N);
     int k = 0;
     for (int i = 1;i < N; i++) {
         auto y = rnorm(1, x[i - 1], sigma)[0];
         if (u[i] <= (exp(-abs(y))/exp(-abs(x[i-1])))) {
             x(i) = y;
         } else {
             x[i] = x[i-1];
             k++;
         }
     }
     return List::create(Named("x") = x, Named("k") = k);
 }')


#define R function
rwMR <- function(sigma, x0, N){
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (exp(-abs(y)) / exp(-abs(x[i-1])))) x[i] = y 
    else {
      x[i] = x[i-1]
      k = k+1
    }
  }
  return(list(x = x, k = k))
}

## -----------------------------------------------------------------------------
out1<-rwMC(0.5,25,2000)
out2<-rwMR(0.5,25,2000)
par(mfrow=c(1,2))
plot(out1$x, type="l", main='rwMC(sigma=0.5)', ylab="x",xlab = "")
plot(out2$x,type="l", main='rwMR(sigma=0.5)', ylab="x",xlab = "")
par(mfrow=c(1,1))
#compare the generated random numbers
qqplot(out1$x,out2$x,xlab = 'rwMR',ylab = 'rwMC',main='qqplot of two functions(sigma=0.5)')
abline(a=0,b=1)
#compare the computation time
ts <- microbenchmark(r=rwMR(0.5,25,2000),cpp=rwMC(0.5,25,2000))
summary(ts)

## -----------------------------------------------------------------------------
out1<-rwMC(2,25,2000)
out2<-rwMR(2,25,2000)
par(mfrow=c(1,2))
plot(out1$x, type="l", main='rwMC(sigma=2)', ylab="x",xlab = "")
plot(out2$x,type="l", main='rwMR(sigma=2)', ylab="x",xlab = "")
par(mfrow=c(1,1))
#compare the generated random numbers
qqplot(out1$x,out2$x,xlab = 'rwMR',ylab = 'rwMC',main='qqplot of two functions(sigma=2)')
abline(a=0,b=1)
#compare the computation time
ts <- microbenchmark(r=rwMR(2,25,2000),cpp=rwMC(2,25,2000))
summary(ts)

## -----------------------------------------------------------------------------
out1<-rwMC(10,25,2000)
out2<-rwMR(10,25,2000)
par(mfrow=c(1,2))
plot(out1$x, type="l", main='rwMC(sigma=10)', ylab="x",xlab = "")
plot(out2$x,type="l", main='rwMR(sigma=10)', ylab="x",xlab = "")
par(mfrow=c(1,1))
#compare the generated random numbers
qqplot(out1$x,out2$x,xlab = 'rwMR',ylab = 'rwMC',main='qqplot of two functions(sigma=10)')
abline(a=0,b=1)
#compare the computation time
ts <- microbenchmark(r=rwMR(10,25,2000),cpp=rwMC(10,25,2000))
summary(ts)

## -----------------------------------------------------------------------------
out1<-rwMC(16,25,2000)
out2<-rwMR(16,25,2000)
par(mfrow=c(1,2))
plot(out1$x, type="l", main='rwMC(sigma=16)', ylab="x",xlab = "")
plot(out2$x,type="l", main='rwMR(sigma=16)', ylab="x",xlab = "")
par(mfrow=c(1,1))
#compare the generated random numbers
qqplot(out1$x,out2$x,xlab = 'rwMR',ylab = 'rwMC',main='qqplot of two functions(sigma=16)')
abline(a=0,b=1)
#compare the computation time
ts <- microbenchmark(r=rwMR(16,25,2000),cpp=rwMC(16,25,2000))
summary(ts)

