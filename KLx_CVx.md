```R
KL_1 <- function(x) 1-log(x^(2))-log(gamma(1/x^(2)))+((digamma(1/x^(2))-1)/x^(2))-digamma(1/x^(2)) #Gamma Distribution
```


```R
KL_3 <- function(x) (1/2)*(log((x^(2)+1)/log(x^(2)+1))+log(exp(1)/(2*pi))) #Lognormal Distribution
```


```R
KL_4 <- function(x) -log(x) #Shifted Exponential Distribution
```


```R
f2 <- function(y, x)  sqrt(1/(2*pi*x^(2)*y^(3)))*exp(-(1/(2*x^(2)))*(((y-1)^(2))/y)) #Inverse Gaussian distribution
```


```R
logf2 <- function(y,x) -((1/2)*log(2*pi*x^(2)))-((3/2)*log(y))-(((y-1)^(2))/(2*(x^(2))*y)) 
#Log of the pdf of the Inverse Gaussian Distribution
```


```R
f3 <- function(y,x) f2(y,x)*logf2(y,x) #Integrand for the entropy function
```


```R
f2int <- function(x) -(integrate(f3, lower=0, upper=Inf, x=x)$value) #Entropy of the Inverse Gaussian function
```


```R
v.h <- Vectorize(f2int) #Vectorize with respect to CV
```


```R
KL_2 <- function(x) 1-v.h(x) #Formula for the KL distance from the exponential distribution with mean =1 
```


```R
x <- seq(0.1,2.5,0.01) #Values of CV-range
```


```R
plot(x, KL_1(x), yaxs='i', xaxs='i', xlim=c(0, 2.5), ylim=c(0,2), ylab="KL distance (ISI distributions)", xlab="CV(ISI)", type="l", col="red")
lines(x, KL_2(x), col="blue", type="l", lty=2)
lines(x, KL_3(x), col="forestgreen", type="l", lty=3)
lines(x, KL_4(x), col="brown", type="l", lty=4)
legend(1.2, 1.2, legend=c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted Exp."),
       col=c("red", "blue","forestgreen", "brown"), lty=c(1,2,3,4), cex=0.65)
```


![png](output_10_0.png)


## From the AIFR perspective:


```R
y <- function(r) (r)/sqrt(1+r^(2))
```


```R
KL_1 <- function(r) 1-log((y(r))^(2))-log(gamma(1/(y(r))^(2)))+((digamma(1/(y(r))^(2))-1)/(y(r))^(2))-digamma(1/(y(r))^(2)) #Gamma Distribution
```


```R
y1 <- function(r) r
```


```R
KL_3 <- function(r) (1/2)*(log(((y1(r))^(2)+1)/log((y1(r))^(2)+1))+log(exp(1)/(2*pi))) #Lognormal Distribution
```


```R
y2 <- function(r) r
```


```R
f2 <- function(z, r)  sqrt(1/(2*pi*(y2(r))^(2)*z^(3)))*exp(-(1/(2*(y2(r))^(2)))*(((z-1)^(2))/z)) #Inverse Gaussian distribution
```


```R
logf2 <- function(z,r) -((1/2)*log(2*pi*(y2(r))^(2)))-((3/2)*log(z))-(((z-1)^(2))/(2*((y2(r))^(2))*z)) 
#Log of the pdf of the Inverse Gaussian Distribution
```


```R
f3 <- function(z,r) f2(z,r)*logf2(z,r) #Integrand for the entropy function
```


```R
f2int <- function(r) -(integrate(f3, lower=0, upper=Inf, r=r)$value) #Entropy of the Inverse Gaussian function
```


```R
v.h <- Vectorize(f2int) #Vectorize with respect to CV|
```


```R
KL_2 <- function(r) 1-v.h(r) #Formula for the KL distance from the exponential distribution with mean =1 
```


```R
library("expint")
```


```R
dat <- data.frame(CV_X= seq(0.227,2.5,0.001))
```


```R
CV_R4 <- function(CV_X){
    sapply(CV_X, function(x) if (x<1) sqrt(((exp((1/x)-1)*expint_E1((1/x)-1, scale = FALSE))/x)-1) else NaN)
}
```


```R
KL_4 <- function(x) -log(x) #Shifted Exponential Distribution
```


```R
dat$CV_R4 <- CV_R4(dat$CV_X)
```


```R
dat$KL_4 <- KL_4(dat$CV_X)
```


```R
r <- seq(0.2,2.5,0.001)  
```


```R
plot(r, KL_1(r), xaxs='i', yaxs='i',xlim=c(0,2.5), ylim=c(0,2), ylab="KL distance (ISI distributions)", xlab="CV(R)", type="l", col="red")
lines(r, KL_3(r), type="l", lty=2, col="blue")
lines(r, KL_2(r), col="green4", lwd=2.5, lty=3)
lines(dat$CV_R4, dat$KL_4, type="l", lty=5, col="brown")
legend(1.2, 1.8, legend=c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted Exp."),
       col=c("red", "blue","green4", "brown"), lty=c(1,2,3,5), cex=0.45)
```


![png](output_30_0.png)



```R

```
