```R
KL_1 <- function(x) 1-log(x^(2))-log(gamma(1/x^(2)))+((digamma(1/x^(2))-1)/x^(2))-digamma(1/x^(2)) #Gamma Distribution
```


```R
KL_2 <- function(x) (1/2)*(log((x^(2)+1)/log(x^(2)+1))+log(exp(1)/(2*pi))) #Lognormal Distribution
```


```R
KL_4 <- function(x) -log(x) #Shifted Exponential Distribution
```


```R
f3 <- function(y, x)  sqrt(1/(2*pi*x^(2)*y^(3)))*exp(-(1/(2*x^(2)))*(((y-1)^(2))/y)) #Inverse Gaussian distribution
```


```R
logf3 <- function(y,x) -((1/2)*log(2*pi*x^(2)))-((3/2)*log(y))-(((y-1)^(2))/(2*(x^(2))*y)) 
#Log of the pdf of the Inverse Gaussian Distribution
```


```R
f4 <- function(y,x) f3(y,x)*logf3(y,x) #Integrand for the entropy function
```


```R
f2int <- function(x) -(integrate(f4, lower=0, upper=Inf, x=x)$value) #Entropy of the Inverse Gaussian function
```


```R
v.h <- Vectorize(f2int) #Vectorize with respect to CV
```


```R
KL_3 <- function(x) 1-v.h(x) #Formula for the KL distance from the exponential distribution with mean =1 
```


```R
x <- seq(0.1,2.5,0.01) #Values of CV-range
```


```R

```

## From the AIFR perspective:


```R
yr <- function(r) (r)/sqrt(1+r^(2))
```


```R
KL_1r <- function(r) 1-log((yr(r))^(2))-log(gamma(1/(yr(r))^(2)))+((digamma(1/(yr(r))^(2))-1)/(yr(r))^(2))-digamma(1/(yr(r))^(2)) #Gamma Distribution
```


```R
y1r <- function(r) r
```


```R
KL_3r <- function(r) (1/2)*(log(((y1r(r))^(2)+1)/log((y1r(r))^(2)+1))+log(exp(1)/(2*pi))) #Lognormal Distribution
```


```R
y2r <- function(r) r
```


```R
f2r <- function(z, r)  sqrt(1/(2*pi*(y2r(r))^(2)*z^(3)))*exp(-(1/(2*(y2r(r))^(2)))*(((z-1)^(2))/z)) #Inverse Gaussian distribution
```


```R
logf2r <- function(z,r) -((1/2)*log(2*pi*(y2r(r))^(2)))-((3/2)*log(z))-(((z-1)^(2))/(2*((y2r(r))^(2))*z)) 
#Log of the pdf of the Inverse Gaussian Distribution
```


```R
f3r <- function(z,r) f2r(z,r)*logf2r(z,r) #Integrand for the entropy function
```


```R
f2intr <- function(r) -(integrate(f3r, lower=0, upper=Inf, r=r)$value) #Entropy of the Inverse Gaussian function
```


```R
v.hr <- Vectorize(f2intr) #Vectorize with respect to CV|
```


```R
KL_2r <- function(r) 1-v.hr(r) #Formula for the KL distance from the exponential distribution with mean =1 
```


```R
library("expint")
```


```R
dat <- data.frame(CV_X= seq(0.2,2.5,0.01))
```


```R
CV_R4r <- function(CV_X){
    sapply(CV_X, function(x) if (x<1) sqrt(((exp((1/x)-1)*expint_E1((1/x)-1, scale = FALSE))/x)-1) else NaN)
}
```


```R
KL_4r <- function(x) -log(x) #Shifted Exponential Distribution
```


```R
dat$CV_R4r <- CV_R4r(dat$CV_X)
```


```R
dat$KL_4r <- KL_4r(dat$CV_X)
```


```R
r <- seq(0.2,2.5,0.01)  
```


```R
KL1 <- KL_1r(r)
```


```R
KL2 <- KL_3r(r)
```


```R
KL3 <- KL_2r(r)
```


```R
plot(r, KL_1r(r), xaxs='i', yaxs='i',xlim=c(0,2.5), ylim=c(0,2), ylab="KL distance", xlab="CV(R)", type="o", pch=20, cex=0.4, col="red")
lines(r, KL_3r(r), type="o", pch=20, cex=0.4, lty=2, col="blue")
lines(r, KL_2r(r), col="green4", lwd=2.5, lty=3, type="o", pch=20, cex=0.4,)
lines(dat$CV_R4r, dat$KL_4r, type="o", pch=20, cex=0.4, lty=5, col="brown")
legend(1.2, 1.8, legend=c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted Exp."),
       col=c("red", "blue","green4", "brown"), lty=c(1,2,3,5), cex=0.65)
```


![png](output_33_0.png)



```R
plot(x, KL_1(x), yaxs='i', xaxs='i', xlim=c(0, 2.5), ylim=c(0,2), ylab="KL distance", xlab="CV(T)", type="o", pch=20, cex=0.4, col="red")
lines(x, KL_2(x), col="blue", type="o", lty=2, pch=20, cex=0.4)
lines(x, KL_3(x), col="forestgreen", type="o", lty=3, pch=20, cex =0.4)
lines(x, KL_4(x), col="brown", type="o", lty=4, pch=20, cex =0.4)
legend(1.2, 1.8, legend=c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted Exp."),
       col=c("red", "blue","forestgreen", "brown"), lty=c(1,2,3,4), cex=0.65)
```


![png](output_34_0.png)

