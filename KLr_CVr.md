# Gamma Distribution


```R
y <- function(r) (r)/sqrt(1+r^(2))
```


```R
f <- function(x,r) ((1/(y(r))^(2))^(1/(y(r))^(2)))*(1/gamma(1/(y(r))^(2)))*(x^((1/(y(r))^(2))-1))*exp(-x/((y(r))^(2)))
```


```R
f1 <- function(x,r) f((1/x),r)/(x^(3))
```


```R

logf1 <- function(x,r) ((1/(y(r))^(2))*log(1/((y(r))^(2))))-log(gamma(1/(y(r))^(2)))-((1/(y(r))^(2))+2)*log(x)-(1/(((y(r))^(2))*x))
```


```R
f2 <- function(x,r) f1(x,r)*logf1(x,r)
```


```R
f3 <- function(r) -(integrate(f2, lower=0, upper=Inf, r=r)$value)
```


```R
h.f <- Vectorize(f3) #Vectorize with respect to CV
```


```R
KL_f <- function(r) 1-h.f(r) #Formula for the KL distance from the exponential distribution with mean =1 
```

## Lognormal Distribution


```R
y1 <- function(r) r
```


```R
g <- function(x,r) (1/(x*sqrt(2*pi*log(1+(y1(r))^(2)))))*exp(-(1/8)*(((log(1+(y1(r))^(2)))+2*log(x))^(2))/log(1+(y1(r))^(2)))
```


```R
g1 <- function(x,r) g((1/x),r)/(x^(3))
```


```R
logg1 <- function(x,r) -(1/2)*log(2*pi*log(1+(y1(r))^(2)))-2*log(x)-(1/8)*(((log(1+(y1(r))^(2))-2*log(x))^(2))/log(1+(y1(r))^(2)))
```


```R
g2 <- function(x,r) g1(x,r)*logg1(x,r)
```


```R
g3 <- function(r) -(integrate(g2, lower=0, upper=Inf, r=r)$value)
```


```R
h.g <- Vectorize(g3) #Vectorize with respect to CV
```


```R
KL_g <- function(r) 1-h.g(r) #Formula for the KL distance from the exponential distribution with mean =1 
```

## Inverse Gaussian Distribution


```R
y2 <- function(r) r
```


```R
j <- function(x,r) sqrt(1/(2*pi*(y2(r))^(2)*x^(3)))*exp(-(1/(2*(y2(r))^(2)))*(((x-1)^(2))/x))
```


```R
j1 <- function(x,r) j((1/x),r)/(x^(3))
```


```R
logj1 <- function(x,r) -(1/2)*log(2*pi*(y2(r))^(2))-(3/2)*log(x)-(1/(2*(y2(r))^(2)))*(((1-x)^(2))/x)
```


```R
j2 <- function(x,r) j1(x,r)*logj1(x,r)
```


```R
j3 <- function(r) -(integrate(j2, lower=0, upper=Inf, r=r)$value)
```


```R
h.j <- Vectorize(j3) #Vectorize with respect to CV
```


```R
KL_j <- function(r) 1-h.j(r) #Formula for the KL distance from the exponential distribution with mean =1 
```

## Shifted Exponential 


```R
dat <- data.frame(CV_X= seq(0.2,200,0.01))
```


```R
CV_R4 <- function(CV_X){
    sapply(CV_X, function(x) if (x<1) sqrt(((exp((1/x)-1)*expint_E1((1/x)-1, scale = FALSE))/x)-1) else NaN)
}
```


```R
library("expint")
```


```R
s <- function(x,y) (1/y)*exp(-(1/y)*(x-1+y))
```


```R
s1 <- function(x,y) s((1/x),y)*(1/(x^{3}))
```


```R
logs1 <- function(x,y){-3*log(x)-log(y)-(1/y)*((1/x)-1+y)}
```


```R
s2 <- function(x,y) s1(x,y)*logs1(x,y)
```


```R
s3 <- function(y) -(integrate(s2, lower=0, upper=(1/(1-y)), y=y)$value)
```


```R
h.s <- Vectorize(s3) 
```


```R
KL_s <- function(y){
    sapply(y, function(x) if (x<1) 1-h.s(x) else NaN)
        }
```


```R
dat$CV_R4 <- CV_R4(dat$CV_X)
```


```R
dat$KL_s <- KL_s(dat$CV_X)
```

## Plotting


```R
r <- seq(0.2,2.5,0.01) #Values of CV-range
```


```R
plot(r, KL_f(r), xaxs='i', yaxs='i',xlim=c(0,2.5), ylim=c(0,2), ylab="KL distance", xlab="CV (AIFR)", type="o", pch=20, cex=0.4, col="red")
lines(r, KL_g(r), type="o", pch=20, cex=0.4, lty=2, col="blue")
lines(r, KL_j(r), type="o", pch=20, cex=0.4, lty=3, col="green4")
lines(dat$CV_R4, dat$KL_s, type="o", pch=20, cex=0.4, lty=5, col="brown")
legend(1.2, 1.8, legend=c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted Exp."),
       col=c("red", "blue","green4", "brown"), lty=c(1,2,3,5), cex=0.65)
```


![png](output_42_0.png)

