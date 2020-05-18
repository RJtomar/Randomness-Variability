## Gamma Distribution


```R
f <- function(x,y) ((1/y^(2))^(1/y^(2)))*(1/gamma(1/y^(2)))*x^((1/y^(2))-1)*exp(-x/y^(2))
```


```R
f1 <- function(x,y) f((1/x),y)/(x^(3))
```


```R
logf1 <- function(x,y) -((1/y^(2))*log(y^(2)))-log(gamma(1/y^(2)))-((1/y^(2))+2)*log(x)-(1/(x*y^(2)))
```


```R
f2 <- function(x,y) f1(x,y)*logf1(x,y)
```


```R
f3 <- function(y) -(integrate(f2, lower=0, upper=Inf, y=y)$value)
```


```R
h.f <- Vectorize(f3) #Vectorize with respect to CV
```


```R
KL_f <- function(y) 1-h.f(y) #Formula for the KL distance from the exponential distribution with mean =1 
```


```R
y <- seq(0,2.5, 0.01)
```


```R
KL_f(y)
```

    Warning message in f((1/x), y):
    “value out of range in 'gammafn'”Warning message in logf1(x, y):
    “value out of range in 'gammafn'”


    Error in integrate(f2, lower = 0, upper = Inf, y = y): non-finite function value
    Traceback:


    1. KL_f(y)

    2. h.f(y)   # at line 1 of file <text>

    3. do.call("mapply", c(FUN = FUN, args[dovec], MoreArgs = list(args[!dovec]), 
     .     SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES))

    4. mapply(FUN = function (y) 
     . -(integrate(f2, lower = 0, upper = Inf, y = y)$value), y = c(0, 
     . 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 
     . 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 
     . 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 
     . 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 
     . 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 
     . 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 
     . 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 
     . 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 
     . 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 
     . 1, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 
     . 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.2, 1.21, 
     . 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3, 1.31, 1.32, 
     . 1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 1.4, 1.41, 1.42, 1.43, 
     . 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53, 1.54, 
     . 1.55, 1.56, 1.57, 1.58, 1.59, 1.6, 1.61, 1.62, 1.63, 1.64, 1.65, 
     . 1.66, 1.67, 1.68, 1.69, 1.7, 1.71, 1.72, 1.73, 1.74, 1.75, 1.76, 
     . 1.77, 1.78, 1.79, 1.8, 1.81, 1.82, 1.83, 1.84, 1.85, 1.86, 1.87, 
     . 1.88, 1.89, 1.9, 1.91, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97, 1.98, 
     . 1.99, 2, 2.01, 2.02, 2.03, 2.04, 2.05, 2.06, 2.07, 2.08, 2.09, 
     . 2.1, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.2, 
     . 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.31, 
     . 2.32, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 
     . 2.43, 2.44, 2.45, 2.46, 2.47, 2.48, 2.49, 2.5), MoreArgs = structure(list(), .Names = character(0)), 
     .     SIMPLIFY = TRUE, USE.NAMES = TRUE)

    5. (function (y) 
     . -(integrate(f2, lower = 0, upper = Inf, y = y)$value))(y = dots[[1L]][[1L]])

    6. integrate(f2, lower = 0, upper = Inf, y = y)   # at line 1 of file <text>


## Lognormal Distribution


```R
g <- function(x,y) (1/(x*sqrt(2*pi*log(1+y^(2)))))*exp(-(1/8)*(((log(1+(y^(2)))+2*log(x))^(2))/log(1+(y^(2)))))
```


```R
g1 <- function(x,y) g((1/x),y)/(x^(3))
```


```R
logg1 <- function(x,y) -2*log(x)-(1/2)*log(2*pi*log(1+(y^(2))))-(1/8)*(((log(1+(y^(2)))-2*log(x))^(2))/log(1+(y^(2))))
```


```R
g2 <- function(x,y) g1(x,y)*logg1(x,y)
```


```R
g3 <- function(y) -(integrate(g2, lower=0, upper=Inf, y=y)$value)
```


```R
h.g <- Vectorize(g3) #Vectorize with respect to CV
```


```R


KL_g <- function(y) 1-h.g(y) #Formula for the KL distance from the exponential distribution with mean =1 
```

## Inverse Gaussian Distribution


```R
j <- function(x,y) sqrt(1/(2*pi*y^(2)*x^(3)))*exp(-(1/(2*y^(2)))*(((x-1)^(2))/x))
```


```R
j1 <- function(x,y) j((1/x),y)/(x^(3))
```


```R
logj1 <- function(x,y) -(1/2)*log(2*pi*y^(2)*x^(3))-((1/(2*y^(2)))*(((1-x)^(2))/x))
```


```R
j2 <- function(x,y) j1(x,y)*logj1(x,y)
```


```R
j3 <- function(y) -(integrate(j2, lower=0, upper=Inf, y=y)$value)
```


```R
h.j <- Vectorize(j3) #Vectorize with respect to CV
```


```R
KL_j <- function(y) 1-h.j(y) #Formula for the KL distance from the exponential distribution with mean =1 
```


```R
#plot(y, KL_f(y), xlim=c(0,2), ylim=c(0,1.2), ylab="KL", xlab="CV (ISI)", type="l", col="red")
#lines(y, KL_g(y), type="l", lty=2, col="blue")
#lines(y, KL_j(y), type="l", lty=3, col="seagreen2")
```

## Shifted Exponential Distribution


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
#If y=1 then evaluate the integral 
#from 0 to infinity

################
#I DON'T KNOW HOW TO DEAL WITH THIS STUPID FUNCTION WITH HALF THE VALUES OF Y SO 
#I AM GIVING IT RANDOM RETURN VALUES FOR WHEN CV IS MORE THAN 1, I AM SORRY
```


```R
h.s <- Vectorize(s3) #Vectorize with respect to CV
```


```R
KL_s <- function(y){
    sapply(y, function(x) if (x<1) 1-h.s(x) else NaN)
        }
```


```R
y <- seq(0.2,2.5, 0.01)
```


```R
y[66]
```


0.85



```R
plot(y, KL_f(y), xaxs='i', yaxs='i', xlim=c(0,2.5), ylim=c(0,2), ylab="KL distance", xlab="CV(T)", type="o", pch=20, cex=0.4, col="red")
lines(y, KL_g(y), type="o", pch=20, cex=0.4, lty=2, col="blue")
lines(y, KL_j(y), type="o", pch=20, cex=0.4, lty=3, col="green4")
lines(y, KL_s(y), type="o", pch=20, cex=0.4, lty=5, col="brown")
legend(1.2, 1.8, legend=c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted Exp."),
       col=c("red", "blue","green4", "brown"), lty=c(1,2,3,5), cex=0.65)
```


![png](output_37_0.png)



```R
plot(y, KL_f(y), xaxs='i', yaxs='i', xlim=c(0,2.5), ylim=c(0,2), ylab="KL distance", xlab="CV(X)", type="o", pch=20, cex=0.4, col="red")
lines(y, KL_g(y), type="o", pch=20, cex=0.4, lty=2, col="blue")
lines(y, KL_j(y), type="o", pch=20, cex=0.4, lty=3, col="green4")
lines(y, KL_s(y), type="o", pch=20, cex=0.4, lty=5, col="brown")
legend(1.2, 1.8, legend=c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted Exp."),
       col=c("red", "blue","green4", "brown"), lty=c(1,2,3,5), cex=0.65)
```


![png](output_38_0.png)



```R

```
