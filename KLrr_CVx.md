## Gamma Distribution


```R
f <- function(x,y) ((1/y^(2))^(1/y^(2)))*(1/gamma(1/y^(2)))*x^((1/y^(2))-1)*exp(-x/y^(2))
```


```R
f1 <- function(x,y) f((1/x),y)/(x^(3))
```


```R
f1byx <- function(x,y) f1(x,y)/x
```


```R
f1logx <- function(x,y) f1(x,y)*log(x)
```


```R
logf1 <- function(x,y) -(1/y^(2))*log(y^(2))-log(gamma(1/y^(2)))+((1/y^(2))+2)*log(1/x)-(1/((y^(2))*x))
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
f4 <- function(y) integrate(f1byx, lower=0, upper=Inf, y=y)$value
```


```R
f4v <- Vectorize(f4)
```


```R
f5 <- function(y) integrate(f1logx, lower=0, upper=Inf, y=y)$value
```


```R
f5v <- Vectorize(f5)
```


```R
KL_f <- function(y) -h.f(y)+f4v(y)+(3*f5v(y))
```

# Lognormal Distribution


```R
g <- function(x,y) (1/(x*sqrt(2*pi*log(1+y^(2)))))*exp(-(1/8)*(((log(1+(y^(2)))+2*log(x))^(2))/log(1+(y^(2)))))
```


```R
g1 <- function(x,y) g((1/x),y)/(x^(3))
```


```R
g1byx <- function(x,y) g1(x,y)/x
```


```R
g1logx <- function(x,y) g1(x,y)*log(x)
```


```R
logg1 <- function(x,y) 2*log(1/x)-(1/2)*log(2*pi*log(1+(y^(2))))-(1/8)*(((log(1+(y^(2)))+2*log(1/x))^(2))/log(1+(y^(2))))
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

```


```R
g4 <- function(y) (integrate(g1byx, lower=0, upper=Inf, y=y)$value)
```


```R
g4v <- Vectorize(g4)
```


```R
g5 <- function(y) integrate(g1logx, lower=0, upper=Inf, y=y)$value
```


```R
g5v <- Vectorize(g5)
```


```R
KL_g <- function(y) -h.g(y)+g4v(y)+3*g5v(y)
```

# Inverse Gamma Distribution


```R
j <- function(x,y) sqrt(1/(2*pi*y^(2)*x^(3)))*exp(-(1/(2*y^(2)))*(((x-1)^(2))/x))
```


```R
j1 <- function(x,y) j((1/x),y)/(x^(3))
```


```R
j1byx <- function(x,y) j1(x,y)/x
```


```R
j1logx <- function(x,y) j1(x,y)*log(x)
```


```R
logj1 <- function(x,y) -(1/2)*log(2*pi*y^(2))+(3/2)*log(1/x)-(1/(2*y^(2)))*(((1-(1/x))^(2))/(1/x))
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

```


```R
j4 <- function(y) (integrate(j1byx, lower=0, upper=Inf, y=y)$value)
```


```R
j4v <- Vectorize(j4)
```


```R
j5 <- function(y) integrate(j1logx, lower=0, upper=Inf, y=y)$value
```


```R
j5v <- Vectorize(j5)
```


```R
KL_j <- function(y) -h.j(y)+j4v(y)+3*j5v(y)
```

# Shifted Exponential Distribution


```R
s <- function(x,y) (1/y)*exp(-(1/y)*(x-1+y))
```


```R
s1 <- function(x,y) s((1/x),y)*(1/(x^{3}))
```


```R
s1byx <- function(x,y) s1(x,y)/x
```


```R
s1logx <- function(x,y) s1(x,y)*log(x)
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
s4 <- function(y) (integrate(s1byx, lower=0, upper=(1/(1-y)), y=y)$value)
```


```R
s4v <- Vectorize(s4)
```


```R
s5 <- function(y) integrate(s1logx, lower=0, upper=(1/(1-y)), y=y)$value
```


```R
s5v <- Vectorize(s5)
```


```R
KL_s <- function(y){
    sapply(y, function(x) if (x<1) -h.s(x)+s4v(x)+3*s5v(x) else NaN)
        }
```

# ISI Duration and General Plots


```R
y <- seq(0.2,2.5,0.01)  

```


```R
y[75]
```


0.94



```R
plot(y, KL_f(y), xaxs='i', yaxs='i',xlim=c(0,2.5), ylim=c(0,2), ylab="KL distance", xlab="CV(T)", type="o", pch=20, cex=0.4, col="red")
lines(y, KL_g(y), type="o", pch=20, cex=0.4, lty=2, col="blue")
lines(y, KL_j(y), col="green4", type="o", pch=20, cex=0.4, lwd=2.5, lty=3)
lines(y, KL_s(y), col="brown", lty=5, type="o", pch=20, cex=0.4)
legend(1.2, 1.8, legend=c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted Exp."),
       col=c("red", "blue","green4", "brown"), lty=c(1,2,3,5), cex=0.65)
```


![png](output_61_0.png)



```R

```
