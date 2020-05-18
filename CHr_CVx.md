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
Ch_f <- function(y) exp(-(1-h.f(y))) #Formula for the KL distance from the exponential distribution with mean =1 
```

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


Ch_g <- function(y) exp(-(1-h.g(y))) #Formula for the KL distance from the exponential distribution with mean =1 
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
Ch_j <- function(y) exp(-(1-h.j(y))) #Formula for the KL distance from the exponential distribution with mean =1 
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
Ch_s <- function(y){
    sapply(y, function(x) if (x<1) exp(-(1-h.s(x))) else NaN)
        }
```


```R
y <- seq(0.2, 200, 0.01)
```


```R
df <- data.frame("x"=y, "y1"=Ch_f(y), "y2"=Ch_g(y), "y3"=Ch_j(y), "y4"=Ch_s(y))
```


```R
library(ggplot2)
library(reshape2)
```


```R
df <- melt(df, id.vars="x", variable="Distributions")
```


```R
library("repr")
options(repr.plot.width=10, repr.plot.height=3.5)
```


```R
library("ggplot2")
g <- ggplot(df, aes(x= x, y = value, col= Distributions, linetype=Distributions)) + geom_line(lwd=1) + 
  scale_linetype_manual(values = c('solid', 'solid', 'dashed', 'solid'), labels= c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted exponential") ) + scale_color_manual(name = "Distributions", 
  labels = c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted exponential"), values= c('red', 'green', 'blue', 'violet')) +
  guides(linetype=FALSE) 
g1 <- g + coord_cartesian(xlim=c(0,200), ylim=c(0,1.2))
g1 + labs ( y=expression(C[h](R)), x=expression(C[V](T))) + theme_classic() 

```

    Warning message:
    “Removed 19901 rows containing missing values (geom_path).”


![png](output_38_1.png)



```R
plot(y, Ch_f(y), xaxs='i', yaxs='i', xlim=c(0,4), ylim=c(0,1.2), ylab="Dispersion Measure (Instantaneous Rate)", xlab="CV(T)", type="o", pch=20, cex=0.4, col="red")
lines(y, Ch_g(y), type="o", pch=20, cex=0.4, lty=2, col="blue")
lines(y, Ch_j(y), type="o", pch=20, cex=0.4, lty=3, col="green4")
lines(y, Ch_s(y), type="o", pch=20, cex=0.4, lty=5, col="brown")
legend(1.2, 1.8, legend=c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted Exp."),
       col=c("red", "blue","green4", "brown"), lty=c(1,2,3,5), cex=0.65)
```


![png](output_39_0.png)



```R
plot(y, KL_f(y), xaxs='i', yaxs='i', xlim=c(0,2.5), ylim=c(0,2), ylab="KL distance", xlab="CV(X)", type="o", pch=20, cex=0.4, col="red")
lines(y, KL_g(y), type="o", pch=20, cex=0.4, lty=2, col="blue")
lines(y, KL_j(y), type="o", pch=20, cex=0.4, lty=3, col="green4")
lines(y, KL_s(y), type="o", pch=20, cex=0.4, lty=5, col="brown")
legend(1.2, 1.8, legend=c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted Exp."),
       col=c("red", "blue","green4", "brown"), lty=c(1,2,3,5), cex=0.65)
```


![png](output_40_0.png)



```R

```
