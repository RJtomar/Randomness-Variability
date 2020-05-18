```R
Cht_1 <- function(x) exp(-(1-log(x^(2))-log(gamma(1/x^(2)))+((digamma(1/x^(2))-1)/x^(2))-digamma(1/x^(2)))) #Gamma Distribution
```


```R
Cht_2 <- function(x) exp(-((1/2)*(log((x^(2)+1)/log(x^(2)+1))+log(exp(1)/(2*pi))))) #Lognormal Distribution
```


```R
Cht_4 <- function(x){
    sapply(x, function(x) if (x<1) exp(-(-log(x))) else NaN)
}
#Shifted Exponential Distribution
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
Cht_3 <- function(x) exp(-(1-v.h(x))) #Formula for the KL distance from the exponential distribution with mean =1 
```


```R
t <- seq(0.1,500,0.01) #Values of CV-range
```


```R
df <- data.frame("CVX"=t, "y1"= Cht_1(t), "y2" = Cht_2(t), "y3" = Cht_3(t), "y4" = Cht_4(t))
```


```R
require(ggplot2)
require(reshape2)
```


```R
df <- melt(df, id.vars="CVX", variable= "Distributions")
```


```R
library("repr")
options(repr.plot.width=10, repr.plot.height=3.5)
```


```R
library("ggplot2")
g <- ggplot(df, aes(x= CVX, y = value, col= Distributions, linetype=Distributions)) + geom_line(lwd=1) + 
  scale_linetype_manual(values = c('solid', 'solid', 'dashed', 'solid'), labels= c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted exponential") ) + scale_color_manual(name = "Distributions", 
  labels = c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted exponential"), values= c('red', 'green', 'blue', 'violet')) +
  guides(linetype=FALSE) 
g1 <- g + coord_cartesian(xlim=c(0,500), ylim=c(0,1.2))
g1 + labs ( y=expression(C[h](T)), x=expression(C[V](T))) + theme_classic() 

```

    Warning message:
    “Removed 49901 rows containing missing values (geom_path).”


![png](output_14_1.png)



```R
plot(t, Cht_1(t), yaxs='i', xaxs='i', xlim=c(0, 4), ylim=c(0,1.2), ylab="Dispersion Measure (ISI distributions)", xlab="CV(T)", type="o", pch=20, cex=0.4, col="red")
lines(t, Cht_2(t), col="blue", type="o", lty=2, pch=20, cex=0.4)
lines(t, Cht_3(t), col="forestgreen", type="o", lty=3, pch=20, cex =0.4)
lines(t, Cht_4(t), col="brown", type="o", lty=4, pch=20, cex =0.4)
legend(1.2, 1.8, legend=c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted Exp."),
       col=c("red", "blue","forestgreen", "brown"), lty=c(1,2,3,4), cex=0.65)
```


![png](output_15_0.png)


## From the AIFR perspective:


```R
yr <- function(r) (r)/sqrt(1+r^(2))
```


```R
Chr_1 <- function(r) exp(-(1-log((yr(r))^(2))-log(gamma(1/(yr(r))^(2)))+((digamma(1/(yr(r))^(2))-1)/(yr(r))^(2))-digamma(1/(yr(r))^(2)))) #Gamma Distribution
```


```R
y1r <- function(r) r
```


```R
Chr_2 <- function(r) exp(-((1/2)*(log(((y1r(r))^(2)+1)/log((y1r(r))^(2)+1))+log(exp(1)/(2*pi))))) #Lognormal Distribution
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
Chr_3 <- function(r) exp(-(1-v.hr(r))) #Formula for the KL distance from the exponential distribution with mean =1 
```


```R
library("expint")
```


```R
dat <- data.frame(CV_X= seq(0.2,500,0.01))
```


```R
CV_R4r <- function(CV_X){
    sapply(CV_X, function(x) if (x<1) sqrt(((exp((1/x)-1)*expint_E1((1/x)-1, scale = FALSE))/x)-1) else NaN)
}
```


```R
Chr_4 <- function(x) exp(-(-log(x))) #Shifted Exponential Distribution
```


```R
dat$CV_R4r <- CV_R4r(dat$CV_X)
```


```R
dat$Chr_4 <- Chr_4(dat$CV_X)
```


```R
r <- seq(0.2,500,0.01)  
```


```R
df1 <- data.frame("CVR"=r, "z1" = Chr_1(r), "z2" = Chr_2(r), "z3" = Chr_3(r))
```


```R
df2 <- data.frame("xr"= dat$CV_R4r, "yr" = dat$Chr_4)
```


```R
library(ggplot2)
library(reshape2)
```


```R
library("repr")
options(repr.plot.width=10, repr.plot.height=3.5)
```


```R
g <- ggplot()+
geom_line(data=df1, aes(x=CVR, y=z1, color="Gamma"), size=1, show.legend = TRUE)+
geom_line(data=df1, aes(x=CVR, y=z2, color="Lognormal"), size=1, show.legend = TRUE)+
geom_line(data=df1, aes(x=CVR, y=z3, color="Inverse Gaussian"), size=1, linetype="dashed", show.legend = TRUE)+
geom_line(data=df2, aes(x=xr, y=yr, color="Shifted exponential"), size=1, show.legend = TRUE) + 
scale_colour_manual(name = "Distributions", 
         values =c(Gamma="red", Lognormal="green", 'Inverse Gaussian'="blue", 'Shifted exponential'="violet"))
g1 <- g + coord_cartesian(xlim=c(0,500), ylim=c(0,1.2))
g1 + labs ( y=expression(C[h](T)), x=expression(C[V](R))) + theme_classic()

```

    Warning message:
    “Removed 49901 rows containing missing values (geom_path).”


![png](output_39_1.png)



```R
plot(r, Chr_1(r), xaxs='i', yaxs='i',xlim=c(0,4), ylim=c(0,1.2), ylab="Dispersion Measure (ISI distributions)", xlab="CV(R)", type="o", pch=20, cex=0.4, col="red")
lines(r, Chr_2(r), type="o", pch=20, cex=0.4, lty=2, col="blue")
lines(r, Chr_3(r), col="green4", lwd=2.5, lty=3, type="o", pch=20, cex=0.4)
lines(dat$CV_R4r, dat$Chr_4, type="o", pch=20, cex=0.4, lty=5, col="brown")
legend(1.2, 1.8, legend=c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted Exp."),
       col=c("red", "blue","green4", "brown"), lty=c(1,2,3,5), cex=0.65)
```


![png](output_40_0.png)



```R

```
