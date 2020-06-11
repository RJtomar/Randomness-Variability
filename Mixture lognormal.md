```R
m1 <- -1
```


```R
m2 <- -0.5
```


```R
s1 <- 0.2
```


```R
s2 <- 1
```


```R
mm <- function(p) p*exp(m1+(s1^{2}/2))+(1-p)*exp(m2+(s2^{2}/2))
```


```R
CV <- function(p) (1/mm(p))*((p*exp(2*m1+2*s1^{2}))+(1-p)*exp(2*m2+2*s2^{2})-(mm(p))^{2})^{1/2}
```


```R
p <- seq(0,1,0.01)
```


```R
Cv <- CV(p)
```


```R
which.max(Cv)
```


52



```R
Cv[52]
```


1.42466339981745



```R
phi <- function(x, m, s) (1/sqrt(2*pi*s^{2}))*exp((-1/(2*s^{2}))*(x-m)^{2})
```


```R
g <- function(x, p) p*phi(x,m1,s1) +(1-p)*phi(x,m2,s2)
```


```R
f <- function(x,p) (1/x)*(g(log(x),p))
```


```R
logf <- function(x,p) -log(x)+log(g(log(x),p))
```


```R
ft <- function(x,p) f(x,p)*logf(x,p)
```


```R
h <- function(p) -(integrate(ft, lower=0, upper=Inf, p=p)$value)
```


```R
he <- Vectorize(h)
```


```R
CH <- function(p) (1/mm(p))*exp(he(p)-1)
```


```R
Ch <- CH(p)
```


```R
df <- data.frame("CV"=Cv, "Ch"=Ch)
```


```R
library("repr")
options(repr.plot.width=4, repr.plot.height=4)
```


```R
 library("ggplot2")
g <- ggplot()+
  geom_path(data=df, aes(x=Cv, y=Ch), size=1.8, arrow= arrow(), color="blue1")
g1 <- g + coord_cartesian(xlim=c(0,1.5), ylim=c(0,1)) + scale_x_continuous(breaks = seq(0, 1.5, by = 0.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5))
g1 + labs ( y=expression(C[h](T)), x=expression(C[V](T))) + theme_classic()
```


![png](output_21_0.png)



```R
fr <- function(x,p) (1/mm(p))*(1/x^{2})*(g(-log(x),p))
```


```R
logfr <- function(x,p) -log(mm(p))-2*log(x)+log(g(-log(x),p))
```


```R
frr <- function(x,p) fr(x,p)*logfr(x,p)
```


```R
hr <- function(p) -(integrate(frr, lower=0, upper=Inf, p=p)$value)
```


```R
hrv <- Vectorize(hr)
```


```R
CHr <- function(p) (1/mm(p))*exp(hrv(p)-1)
```


```R
Chr <- CHr(p)
```


```R
which.max(Chr)
```


88



```R
p[88]
```


0.87



```R
df1 <- data.frame("CV"=Cv, "Ch2"=Chr)
```


```R
library("repr")
options(repr.plot.width=4, repr.plot.height=4)
```


```R

g <- ggplot()+
  geom_path(data=df, aes(x=Cv, y=Chr), size=1.8, arrow= arrow(), color="maroon")
g1 <- g + coord_cartesian(xlim=c(0,1.5), ylim=c(0.5,3)) + scale_x_continuous(breaks = seq(0, 1.5, by = 0.5)) +
  scale_y_continuous(breaks = seq(0, 3, by = 1))
g1 + labs ( y=expression(C[h](R)), x=expression(C[V](T))) + theme_classic()
```


![png](output_33_0.png)



```R

```


```R
df2 <- data.frame("Ch1"=Ch, "Ch2"=Chr)
```


```R
library("repr")
options(repr.plot.width=4, repr.plot.height=4)
```


```R
 library("ggplot2")
g <- ggplot()+
  geom_path(data=df2, aes(x=Ch1, y=Ch2), size=1.8, arrow= arrow(), color="forestgreen")
g1 <- g + coord_cartesian(xlim=c(0.2,0.9), ylim=c(0.5,3)) + scale_x_continuous(breaks = seq(0.2, 1, by = 0.2)) +
  scale_y_continuous(breaks = seq(0, 3, by = 1))
g1 + labs ( y=expression(C[h](R)), x=expression(C[h](T))) + theme_classic()
```


![png](output_37_0.png)



```R
max(Ch)
```


2.50662828707811



```R
which.max(Cv)
```


52



```R
p[52]
```


0.51



```R
max(Cv)
```


1.42466339981745



```R
max(Ch)
```


2.50662828707811



```R
which.max(Ch)
```


1



```R
which.max(Chr)
```


88



```R
max(Chr)
```


7.66404691213608



```R
p[88]
```


0.87



```R
CHr(0.87)
```


7.66404691213608



```R
max(mm(p))
```


1



```R

```


```R

```


```R

```


```R

```


```R

```

## Mixture Normal with different variabilities




```R
m1 <- +1
```


```R
m2 <- +2
```


```R
s1 <- 1.85
```


```R
s2 <- 0.85
```


```R
p <- 0.4
```


```R
mm <- p*exp(m1+(s1^{2}/2))+(1-p)*exp(m2+(s2^{2}/2))
```


```R
CV <- (1/mm)*(p*exp(2*m1+2*s1^{2})+(1-p)*exp(2*m2+2*s2^{2})-(mm)^{2})^{1/2}
```


```R
phi <- function(x, m, s) (1/sqrt(2*pi*s^{2}))*exp((-1/(2*s^{2}))*(x-m)^{2})
```


```R
g <- function(x) p*phi(x,m1,s1) +(1-p)*phi(x,m2,s2)
```


```R
f <- function(x) ifelse(x==0, 0, (1/x)*(g(log(x))))
```


```R
f1 <- function(x) ifelse(x==0, 0, (1/x)*phi(log(x), m1, s1))
  
```


```R
f2 <- function(x) ifelse(x==0, 0, (1/x)*phi(log(x), m2, s2))
```


```R
x <- seq(0, 20, 0.001)
```


```R
f1x <- f1(x)
```


```R
f2x <- f2(x)
```


```R
fx <- f(x)
```


```R
dff <- data.frame("Xval"=x, "Fval"=fx)
```


```R
dff1 <- data.frame("X"=x, "Y1"=f1x)
```


```R
dff2 <- data.frame("X"=x, "Y2"=f2x)
```


```R
library("ggplot2")
g <- ggplot()+ 
geom_path(data=dff1, aes(x=X, y=Y1, col="red"))
g1 <- g + coord_cartesian(xlim=c(0,20), ylim=c(0, 5))
g1 + theme_classic()
```


![png](output_74_0.png)



```R
library("ggplot2")
g <- ggplot()+ 
geom_path(data=dff2, aes(x=X, y=Y2, col="red"))
g1 <- g + coord_cartesian(xlim=c(0,20), ylim=c(0, 1))
g1 + theme_classic()
```


![png](output_75_0.png)



```R
library("ggplot2")
g <- ggplot()+ 
geom_path(data=dff, aes(x=Xval, y=Fval, col="red"))
g1 <- g + coord_cartesian(xlim=c(0,20), ylim=c(0, 0.20))
g1 + theme_classic()
```


![png](output_76_0.png)



```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```
