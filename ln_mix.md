```R
p <- 0.4
```


```R
m1 <- +1
```


```R
m2 <- +2
```


```R
s1 <- function(v) v
```


```R
s2 <- function(v) v-1
```


```R
mm <- function(v) p*exp(m1+((s1(v))^{2}/2))+(1-p)*exp(m2+((s2(v))^{2}/2))
```


```R
CV <- function(v) (1/mm(v))*(p*exp(2*m1+2*s1(v)^{2})+(1-p)*exp(2*m2+2*s2(v)^{2})-(mm(v))^{2})^{1/2}
```


```R
phi <- function(x, m, s) (1/sqrt(2*pi*s^{2}))*exp((-1/(2*s^{2}))*(x-m)^{2})
```


```R
g <- function(x,v) p*phi(x,m1,s1(v)) +(1-p)*phi(x,m2,s2(v))
```


```R
f <- function(x,v) ifelse(x==0, 0, (1/x)*(g(log(x),v)))
```


```R
logf <- function(x,v) -log(x)+log(g(log(x),v))
```


```R
ft <- function(x,v) f(x,v)*logf(x,v)
```


```R
h <- function(v) -(integrate(ft, lower=0, upper=Inf, v=v)$value)
```


```R
he <- Vectorize(h)
```


```R
CH <- function(v) (1/mm(v))*exp(he(v))
```


```R
v <- seq(1.01, 2, 0.01)
```


```R
Ch <- CH(v)
```


```R
Cv <- CV(v)
```


```R
df <- data.frame("CV"=v, "Ch"=Cv)
```


```R
library("repr")
options(repr.plot.width=5, repr.plot.height=5)
```


```R
 library("ggplot2")
g <- ggplot()+
  geom_path(data=df, aes(x=Cv, y=Ch), size=1, arrow= arrow(), color="skyblue")
g1 <- g + coord_cartesian(xlim=c(1.4,8), ylim=c(0,3)) 
g1 + labs ( y=expression(C[h](T)), x=expression(C[V](T))) + theme_classic()
#For normal std dev ranging from (1,2)
```


![png](output_20_0.png)



```R
 library("ggplot2")
g <- ggplot()+
  geom_path(data=df, aes(x=Cv, y=Ch), size=1, arrow= arrow(), color="skyblue")
g1 <- g + coord_cartesian(xlim=c(0.4,2), ylim=c(1.4,2.6)) 
g1 + labs ( y=expression(C[h](T)), x=expression(C[V](T))) + theme_classic()
#For normal std dev ranging from (1.2,2.2)
```


![png](output_21_0.png)



```R
 library("ggplot2")
g <- ggplot()+
  geom_path(data=df, aes(x=v, y=Cv), size=1, arrow= arrow(), color="skyblue")
g1 <- g + coord_cartesian(xlim=c(0,4), ylim=c(0,7)) 
g1 + labs ( y=expression(C[V](T)), x="s") + theme_classic()
#For normal std dev ranging from (1.01,2) and (0.01,1)
```


![png](output_22_0.png)



```R
v[57]
```


1.57



```R
fr <- function(x,v) (1/mm(v))*(1/x^{2})*(g(-log(x),v))
```


```R
logfr <- function(x,v) -log(mm(v))-2*log(x)+log(g(-log(x),v))
```


```R
frr <- function(x,v) fr(x,v)*logfr(x,v)
```


```R
hr <- function(v) -(integrate(frr, lower=0, upper=Inf, v=v)$value)
```


```R
hrv <- Vectorize(hr)
```


```R
CHr <- function(v) (1/mm(v))*exp(hrv(v))
```


```R
v <- seq(1,2, 0.01)
```


```R
Chr <- CHr(v)
```


```R
which.max(CHr(v))
```


33



```R
Cv[33]
```


1.28883012952981



```R
df1 <- data.frame("CV"=Cv, "Ch2"=Chr)
```


```R
library("repr")
options(repr.plot.width=5, repr.plot.height=5)
```


```R

g <- ggplot()+
  geom_path(data=df, aes(x=Cv, y=Chr), size=1, arrow= arrow(), color="maroon")
g1 <- g + coord_cartesian(xlim=c(1.4,8), ylim=c(0,.03)) 
g1 + labs ( y=expression(C[h](R)), x=expression(C[V](T))) + theme_classic()
#For normal std dev ranging from (-1,2)
```


![png](output_36_0.png)



```R

g <- ggplot()+
  geom_path(data=df, aes(x=Cv, y=Chr), size=1, arrow= arrow(), color="maroon")
g1 <- g + coord_cartesian(xlim=c(0.38,2), ylim=c(0.018,0.055)) 
g1 + labs ( y=expression(C[h](R)), x=expression(C[V](T))) + theme_classic()
#For normal std dev ranging from (0.2,1.2)
```


![png](output_37_0.png)



```R

g <- ggplot()+
  geom_path(data=df, aes(x=Cv, y=Chr), size=1, arrow= arrow(), color="maroon")
g1 <- g + coord_cartesian(xlim=c(0.5,6), ylim=c(00.005,0.04)) 
g1 + labs ( y=expression(C[h](R)), x=expression(C[V](T))) + theme_classic()
#For normal std dev ranging from (1.01,2) and (0.01,1)
```


![png](output_38_0.png)



```R
df2 <- data.frame("Ch1"=Ch, "Ch2"=Chr)
```


```R
library("repr")
options(repr.plot.width=5, repr.plot.height=5)
```


```R
 library("ggplot2")
g <- ggplot()+
  geom_path(data=df2, aes(x=Ch1, y=Ch2), size=1, arrow= arrow(), color="forestgreen")
g1 <- g + coord_cartesian(xlim=c(1,2.5), ylim=c(0,0.032)) 
g1 + labs ( y=expression(C[h](R)), x=expression(C[h](T))) + theme_classic()
#For normal std dev ranging from (-1,2)
```


![png](output_41_0.png)



```R
 library("ggplot2")
g <- ggplot()+
  geom_path(data=df2, aes(x=Ch1, y=Ch2), size=1, arrow= arrow(), color="forestgreen")
g1 <- g + coord_cartesian(xlim=c(1.4,2.5), ylim=c(0.015,0.06)) 
g1 + labs ( y=expression(C[h](R)), x=expression(C[h](T))) + theme_classic()
#For normal std dev ranging from (-1,2)
```


![png](output_42_0.png)



```R
 library("ggplot2")
g <- ggplot()+
  geom_path(data=df2, aes(x=Ch1, y=Ch2), size=1, arrow= arrow(), color="forestgreen")
g1 <- g + coord_cartesian(xlim=c(0.4,2.5), ylim=c(0.005,0.04)) 
g1 + labs ( y=expression(C[h](R)), x=expression(C[h](T))) + theme_classic()
#For normal std dev ranging from (1.01,2) and (0.01, 1)
```


![png](output_43_0.png)



```R
td <- data.frame("v"=v, "m"=mm(v))
```


```R
 library("ggplot2")
g <- ggplot()+
  geom_path(data=td, aes(x=v, y=m), size=1, arrow= arrow(), color="forestgreen")
g1 <- g + coord_cartesian(xlim=c(0,10), ylim=c(0,20)) 
g1 + labs ( y=expression(C[h](R)), x=expression(C[h](T))) + theme_classic()
#For normal std dev ranging from (1.01,2) and (0.01, 1)
```


![png](output_45_0.png)



```R

```
