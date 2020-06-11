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
mm <- function(v) p*exp(m1+(((s1(v))^{2})/2))+(1-p)*exp(m2+(((s2(v))^{2})/2))
```


```R
CV <- function(v) (1/mm(v))*((p*exp(2*m1+2*(s1(v))^{2}))+(1-p)*exp(2*m2+2*(s2(v))^{2})-(mm(v))^{2})^{1/2}
```


```R
phi <- function(x, m, s) (1/sqrt(2*pi*s^{2}))*exp((-1/(2*s^{2}))*(x-m)^{2})
```


```R

g <- function(x,v) p*phi(x,m1,s1(v)) +(1-p)*phi(x,m2,s2(v))
```


```R
f <- function(x,v) (1/x)*(g(log(x),v))
#f <- function(x,v) ifelse(x==0, 0, (1/x)*(g(log(x),v)))
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
CH <- function(v) (1/mm(v))*exp(he(v)-1)
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
options(repr.plot.width=4, repr.plot.height=4)
```


```R
 library("ggplot2")
g <- ggplot()+
  geom_path(data=df, aes(x=Cv, y=Ch), size=1.8, arrow= arrow(), color="blue1")
g1 <- g + coord_cartesian(xlim=c(0,6.5), ylim=c(0,6.5)) + scale_x_continuous(breaks = seq(0, 6.5, by = 2)) +
  scale_y_continuous(breaks = seq(0, 6.5, by = 3))
g1 + labs ( y=expression(C[h](T)), x=expression(C[V](T))) + theme_classic()
```


![png](output_20_0.png)



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
CHr <- function(v) (1/mm(v))*exp(hrv(v)-1)
```


```R
v <- seq(1.01,2, 0.01)
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
g <- ggplot()+
  geom_path(data=df, aes(x=Cv, y=Chr), size=1.8, arrow= arrow(), color="maroon")
g1 <- g + coord_cartesian(xlim=c(0,6), ylim=c(0,0.015)) + scale_x_continuous(breaks = seq(0, 6, by = 2)) +
  scale_y_continuous(breaks = seq(0, 0.015, by = 0.007))
g1 + labs ( y=expression(C[h](R)), x=expression(C[V](T))) + theme_classic()
```


![png](output_32_0.png)



```R
df2 <- data.frame("Ch1"=Ch, "Ch2"=Chr)
```


```R
 library("ggplot2")
g <- ggplot()+
  geom_path(data=df2, aes(x=Ch1, y=Ch2), size=1.8, arrow= arrow(), color="forestgreen")
g1 <- g + coord_cartesian(xlim=c(0,1), ylim=c(0,0.015)) + scale_x_continuous(breaks = seq(0, 1, by = 0.3)) +
  scale_y_continuous(breaks = seq(0, 0.015, by = 0.007))
g1 + labs ( y=expression(C[h](R)), x=expression(C[h](T))) + theme_classic()
```


![png](output_34_0.png)



```R

```
