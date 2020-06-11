# Varying the probability


```R
a=5
```


```R
b=1.25
```


```R
f <- function(t,p) p*a*exp(-a*t)+(1-p)*b*exp(-b*t)
```


```R
logf <- function(t,p) ifelse( t > 500, 0, log(f(t,p)))
```


```R
f1 <- function(t,p) f(t,p)*logf(t,p)
```


```R
h <- function(p)  -integrate(Vectorize(f1), 0, Inf, p=p)$value
```


```R
h_f <- Vectorize(h)
```


```R
m <- function(p) (p*b+(1-p)*a)/(a*b)
```


```R
plot(p, m(p), col="red", type="p")
```


![png](output_9_0.png)



```R
Cv <- function(p) sqrt(((2*p*b^(2)+2*(1-p)*a^(2))/(p*b+(1-p)*a)^(2))-1)
```


```R
p[80]
```


0.8



```R
Cv(0.8)
```


1.45773797371133



```R
Ch <- function(p) (1/m(p))*exp(h_f(p)-1)
```


```R
p <- seq(0.01, 0.99, 0.01)
```


```R
CV <- Cv(p)
```


```R
CH <- Ch(p)
```


```R
df <- data.frame("en1"=CV, "en2"=CH)
```


```R
library("repr")
options(repr.plot.width=4, repr.plot.height=4)
```


```R
 library("ggplot2")
g <- ggplot()+
  geom_path(data=df, aes(x=en1, y=en2), size=1.8, arrow= arrow(), col="blue1")
g1 <- g + coord_cartesian(xlim=c(1,1.5), ylim=c(0.92,1.02)) + scale_x_continuous(breaks = seq(1, 1.5, by = 0.15)) +
  scale_y_continuous(breaks = seq(0.92, 1.02, by = 0.05))
g1 + labs ( y=expression(C[h](T)), x=expression(C[V](T))) + theme_classic()
```


![png](output_19_0.png)



```R
fR <- function(t,p) (1/m(p))*(1/t^{3})*f(1/t,p)
```


```R
logfR <- function(t, p) ifelse((1/t) > 500, 0, -log(m(p))-3*log(t)+log(f(1/t,p)))
```


```R
FR <- function(t,p) fR(t,p)*logfR(t,p)
```


```R
H <- function(p) -(integrate(FR, lower=0, upper=Inf, p=p)$value)
```


```R
HH <- Vectorize(H)
```


```R
ChR <- function(p) (1/m(p))*(exp(HH(p)-1))
```


```R
CHR <- ChR(p)
```


```R
df1 <- data.frame("df11"=CV, "df12"=CHR)
```


```R
g <- ggplot()+
  geom_path(data=df1, aes(x=df11, y=df12), size=1.8, arrow= arrow(), color="maroon")
g1 <- g + coord_cartesian(xlim=c(1,1.6), ylim=c(1,20)) + scale_x_continuous(breaks = seq(1, 1.6, by = 0.2)) +
  scale_y_continuous(breaks = seq(1, 20, by = 9))
g1 + labs ( y=expression(C[h](R)), x=expression(C[V](T))) + theme_classic()
```


![png](output_28_0.png)



```R
df2 <- data.frame("Ch1"=CH, "Ch2"=CHR)
```


```R
 library("ggplot2")
g <- ggplot()+
  geom_path(data=df2, aes(x=Ch1, y=Ch2), size=1.8, arrow= arrow(), color="forestgreen")
g1 <- g + coord_cartesian(xlim=c(0.5,1.5), ylim=c(1,20))+ scale_x_continuous(breaks = seq(0.5, 1.5, by = 0.3)) +
  scale_y_continuous(breaks = seq(1, 20, by = 9))
g1 + labs ( y=expression(C[h](R)), x=expression(C[h](T))) + theme_classic()
```


![png](output_30_0.png)



```R

```
