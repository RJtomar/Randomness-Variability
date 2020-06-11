## Varying the shape parameters


```R
p= 0.5
```


```R
b= function(a) ((1-p)*a)/(m*a-p)
```


```R
m <- 100
```


```R
f <- function(t,a) p*a*exp(-a*t)+(1-p)*b(a)*exp(-b(a)*t)
```


```R
logf <- function(t,a) ifelse( t > 100, 0, log(f(t,a)))
```


```R
f1 <- function(t,a) f(t,a)*logf(t,a)
```


```R
h <- function(a)  -integrate(f1, 0, 190, a=a)$value
```


```R
h_f <- Vectorize(h)
```


```R
Cv <- function(a) sqrt(((2*p*(b(a))^(2)+2*(1-p)*a^(2))/(p*b(a)+(1-p)*a)^(2))-1)
```


```R
plot(a, b(a), col="red", type="p")
```


![png](output_10_0.png)



```R
Ch <- function(a) (1/m)*exp(h_f(a)-1)
```


```R
a <- seq(0.01, 2, by=0.001)
```


```R
CH <- Ch(a)
which.max(CH)
```


20



```R
a[20]
```


0.029



```R
Ch(0.029)
```


0.105758280095732



```R
CV <- Cv(a)
```


```R
Cv(0.029)
```


1.36326878626005



```R
CV[20]
```


1.36326878626005



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
g1 <- g + coord_cartesian(xlim=c(0.9, 1.9), ylim=c(0, 0.12)) + scale_x_continuous(breaks = seq(0.9, 1.9, by = 0.3)) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.05))
g1 + labs ( y=expression(C[h](T)), x=expression(C[V](T))) + theme_classic()
```


![png](output_21_0.png)



```R
fR <- function(t,a) (1/m)*(1/t^{3})*f(1/t,a)
```


```R
logfR <- function(t, a) ifelse((1/t) > 500, 0, -log(m)-3*log(t)+log(f(1/t,a)))
```


```R
FR <- function(t,a) fR(t,a)*logfR(t,a)
```


```R
H <- function(a) -(integrate(FR, lower=0, upper=Inf, a=a)$value)
```


```R
HH <- Vectorize(H)
```


```R
ChR <- function(a) (1/m)*(exp(HH(a)-1))
```


```R
CHR <- ChR(a)
```


```R

```


1.6559656464499



```R
df1 <- data.frame("df11"=CV, "df12"=CHR)
```


```R
g <- ggplot()+
  geom_path(data=df1, aes(x=df11, y=df12), size=1.8, arrow= arrow(), color="maroon")
g1 <- g + coord_cartesian(xlim=c(0.9,1.8), ylim=c(5e-05,22e-05)) + scale_x_continuous(breaks = seq(0.9, 1.8, by = 0.3)) +
  scale_y_continuous(breaks = seq(5e-05, 22e-5, by = 8e-5))
g1 + labs ( y=expression(C[h](R)), x=expression(C[V](T))) + theme_classic() + theme(axis.text.y=element_text(angle=50, size=6, vjust=.5))
```


![png](output_31_0.png)



```R
df2 <- data.frame("Ch1"=CH, "Ch2"=CHR)
```


```R
 library("ggplot2")
g <- ggplot()+
  geom_path(data=df2, aes(x=Ch1, y=Ch2), size=1.8, arrow= arrow(), color="forestgreen")
g1 <- g + coord_cartesian(xlim=c(0,0.12), ylim=c(0.00007,0.00022))+ scale_x_continuous(breaks = seq(0, 0.12, by = 0.04)) +
  scale_y_continuous(breaks = seq(0.00007, 0.00022, by = 0.00007))
g1 + labs ( y=expression(C[h](R)), x=expression(C[h](T))) + theme_classic() + theme(axis.text.y=element_text(angle=50, size=6, vjust=.5))
```


![png](output_33_0.png)



```R

```
