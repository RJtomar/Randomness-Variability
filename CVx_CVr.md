```R
library("pracma")
library("expint")
```

    
    Attaching package: ‘expint’
    
    The following objects are masked from ‘package:pracma’:
    
        expint, expint_E1, expint_Ei, gammainc
    



```R
CV_R1p <- function(X) X/sqrt(1-(X^{2}))
```


```R
CV_R2p <- function(X) X
```


```R
CV_R3p <- function(X) X
```


```R

```


```R
CV_R4p <- function(X){
    sapply(X, function(x) if (x<1) sqrt(((exp((1/x)-1)*expint_E1((1/x)-1, scale = FALSE))/x)-1) else NaN)
}
```


```R
CV_R5p <- function(X) sqrt((2/(1+X^(2)))-1)
```


```R
CV_X <- seq(0.01,2.5, by=0.01)
```


```R
CV_X2 <- seq(0.01, 2.5, by=0.1)
```


```R
CV_R1 <- CV_R1p(CV_X)
```

    Warning message in sqrt(1 - (X^{:
    “NaNs produced”


```R
CV_R2 <- CV_R2p(CV_X)
```


```R
CV_R3 <- CV_R3p(CV_X)
```


```R
CV_R4 <- CV_R4p(CV_X)
```


```R
CV_R5 <- CV_R5p(CV_X)
```

    Warning message in sqrt((2/(1 + X^(2))) - 1):
    “NaNs produced”


```R
P <- data.frame("ISI"=CV_X, "Gamma"=CV_R1, "LogN"= CV_R2, "IGau"= CV_R3, "SExp" = CV_R4, "new"=CV_R5)
```


```R
require(ggplot2)
require(reshape2)
```


```R
df <- melt(P ,  id.vars = 'ISI', variable='Distributions')
```


```R
library("repr")
options(repr.plot.width=5, repr.plot.height=3.5)
```


```R
library("ggplot2")
g <- ggplot(df, aes(x= ISI, y = value, col= Distributions, linetype=Distributions)) + geom_line(lwd=1) + 
  scale_linetype_manual(values = c('solid', 'solid', 'dashed', 'solid', 'solid'), labels= c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted exponential", "NEW") ) + scale_color_manual(name = "Distributions", 
  labels = c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted exponential", "NEW"), values= c('red', 'green', 'blue', 'violet', 'brown')) +
  guides(linetype=FALSE) 
g1 <- g + coord_cartesian(xlim=c(0,1.5), ylim=c(0,1.5))
g1 + labs ( y=expression(C[V](R)), x=expression(C[V](T))) + theme_classic() 

```

    Warning message:
    “Removed 451 rows containing missing values (geom_path).”


![png](output_18_1.png)



```R
library("ggplot2")
ggplot(df, aes(x= ISI, y = value, col= Distributions, linetype=Distributions)) + geom_line(lwd=1)  + scale_color_hue(name = "Distributions", labels = c("T999", "T888", "yy", "cc")) + 
scale_shape_discrete(name="", labels= c("T999", "T888", "yy", "cc"))  + 
scale_linetype_manual(values= c('solid', 'solid', 'dashed', 'solid')) + guides(linetype = FALSE)
```

    Warning message:
    “Removed 301 rows containing missing values (geom_path).”


![png](output_19_1.png)



```R
library("ggplot2")
g <- ggplot(df, aes(ISI,value)) + geom_line(aes(colour = Distributions)) + scale_linetype_manual(values = c(rep("solid", 2), rep("dashed", 2)))
g1 <- g + coord_cartesian(xlim=c(0,1.5), ylim=c(0,1.5))
g1 + labs ( y=expression(C[V](R)), x=expression(C[V](T))) + theme_classic() + scale_color_manual(name="distributions",
                                                                                                labels=c("Gamma",
                                                                                                        "Lognormal",
                                                                                                        "Inverse Gaussian",
                                                                                                        "Shifted exponential"),
                                                                                                values = c("Gamma"="blue",
                                                                                                           "LogN"="red",
                                                                                                           "IGau"="yellow",
                                                                                                           "SExp"="green"))

```

    Warning message:
    “Removed 301 rows containing missing values (geom_path).”


![png](output_20_1.png)



```R
#g1 + labs ( y=expression(C[V](R)), x=expression(C[V](T))) + theme_classic() + scale_color_manual(name="distributions",
                                                                                                labels=c("Gamma",
                                                                                                        "Lognormal",
                                                                                                        "Inverse Gaussian",
                                                                                                        "Shifted exponential"),
                                                                                                values = c("Gamma"="blue",
                                                                                                           "LogN"="red",
                                                                                                           "IGau"="yellow",
                                                                                                           "SExp"="green"))
```


```R

```


    Error: Cannot add ggproto objects together. Did you forget to add this object to a ggplot object?
    Traceback:


    1. coord_cartesian(xlim(c(0, 2.5)) + ylim(c(0, 2)))

    2. ggproto(NULL, CoordCartesian, limits = list(x = xlim, y = ylim), 
     .     expand = expand, default = default, clip = clip)

    3. `+.gg`(xlim(c(0, 2.5)), ylim(c(0, 2)))

    4. stop("Cannot add ggproto objects together.", " Did you forget to add this object to a ggplot object?", 
     .     call. = FALSE)



```R


plot(CV_X, CV_R1(CV_X), xlim=c(0,2.5), ylim=c(0,2.0),  xaxs='i', yaxs='i', type="o", pch=20, cex=0.4, col="red", xlab="CV(T)", ylab="CV(R)")
lines(CV_X, CV_R2(CV_X), type="o", pch=20, cex=0.4, lty=3, col="blue3")
lines(CV_X, CV_R3(CV_X), type="o", pch=20, cex=0.4, lty=3, col="green4")
lines(CV_X, CV_R4(CV_X), type="o", pch=20, cex=0.4, lty=3, col="brown")
legend(1.9, 1.8, legend=c("Gamma", "Lognormal", "Inverse Gaussian", "Shifted Exp."),
       col=c("red", "blue3","green4", "brown"), lty=c(1,2,5,3), cex=0.65)
```

    Warning message in sqrt(1 - (CV_X^{:
    “NaNs produced”


![png](output_23_1.png)



```R

```


```R

```
