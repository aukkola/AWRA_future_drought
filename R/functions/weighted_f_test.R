ratio <- 1

var_x <- wtd.var(x)
var_y <- wtd.var(y)


x.root=1
y.root=1


ESTIMATE <- var_x / var_y
STATISTIC <- ESTIMATE / ratio

DF.x <- 

PVAL <- pf(STATISTIC, DF.x, DF.y)





x <- rnorm(50)
y <- rnorm(50)

res.x <- wle.normal(x,group=5)
res.y <- wle.normal(y,group=5)

wle.var.test(res.x, res.y) # Do x and y have the same variance?
