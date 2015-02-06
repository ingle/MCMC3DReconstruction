# obj function value for 2D case
# Jan 30, 2015
optfun2D <- function(x)
{
    norm(AB$A%*%x-ellipsdat$fdatn,type='f')^2 + lam *norm(AB$B%*%x,type='f')^2
}
