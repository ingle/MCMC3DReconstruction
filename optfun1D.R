# obj function value for 2D case
# Feb 3, 2015
optfun1D <- function(x)
{
    norm(AB$A%*%x-sqrdat$fdatn,type='f')^2 + lam *norm(AB$B%*%x,type='f')^2
}
