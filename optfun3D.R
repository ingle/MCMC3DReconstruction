# obj function value for 3D case
# Feb 17, 2015
optfun3D <- function(x)
{
    norm(AB$A%*%x-ellipsoiddat$fdatn,type='f')^2 + lam *norm(AB$B%*%x,type='f')^2
}
