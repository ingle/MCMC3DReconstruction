# obj function value for 1D case
# Jan 23, 2015
optfun <- function(x)
{
    1/sig^2*norm(AB$A%*%x-sqrdat$fdatn,type='f')^2
            lam    *norm(AB$B%*%x,type='f')^2
}
