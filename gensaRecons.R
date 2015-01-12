# make sure data from mcmcRecons is in workspace
#
library(GenSA)

optfun <- function(x)
{
    1/sig^2*norm(AB$A%*%x-ellipsdat$fdatn,type='f')^2
            lam    *norm(AB$B%*%x,type='f')^2
}

set.seed(10)
lower <- rep(0,Ny*Nx)
upper <- rep(10,Ny*Nx)

#xoptgsa <- GenSA( par=as.vector(x0init), lower=lower, upper=upper, fn=optfun,
#               control=list(verbose=TRUE,max.time=10) ) 
# par = matrix(0, nrow=Ny, ncol=Nx)
# par[seq(Ny/2-30,Ny/2+30), seq(Nx/2-20,Nx/2+20)]=2
# par = as.vector(par)
# 
xoptsa <- optim( par = as.vector(x0init), fn=optfun,
                 method="SANN", control=list(REPORT=1,maxit=50000) )

xopt <- t(matrix(xoptgsa$par, nrow=Ny))
image(x=vec$xvec, y=vec$yvec, z=xopt)  # note that image() plots data transposed

