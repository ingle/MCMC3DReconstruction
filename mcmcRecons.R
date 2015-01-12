set.seed(10)
source('makeAB.R')
source('getEllipseData.R')

a=1
b=2
limX=2
limY=3
N=1000
Nx=100
Ny=200

ellipsdat <- getEllipseData(a=a,b=b,limX=limX,limY=limY,N=N)
vec = list(xvec=seq(-limX,limX,,Nx), yvec=seq(-limY,limY,,Ny))
AB <- makeAB(vec=vec, dat=ellipsdat)
lam=0.1
sig=0.1
ellipsdat$fdatn <- ellipsdat$fdat + rnorm(ellipsdat$fdat,sd=sig)

library('Matrix')
vec$fvec <- solve( t(AB$A)%*%AB$A + lam*(t(AB$B)%*%AB$B), t(AB$A)%*%ellipsdat$fdatn )
vec$fvec <- t(matrix(vec$fvec, nrow=Ny))


# image(x=vec$xvec,y=vec$yvec,z=vec$fvec, asp=1.5)

# mcmc loop
# x0 <- runif(n=Ny*Nx,min=0, max=2)
# MAXITER <- 1000
# updtcnt = 0
# for ( i in seq(1,MAXITER) )
# {
#     xr <- rnorm(n=Ny*Nx, mean=x0, sd=1.5)
#     randind <- sample.int(Ny*Nx, size=100)
#     x1 <- x0
#     x1[randind] <- xr[randind]
#     expo <- -1/sig^2*(norm(AB$A%*%x1-ellipsdat$fdatn,type='f')^2
#                       -norm(AB$A%*%x0-ellipsdat$fdatn,type='f')^2)
#             -lam    *(norm(AB$B%*%x1,type='f')^2
#                       -norm(AB$B%*%x1,type='f')^2)
#     logr = min(0, expo)
#     u = runif(1)
#     if (log(u)<logr)
#     {
#         x0<-x1
#         updtcnt<-updtcnt+1
#         cat('updated ',updtcnt,'\n')
#     }
# }

# x0 <- t(matrix(x0, nrow=Ny))
# image(x=vec$xvec, y=vec$yvec, z=x0)
