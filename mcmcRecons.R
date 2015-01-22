# set.seed(10)
# source('makeAB.R')
# source('getEllipseData.R')
#
library('VGAM')
# 
# a=1
# b=2
# limX=2
# limY=3
# N=1000
# Nx=100
# Ny=200
# 
# ellipsdat <- getEllipseData(a=a,b=b,limX=limX,limY=limY,N=N)
# vec = list(xvec=seq(-limX,limX,,Nx), yvec=seq(-limY,limY,,Ny))
# AB <- makeAB(vec=vec, dat=ellipsdat)
# lam=10
# sig=0.1
# ellipsdat$fdatn <- ellipsdat$fdat + rnorm(ellipsdat$fdat,sd=sig)
# 
# library('Matrix')
# vec$fvec <- solve( t(AB$A)%*%AB$A + lam*(t(AB$B)%*%AB$B), t(AB$A)%*%ellipsdat$fdatn )
# vec$fvec <- matrix(vec$fvec, nrow=Ny)
# 
# cat('now creating initial state\n')
# x0 <- matrix(0, nrow=Ny, ncol=Nx)
# for ( y in seq(1,length(vec$yvec)) )
# {
#     for ( x in seq(1,length(vec$xvec)) ) 
#     {
#         curdist <- Inf
#         for ( dpt in seq(1,N) )
#         {
#             dist <- (ellipsdat$xdat[dpt] - vec$xvec[x])^2 + (ellipsdat$ydat[dpt] - vec$yvec[y])^2
#             if (dist < curdist)
#             {
#                 x0[y,x] <- ellipsdat$fdatn[dpt] 
#                 curdist <- dist
#             }
#         }
#     }
#     cat('row ',y,' done\n')
# }
# cat('initial state created\n')

x0 = as.vector(x0init)

# image(x=vec$xvec,y=vec$yvec,z=vec$fvec, asp=1.5)

# mcmc loop
x0 <- runif(n=Ny*Nx,min=0, max=2)
MAXITER <- 100
updtcnt = 0
for ( i in seq(1,MAXITER) )
{
    for (randind in sample(seq(1,Ny*Nx)))
    {
        #xr <- rnorm(n=Ny*Nx, mean=x0, sd=1.5)
        #randind <- sample.int(Ny*Nx, size=2)
        
        nbd = 1000
        indvec <- seq( max(1,randind-nbd), min(Ny*Nx,randind+nbd) )
        subvec <- solve( M[indvec,indvec], c[indvec] )
        x1 <- x0
        # only replace a subset around the center in indvec
        if ( randind-nbd/2 >=1 && randind+nbd/2 <= Ny*Nx )
        {
            subindvec <- seq( randind-nbd/2, randind+nbd/2 )
            #x1[subindvec] <- subvec[seq(length(subvec)/2-nbd/2, length(subvec)/2+nbd/2)]
            x1[indvec] <- rlaplace(n=length(subindvec), 
                        loc = as.vector(subvec[seq(length(subvec)/2-nbd/2, length(subvec)/2+nbd/2)]),
                        scale=0.1/log(exp(i)) )
        }
        else
        {
            #x1[indvec] <- as.vector(subvec)
            x1[indvec] <- rlaplace(n=length(indvec), loc=as.vector(subvec), scale=0.1/log(exp(i)) )
        }
        expo <- -1/sig^2*(norm(AB$A%*%x1-ellipsdat$fdatn,type='f')^2
                          -norm(AB$A%*%x0-ellipsdat$fdatn,type='f')^2)
                -lam    *(norm(AB$B%*%x1,type='f')^2
                          -norm(AB$B%*%x1,type='f')^2)
        logr = min(0, expo)
        u = runif(1)

        #if (log(u)<logr & optfun(x1) < optfun(x0))
        if ( log(u) < logr )
        {
            x0<-x1
            updtcnt<-updtcnt+1
            cat('updated ',updtcnt,' new fnval = ',optfun(x0),'\n')
        }
        if( runif(1)<1/2000 )
        {
            image(x=vec$xvec, y=vec$yvec, z=t(matrix(x0, nrow=Ny)))
            title(main=paste('iter ',i,'updated ',updtcnt,' new fnval = ',optfun(x0),'\n'))
        }
    }
}

x0 <- t(matrix(x0, nrow=Ny))
image(x=vec$xvec, y=vec$yvec, z=x0)
