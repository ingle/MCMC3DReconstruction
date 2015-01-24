# We will use the ICM method shown here:
# http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/AV0809/ORCHARD/restore_image.html
set.seed(10)
#set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )
source('makeAB.R')
source('getSquareWaveData.R')

library('VGAM')

a=1
limX=2
N=50
Nx=500

sqrdat <- getSquareWaveData(a=a,limX=limX,N=N)
vec = list( xvec=seq(-limX,limX,,Nx) )
AB <- makeAB_1D(vec=vec, dat=sqrdat)
# library('SparseM')
# SparseM::image(t(AB$A)%*%AB$A, lwd=2)
# SparseM::image( AB$B, lwd=2 )


lam=1000
sig=0.01
sqrdat$fdatn <- sqrdat$fdat + rnorm(sqrdat$fdat,sd=sig)

library('Matrix')
M <- t(AB$A)%*%AB$A + lam*(t(AB$B)%*%AB$B)
c <- t(AB$A)%*%sqrdat$fdatn
vec$fvec <- solve( M, c )
vecdf = data.frame( xvec=vec$xvec, fvec=as.vector(vec$fvec) )

library('ggplot2')

gp <- ggplot( sqrdat, aes(x=xdat, y=fdat) ) + geom_point()
gp <- gp + geom_point( aes(x=xdat, y=fdatn), color='red' )
gp <- gp + geom_line ( data=vecdf, aes(x=xvec, y=fvec ), color='green' )

cat('now creating initial state\n')
x0 <- rep(0, Nx) 
for ( x in seq(1,length(vec$xvec)) ) 
{
    curdist <- Inf
    for ( dpt in seq(1,N) )
    {
        dist <- abs(sqrdat$xdat[dpt] - vec$xvec[x])
        if (dist < curdist)
        {
            x0[x] <- sqrdat$fdatn[dpt] 
            curdist <- dist
        }
    }
}

cat('initial state created\n')

x0init = as.vector(x0)

gp <- gp + geom_line( data=vecdf, aes(x=xvec, y=x0init), color='blue' )
print(gp)

# 
# mcmc loop
#
source('optfun.R')

#x0 <- runif(n=Nx,min=0, max=2)
x0 <- x0init
MAXITER <- 100
updtcnt = 0
for ( i in seq(1,MAXITER) )
{
    for (randind in sample(seq(1,Nx)))
    {
        #xr <- rnorm(n=Ny*Nx, mean=x0, sd=1.5)
        #randind <- sample.int(Ny*Nx, size=2)
        
        nbd = 200
        indvec <- seq( max(1,randind-nbd), min(Nx,randind+nbd) )
        subvec <- solve( M[indvec,indvec], c[indvec] )
        x1 <- x0

        nbdwrite = 20 
        # only replace a subset around the center in indvec
        if ( randind-nbdwrite>=1 && randind+nbdwrite <= Nx )
        {
            subindvec <- seq( randind-nbdwrite, randind+nbdwrite )
            x1[subindvec] <- subvec[seq(length(subvec)/2-nbdwrite, length(subvec)/2+nbdwrite)]
            #x1[indvec] <- rlaplace(n=length(subindvec), 
            #            loc = as.vector(subvec[seq(length(subvec)/2-nbd/2, length(subvec)/2+nbd/2)]),
            #            scale=0.1/log(exp(i)) )
        }
        else
        {
            x1[indvec] <- as.vector(subvec)
            #x1[indvec] <- rlaplace(n=length(indvec), loc=as.vector(subvec), scale=0.1/log(exp(i)) )
        }
        expo <- -1/sig^2*(norm(AB$A%*%x1-sqrdat$fdatn,type='f')^2
                          -norm(AB$A%*%x0-sqrdat$fdatn,type='f')^2)
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
        if( runif(1)<1/200 )
        {
            plot( x=vec$xvec, y=x0 )
            #image(x=vec$xvec, y=vec$yvec, z=t(matrix(x0, nrow=Ny)))
            title(main=paste('iter ',i,'updated ',updtcnt,' new fnval = ',optfun(x0),'\n'))
        }
    }
}
# 
# x0 <- t(matrix(x0, nrow=Ny))
# image(x=vec$xvec, y=vec$yvec, z=x0)
