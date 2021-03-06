# We will use the ICM method shown here:
# http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/AV0809/ORCHARD/restore_image.html

# In addition to an L2 penalty on the 2nd derivative, we will impose L1 penalty
# on the first derivaitve
# Feb 25, 2015

set.seed(10)
#set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )
source('makeAB1D.R')
source('getSquareWaveData.R')

library('VGAM')

a=1
limX=2
N=50
Nx=500

sqrdat <- getSquareWaveData(a=a,limX=limX,N=N)
vec = list( xvec=seq(-limX,limX,,Nx) )
AB <- makeAB1D(vec=vec, dat=sqrdat)
# library('SparseM')
# SparseM::image(t(AB$A)%*%AB$A, lwd=2)
# SparseM::image( AB$B, lwd=2 )


lam=1000
sig=0.1
sqrdat$fdatn <- sqrdat$fdat + rlaplace(sqrdat$fdat,scale=sig)
# sqrdat$fdatn <- sqrdat$fdat + rnorm(sqrdat$fdat,sd=sig)

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
source('optfun1D.R')
source('truefun1D.R')

bestval <- optfun1D( as.vector(vec$fvec) )
nnbval  <- optfun1D( as.vector(x0init) )

#x0 <- runif(n=Nx,min=0, max=2) <<< this is stupid, the first few loops are wasted converging to nnb soln
x0 <- x0init    # smarter to init using nnb solution so we can directly start from big wt
MAXITER <- 5000
updtcnt = 0

# Let's form the vector update step matrix, which has two off diagonals of all 1's and diagonal of 0
ri = c( seq(2,length(x0)  ),  seq(1,length(x0)-1) )
ci = c( seq(1,length(x0)-1),  seq(2,length(x0)  ) )
vali = rep( 1,length(ri) )

UpdtMat <- sparseMatrix( i=ri, j=ci, x=vali, dims=c(length(x0), length(x0)) )

for ( i in seq(1,MAXITER) )
{
    #readline('');
    #wt <- 0.0001 * (i<30) + 100000 * (i>=30)
    wt1 <- 0.0000
    wt2 <- 1e-6

    dx <- mean( diff( vecdf$xvec ) )

    aa <- 2 + 8*wt1/dx^2
    xjav <- 0.5 * UpdtMat %*% x0
    bb <- 8 * wt1/dx^2 * xjav + 2*as.vector(x0init)

    #xnew1 <- 1/aa * (bb + 2*wt2/dx^2)
    #xnew2 <- 1/aa * (bb - 2*wt2/dx^2)
    #xnew3 <- 1/aa * (bb)
    
    #x1 <- xnew1 * (bb < aa*xjav - 2*wt2/dx^2) + xnew2 * (bb >= aa*xjav + 2*wt2/dx^2) +
    #      xnew3 * ( (aa*xjav-2*wt2/dx^2 <= bb) & (bb < aa*xjav+2*wt2/dx^2) )

    # Try simplified case with wt1=0. We get the shrinkage estimator.

    # xnew1 <- x0init + wt2/dx^2
    # xnew2 <- x0init - wt2/dx^2
    xnew1 <- x0 + wt2/dx^2
    xnew2 <- x0 - wt2/dx^2

    x1 <- xnew1 * (x0 < 0.5*(UpdtMat%*%x0) - wt2/dx^2) +
          xnew2 * (x0 > 0.5*(UpdtMat%*%x0) + wt2/dx^2) +
          x0 * ((x0 > 0.5*(UpdtMat%*%x0) - wt2/dx^2) & (x0 < 0.5*(UpdtMat%*%x0) + wt2/dx^2))

 
    # x1 <- xnew1 * (x0init < 0.5*(UpdtMat%*%x0) - wt2/dx^2) +
    #       xnew2 * (x0init > 0.5*(UpdtMat%*%x0) + wt2/dx^2) +
    #       x0init * ((x0init > 0.5*(UpdtMat%*%x0) - wt2/dx^2) & (x0init < 0.5*(UpdtMat%*%x0) + wt2/dx^2))

    # x1 <- (dx^2 * x0init + 2*wt * UpdtMat %*% x0) / (4*wt+dx^2) 

    # x1 <- x0
    # dif <- abs(vecdf$xvec[randind] - sqrdat$xdat)
    # nearDat <- sqrdat$fdatn[ which(dif==min(dif)) ]
    # stopifnot( length(nearDat)==1 )

    # x1[randind] = ( dx^2 * nearDat + 2*wt * ( x0[min(Nx,randind+1)]+x0[max(1,randind-1)] ))/(4*wt + dx^2)

    expo <- - (norm(AB$A%*%x1-sqrdat$fdatn,type='f')^2
                      -norm(AB$A%*%x0-sqrdat$fdatn,type='f')^2) -
                      lam * (norm(AB$B%*%x1,type='f')^2
                      -norm(AB$B%*%x1,type='f')^2)
    logr = min(0, expo)
    u = runif(1)

    #if (log(u)<logr & optfun(x1) < optfun(x0))
    #if ( log(u) < logr )
    if( 1 )
    #if ( log(u) < logr )
    {
        x0<-x1
        updtcnt<-updtcnt+1
        cat('updated ',updtcnt,' new fnval = ',optfun1D(x0),'\n')
    }
    if( runif(1)<=1 )
    {
        plot( x=vec$xvec, y=x0, ylim=c(-1e-2,1.1) )
        #image(x=vec$xvec, y=vec$yvec, z=t(matrix(x0, nrow=Ny)))
        err <- 1/Nx * sum((x0-truefun1D(vec$xvec))^2)
        errnnb <- 1/Nx * sum((as.vector(x0init)-truefun1D(vec$xvec))^2)
        title(main=paste('iter ',i,'updated ',updtcnt,'wt=', wt,'\n',
                         'new fnval = ',optfun1D(x0),'\n',
                         'err = ', err, 'errnnb = ', errnnb))
                         #'bestval = ', bestval, 'nnbval = ', nnbval,
    }

    # gp <- ggplot( sqrdat, aes(x=xdat, y=fdat) ) + geom_point()
    # gp <- gp + geom_line( aes(x=vecdf$xvec, y=x0init), color='black' )
    # #gp <- gp + geom_line( aes(x=vecdf$xvec, y=x0init+wt2/dx^2), color='red' )
    # #gp <- gp + geom_line( aes(x=vecdf$xvec, y=x0init-wt2/dx^2), color='blue' )
    # gp2 <- gp + geom_point( aes(x=vecdf$xvec, y=1*as.vector(x0init > 0.5*(UpdtMat%*%x0)+wt2/dx^2)), color='red' )
    # gp2 <- gp2 + geom_point( aes(x=vecdf$xvec, y=1*as.vector(x0init < 0.5*(UpdtMat%*%x0)-wt2/dx^2)), color='blue' )
    # print(gp2)
}
# 
# x0 <- t(matrix(x0, nrow=Ny))
# image(x=vec$xvec, y=vec$yvec, z=x0)
