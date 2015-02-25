# Here we will use the ICM method shown here:
# http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/AV0809/ORCHARD/restore_image.html
# for full 3D reconstruction
set.seed(10)
source('makeAB3D.R')
source('getEllipsoidData.R')

library('VGAM')
library('Matrix')
 
a=1
b=2
c=3
limX=2
limY=3
limZ=4
N=5000
Nx=50
Ny=75
Nz=100

# ellipsoiddat <- getEllipsoidData(a=a,b=b,c=c,limX=limX,limY=limY,limZ=limZ,N=N)
# vec = list(xvec=seq(-limX,limX,,Nx), yvec=seq(-limY,limY,,Ny), zvec=seq(-limZ,limZ,,Nz))
# AB <- makeAB3D(vec=vec, dat=ellipsoiddat)
# lam=10
# sig=0.1
# ellipsoiddat$fdatn <- ellipsoiddat$fdat + rnorm(ellipsoiddat$fdat,sd=sig)
# 
# vec$fvec <- solve( t(AB$A)%*%AB$A + lam*(t(AB$B)%*%AB$B), t(AB$A)%*%ellipsdat$fdatn )
# vec$fvec <- matrix(vec$fvec, nrow=Ny)
# 
# cat('now creating initial state\n')
# x0 <- matrix(0, nrow=Ny, ncol=Nx)
# for ( z in seq(1,length(vec$zvec)) )
# {
#     for ( y in seq(1,length(vec$yvec)) )
#     {
#         for ( x in seq(1,length(vec$xvec)) ) 
#         {
#             curdist <- Inf
#             for ( dpt in seq(1,N) )
#             {
#                 dist <- (ellipsoiddat$xdat[dpt] - vec$xvec[x])^2 + 
#                         (ellipsoiddat$ydat[dpt] - vec$yvec[y])^2 +
#                         (ellipsoiddat$zdat[dpt] - vec$zvec[x])^2
#                 if (dist < curdist)
#                 {
#                     x0[y,x,z] <- ellipsoiddat$fdatn[dpt] 
#                     curdist <- dist
#                 }
#             }
#         }
#     }
#     cat('slice',z,' done\n')
# }
# cat('initial state created\n')

# x0init <- x0
x0 <- as.vector(x0init)

bestvec <- solve( t(AB$A)%*%AB$A + lam*(t(AB$B)%*%AB$B), t(AB$A)%*%ellipsdat$fdatn )
bestval <- optfun2D(bestvec)

# image(x=vec$xvec,y=vec$yvec,z=vec$fvec, asp=1.5)

source('optfun3D.R')
# mcmc loop
# x0 <- x0init
MAXITER <- 5000
updtcnt = 0

dx = mean( diff(vec$xvec) )
dy = mean( diff(vec$yvec) )

wt <- 0.01
    
lamdx2 <- wt/dx/dx
lamdy2 <- wt/dy/dy

sub2ind <- function( siz, r, c )
{
    (c-1)*siz[1] + r
}

sub2ind3 <- function( siz, r, c, s ) # 3dsize, row, col, slice
{
    r + (c-1) * siz[1] + (s-1) * siz[1]*siz[2]
}

alljj <- seq(1,Ny*Nx)
tmpy <- (alljj-1)%%Ny + 1
tmpx <- floor((alljj-1)/Ny)+1

ik <- as.vector( c( 
                    alljj[tmpy-1>=1],
                    alljj[tmpy+1<=Ny],
                    alljj[tmpx-1>=1],
                    alljj[tmpx+1<=Nx]
                  )
               )

siz <- c(Ny,Nx)
jk <- as.vector( c(
                       sub2ind(siz, tmpy[tmpy-1>=1]-1, tmpx[tmpy-1>=1]),
                       sub2ind(siz, tmpy[tmpy+1<=Ny]+1,tmpx[tmpy+1<=Ny]),
                       sub2ind(siz, tmpy[tmpx-1>=1], tmpx[tmpx-1>=1]-1),
                       sub2ind(siz, tmpy[tmpx+1<=Nx], tmpx[tmpx+1<=Nx]+1)
                  )
               )
xk <- as.vector( c( 
                       rep(2*lamdy2 ,sum(tmpy-1>=1)),
                       rep(2*lamdy2, sum(tmpy+1<=Ny)),
                       rep(2*lamdx2, sum(tmpx-1>=1)),
                       rep(2*lamdx2, sum(tmpx+1<=Nx))
                   )
               )

UpdtMat <- sparseMatrix( i=ik, j=jk, x=xk, dims=c(Ny*Nx,Ny*Nx) )

for ( i in seq(1,MAXITER) )
{
        
    x1 <- ( as.vector(x0init) + UpdtMat %*% x0 ) / (1 + 4*(lamdx2+lamdy2) )

    expo <- (norm(AB$A%*%x1-ellipsdat$fdatn,type='f')^2
                      -norm(AB$A%*%x0-ellipsdat$fdatn,type='f')^2)
            -lam    *(norm(AB$B%*%x1,type='f')^2
                      -norm(AB$B%*%x1,type='f')^2)
    logr = min(0, expo)
    u = runif(1)

    if ( log(u) < logr )
    {
        x0 <- x1
        updtcnt<-updtcnt+1
        #cat('updated ',updtcnt,' new fnval = ',optfun2D(x0),'\n')
    }
    #if( runif(1)<1/2000 )
    #{
    #    image(x=vec$xvec, y=vec$yvec, z=t(matrix(x0, nrow=Ny)))
    #    title(main=paste('iter ',i,'updated ',updtcnt,' new fnval = ',optfun2D(x0),'\n'))
    #}
    
    cat('image displayed ',i,' \n')
    image(x=vec$xvec, y=vec$yvec, z=t(matrix(x0, nrow=Ny)))
    title(main=paste('iter ',i,'updated ',updtcnt,'wt=', wt,'\n',
                         'new fnval = ',optfun2D(x0),'\n',
                         'bestval = ', bestval, 'nnbval = ', optfun2D(as.vector(x0init)) ))

}

