# Here we will use the ICM method shown here:
# http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/AV0809/ORCHARD/restore_image.html
# for full 3D reconstruction
set.seed(10)
# source('makeAB3D.R')
source('getEllipsoidData.R')

library('VGAM')
library('Matrix')
library('RANN')
 
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

ellipsoiddat <- getEllipsoidData(a=a,b=b,c=c,limX=limX,limY=limY,limZ=limZ,N=N)
vec = list(xvec=seq(-limX,limX,,Nx), yvec=seq(-limY,limY,,Ny), zvec=seq(-limZ,limZ,,Nz))
# AB <- makeAB3D(vec=vec, dat=ellipsoiddat)
# lam=10
sig=0.1
ellipsoiddat$fdatn <- ellipsoiddat$fdat + rnorm(ellipsoiddat$fdat,sd=sig)
# 
# vec$fvec <- solve( t(AB$A)%*%AB$A + lam*(t(AB$B)%*%AB$B), t(AB$A)%*%ellipsdat$fdatn )
# vec$fvec <- matrix(vec$fvec, nrow=Ny)
# 
cat('now creating initial state\n')
grid <- expand.grid( vec$yvec, vec$xvec, vec$zvec )
cloud <- cbind( ellipsoiddat$ydat, ellipsoiddat$xdat, ellipsoiddat$zdat )
nnres <- nn2( data = cloud, query = grid, k=1 )
x0init = ellipsoiddat$fdatn[nnres$nn.idx]
cat('initial state created\n')

# x0init <- x0
x0 <- as.vector(x0init)

# image(x=vec$xvec,y=vec$yvec,z=vec$fvec, asp=1.5)

# mcmc loop
# x0 <- x0init
MAXITER <- 5000
updtcnt = 0

dy = mean( diff(vec$yvec) ) # row
dx = mean( diff(vec$xvec) ) # col
dz = mean( diff(vec$zvec) ) # sli

wt <- 0.1
    
lamdy2 <- wt/dy/dy
lamdx2 <- wt/dx/dx
lamdz2 <- wt/dz/dz

sub2ind2d <- function( siz, r, c )
{
    (c-1)*siz[1] + r
}

sub2ind3d <- function( siz, r, c, s ) # 3dsize, row, col, slice
{
    r + (c-1) * siz[1] + (s-1) * siz[1]*siz[2]
}

alljj <- seq(1,Ny*Nx*Nz)
tmpy <- ((alljj-1)%%(Ny*Nx))%%Ny + 1
tmpx <- floor(((alljj-1)%%(Ny*Nx))/Ny) + 1
tmpz <- floor(floor((alljj-1)/Ny)/Nx) + 1

ik <- as.vector( c( 
                    alljj[tmpy-1>=1 ],
                    alljj[tmpy+1<=Ny],
                    alljj[tmpx-1>=1 ],
                    alljj[tmpx+1<=Nx],
                    alljj[tmpz-1>=1 ],
                    alljj[tmpz+1<=Nz]
                  )
               )

siz <- c(Ny,Nx,Nz)
jk <- as.vector( c(
                       sub2ind3d(siz, tmpy[tmpy-1>=1 ]-1 , tmpx[tmpy-1>=1 ]   , tmpz[tmpy-1>=1 ]   ),
                       sub2ind3d(siz, tmpy[tmpy+1<=Ny]+1 , tmpx[tmpy+1<=Ny]   , tmpz[tmpy+1<=Ny]   ),
                       sub2ind3d(siz, tmpy[tmpx-1>=1 ]   , tmpx[tmpx-1>=1 ]-1 , tmpz[tmpx-1>=1 ]   ),
                       sub2ind3d(siz, tmpy[tmpx+1<=Nx]   , tmpx[tmpx+1<=Nx]+1 , tmpz[tmpx+1<=Nx]   ),
                       sub2ind3d(siz, tmpy[tmpz-1>=1 ]   , tmpx[tmpz-1>=1 ]   , tmpz[tmpz-1>=1]-1 ),
                       sub2ind3d(siz, tmpy[tmpz+1<=Nz]   , tmpx[tmpz+1<=Nz]   , tmpz[tmpz+1<=Nz]+1 )
                   )
               )
xk <- as.vector( c( 
                       rep( 2*lamdy2 ,sum(tmpy-1>=1 ) ),
                       rep( 2*lamdy2, sum(tmpy+1<=Ny) ),
                       rep( 2*lamdx2, sum(tmpx-1>=1 ) ),
                       rep( 2*lamdx2, sum(tmpx+1<=Nx) ),
                       rep( 2*lamdz2, sum(tmpz-1>=1 ) ),
                       rep( 2*lamdz2, sum(tmpz+1<=Nz) )
                   )
               )

UpdtMat <- sparseMatrix( i=ik, j=jk, x=xk, dims=c(Ny*Nx*Nz, Ny*Nx*Nz) )

cat('updmat created\n')

x0 <- as.vector(x0init)

for ( i in seq(1,MAXITER) )
{
        
    x1 <- ( as.vector(x0init) + UpdtMat %*% x0 ) / (1 + 4*(lamdx2+lamdy2+lamdz2) )
    
    x0 <- x1
    updtcnt<-updtcnt+1
    #if( runif(1)<1/2000 )
    #{
    #    image(x=vec$xvec, y=vec$yvec, z=t(matrix(x0, nrow=Ny)))
    #    title(main=paste('iter ',i,'updated ',updtcnt,' new fnval = ',optfun2D(x0),'\n'))
    #}
    
    Arr <- array( x1, c(Ny, Nx, Nz) )
    cat('image displayed ',i,' \n')
    image(x=vec$xvec, y=vec$yvec, z=t(matrix(Arr[,,50], nrow=Ny)))
    title(main=paste('iter ',i,'updated ',updtcnt ))

}

