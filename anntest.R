set.seed(10)
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

# data = c(0,0,0,1,1,1,2,2,2,3,3,3)
# data = t(Matrix( data, nrow=3))
# query = c(0.1,0.2,0.1,2.55,2.4,10)
# query = t(Matrix( query, nrow=3))

grid <- expand.grid( vec$yvec, vec$xvec, vec$zvec )
cloud <- cbind( ellipsoiddat$ydat, ellipsoiddat$xdat, ellipsoiddat$zdat )

nnres <- nn2( data=cloud, query=grid, k=1 )

x0 = ellipsoiddat$fdatn[nnres$nn.idx]
Arr = array( x0, c(Ny, Nx, Nz) )

for (i in 1:Nz)
{
    z = t(matrix(Arr[,,i], nrow=Ny))
    image(x=vec$xvec, y=vec$yvec, z=z, zlim=c(-1.1,2.1), col=heat.colors(100))
    title(main=paste('z = ', vec$zvec[i]))
    Sys.sleep(0.2)
}




