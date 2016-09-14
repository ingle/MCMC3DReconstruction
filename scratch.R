library(plot3D)
source('full3drecons.R')
source('getSmoothEllipsoidSheafData.R')

Nx = 100
Ny = 100
Nz = 100

vec <- list(xvec=seq(-2,2, length.out=Nx),
                yvec=seq(-2,2, length.out=Ny),
                zvec=seq(0,4.5,length.out=Nz))

dat <- getSmoothEllipsoidSheafDataCylBore(P=16)
datn <- dat
datn$fdat <- dat$fdat + rnorm( dat$fdat, mean=0, sd=3*exp(-10/20) )
testrecon <- full3drecons(dat=datn, vec=vec, wt=0.01, maxiter=100)
frecons3d <- array( testrecon, c(Ny,Nx,Nz) )

zlabcm = seq(1, 4.5, length.out=length(vec$zvec))
slice3D( vec$yvec, vec$xvec, zlabcm, 
         colvar=frecons3d[,,seq(dim(frecons3d)[3],1,by=-1)],
         theta=35, phi=20, facets=T, xs=0, ys=0, zs=zlabcm[40],
         ticktype='detailed', clab='SWV (m/s)', xlab='y (cm)',
         ylab='x (cm)', zlab='z (cm)', cex.axis=1.3, 
         cex.lab=1.5, cex.sub=1.5 )

