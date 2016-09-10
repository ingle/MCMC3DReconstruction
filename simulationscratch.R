source('getEllipsoidSheafData.R')
source('full3drecons.R')

library(Matrix)
library(VGAM)

set.seed(10)

NSIMS = 100

P <- c(4,6,12,16)
SNR <- c(5, 10, 20)

RhoMax=2; NRho=50; 
ZMax=4.5; NZ=100; a=1; b=1; c=1.5; 
zshift=2.25

Nx <- 100
Ny <- 100
Nz <- 100

vec <- list(xvec=seq(-2,2, length.out=Nx),
            yvec=seq(-2,2, length.out=Ny),
            zvec=seq(0,4.5,length.out=Nz))
grid <- expand.grid( vec$yvec, vec$xvec, vec$zvec )
colnames(grid) <- c("ygrid", "xgrid", "zgrid")
ftrue <- rep(1, length(grid$xgrid))
ftrue[grid$xgrid^2/a^2 + grid$ygrid^2/b^2 + (grid$zgrid-zshift)^2/c^2 <=1] <- 4
ftrue3d <- array( ftrue, c(Ny,Nx,Nz) )

mse = array(0, dim=c(2,NSIMS,length(P),length(SNR)))

for(snr in SNR)
{
    for( p in P )
    {
        cat('--------------- processing P=',p,'--------------------\n', sep='')

        dat <- getEllipsoidSheafData(P=p)

        for( nsim in 1:NSIMS )
        {
            cat('nsim = ', nsim,'\n', sep='')

            datn <- dat;
            datn$fdat <- dat$fdat + rnorm( dat$fdat, mean=0,
                                          sd=3*exp(-snr/20) )

            testrecon <- full3drecons(dat=datn, vec=vec, wt=0.01,
                                      maxiter=1000)
            testrecon3d <- array(testrecon, c(Ny,Nx,Nz))
            
            testreconnnb <- full3drecons( dat=datn, vec=vec, wt=0, maxiter=1 )
            testreconnnb3d <- array(testreconnnb, c(Ny,Nx,Nz))

            mse[1,nsim,p==P] = 1/(Nx*Ny*Nz)*sum( (as.vector(ftrue3d)-as.vector(testrecon3d))^2 )
            mse[2,nsim,p==P] = 1/(Nx*Ny*Nz)*sum( (as.vector(ftrue3d)-as.vector(testreconnnb3d))^2 )
        }
    }
}

