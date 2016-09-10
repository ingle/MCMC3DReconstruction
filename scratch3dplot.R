library(RANN)
library(Matrix)
library(VGAM)
source('full3drecons.R')

# dat <- read.table('xyzfdata_Oct26_2013_1_6.dat', header=T)



# datn <- dat;
# datn$fdat <- dat$fdat + rnorm( dat$fdat, mean=0, sd=3*exp(-15/20) )

# testrecon <- full3drecons(dat=datn, vec=vec, wt=0.01, maxiter=1000)
# testrecon3d <- array(testrecon, c(100,100,100))
# library(plot3D)
# slice3D( vec$yvec, vec$xvec, vec$zvec, 
#         colvar=testrecon3d[,,seq(dim(testrecon3d)[3],1,by=-1)],
#           theta=35, phi=20, facets=T, xs=0, ys=0, zs=vec$zvec[40] )

#scatter3D( datn$xdat, datn$ydat, datn$zdat, colvar=datn$fdat )

bkgx = c(73, 83)
bkgy = c(67, 88)
bkgz = c(30, 50)

incx = c(31, 41)
incy = c(40, 55)
incz = c(30, 50)

allstats <- expand.grid( c(4,6,12,16), 1:5 )
colnames(allstats) <- c("nsli", "nset");
allstats$bkgmean  <- rep(0,times=20)
allstats$bkgstdv  <- rep(0,times=20)
allstats$incmean  <- rep(0,times=20)
allstats$incstdv  <- rep(0,times=20)
allstats$bkgsnr   <- rep(0,times=20)
allstats$incsnr   <- rep(0,times=20)
allstats$c        <- rep(0,times=20)
allstats$cnr      <- rep(0,times=20)

for( idx in 1:nrow(allstats) )
{
    load(paste('processed_',allstats$nset[idx],'_',allstats$nsli[idx],'.RData',sep=''))
    bkgroi <- frecons3d[bkgx[1]:bkgx[2], bkgy[1]:bkgy[2], bkgz[1]:bkgz[2]]
    incroi <- frecons3d[incx[1]:incx[2], incy[1]:incy[2], incz[1]:incz[2]]

    allstats$bkgmean[idx] = mean(bkgroi)
    allstats$bkgstdv[idx] = sd  (bkgroi)
    allstats$incmean[idx] = mean(incroi)
    allstats$incstdv[idx] = sd  (incroi)

    allstats$bkgsnr[idx]  = mean(bkgroi) / sd(bkgroi)
    allstats$incsnr[idx]  = mean(incroi) / sd(incroi)

    allstats$c[idx] = mean(incroi) / mean(bkgroi)
    allstats$cnr[idx] = ( mean(incroi)-mean(bkgroi) ) / (sd(bkgroi) + sd(incroi))
}

allstatsnnb <- expand.grid( c(4,6,12,16), 1:5 )
colnames(allstatsnnb) <- c("nsli", "nset");
allstatsnnb$bkgmean  <- rep(0,times=20)
allstatsnnb$bkgstdv  <- rep(0,times=20)
allstatsnnb$incmean  <- rep(0,times=20)
allstatsnnb$incstdv  <- rep(0,times=20)
allstatsnnb$bkgsnr   <- rep(0,times=20)
allstatsnnb$incsnr   <- rep(0,times=20)
allstatsnnb$c        <- rep(0,times=20)
allstatsnnb$cnr      <- rep(0,times=20)

for( idx in 1:nrow(allstatsnnb) )
{
    load(paste('processed_nnb_',allstatsnnb$nset[idx],'_',allstatsnnb$nsli[idx],'.RData',sep=''))
    bkgroi <- frecons3d[bkgx[1]:bkgx[2], bkgy[1]:bkgy[2], bkgz[1]:bkgz[2]]
    incroi <- frecons3d[incx[1]:incx[2], incy[1]:incy[2], incz[1]:incz[2]]

    allstatsnnb$bkgmean[idx] = mean(bkgroi)
    allstatsnnb$bkgstdv[idx] = sd  (bkgroi)
    allstatsnnb$incmean[idx] = mean(incroi)
    allstatsnnb$incstdv[idx] = sd  (incroi)

    allstatsnnb$bkgsnr[idx]  = mean(bkgroi) / sd(bkgroi)
    allstatsnnb$incsnr[idx]  = mean(incroi) / sd(incroi)

    allstatsnnb$c[idx] = mean(incroi) / mean(bkgroi)
    allstatsnnb$cnr[idx] = ( mean(incroi)-mean(bkgroi) ) / (sd(bkgroi) + sd(incroi))
}

load( 'processed_1_16.RData' )

for( nPln in 1:90 )
{
    image(1:100, 1:100, frecons3d[,,nPln] )
    title( paste("nPln=",nPln,sep='') )
    points( c(73,83,73,83), c(67, 67, 88, 88) )
    points( c(31,41,31,41), c(40,40,55,55) )
    tmp=readline()
}

