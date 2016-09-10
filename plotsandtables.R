# makes various plots for the full3d paper
# May 11, 2015

#folderpath = '/home/ingle/Documents/Papers Ultrasound/Full3D Reconstruction/'
folderpath = '~/Insync/Papers Ultrasound/Full3D Reconstruction/'

# 1. sampling pattern plot
library(plot3D)
source('getEllipsoidSheafData.R')

samplingpatt <- getEllipsoidSheafData(P=4, 2, 20, 4.5, 15)
scatter3D(x=samplingpatt$xdat, y=samplingpatt$ydat, z=samplingpatt$zdat, 
          colvar=samplingpatt$fdat, theta=25, phi=50)

dev.copy2eps(file=paste(folderpath,'samplingpattern.eps',sep=''))
dev.off()

# 2. image statistics (snr/c/etc.)
source('calcstats.R')
source('calcstatsFea.R')

# 3. volume reconstruction
library(plot3D)
load('processed_1_4.RData')
zlabcm = seq(1, 4.5, length.out=length(vec$zvec))
slice3D( vec$yvec, vec$xvec, zlabcm, 
         colvar=frecons3d[,,seq(dim(frecons3d)[3],1,by=-1)],
         theta=35, phi=20, facets=T, xs=0, ys=0, zs=zlabcm[40],
         ticktype='detailed', clab='SWV (m/s)', xlab='y (cm)',
         ylab='x (cm)', zlab='z (cm)', cex.axis=1.3, 
         cex.lab=1.5, cex.sub=1.5 )
dev.copy2eps(file=paste(folderpath,'mrf4.eps',sep=''))
dev.off()

load('processed_nnb_1_4.RData')
zlabcm = seq(1, 4.5, length.out=length(vec$zvec))
slice3D( vec$yvec, vec$xvec, zlabcm, 
         colvar=frecons3d[,,seq(dim(frecons3d)[3],1,by=-1)],
         theta=35, phi=20, facets=T, xs=0, ys=0, zs=zlabcm[40],
         ticktype='detailed', clab='SWV (m/s)', xlab='y (cm)',
         ylab='x (cm)', zlab='z (cm)', cex.axis=1.3, 
         cex.lab=1.5, cex.sub=1.5 )
dev.copy2eps(file=paste(folderpath,'nnb4.eps',sep=''))
dev.off()

tab = read.table('processed_linearnb_1_4.csv', header = T)
frecons3d = array(tab$f, c(100,100,90))
frecons3d[frecons3d<0] = 0
frecons3d[frecons3d>4.5] = 4.5
slice3D( vec$yvec, vec$xvec, zlabcm, 
         colvar=frecons3d[,,seq(dim(frecons3d)[3],1,by=-1)],
         theta=35+270, phi=20, facets=T, xs=0.1, ys=0.1, zs=zlabcm[40],
         ticktype='detailed', clab='SWV (m/s)', xlab='y (cm)',
         ylab='x (cm)', zlab='z (cm)', cex.axis=1.3, 
         cex.lab=1.5, cex.sub=1.5 )
dev.copy2eps(file=paste(folderpath,'lnb4.eps',sep=''))
dev.off()


load('processed_1_16.RData')
zlabcm = seq(1, 4.5, length.out=length(vec$zvec))
slice3D( vec$yvec, vec$xvec, zlabcm, 
         colvar=frecons3d[,,seq(dim(frecons3d)[3],1,by=-1)],
         theta=35, phi=20, facets=T, xs=0, ys=0, zs=zlabcm[40],
         ticktype='detailed', clab='SWV (m/s)', xlab='y (cm)',
         ylab='x (cm)', zlab='z (cm)', cex.axis=1.3, 
         cex.lab=1.5, cex.sub=1.5 )
dev.copy2eps(file=paste(folderpath,'mrf16.eps',sep=''))
dev.off()

load('processed_nnb_1_16.RData')
zlabcm = seq(1, 4.5, length.out=length(vec$zvec))
slice3D( vec$yvec, vec$xvec, zlabcm, 
         colvar=frecons3d[,,seq(dim(frecons3d)[3],1,by=-1)],
         theta=35, phi=20, facets=T, xs=0, ys=0, zs=zlabcm[40],
         ticktype='detailed', clab='SWV (m/s)', xlab='y (cm)',
         ylab='x (cm)', zlab='z (cm)', cex.axis=1.3, 
         cex.lab=1.5, cex.sub=1.5 )
dev.copy2eps(file=paste(folderpath,'nnb16.eps',sep=''))
dev.off()

tab = read.table('processed_linearnb_1_16.csv', header = T)
frecons3d = array(tab$f, c(100,100,90))
frecons3d[frecons3d<0] = 0
frecons3d[frecons3d>4.5] = 4.5
slice3D( vec$yvec, vec$xvec, zlabcm, 
         colvar=frecons3d[,,seq(dim(frecons3d)[3],1,by=-1)],
         theta=35+270, phi=20, facets=T, xs=0.1, ys=0.1, zs=zlabcm[40],
         ticktype='detailed', clab='SWV (m/s)', xlab='y (cm)',
         ylab='x (cm)', zlab='z (cm)', cex.axis=1.3, 
         cex.lab=1.5, cex.sub=1.5 )
dev.copy2eps(file=paste(folderpath,'lnb16.eps',sep=''))
dev.off()

#4. Fea data volume reconstruction
load('processed_fea_1_4.RData')
zlabcm = seq(1, 4.5, length.out=length(vec$zvec))
slice3D( vec$yvec, vec$xvec, zlabcm, 
         colvar=frecons3d[,,seq(dim(frecons3d)[3],1,by=-1)],
         theta=35, phi=20, facets=T, xs=0, ys=0, zs=zlabcm[35],
         ticktype='detailed', clab='SWV (m/s)', xlab='y (cm)',
         ylab='x (cm)', zlab='z (cm)', cex.axis=1.3, 
         cex.lab=1.5, cex.sub=1.5 )
dev.copy2eps(file=paste(folderpath,'mrf_fea4.eps',sep=''))
dev.off()

load('processed_fea_nnb_1_4.RData')
zlabcm = seq(1, 4.5, length.out=length(vec$zvec))
slice3D( vec$yvec, vec$xvec, zlabcm, 
         colvar=frecons3d[,,seq(dim(frecons3d)[3],1,by=-1)],
         theta=35, phi=20, facets=T, xs=0, ys=0, zs=zlabcm[35],
         ticktype='detailed', clab='SWV (m/s)', xlab='y (cm)',
         ylab='x (cm)', zlab='z (cm)', cex.axis=1.3, 
         cex.lab=1.5, cex.sub=1.5 )
dev.copy2eps(file=paste(folderpath,'nnb_fea4.eps',sep=''))
dev.off()

tab = read.table('processed_fea_linearnb_1_4.csv', header = T)
frecons3d = array(tab$f, c(100,100,90))
frecons3d[frecons3d<0] = 0
frecons3d[frecons3d>4.5] = 4.5
slice3D( vec$yvec, vec$xvec, zlabcm, 
         colvar=frecons3d[,,seq(dim(frecons3d)[3],1,by=-1)],
         theta=35+270, phi=20, facets=T, xs=0, ys=0, zs=zlabcm[35],
         ticktype='detailed', clab='SWV (m/s)', xlab='y (cm)',
         ylab='x (cm)', zlab='z (cm)', cex.axis=1.3, 
         cex.lab=1.5, cex.sub=1.5 )
dev.copy2eps(file=paste(folderpath,'lnb_fea4.eps',sep=''))
dev.off()

load('processed_fea_1_16.RData')
zlabcm = seq(1, 4.5, length.out=length(vec$zvec))
slice3D( vec$yvec, vec$xvec, zlabcm, 
         colvar=frecons3d[,,seq(dim(frecons3d)[3],1,by=-1)],
         theta=35, phi=20, facets=T, xs=0, ys=0, zs=zlabcm[35],
         ticktype='detailed', clab='SWV (m/s)', xlab='y (cm)',
         ylab='x (cm)', zlab='z (cm)', cex.axis=1.3, 
         cex.lab=1.5, cex.sub=1.5 )
dev.copy2eps(file=paste(folderpath,'mrf_fea16.eps',sep=''))
dev.off()

load('processed_fea_nnb_1_16.RData')
zlabcm = seq(1, 4.5, length.out=length(vec$zvec))
slice3D( vec$yvec, vec$xvec, zlabcm, 
         colvar=frecons3d[,,seq(dim(frecons3d)[3],1,by=-1)],
         theta=35, phi=20, facets=T, xs=0, ys=0, zs=zlabcm[35],
         ticktype='detailed', clab='SWV (m/s)', xlab='y (cm)',
         ylab='x (cm)', zlab='z (cm)', cex.axis=1.3, 
         cex.lab=1.5, cex.sub=1.5 )
dev.copy2eps(file=paste(folderpath,'nnb_fea16.eps',sep=''))
dev.off()

tab = read.table('processed_fea_linearnb_1_16.csv', header = T)
frecons3d = array(tab$f, c(100,100,90))
frecons3d[frecons3d<0] = 0
frecons3d[frecons3d>4.5] = 4.5
slice3D( vec$yvec, vec$xvec, zlabcm, 
         colvar=frecons3d[,,seq(dim(frecons3d)[3],1,by=-1)],
         theta=35, phi=20, facets=T, xs=0.1, ys=0.1, zs=zlabcm[35],
         ticktype='detailed', clab='SWV (m/s)', xlab='y (cm)',
         ylab='x (cm)', zlab='z (cm)', cex.axis=1.3, 
         cex.lab=1.5, cex.sub=1.5 )
dev.copy2eps(file=paste(folderpath,'lnb_fea16.eps',sep=''))
dev.off()
