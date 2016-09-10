library(RANN)
library(Matrix)
library(VGAM)

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

#bkgx = c(73, 83)
#bkgy = c(67, 88)
#bkgz = c(30, 50)

#incx = c(31, 41)
#incy = c(40, 55)
#incz = c(30, 50)

# we use slightly different ROI box limits for fea dataset
bkgx = c(73, 83)
bkgy = c(67, 88)
bkgz = c(45, 65)

incx = c(31, 41)
incy = c(40, 55) 
incz = c(45, 65)

# load ground truth FEA SWV map
ground_truth = read.table('groundtruth_fea.csv',sep=' ', header=T)
ground_truth_3d = array(ground_truth$f, c(100,100,90))
mask_3d = array(ground_truth$mask, c(100,100,90))

MSE_mat = matrix(NA, 3,4)
rownames(MSE_mat) <- c('mrf', 'nnb', 'lnb')
colnames(MSE_mat) <- c(4,6,12,16)

# ------------- analyze mrf reconstructions -----------------------------
allstats <- expand.grid( c(4,6,12,16), 1:1 )
colnames(allstats) <- c("nsli", "nset");
allstats$bkgmean   <- rep(0,times=4)
allstats$bkgstdv   <- rep(0,times=4)
allstats$incmean   <- rep(0,times=4)
allstats$incstdv   <- rep(0,times=4)
allstats$bkgsnr    <- rep(0,times=4)
allstats$incsnr    <- rep(0,times=4)
allstats$c         <- rep(0,times=4)
allstats$cnr       <- rep(0,times=4)

for( idx in 1:nrow(allstats) )
{
  load(paste('processed_fea_',allstats$nset[idx],'_',allstats$nsli[idx],'.RData',sep=''))
  
  MSE_mat['mrf',toString(allstats$nsli[idx])] = (sum(((frecons3d - ground_truth_3d)*mask_3d)**2.0)/sum(mask_3d))
  
  bkgroi <- frecons3d[bkgx[1]:bkgx[2], bkgy[1]:bkgy[2], bkgz[1]:bkgz[2]]
  incroi <- frecons3d[incx[1]:incx[2], incy[1]:incy[2], incz[1]:incz[2]]
  
  allstats$bkgmean[idx] = mean(bkgroi)
  allstats$bkgstdv[idx] = sd  (bkgroi)
  allstats$incmean[idx] = mean(incroi)
  allstats$incstdv[idx] = sd  (incroi)
  
  allstats$bkgsnr[idx]  = mean(bkgroi) / sd(bkgroi)
  allstats$incsnr[idx]  = mean(incroi) / sd(incroi)
  
  allstats$c[idx] = mean(incroi) / mean(bkgroi)
  allstats$cnr[idx] = ( mean(incroi)-mean(bkgroi) ) / sqrt(sd(bkgroi)^2 + sd(incroi)^2)
}

# get overall means and stdevs over the ROIs
fourbkgmeans <- numeric(0)
fourincmeans <- numeric(0)
fourbkgstdvs <- numeric(0)
fourincstdvs <- numeric(0)
for( nsli in c(4,6,12,16) )
{
  bkgroi <- numeric(0)
  incroi <- numeric(0)
  for( nset in 1:1 )
  {
    load(paste('processed_fea_',nset,'_',nsli,'.RData',sep=''))
    bkgroi <- c(bkgroi, as.vector(frecons3d[bkgx[1]:bkgx[2], bkgy[1]:bkgy[2], bkgz[1]:bkgz[2]]))
    incroi <- c(incroi, as.vector(frecons3d[incx[1]:incx[2], incy[1]:incy[2], incz[1]:incz[2]]))
  }
  fourbkgmeans[nsli==c(4,6,12,16)] <- mean(bkgroi)
  fourbkgstdvs[nsli==c(4,6,12,16)] <- sd(  bkgroi)
  fourincmeans[nsli==c(4,6,12,16)] <- mean(incroi)
  fourincstdvs[nsli==c(4,6,12,16)] <- sd(  incroi)
}

# ------------- analyze nnb reconstructions -----------------------------
allstatsnnb <- expand.grid( c(4,6,12,16), 1:1 )
colnames(allstatsnnb) <- c("nsli", "nset")
allstatsnnb$bkgmean   <- rep(0,times=4)
allstatsnnb$bkgstdv   <- rep(0,times=4)
allstatsnnb$incmean   <- rep(0,times=4)
allstatsnnb$incstdv   <- rep(0,times=4)
allstatsnnb$bkgsnr    <- rep(0,times=4)
allstatsnnb$incsnr    <- rep(0,times=4)
allstatsnnb$c         <- rep(0,times=4)
allstatsnnb$cnr       <- rep(0,times=4)

for( idx in 1:nrow(allstatsnnb) )
{
  load(paste('processed_fea_nnb_',allstatsnnb$nset[idx],'_',allstatsnnb$nsli[idx],'.RData',sep=''))
  
  MSE_mat['nnb',toString(allstats$nsli[idx])] = (sum(((frecons3d - ground_truth_3d)*mask_3d)**2.0)/sum(mask_3d))
  
  bkgroi <- frecons3d[bkgx[1]:bkgx[2], bkgy[1]:bkgy[2], bkgz[1]:bkgz[2]]
  incroi <- frecons3d[incx[1]:incx[2], incy[1]:incy[2], incz[1]:incz[2]]
  
  allstatsnnb$bkgmean[idx] = mean(bkgroi)
  allstatsnnb$bkgstdv[idx] = sd  (bkgroi)
  allstatsnnb$incmean[idx] = mean(incroi)
  allstatsnnb$incstdv[idx] = sd  (incroi)
  
  allstatsnnb$bkgsnr[idx]  = mean(bkgroi) / sd(bkgroi)
  allstatsnnb$incsnr[idx]  = mean(incroi) / sd(incroi)
  
  allstatsnnb$c[idx] = mean(incroi) / mean(bkgroi)
  allstatsnnb$cnr[idx] = ( mean(incroi)-mean(bkgroi) ) / sqrt(sd(bkgroi)^2 + sd(incroi)^2)
}

# get overall means and stdevs over the ROIs
fourbkgmeansnnb <- numeric(0)
fourincmeansnnb <- numeric(0)
fourbkgstdvsnnb <- numeric(0)
fourincstdvsnnb <- numeric(0)
for( nsli in c(4,6,12,16) )
{
  bkgroi <- numeric(0)
  incroi <- numeric(0)
  for( nset in 1:1 )
  {
    load(paste('processed_fea_nnb_',nset,'_',nsli,'.RData',sep=''))
    bkgroi <- c(bkgroi, as.vector(frecons3d[bkgx[1]:bkgx[2], bkgy[1]:bkgy[2], bkgz[1]:bkgz[2]]))
    incroi <- c(incroi, as.vector(frecons3d[incx[1]:incx[2], incy[1]:incy[2], incz[1]:incz[2]]))
  }
  fourbkgmeansnnb[nsli==c(4,6,12,16)] <- mean(bkgroi)
  fourbkgstdvsnnb[nsli==c(4,6,12,16)] <- sd(  bkgroi)
  fourincmeansnnb[nsli==c(4,6,12,16)] <- mean(incroi)
  fourincstdvsnnb[nsli==c(4,6,12,16)] <- sd(  incroi)
}

# ---------------------- analyze linear interpn reconstructions -----------------------------
allstatslnb <- expand.grid( c(4,6,12,16), 1:1 )
colnames(allstatslnb) <- c("nsli", "nset")
allstatslnb$bkgmean   <- rep(0,times=4)
allstatslnb$bkgstdv   <- rep(0,times=4)
allstatslnb$incmean   <- rep(0,times=4)
allstatslnb$incstdv   <- rep(0,times=4)
allstatslnb$bkgsnr    <- rep(0,times=4)
allstatslnb$incsnr    <- rep(0,times=4)
allstatslnb$c         <- rep(0,times=4)
allstatslnb$cnr       <- rep(0,times=4)

for( idx in 1:nrow(allstatslnb) )
{
  #load(paste('processed_nnb_',allstatsnnb$nset[idx],'_',allstatsnnb$nsli[idx],'.RData',sep=''))
  lnbcsv = read.table(paste('processed_fea_linearnb_',allstatslnb$nset[idx],'_',allstatslnb$nsli[idx],'.csv',sep=''),sep=' ', header=T)
  frecons3d = array(lnbcsv$f, c(100,100,90))
  frecons3d[frecons3d<0] = 0
  frecons3d[frecons3d>4] = 4
  
  MSE_mat['lnb',toString(allstats$nsli[idx])] = (sum(((frecons3d - ground_truth_3d)*mask_3d)**2.0)/sum(mask_3d))
  
  bkgroi <- frecons3d[bkgx[1]:bkgx[2], bkgy[1]:bkgy[2], bkgz[1]:bkgz[2]]
  incroi <- frecons3d[incx[1]:incx[2], incy[1]:incy[2], incz[1]:incz[2]]
  
  allstatslnb$bkgmean[idx] = mean(bkgroi)
  allstatslnb$bkgstdv[idx] = sd  (bkgroi)
  allstatslnb$incmean[idx] = mean(incroi)
  allstatslnb$incstdv[idx] = sd  (incroi)
  
  allstatslnb$bkgsnr[idx]  = mean(bkgroi) / sd(bkgroi)
  allstatslnb$incsnr[idx]  = mean(incroi) / sd(incroi)
  
  allstatslnb$c[idx] = mean(incroi) / mean(bkgroi)
  allstatslnb$cnr[idx] = ( mean(incroi)-mean(bkgroi) ) / sqrt(sd(bkgroi)^2 + sd(incroi)^2)
}

# get overall means and stdevs over the ROIs
fourbkgmeanslnb <- numeric(0)
fourincmeanslnb <- numeric(0)
fourbkgstdvslnb <- numeric(0)
fourincstdvslnb <- numeric(0)
for( nsli in c(4,6,12,16) )
{
  bkgroi <- numeric(0)
  incroi <- numeric(0)
  for( nset in 1:1 )
  {
    #load(paste('processed_nnb_',nset,'_',nsli,'.RData',sep=''))
    lnbcsv = read.table(paste('processed_fea_linearnb_',nset,'_',nsli,'.csv',sep=''), sep=' ', header=T)
    frecons3d = array(lnbcsv$f, c(100,100,90))    
    frecons3d[frecons3d<0] = 0
    frecons3d[frecons3d>4] = 4
    bkgroi <- c(bkgroi, as.vector(frecons3d[bkgx[1]:bkgx[2], bkgy[1]:bkgy[2], bkgz[1]:bkgz[2]]))
    incroi <- c(incroi, as.vector(frecons3d[incx[1]:incx[2], incy[1]:incy[2], incz[1]:incz[2]]))
  }
  fourbkgmeanslnb[nsli==c(4,6,12,16)] <- mean(bkgroi)
  fourbkgstdvslnb[nsli==c(4,6,12,16)] <- sd(  bkgroi)
  fourincmeanslnb[nsli==c(4,6,12,16)] <- mean(incroi)
  fourincstdvslnb[nsli==c(4,6,12,16)] <- sd(  incroi)
}

#--------------- print results in LaTeX code format ---------------------------
cat('printing mean and std SWV (m/s) values MRF alg:\n')
cat(sprintf('%.2f', c(rbind( fourbkgmeans, fourbkgstdvs ) )), '\n')
cat(sprintf('%.2f', c(rbind( fourincmeans, fourincstdvs ) )), '\n\n\n')

cat('printing mean and std SWV (m/s) values NNB alg:\n')
cat(sprintf('%.2f', c(rbind( fourbkgmeansnnb, fourbkgstdvsnnb ) )), '\n')
cat(sprintf('%.2f', c(rbind( fourincmeansnnb, fourincstdvsnnb ) )), '\n\n\n')

cat('printing mean and std SWV (m/s) values LNB alg:\n')
cat(sprintf('%.2f', c(rbind( fourbkgmeanslnb, fourbkgstdvslnb ) )), '\n')
cat(sprintf('%.2f', c(rbind( fourincmeanslnb, fourincstdvslnb ) )), '\n\n\n')

cat('printing snr values (mean calculated first, then converted to dB):\n')
cat( sprintf('%.2f',
             c(20*log10( mean(allstats$bkgsnr[allstats$nsli==4 ]) ), 
               20*log10( mean(allstats$bkgsnr[allstats$nsli==6 ]) ),
               20*log10( mean(allstats$bkgsnr[allstats$nsli==12]) ),
               20*log10( mean(allstats$bkgsnr[allstats$nsli==16]) ))
), '\n'
)
cat( sprintf('%.2f',
             c(20*log10( mean(allstatsnnb$bkgsnr[allstatsnnb$nsli==4 ]) ), 
               20*log10( mean(allstatsnnb$bkgsnr[allstatsnnb$nsli==6 ]) ),
               20*log10( mean(allstatsnnb$bkgsnr[allstatsnnb$nsli==12]) ),
               20*log10( mean(allstatsnnb$bkgsnr[allstatsnnb$nsli==16]) ))
), '\n' 
)
cat( sprintf('%.2f',
             c(20*log10( mean(allstatslnb$bkgsnr[allstatslnb$nsli==4 ]) ), 
               20*log10( mean(allstatslnb$bkgsnr[allstatslnb$nsli==6 ]) ),
               20*log10( mean(allstatslnb$bkgsnr[allstatslnb$nsli==12]) ),
               20*log10( mean(allstatslnb$bkgsnr[allstatslnb$nsli==16]) ))
), '\n\n\n' 
)
cat( sprintf('%.2f',
             c(20*log10( mean(allstats$incsnr[allstats$nsli==4 ]) ), 
               20*log10( mean(allstats$incsnr[allstats$nsli==6 ]) ),
               20*log10( mean(allstats$incsnr[allstats$nsli==12]) ),
               20*log10( mean(allstats$incsnr[allstats$nsli==16]) ))
), '\n'
)
cat( sprintf('%.2f',
             c(20*log10( mean(allstatsnnb$incsnr[allstatsnnb$nsli==4 ]) ), 
               20*log10( mean(allstatsnnb$incsnr[allstatsnnb$nsli==6 ]) ),
               20*log10( mean(allstatsnnb$incsnr[allstatsnnb$nsli==12]) ),
               20*log10( mean(allstatsnnb$incsnr[allstatsnnb$nsli==16]) ))
), '\n' 
)
cat( sprintf('%.2f',
             c(20*log10( mean(allstatslnb$incsnr[allstatslnb$nsli==4 ]) ), 
               20*log10( mean(allstatslnb$incsnr[allstatslnb$nsli==6 ]) ),
               20*log10( mean(allstatslnb$incsnr[allstatslnb$nsli==12]) ),
               20*log10( mean(allstatslnb$incsnr[allstatslnb$nsli==16]) ))
), '\n\n\n' 
)

cat('printing contrast values (mean calculated first, then converted to dB):\n')
cat( sprintf('%.2f',
             c(20*log10( mean(allstats$c[allstats$nsli==4 ]) ), 
               20*log10( mean(allstats$c[allstats$nsli==6 ]) ),
               20*log10( mean(allstats$c[allstats$nsli==12]) ),
               20*log10( mean(allstats$c[allstats$nsli==16]) ))
), '\n'
)
cat( sprintf('%.2f',
             c(20*log10( mean(allstatsnnb$c[allstatsnnb$nsli==4 ]) ), 
               20*log10( mean(allstatsnnb$c[allstatsnnb$nsli==6 ]) ),
               20*log10( mean(allstatsnnb$c[allstatsnnb$nsli==12]) ),
               20*log10( mean(allstatsnnb$c[allstatsnnb$nsli==16]) ))
), '\n' 
)
cat( sprintf('%.2f',
             c(20*log10( mean(allstatslnb$c[allstatslnb$nsli==4 ]) ), 
               20*log10( mean(allstatslnb$c[allstatslnb$nsli==6 ]) ),
               20*log10( mean(allstatslnb$c[allstatslnb$nsli==12]) ),
               20*log10( mean(allstatslnb$c[allstatslnb$nsli==16]) ))
), '\n\n\n' 
)

cat('printing cnr values (mean calculated first, then converted to dB):\n')
cat( sprintf('%.2f',
             c(20*log10( mean(allstats$cnr[allstats$nsli==4 ]) ), 
               20*log10( mean(allstats$cnr[allstats$nsli==6 ]) ),
               20*log10( mean(allstats$cnr[allstats$nsli==12]) ),
               20*log10( mean(allstats$cnr[allstats$nsli==16]) ))
), '\n'
)
cat( sprintf('%.2f',
             c(20*log10( mean(allstatsnnb$cnr[allstatsnnb$nsli==4 ]) ), 
               20*log10( mean(allstatsnnb$cnr[allstatsnnb$nsli==6 ]) ),
               20*log10( mean(allstatsnnb$cnr[allstatsnnb$nsli==12]) ),
               20*log10( mean(allstatsnnb$cnr[allstatsnnb$nsli==16]) ))
), '\n' 
)
cat( sprintf('%.2f',
             c(20*log10( mean(allstatslnb$cnr[allstatslnb$nsli==4 ]) ), 
               20*log10( mean(allstatslnb$cnr[allstatslnb$nsli==6 ]) ),
               20*log10( mean(allstatslnb$cnr[allstatslnb$nsli==12]) ),
               20*log10( mean(allstatslnb$cnr[allstatslnb$nsli==16]) ))
), '\n\n\n' 
)

print(MSE_mat)