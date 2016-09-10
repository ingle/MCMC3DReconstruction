saveplots = T

if (!file.exists('simulatedmseresults.RData'))
{
    cat('will simulate for many hours. go home take a nap.\n')
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
            cat('--------------- processing P=',p,' SNR=',snr,'--------------------\n', sep='')

            dat <- getEllipsoidSheafData(P=p)

            for( nsim in 1:NSIMS )
            {
                cat('nsim = ', nsim,'\n', sep='')

                datn <- dat
                datn$fdat <- dat$fdat + rnorm( dat$fdat, mean=0,
                                              sd=3*exp(-snr/20) )

                testrecon <- full3drecons(dat=datn, vec=vec, wt=0.01,
                                          maxiter=1000)
                testrecon3d <- array(testrecon, c(Ny,Nx,Nz))
                
                testreconnnb <- full3drecons( dat=datn, vec=vec, wt=0, maxiter=1 )
                testreconnnb3d <- array(testreconnnb, c(Ny,Nx,Nz))

                mse[1,nsim,p==P,snr==SNR] = 1/(Nx*Ny*Nz)*sum( (as.vector(ftrue3d)-as.vector(testrecon3d))^2 )
                mse[2,nsim,p==P,snr==SNR] = 1/(Nx*Ny*Nz)*sum( (as.vector(ftrue3d)-as.vector(testreconnnb3d))^2 )
            }
        }
    }

    save.image("simulatedmseresults.RData")
    cat('results saved. done.\n')
} else
{
    cat('data file exists. loading from disk.\n')
    load('simulatedmseresults.RData')
    #    msesim <- as.data.frame( list(   mse=as.vector(mse), 
    #                     alg=rep(rbind(rep('full3d',NSIMS),rep('nnb',NSIMS)),length(P)),
    #                     P=rep(P,each=NSIMS*2)
    #                  ) )
    msesimfull3d <-  rbind(data.frame( snr='5dB', slices=P, 
                                      meanerr=apply(mse4lev[,,,1],c(1,3),mean)[1,], 
                                      sdev=apply(mse4lev[,,,1],c(1,3),sd)[1,] ),
                         data.frame( snr='10dB', slices=P, 
                                    meanerr=apply(mse4lev[,,,2],c(1,3),mean)[1,], 
                                    sdev=apply(mse4lev[,,,2],c(1,3),sd)[1,] ),
                         data.frame( snr='15dB', slices=P, 
                                    meanerr=apply(mse4lev[,,,3],c(1,3),mean)[1,], 
                                    sdev=apply(mse4lev[,,,3],c(1,3),sd)[1,] ),
                         data.frame( snr='20dB', slices=P, 
                                    meanerr=apply(mse4lev[,,,4],c(1,3),mean)[1,], 
                                    sdev=apply(mse4lev[,,,4],c(1,3),sd)[1,] )
                         )
    msesimnnb    <-  rbind(data.frame( snr='5dB', slices=P, 
                                      meanerr=apply(mse4lev[,,,1],c(1,3),mean)[2,], 
                                      sdev=apply(mse4lev[,,,1],c(1,3),sd)[2,] ),
                         data.frame( snr='10dB', slices=P, 
                                    meanerr=apply(mse4lev[,,,2],c(1,3),mean)[2,], 
                                    sdev=apply(mse4lev[,,,2],c(1,3),sd)[2,] ),
                         data.frame( snr='15dB', slices=P, 
                                    meanerr=apply(mse4lev[,,,3],c(1,3),mean)[2,], 
                                    sdev=apply(mse4lev[,,,3],c(1,3),sd)[2,] ),
                         data.frame( snr='20dB', slices=P, 
                                    meanerr=apply(mse4lev[,,,4],c(1,3),mean)[2,], 
                                    sdev=apply(mse4lev[,,,4],c(1,3),sd)[2,] )
                         )
    library('ggplot2')
    graphics.off()
    
    f1 = ggplot(data = msesimfull3d, aes(x = slices, y = meanerr, group = snr) ) + 
        geom_errorbar(aes(ymin = meanerr-sdev, ymax = meanerr + sdev), width=0.3) + 
        geom_line() + 
        geom_point(aes(shape=snr, fill=snr), size=5) + 
        scale_x_continuous("Number of slices") + 
        scale_y_continuous(bquote('MSE ' ~m^2 / s^2)) + 
        scale_shape_manual(values=c(24,23,22,21)) + 
        scale_fill_manual(values=c("white","black", "white", "black")) + 
        theme_bw() + 
        theme(
             axis.title.x = element_text(face="bold", size=20),
             axis.title.y = element_text(face="bold", size=20, angle=90),
             axis.text.x  = element_text(size=20),
             axis.text.y  = element_text(size=20),
             panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(), 
             legend.position = c(0.8,0.8), 
             legend.title = element_blank(), 
             legend.text = element_text(size=20),
             legend.key.size = unit(1.5, "lines"),
             legend.key = element_blank() 
             )
    print(f1)
    if(saveplots)
    {
      dev.copy2eps(file='~/Insync/Papers Ultrasound/Full3D Reconstruction/full3dmse.eps')
    }
    #dev.copy2eps(file='/home/ingle/Documents/Papers Ultrasound/Full3D Reconstruction/full3dmse.eps')

    dev.new()
    f2 = ggplot(data = msesimnnb, aes(x = slices, y = meanerr, group = snr) ) + 
        geom_errorbar(aes(ymin = meanerr-10*sdev, ymax = meanerr + 10*sdev), width=0.3) + 
        geom_line() + 
        geom_point(aes(shape=snr, fill=snr), size=5) + 
        scale_x_continuous("Number of slices") + 
        scale_y_continuous(bquote('MSE ' ~m^2 / s^2)) + 
        scale_shape_manual(values=c(24,23,22,21)) + 
        scale_fill_manual(values=c("white","black", "white", "black")) + 
        theme_bw() + 
        theme(
             axis.title.x = element_text(face="bold", size=20),
             axis.title.y = element_text(face="bold", size=20, angle=90),
             axis.text.x  = element_text(size=20),
             axis.text.y  = element_text(size=20),
             panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(), 
             legend.position = c(0.8,0.65), 
             legend.title = element_blank(), 
             legend.text = element_text(size=20),
             legend.key.size = unit(1.5, "lines"),
             legend.key = element_blank() 
             )
    print(f2)
    if(saveplots)
    {
      dev.copy2eps(file='~/Insync/Papers Ultrasound/Full3D Reconstruction/nnbmse.eps')
    }
    
    dev.new()
    f3 = ggplot(data = NULL, aes(x = slices, y = 20*log10(meanerr), group = snr) ) + 
      geom_line(data=msesimfull3d) + 
      geom_line(data=msesimnnb, linetype="dotted") + 
      geom_point(data=msesimfull3d, aes(shape=snr, fill=snr), size=3) + 
      geom_point(data=msesimnnb, aes(shape=snr, fill=snr), size=3) + 
      scale_x_continuous("Number of slices") + 
      scale_y_continuous(bquote('20 log(MSE)')) + 
      scale_shape_manual(values=c(24,23,22,21)) + 
      scale_fill_manual(values=c("white","black", "white", "black")) +
      theme_bw() + 
      theme( 
             axis.title.x = element_text(face="bold", size=20),
             axis.title.y = element_text(face="bold", size=20, angle=90),
             axis.text.x  = element_text(size=20),
             axis.text.y  = element_text(size=20),
             panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(), 
             legend.position = c(0.8,0.25), 
             legend.title = element_blank(), 
             legend.text = element_text(size=20),
             legend.key.size = unit(1.5, "lines"),
             legend.key = element_blank() 
           ) +
    annotate("text", x=8,y=9, label="NNB", size=5) +
    annotate("text", x=8,y=-12, label="MRF", size=5) +
    annotate("text", x=14,y=-4, label="SNR", size=8)
    
    print(f3)
    if(saveplots)
    {
      dev.copy2eps(file='~/Insync/Papers Ultrasound/Full3D Reconstruction/nnbmrf_mse_db.eps')
    }
    if(saveplots)
    {
      graphics.off()
    }

    # gp<-qplot( slices, meanerr, shape=snr, linetype=snr, data=msesimfull3d ) + geom_line() +
    #                    geom_errorbar(aes(ymin=meanerr-sdev, ymax=meanerr+sdev)) +
    #                    xlab('Number of slices (P)') +
    #                    ylab('MSE') +
    #                    ggtitle('Full 3D reconstruction')
    # print(gp)

    # dev.off()

    # gp<-qplot( slices, meanerr, shape=snr, linetype=snr, data=msesimnnb ) + geom_line() +
    #                    geom_errorbar(aes(ymin=meanerr-sdev, ymax=meanerr+sdev)) +
    #                    xlab('Number of slices (P)') +
    #                    ylab('MSE') +
    #                    ggtitle('Nearest neighbor reconstruction')
    # print(gp)

    if(saveplots)
    {
      cat('plots saved\n')
    }
    # dev.off()
}

