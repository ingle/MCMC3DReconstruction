#library(plot3D) # for 3d slice plots
source('full3drecons.R') # for 3d recons function

Nx <- 100
Ny <- 100
Nz <- 90

#nset <- 1 # data set number 1-5
#nsli <- 12 # number of slices in sheaf


for (nset in 1:1)
{
    for(nsli in c(4,6,12,16))
    {
        rm(list=ls()[! ls() %in% c("nset","nsli","Nx","Ny","Nz","full3drecons")])
        cat('cleared wspc. processing nset =',nset,'nslice =',nsli,'\n')
        filename <- paste('fea_sheaf_xyzfdata_',nsli,'.csv', sep='')
        tab <- read.table(filename, header=T, sep=' ')

        vec <- list( xvec=seq(-1.9,1.9,length.out=Nx),
                    yvec=seq(-1.9,1.9,length.out=Ny),
                    zvec=seq(1,4.5,length.out=Nz) )

        dat <- list( xdat=tab$x, ydat=tab$y, zdat=tab$z, fdat=tab$f )

        frecons <- full3drecons( dat, vec, wt=0.0, maxiter=2 )

        frecons3d <- array(frecons, c(Ny,Nx,Nz))
        frecons3d[frecons3d<0] <- 0
        frecons3d[frecons3d>4.5] <- 4.5

        save(list=ls(), file=paste('processed_fea_nnb_',nset,'_',nsli,'.RData',sep=''))

        cat('processed_fea_nnb_',nset,'_',nsli,'\n',sep='')

        # slice3D( vec$yvec,
        #         vec$xvec,
        #         vec$zvec,
        #         colvar=frecons3d[,,seq(dim(frecons3d)[3],1,by=-1)],
        #         theta=35, phi=20,
        #         facets=T,
        #         xs=0, ys=0, zs=vec$zvec[40])
    }
}
