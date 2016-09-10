full3drecons <- function( dat, vec, wt=0.01, maxiter=1000, minchange=1e-5, verbose=FALSE )
{
    # full3drecons( dat, vec, wt=0.01, maxiter=1000, minchange=1e-5 )
    # 
    # dat contains 4 fields $xdat, $ydat, $zdat, $fdat for
    # noisy data fdat = f(xdat, ydat, zdat)+noise
    #
    # vec contains 3 fields $xvec, $yvec, $zvec which defines
    # the rectangular 3d grid where f is to be reconstructed
    # 
    # Note: x refers to cols
    #       y refers to rows
    #       z refers to slices
    # Vectorization in this code is done by reading columns top
    # to bottom for each slice.
    #
    # Dependencies: Matrix, VGAM, RANN


    library('VGAM')
    library('Matrix')
    library('RANN')

    Ny = length(vec$yvec)
    Nx = length(vec$xvec)
    Nz = length(vec$zvec)

    if(verbose)
    {
        cat('creating initial state using nnb\n')
    }
    grid <- expand.grid( vec$yvec, vec$xvec, vec$zvec )
    cloud <- cbind( dat$ydat, dat$xdat, dat$zdat )
    nnres <- nn2( data = cloud, query = grid, k=1 )
    x0init = dat$fdat[nnres$nn.idx]
    if(verbose)
    {
        cat('initial state created\n')
    }

    dy = mean( diff(vec$yvec) ) # row
    dx = mean( diff(vec$xvec) ) # col
    dz = mean( diff(vec$zvec) ) # sli
    lamdy2 <- wt/dy/dy
    lamdx2 <- wt/dx/dx
    lamdz2 <- wt/dz/dz

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

    if(verbose)
    {
        cat('updmat created\n')
    }

    x0 <- as.vector(x0init)

    for ( i in seq(1,maxiter) )
    {
            
        x1 <- ( as.vector(x0init) + UpdtMat %*% x0 ) / (1 + 4*(lamdx2+lamdy2+lamdz2) )

        if( i%%50==0 && verbose )
            cat('[iter ',i,']','|x-xprev|=',norm(x1-x0), '\n', sep='')

        if( norm(x1-x0) < minchange && verbose )
        {
            cat('convergence criterion met\n')
            break
        }
        
        x0 <- x1
        
        Arr <- array( x1, c(Ny, Nx, Nz) )
        if(verbose)
        {
            cat('image displayed ',i,' \n')
            image(x=vec$xvec, y=vec$yvec, z=t(matrix(Arr[,,50], nrow=Ny)))
            title(main=paste('iter ',i))
        }

    }

    if (verbose) {
      cat('max iterations complete\n')
    }
    
    x1
}

