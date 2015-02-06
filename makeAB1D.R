# creates A and B matrices
# A is the interpolation matrix
# B is a smoothness (Laplacian operator) matrix
#
# Tip: this is very similar to the 2D version of
# the code in makeAB() except, we have gotten
# rid of Ny, and X plays the role of Y.

makeAB1D <- function( vec, dat, verbose=FALSE )
{
    if(verbose)
    {
        print('in func makeAB_1D')
    }
    Nx <- length( vec$xvec )
    N  <- length( dat$xdat )

    tmpX <- hist(dat$xdat, breaks=vec$xvec, plot=F)
    indX <- findInterval(dat$xdat, tmpX$breaks)
    indov <- (indX==Nx)
    indX[indov] <- indX[indov]-1
    
    ind <- indX 
    delX <- diff(vec$xvec)
    tx = (dat$xdat - vec$xvec[indX])/delX[indX];

    library('Matrix')

    ik <- as.vector( rep(seq(1,N),2) )
    jk <- as.vector( c(ind, ind+1))
    xk <- as.vector( c( (1-tx), tx ) )

    A <- sparseMatrix( i=ik, j=jk, x=xk, dims=c(N,Nx) )
    
    if(verbose)
    {
        print('makeAB_1D: A created')
    }

    alljj <- seq(1,Nx)
    tmpx <- (alljj-1)%%Nx + 1

    ik <- as.vector( c( seq(1,Nx), 
                     alljj[tmpx-1>=1],
                     alljj[tmpx+1<=Nx] )
                   )

    if(verbose)
    {
        print('makeAB_1D: ik created for B')
    }
    siz <- c(Nx)
    mindel <- mean( delX )
    jk <- as.vector(     c(
                           seq(1,Nx),
                           tmpx[tmpx-1>=1]-1,
                           tmpx[tmpx+1<=Nx]+1
                          )
                   )

    if(verbose)
    {
        print('makeAB_1D: jk created for B')
    }
    
    xk <- as.vector(     c( 
                            rep(-2*(mindel/delX[1])^2 , Nx),
                            rep((mindel/delX[1])^2 ,sum(tmpx-1>=1)),
                            rep((mindel/delX[1])^2, sum(tmpx+1<=Nx))
                          )
                   )

    if(verbose)
    {
        print('makeAB_1D: xk created for B')
    }

    B <- sparseMatrix( i=ik, j=jk, x=xk, dims=c(Nx,Nx) )
    if(verbose)
    {
        print('makeAB_1D: B created')
    }
    AB <- list( A=A, B=B )
}

# sub2ind <- function( siz, r, c )
# {
#     (c-1)*siz[1] + r
# }
