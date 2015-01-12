makeAB <- function( vec, dat, verbose=FALSE )
{
    if(verbose)
    {
        print('in func makeAB')
    }
    stopifnot( length(dat$xdat)==length(dat$ydat) )
    Nx <- length( vec$xvec )
    Ny <- length( vec$yvec )
    N  <- length( dat$xdat )

    tmpX <- hist(dat$xdat, breaks=vec$xvec, plot=F)
    tmpY <- hist(dat$ydat, breaks=vec$yvec, plot=F)
    indX <- findInterval(dat$xdat, tmpX$breaks)
    indY <- findInterval(dat$ydat, tmpY$breaks)
    indov <- (indX==Nx)
    indX[indov] <- indX[indov]-1
    indov <- (indY==Ny)
    indY[indov] <- indY[indov]-1
    
    ind <- indY + Ny * (indX-1)
    delX <- diff(vec$xvec)
    delY <- diff(vec$yvec)
    tx = (dat$xdat - vec$xvec[indX])/delX[indX];
    ty = (dat$ydat - vec$yvec[indY])/delY[indY];

    library('Matrix')

    ik <- as.vector( rep(seq(1,N),4) )
    jk <- as.vector( c(ind, ind+1, ind+Ny, ind+Ny+1))
    xk <- as.vector( c( (1-tx)*(1-ty), (1-tx)*ty, tx*(1-ty), tx*ty ) )

    A <- sparseMatrix( i=ik, j=jk, x=xk, dims=c(N,Nx*Ny) )
    
    if(verbose)
    {
        print('makeAB: A created')
    }

    alljj <- seq(1,Ny*Nx)
    tmpy <- (alljj-1)%%Ny + 1
    tmpx <- floor((alljj-1)/Ny)+1

    ik <- as.vector( c(seq(1,Ny*Nx), 
                     alljj[tmpy-1>=1],
                     alljj[tmpy+1<=Ny],
                     alljj[tmpx-1>=1],
                     alljj[tmpx+1<=Nx])
                   )

    siz <- c(Ny,Nx)
    mindel <- min( delX[1], delY[1] )
    jk <- as.vector( c(seq(1,Ny*Nx),
                           sub2ind(siz, tmpy[tmpy-1>=1]-1, tmpx[tmpy-1>=1]),
                           sub2ind(siz, tmpy[tmpy+1<=Ny]+1,tmpx[tmpy+1<=Ny]),
                           sub2ind(siz, tmpy[tmpx-1>=1], tmpx[tmpx-1>=1]-1),
                           sub2ind(siz, tmpy[tmpx+1<=Nx], tmpx[tmpx+1<=Nx]+1)
                          )
                   )
    xk <- as.vector( c( rep(-2*(mindel/delX[1])^2 -2*(mindel/delY[1])^2, Ny*Nx),
                            rep((mindel/delY[1])^2 ,sum(tmpy-1>=1)),
                            rep((mindel/delY[1])^2, sum(tmpy+1<=Ny)),
                            rep((mindel/delX[1])^2, sum(tmpx-1>=1)),
                            rep((mindel/delX[1])^2, sum(tmpx+1<=Nx))
                          )
                   )

    B <- sparseMatrix( i=ik, j=jk, x=xk, dims=c(Ny*Nx,Ny*Nx) )
    if(verbose)
    {
        print('makeAB: B created')
    }
    AB <- list(A=A,B=B)
}

sub2ind <- function( siz, r, c )
{
    (c-1)*siz[1] + r
}
