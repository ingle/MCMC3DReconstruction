getSmoothEllipsoidSheafData <- function(P, RhoMax=2, NRho=50, ZMax=4.5, NZ=100, a=1, b=1, c=1.5, zshift=2.25)
{
    Rho <- seq(0, RhoMax, length.out=NRho)
    Theta <- seq(0, 2*P-1)*pi/P
    Zs = seq(0, ZMax, length.out=NZ)

    Msh <- list(RhoM=outer(Theta*0, Rho, FUN='+'), ThetaM=outer(Theta,Rho*0,FUN='+'))
    RhoMVec <- as.vector(Msh$RhoM)
    ThetaMVec <- as.vector(Msh$ThetaM)

    XMVec <- RhoMVec * cos(ThetaMVec)
    YMVec <- RhoMVec * sin(ThetaMVec)
    ZMVec <- rep(seq(0,ZMax,length.out=NZ))
    ZMVec <- rep(ZMVec, each=length(XMVec))
    XMVec <- rep(XMVec, times=NZ)
    YMVec <- rep(YMVec, times=NZ)

    fdat <- rep(1, length(ZMVec))
    #fdat[XMVec^2/a^2 + YMVec^2/b^2 + (ZMVec-zshift)^2/c^2 <=1] <- 4

    cpt = 1.0 # change point at 1 because that's how you define the equation of an ellipsoid
    alph = -1/0.25*log(1/0.99-1) # 99% transition region in +/-0.25cm of the change point
    scale = 3.0  
    offset = 1.0 # ensure inc=4, bkg=1
    for( i in seq(length(fdat)) )
    {
        scaled_square_dist <- XMVec[i]^2/a^2 + YMVec[i]^2/b^2 + (ZMVec[i]-zshift)^2/c^2
        fdat[i] <- sigmoid( scaled_square_dist, alph, cpt, scale, offset )
    }

    dat <- data.frame( xdat=XMVec, ydat=YMVec, zdat=ZMVec, fdat=fdat )
}

sigmoid <- function( x, alph=1.0, c=0.0, scale=1.0, offset=0.0 )
{
    val <- offset + scale * (1.0 - 1.0/(1.0 + exp(-alph*(x-c))))
}
