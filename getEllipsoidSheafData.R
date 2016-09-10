getEllipsoidSheafData <- function(P, RhoMax=2, NRho=50, ZMax=4.5, NZ=100, a=1, b=1, c=1.5, zshift=2.25)
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
    fdat[XMVec^2/a^2 + YMVec^2/b^2 + (ZMVec-zshift)^2/c^2 <=1] <- 4

    dat <- data.frame( xdat=XMVec, ydat=YMVec, zdat=ZMVec, fdat=fdat )
}
