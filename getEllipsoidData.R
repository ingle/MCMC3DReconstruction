getEllipsoidData <- function(a,b,c,limX,limY,limZ, N)
# Samples an ellipse function given by:
# fdat = 1 if xdat.^2/a^2 + ydat.^2/b^2 + zdat.^2/c^2 <= 1
# else fdat = 0.
#
# Args:
# limX, limY, limZ are box limits along x,y,z
# a,b,c are the three semiprincipal axis lengths along x,y,z
# N number of data points to produce
#
# Returns:
# xdat: xcoor of data points
# ydat: ycoor of data points
# zdat: zcoor of data points
# fdat: ellipsoid data f(xdat, ydat, zdat)
#
# v1 Feb 17, 2015
#
{
    set.seed(10);

    stopifnot(a<limX, b<limY, c<limZ)
    xdat <- runif( N, -limX, limX )
    ydat <- runif( N, -limY, limY )
    zdat <- runif( N, -limZ, limZ )
    fdat <- rep(0, N)

    fdat[ xdat^2/a^2 + ydat^2/b^2 +zdat^2/c^2 <= 1 ] <- 1

    dat <- data.frame( xdat=xdat, ydat=ydat, zdat=zdat, fdat=fdat );
}


