getEllipseData <- function(a,b,limX, limY, N)
# Samples an ellipse function given by:
# fdat = 1 if xdat.^2/a^2 + ydat.^2/b^2 <= 1
# else fdat = 0.
#
# Args:
# limX, limY are box limits along x,y
# a,b are the three semiprincipal axis lengths along x,y
# N number of data points to produce
#
# Returns:
# xdat: xcoor of data points
# ydat: ycoor of data points
# fdat: ellipsoid data f(xdat, ydat)
#
# v1 Jan 10, 2015
#
{
    set.seed(10);

    stopifnot(a<limX, b<limY)
    xdat <- runif( N, -limX, limX )
    ydat <- runif( N, -limY, limY )
    fdat <- rep(0, N)

    fdat[ xdat^2/a^2 + ydat^2/b^2 <= 1 ] <- 1

    dat <- data.frame( xdat=xdat, ydat=ydat, fdat=fdat );
}


