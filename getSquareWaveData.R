getSquareWaveData <- function(a,limX,N)
# Samples a square wave function given by:
# fdat = 1 if -a < xdat < a 
# else fdat = 0.
#
# Args:
# limX, is box limit along x
# a, is half the box width around 0 
# N number of data points to produce
#
# Returns:
# xdat: xcoor of data points
# fdat: square pulse data = f(xdat)
#
# v1 Jan 23, 2015
#
{
#    set.seed(10);

    stopifnot(a<limX)
    xdat <- runif( N, -limX, limX )
    fdat <- rep(0, N)

    fdat[ xdat^2/a^2 <= 1 ] <- 1

    dat <- data.frame( xdat=xdat, fdat=fdat );
}


