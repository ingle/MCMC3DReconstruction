# full3drecon( dat, vec, wt=0.01, maxiter=1000 )
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

fvec <- function( dat, vec, wt=0.01, maxiter=1000 )
{

}
