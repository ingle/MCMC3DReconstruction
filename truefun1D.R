
# obj function value for 2D case
# Feb 3, 2015
truefun1D <- function(xvec)
{
    ftrue <- rep( 0, length(xvec) )
    ftrue[ xvec^2/a^2 <= 1 ] <- 1
    as.vector(ftrue)
}
