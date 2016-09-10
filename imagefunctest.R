# this is an educational script on how image plots a matrix

Nx = 30 # ncols
Ny = 20 # nrows

# for now forget about columns and rows. Just do x first, y next.
ftrue = array(0, c(Nx,Ny)) 
xvec <- seq(-3,3,,Nx)
yvec <- seq(-2,2,,Ny)
xygrid <- expand.grid(xvec, yvec)
colnames(xygrid) <- c('xgrid', 'ygrid')

ftrue[ xygrid$xgrid^2/2^2 + xygrid$ygrid^2/1^2 <=1 ] = 1

image(xvec, yvec, ftrue) #this will do the right thing: short fat ellipse

# ftrue is 30 rows 20 cols
# but image() will render the transpose of ftrue
# note that this is different from how imagesc() behaves
# in Matlab

#  Now treat Nx as ncol and Ny as nrows
ftrue = array(0, c(Ny,Nx))
xvec <- seq(-3,3,,Nx)
yvec <- seq(-2,2,,Ny)
xygrid <- expand.grid(yvec, xvec)
colnames(xygrid) <- c('ygrid', 'xgrid')

ftrue[ xygrid$xgrid^2/2^2 + xygrid$ygrid^2/1^2 <=1 ] = 1

dev.new()
image(yvec, xvec, ftrue) #this will show the transpose, tall skinny ellipse

dev.new()
image(xvec, yvec, t(ftrue)) # fix


