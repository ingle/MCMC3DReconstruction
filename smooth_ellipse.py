import numpy as np
import matplotlib.pyplot as plt

def sigmoid( x, alph=1.0, c=0.0, scale=1.0, offset=0.0 ) :
  return offset + scale * (1 - 1.0/(1.0 + np.exp(-alph*(x-c))))

if __name__ == '__main__' :
  c = 1.0 # change point
  alph = -1/0.25 * np.log( 1/0.99 - 1 ) # 99% transition region in +/-2.5 units of the change point
  scale=3.0
  offset=1.0
  
  x = np.linspace(-5,5,500)
  y = np.linspace(-3,3,250)
  xg, yg = np.meshgrid( x, y, sparse=False, indexing='xy' )
  fg = np.zeros(yg.shape)
  a=2.0
  b=1.0
  for i in range(xg.shape[0]) :
    for j in range(xg.shape[1]) :
      scaled_dist = ( xg[i,j]**2.0/a**2.0 + yg[i,j]**2.0/b**2.0 )
      fg[i,j] = sigmoid( scaled_dist, alph, c, scale, offset )

  



