'''
aingle Sat March 16, 2016
Calculate linear nnb interpolated volumes for
all datasets and save as csv files
'''

from scipy.interpolate import griddata
from scipy.interpolate import interpn
import pandas as pd
import numpy as np
import itertools
import time

def LinearInterpolateDataSet( set_num, slice_num ) :
  Nx = 100/2
  Ny = 100/2
  Nz = 90/2

  #for set_num in range(1,5+1) :
  #for slice_num in [4,6,12,16] :
  fname = 'xyzfdata_Oct26_2013_'+str(set_num)+'_'+str(slice_num)+'.dat'

  print '---------------------------------------------------'
  print time.strftime('%H:%M:%S') + ': Processing '+fname

  # rest of this mimicks reconstructAllSetsNnb.R...
  tt = pd.read_csv(fname, sep=' ')
  print 'raw data read from .dat file'

  x_scaled = np.array(tt.x) * 3.0/100;
  y_scaled = np.array(tt.y) * 3.0/100;
  f_scaled = np.array(tt.f)/0.8678

  print 'x,y,f rescaled'

  points = np.array( [x_scaled, y_scaled, np.array(tt.z)] ).T

  gx, gy, gz = np.mgrid[-2:2:Nx*1j, -2:2:Ny*1j, 1:78:Nz*1j]

  gf_nearest         = griddata( points, f_scaled, (gx,gy,gz), method='nearest' )
  print 'nnb interpolation done'

  gf_linear_unfilled = griddata( points, f_scaled, (gx,gy,gz), method='linear' )
  print 'linear interpolation done' #note this takes the bulk of processing time >10mins for 10^6 points

  # outside the convex hull replace NaN's by nearest nbr interpolated values
  gf_linear = gf_linear_unfilled
  gf_linear[ np.isnan(gf_linear) ] = gf_nearest[ np.isnan(gf_linear) ]
  print 'nans fixed'

  gx_up, gy_up, gz_up = np.mgrid[-2:2:(2*Nx)*1j, -2:2:(2*Ny)*1j, 1:78:(2*Nz)*1j]
  gf_linear_up = interpn( (np.linspace(-2,2,Nx), np.linspace(-2,2,Ny), np.linspace(1,78,Nz)),
      gf_linear, (gx_up, gy_up, gz_up) )
  print 'upsampled by 2x to get to 100x100x90 gridsize'

  recons_df = pd.DataFrame( columns=('x','y','z','f') )
  recons_df['x'] = np.ravel(gx_up, order='F') # 'Fortran' ordering first index changes fastest
  recons_df['y'] = np.ravel(gy_up, order='F')
  recons_df['z'] = np.ravel(gz_up, order='F')
  recons_df['f'] = np.ravel(gf_linear_up, order='F')
  recons_df.to_csv( 'processed_linearnb_'+str(set_num)+'_'+str(slice_num)+'.csv', sep = ' ', index=False)

  print 'processed_linearnb_'+str(set_num)+'_'+str(slice_num)+'.csv written.'

  print time.strftime('%H:%M:%S')+' Done.'

if __name__ == '__main__' :
  for set_num in range(1,5+1) :
    for slice_num in [4,6,12,16] :
      LinearInterpolateDataSet( set_num, slice_num )
