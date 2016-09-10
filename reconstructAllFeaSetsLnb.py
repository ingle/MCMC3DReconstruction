'''
aingle Sat March 14, 2016
Calculate linear nnb interpolated volumes for all FEA volumes and save as csv
files FEA volumes should be first generated using fea_3d_generate.py
'''

from scipy.interpolate import griddata
from scipy.interpolate import interpn
import pandas as pd
import numpy as np
import itertools
import time

def LinearInterpolateFeaDataSet( slice_num ) :
  '''
  LinearInterpolateFeaDataSet( number_of_slices )
  Reads data from fea_sheaf_xyzfdata_#.csv datafile where # is the number_of_slices
  Saves resulting linear griddata interpolated volume to processed_fea_linearnb_ files
  '''
  Nx = 100/2
  Ny = 100/2
  Nz = 90/2

  #for slice_num in [4,6,12,16] :
  fname = 'fea_sheaf_xyzfdata_'+str(slice_num)+'.csv'

  print '---------------------------------------------------'
  print time.strftime('%H:%M:%S') + ': Processing '+fname

  # rest of this mimicks reconstructAllSetsNnb.R...
  tt = pd.read_csv(fname, sep=' ')
  print 'raw data read from .csv file'

  points = np.array( [np.array(tt.x), np.array(tt.y), np.array(tt.z)] ).T

  gx, gy, gz = np.mgrid[-1.9 : 1.9 : Nx*1j, -1.9 : 1.9 : Ny*1j, 1 : 4.5 : Nz*1j]

  gf_nearest         = griddata( points, np.array(tt.f), (gx,gy,gz), method='nearest' )
  print 'nnb interpolation done'

  gf_linear_unfilled = griddata( points, np.array(tt.f), (gx,gy,gz), method='linear' )
  print 'linear interpolation done' #note this takes the bulk of processing time >10mins for 10^6 points

  # outside the convex hull replace NaN's by nearest nbr interpolated values
  gf_linear = gf_linear_unfilled
  gf_linear[ np.isnan(gf_linear) ] = gf_nearest[ np.isnan(gf_linear) ]
  print 'nans fixed'

  gx_up, gy_up, gz_up = np.mgrid[-1.9 :1.9 :(2*Nx)*1j, -1.9 : 1.9 :(2*Ny)*1j, 1:4.5 : (2*Nz)*1j]
  gf_linear_up = interpn( (np.linspace(-1.9,1.9,Nx), np.linspace(-1.9,1.9,Ny), np.linspace(1,4.5,Nz)),
      gf_linear, (gx_up, gy_up, gz_up) )
  print 'upsampled by 2x to get to 100x100x90 gridsize'

  recons_df = pd.DataFrame( columns=('x','y','z','f') )
  recons_df['x'] = np.ravel(gx_up, order='F') # 'Fortran' ordering first index changes fastest
  recons_df['y'] = np.ravel(gy_up, order='F')
  recons_df['z'] = np.ravel(gz_up, order='F')
  recons_df['f'] = np.ravel(gf_linear_up, order='F')
  recons_df.to_csv( 'processed_fea_linearnb_'+str(slice_num)+'.csv', sep = ' ', index=False)

  print 'processed_fea_linearnb_'+str(slice_num)+'.csv written.'

  print time.strftime('%H:%M:%S')+' Done.'

if __name__ == '__main__' :
  for slice_num in [4,6,12,16] :
    LinearInterpolateFeaDataSet( slice_num )
