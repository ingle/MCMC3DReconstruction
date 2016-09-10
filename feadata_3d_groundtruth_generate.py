'''
this script generates the ground truth fea data
for MSE calculations
'''

import scipy.io as sio
import numpy as np
import pandas as pd
import itertools
import scipy.signal as sp
from scipy.interpolate import griddata

if __name__=='__main__' :
  feaswv = sio.loadmat('/home/ingle/Insync/LabComputerMATLAB/feadatamedphyswv_all.mat')
  feaswv = feaswv['swvtrue']
  feaswv = sp.medfilt(feaswv)

  zvec = np.linspace(1,4.5,feaswv.shape[0])
  zvec_rep = np.array(np.repeat(np.matrix(np.array(zvec)).T, feaswv.shape[1], axis=1))

  allx = []
  ally = []
  allz = []
  allf = []
  sheaf_df = pd.DataFrame( columns=('x', 'y', 'z', 'f') )
  horizvec = np.linspace(-3.8/2, 3.8/2, feaswv.shape[1])
  nsli = 1024 
  for ang in np.arange(0, nsli-1+1)*np.pi/float(nsli) :
    xvec = horizvec * np.cos(ang) 
    yvec = horizvec * np.sin(ang)
    xvec_rep = np.array(np.matrix(xvec).repeat( feaswv.shape[0], axis=0 ))
    yvec_rep = np.array(np.matrix(yvec).repeat( feaswv.shape[0], axis=0 ))
    allx.append( np.ravel(xvec_rep) )
    ally.append( np.ravel(yvec_rep) )
    allz.append( np.ravel(zvec_rep) )
    allf.append( np.ravel(feaswv) )

  sheaf_df['x'] = np.array( list(itertools.chain(*allx)) )
  sheaf_df['y'] = np.array( list(itertools.chain(*ally)) )
  sheaf_df['z'] = np.array( list(itertools.chain(*allz)) )
  sheaf_df['f'] = np.array( list(itertools.chain(*allf)) )

  Nx, Ny, Nz = 100, 100, 90
  gx, gy, gz = np.mgrid[-1.9 :1.9 :(Nx)*1j, -1.9 : 1.9 :(Ny)*1j, 1:4.5 : (Nz)*1j]
  points = np.array( [np.array(sheaf_df.x), np.array(sheaf_df.y), np.array(sheaf_df.z)] ).T
  gf = griddata( points, np.array(sheaf_df.f), (gx,gy,gz), method='nearest')

  mask = np.ones(gf.shape)
  # remove needle
  mask[48:51,48:51,:] = 0
  mask[:,:,0:28] = 0
  mask[:,:,74:]  = 0


  ground_truth_df = pd.DataFrame( columns=('x','y','z','f', 'mask') )
  ground_truth_df['x'] = np.ravel(gx, order='F')
  ground_truth_df['y'] = np.ravel(gy, order='F')
  ground_truth_df['z'] = np.ravel(gz, order='F')
  ground_truth_df['f'] = np.ravel(gf, order='F')
  ground_truth_df['mask'] = np.ravel(mask, order='F')
  ground_truth_df.to_csv( 'groundtruth_fea.csv', sep=' ', index=False )

  print 'ground truth data saved'

