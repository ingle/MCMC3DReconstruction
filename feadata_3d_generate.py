'''
In this script we will simulate FEA 3D data using single plane
recons from feadataswv_all and replicating it with noise in
4,6,12,16 planes

 /y
/
---- x
|
|z

'''

import scipy.io as sio
import numpy as np
import pandas as pd
import itertools

feaswv = sio.loadmat('/Users/aingle/Insync/LabComputerMATLAB/fea_data/feadatamedphyswv_all.mat')
feaswv = feaswv['swvpsm']
(r,c) = feaswv.shape
# interp nan's with nearby values
feaswv = np.ravel(feaswv)
mask = np.isnan(feaswv)
feaswv[mask] = np.interp( np.flatnonzero(mask), np.flatnonzero(~mask), feaswv[~mask] )
feaswv = np.reshape( feaswv, (r,c) )

zvec = np.linspace(1, 4.5, feaswv.shape[0])
zvec_rep = np.array(np.repeat(np.matrix(np.array(zvec)).T, feaswv.shape[1], axis=1))
horizvec = np.linspace(-3.8/2, 3.8/2, feaswv.shape[1])

for nsli in [4,6,12,16] :
  allx = []
  ally = []
  allz = []
  allf = []
  sheaf_df = pd.DataFrame( columns=('x', 'y', 'z', 'f') )
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
  sheaf_df.to_csv( 'fea_sheaf_xyzfdata_' + str(nsli) + '.csv', sep=' ', index=False )
  print( 'fea_sheaf_xyzfdata_' + str(nsli) + '.csv saved.')


