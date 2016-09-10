* 4/14/2015

Processing chain:
shear wave velocity images are already available in Oct13 datasets called
multislicedata in the SheafData folder.

writexyzfdata() in Matlab
     |
     |
     V
read xyzfdata .dat file
form vec grid
     |
     |
     V
call full3drecons in R


* 4/21/2015
Simulated data from fabricated ellipsoid inclusion 

Simulation results are generated at different values of slice numbers P,
different SNR levels and using our method aand a nearest neighbor interpolator.
Data is stored in the following format:

mse 2     x   100   x     4        x    3
   alg       NSIMS    P(nslices)       SNR levels (5,10,20)
  full3d
   nnb

Note that SNR level of 15 was separately simulated with results in its own
RData file.


* Apr 22, 2015
I assembled all snr data into one RData file called simulatedmseresults.RData

Note that I had run two separate simulations previously.
The 15dB SNR data is in simulatedmseresults15dB.RData
5, 10, 20dB SNR data is in simulatedmseresults51020dB.RData

* May 8, 2015
xyzfdata files are already processed and results saved in RData files.
The relevant script is reconstructAllSets.R

* May 29, 2015
Plots are created and saved by plotsandtables.R 
Image statistics are calcualted using calstats.R (called as part of
plotsandtables.R)

* Mar 12, 2016
New python script: griddata_test.py calculates griddata linear interpolated
volume.  It plays a small trick for processing speed reasons. Instead of
interpolating to 100x100x90 grid, the griddata function operates on a smaller
grid of 1/2^3 the size and then we use interpn to upsample to the 100x100x90
grid size. Volume recons are stored as .csv files called
"processed_linearnb_1_4.csv" etc. Also note that the 3D indices must be
vectorized using the numpy ravel Fortran option so that the first index goes
fastest (instead of the C option which does the opposite).  The .csv files can
be read into R using read.table for further processing.

Renamed griddata_test.py to reconstructAllSetsLnb.py.

* Mar 14, 2016
feadata_3d_generate.py creates four csv files with simulated sheaves using
single plane of FEA reconstructed data. These csv files have the same format as
the xyzf dat files (four cols).
Two new R scripts reconstructAllFeaSets.R and reconstructAllFeaSetsNnb.R.
New python script reconstructAllFeaSetsLnb.py

* Jul 19. 2016
In order to generate MSE for FEA datasets we need the ground truth FEA in 3D.
We generate it using feadata_3d_groundtruth_generate.py. Basically read in the
true SWV map from the mat file, then use many radial slices (1024) and take NNB
on a grid 100x100x90. Treat this as ground truth. A mask is also included to
remove regions above and below the inclusion and around the needle.

* Sep 5, 2016
Now simulating smooth transition boundary instead of an abrupt boundary. We do
this by modeling the transition as a smooth (sigmoid) function that decreases
from 4 to 1 going from the inclusion into the background. The transition region
is 5 mm wide (from R. DeWall's nano-indentation paper) and it reaches 99% of
the saturation values within the transition width. (See sigmoid function in
getSmoothEllipsoidSheafData.R)

Correspondingly, we also created a new script simulated_smooth_mse.R which is
similar to simulatedmse.R except that it uses smooth transition boundary
ellipsoid instead of the abrupt change model.

CRAZY CRAZY ERROR FIXED. All my SWV values are off by a factor because instead
of multiplying by 0.868 I divided them by 0.868. So I have now corrected the
calcstats.R script to multiply ROI values by 0.868*0.868.
