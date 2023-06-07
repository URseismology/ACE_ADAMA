# ACE_ADAMA

1. ADAMA_Models (On BlueHive) [/scratch/tolugboj_lab/Prj5_HarnomicRFTraces/2_Data/ADAMA_Models] <br>
INCLUDES: 
* ADAMA_1D_fixed1p7.mat: the ADAMA 1D models using a fixed vp/vs ratio = 1.7 (note some of the field is empty because there is no update).
* ADAMA_1D_fixedpoisson.mat: the ADAMA 1D models using poission ratio fixed (note some of the field is empty because there is no update).
* Litho_1D.mat: the Litho1.0 1D models (same order as the ADAMA_1D and the AfrDownSampleDots)
* AfrDownSampleDots.mat: location of the downsampled points in Africa, and the corresponding index in the map (afrmask; set how to use it in the plotting code fig7_VelMesmQlt.m) 

2. ADAMA_Maps  (On BlueHive) [/scratch/tolugboj_lab/Prj5_HarnomicRFTraces/2_Data/ADAMA_Maps] <br>
INCLUDES: 
* RCellVelDist: Rayleigh maps cell and velocity distribution (Fig. 3)
* LCellVelDist: Love maps cell and velocity distribution (Fig. 3)
* RayAvgMap: Rayleigh velocity maps, group and phase
* LoveAvgMap: Love velocity maps, group and phase
* AfrDownSampleDots.mat: location of the downsampled points in Africa, and the corresponding index in the map (afrmask; set how to use it in the plotting code fig7_VelMesmQlt.m) 

3. ADAMA_MCMC  (On Bluehive) [/scratch/tolugboj_lab/Prj5_HarnomicRFTraces/2_Data/ADAMA_MCMC] <br>
INCLUDES:
* (plotpara.mat)
*
*
*

4. Plotting Code <br>
INFO IN THE DIRECTORY

5. forwardADAMA.m: compute the dispersion velocity using ADAMA_1D models.
6. inverseADAMA.m: inverse the ADAMA_1D models (from dispersion curve measurements) using Lovee-Raylee code.
