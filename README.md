# ACE_ADAMA

1. ADAMA_Models [BH: /scratch/tolugboj_lab/Prj5_HarnomicRFTraces/2_Data/ADAMA_Models] <br>
INCLUDES: 
* ADAMA_1D_fixed1p7.mat: the ADAMA 1D models using a fixed vp/vs ratio = 1.7 
* ADAMA_1D_fixedpoisson.mat: the ADAMA 1D models using poission ratio fixed (note: both ADAMA_1D has empty fields because there are no model updates for some locations; both inversion processes use chi-square range = [0.3, 0.5] and msigmaf = 2).
* Litho_1D.mat: the Litho1.0 1D models (same order as the ADAMA_1D and the AfrDownSampleDots)
* AfrDownSampleDots.mat: location of the downsampled points in Africa, and the corresponding index in the map (afrmask; set how to use it in the plotting code fig7_VelMesmQlt.m) 
* ADAMA_Uncer.mat: XXXXXXX (not ready yet!)

2. ADAMA_Maps  (On BlueHive) [BH: /scratch/tolugboj_lab/Prj5_HarnomicRFTraces/2_Data/ADAMA_Maps] <br>
INCLUDES: 
* RCellVelDist: Rayleigh maps cell and velocity distribution (Fig. 3)
* LCellVelDist: Love maps cell and velocity distribution (Fig. 3)
* RayAvgMap: Rayleigh velocity maps, group and phase
* LoveAvgMap: Love velocity maps, group and phase
* AfrDownSampleDots.mat: location of the downsampled points in Africa, and the corresponding index in the map (afrmask; set how to use it in the plotting code fig7_VelMesmQlt.m) 

3. ADAMA_MCMC  (On Bluehive) [BH: /scratch/tolugboj_lab/Prj5_HarnomicRFTraces/2_Data/ADAMA_MCMC] <br>
INCLUDES:
* RayleighPData: Rayleigh phase velocity map (contains cell locations and cell velocity) ensembles for 11 periods(plotpara.mat contains the lat & lon grid of the velocity map).
* RayleighGData: Rayleigh group velocity map ensembles for 11 periods.
* LovePData: Love phase velocity map ensembles for 11 periods.
* LoveGData: Love group velocity map ensembles for 11 periods.

4. Plotting Code <br>
INFO IN THE DIRECTORY

5. forwardADAMA.m: compute the dispersion velocity using ADAMA_1D models.
6. inverseADAMA.m: inverse the ADAMA_1D models (from dispersion curve measurements) using Lovee-Raylee code (note it uses the function raylee_invert.m).


7. Table 4: Calculation of values in Table 4
   
