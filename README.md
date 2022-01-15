# LCLS_Analysis
Data analysis for X-ray emission spectroscopy at LCLS in 2016.  
Runs on LCLS high performance computing clusters.  
See https://pubs.acs.org/doi/10.1021/acs.jpclett.8b03595 for more information.  
Note: This may not be the full final codebase for the analysis in the paper above, this is what I still have many years after publishing.  
  
## Requirements:
Python 2.7.11  
Psana https://confluence.slac.stanford.edu/display/PSDMInternal/psana+-+Reference+Manual ana-0.17.31  
Numpy 1.11  

## Description  
The data are processed in 4 steps:  
1. Extracts the raw detector information:  
    a. X-ray intensity  
    b. The x-ray diffraction data (mainly for determining x-ray intensity and liquid jet interaction volume)  
    c. The x-ray area detector (csPad or epix) that records the x-ray emission spectra  
        - Both raw detector units and calibrate photon number per pixel are recorded per each shot
    d. XTCAV its a new machine that will take the elecron bunch from the linear accelerator and deflect them upward and record an EM intensity of sorts  
        - This gives indirect information on the x-ray pulse duration, structure, x-ray photon energy, and magnitude  
    e. Other machine variables (time of the shot, realtive intensity etc)  
    
2.

## Experimental Setup in Brief
1. A ~10 femtosecond x-ray pulse would interact with a liquid jet.
2. The atoms would ionize the atoms (Manganese atoms in our case) and re-emit x-rays through a 3p electron filling the vacancy in the 1s (K-beta emission lines)
3. The emission line is dispersed through a spectrometer onto a cs-pad detector so that the spatial position determines the x-ray energy
4. The detector records the energy absorbed during emission (very poor energy resolution and is only used to classify the number of photons absorbed)
    4a. This prevents electronic noise being read as a photon (x-ray) and allows for multiphoton detection (ie identify how many photons were absorbed before readout)
    4b. The distribution of the "energy" read by each pixel results in minimally overlaping gaussians to isolate each peak for the 0,1 etc photon cases
5. The data are stored for each shot which occured on a 120Hz time scale
