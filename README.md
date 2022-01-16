# LCLS_Analysis
Data analysis for X-ray emission spectroscopy at LCLS in 2016.  
Runs on LCLS high performance computing clusters.  
See https://pubs.acs.org/doi/10.1021/acs.jpclett.8b03595 for more information.  
Note: This is not be the full final codebase for the analysis in the paper above, but the code for analyzing the data to generate the graphs.  
  
## Requirements:
Python 2.7.11  
Psana https://confluence.slac.stanford.edu/display/PSDMInternal/psana+-+Reference+Manual ana-0.17.31  
Numpy 1.11  

## Description  
The data are processed in 4 steps:  
1. Extracts the raw detector information:
   - X-ray intensity  
   - The x-ray diffraction data (mainly for determining x-ray intensity and liquid jet interaction volume)  
   - The x-ray area detector (csPad or epix) that records the x-ray emission spectra (Both raw ADU counts and photon counts) 
   - XTCAV which gives indirect information on the x-ray pulse duration, structure, x-ray photon energy, and magnitude  
   - Other machine variables (time of the shot, realtive intensity etc)  
    
2. Extracts the emission spectra from the detector images (Calibrating, handling bad pixels, and integrating signal)

3. Bins the emission spectra by experimental conditions including: x-ray pulse width, duration, and photon energy.

4. Analyzes the spectra

## Experimental Setup in Brief
1. A ~10 femtosecond x-ray pulse would interact with a liquid jet.  
2. The atoms would ionize the atoms and re-emit x-rays from Manganese atoms via 3p to 1s electron transition (K-beta emission lines).  
3. The emission line is dispersed through a spectrometer onto a cs-pad detector so that the spatial position determines the x-ray energy.  
4. The detector records the emission spectrum from x-ray liquid jet interaction.   
5. The data are stored for each shot which occured on a 120Hz time scale.  
