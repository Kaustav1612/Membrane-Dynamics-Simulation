# Membrane-Dynamics-Simulation

These scripts are written during my masters in IISER K at Cell Biophysics Lab, primarily related to cell imaging object detection and object tracking.

These scripts are written for a particular kind of analysis work, however these scripts can be tweaked to be used for analysis of biological images.

For any assistance feel free to contact me 

Email: kaustavpal1612@gmail.com
X: https://x.com/KaustavPal1612
Bluesky: https://bsky.app/profile/kaustavpal1612.bsky.social
LinkedIn: [https://www.linkedin.com/in/kaustav-pal-b524a5223/](https://www.linkedin.com/in/kaustavpal1612/)

Keep all the constituent functions in an individual folder all in the same directory 

Folder 1: IRM_Fitting_Plots

RICM data after calibration of intensity to heights, can be further processed here

spectra_fit.m : MATLAB program to fit the experimental Power Spectra to an existing model, model can be altered, here we use only two parameter namely the active timescale and the effective viscosity to fit the power spectra
in ultra low to mid frequency. The excess area from the mode dependent activity provides the active amplitude and tension which is used to fit. (Needs FBR Map coordinates as inputs along with Imin and conversion factor)

Mode_Dependent_PSD.m : Look at the static mode dependent PSD to look at the dominant active force polarity
needs modewise_psd_analysis.m function to work keep them in the same directory

FBR_map.m : Used to plot color plot of the different mechanical properties of the FBR overlaid on the interference image, needs AA files as well as the mechanical parameters outputs from the fits.

main_simulation.m : The main simulation MATLAB file to run active as well as passive membrane sheets of definitve size and resolution with known mechanical parameters ranging from tension to pinning density, also can be used to 
test out fluctuations spectra model to try fitting on. 
