# Swimming analysis in agar

This collection of Python and Matlab code for single-cell detection, cell-tracking and the extraction of movement statistics (diffusion and drift). This code was used to analyze the movement characteristics of cells in expanding poputlations (within soft-agar), as published in our manuscript "Chemotaxis as a navigation strategy to thrive in nutrient replete environments”.

Ying Tang and Jonas Cremer, together with Tomoya Honda, Jerome Wong-Ng, Massimo Vergassola, Terence Hwa, 2019.

## Required packages

The code combines Python and MATLAB scripts.

For cell tracking the u-track package from the Danuser Lab was used:

- utrack 2.0 - Copyright (C) 2007 Danuser Lab. See https://github.com/DanuserLab/u-track.

In addition, the Python code requires the modules NumPy an SciPy. 

The code was developed and tested in Python2.7 and MATLAB R2016b and MATLAB R2018b.

## Running analysis

Analysis works in 3 steps:
1. single-cell detection, 
2. cell tracking (u-track or trackpy)
3. estimation of diffusion and drift. 

The trajectory generation via u-track reads in .mat files and thus a Python script first digests the raw imaging files from the microscopes and detects cells. Positions of detected cells are saved in a .mat file and analyzed using u-track. 

For the analysis of swimming behavior in liquid culture we used the python module trackpy. This can also be used as an alternative to track cells in soft-agar. 
- trackpy 0.3 0  See https://github.com/soft-matter/trackpy
See the alogrithms in “gp_swimming.py”.

To run this code, execute python file readindata_PIP.py. Before running, settings in the readindata_PIP.py and DriftDiffusion.m have to be adjusted.


## Adjust readindata_PIP.py:

1. filename: folder names that contain images.

2. SquenceName: label of experimental runs.

3. thresholdv: the lower cutoff of signal value in each pixel to detect cell. 

4. generate_cluster: if True, it will detect single cells from GFP images of 
confocal microscope. The output is a matrix file of cells' 2 dimensional 
positions for each image. 

5. do_tracking: if True, it will run Matlab package utrack 2.0 to connect 
the cells between images, track the cells and generate the trajectoreis.

6. CalculateDiffusion: if True, run Matlab file DriftDiffusion.m to calculate 
diffusion, drift, etc of cell population dynamics. It estimate mean square 
displacment versus time, in X and Y direction, and do a linear fitting to the it 
for calculating the diffusion and drift coefficients.

## Adjust DriftDiffusion.m:

- PlotDetail: if 1, it will plot the distribution of single-cell statistics of mean 
displacement in both x and y direction.

- PlotDriftDiffusionFit: if 1, it will plot the linear fit between mean square 
displacement and mean displacement versus time.

- ThrowNonMotile: if 1, it will delete non-motile cells. 

- PhysicalLength: length size of images in the unit of mum.

- pixel: pixel number of images.

- tLength: time step between images in the unit of seconds.


## Output

The main output includes:

1. Figures:

- The distribution of mean displacement and mean square displacement 
of single-cells, in both x and y directions.

- The average of mean displacement and mean square displacement versus
time, for both x and y directions. The linear fitting to the latter gives the slope 
corresponding to the diffusion coefficients;

- The auto-correlation function of position variables, in both x and y directions. 

- The distribution of angular change. 

2. Mat file with the statistics:

- The mean displacement and mean square displacement versus time.

- The diffusion coefficient as the slope of linear fitting on mean square 
displacement versus time.

