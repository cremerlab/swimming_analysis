# swimming_analysis
Analysis of movies to extract bacterial swimming behavior and statistics

We used this code to analyze the swimming behavior of cells within liquid. Method details and the biological context are provided in our  manuscript:
- J.Cremer, T.Honda, Y.Tang, J.Wong-Ng, M.Vergassola, T.Hwa. Chemotaxis as a navigation strategy to thrive in nutrient-replete environments

August 2019, Jonas Cremer and all coauthors.

## Required packages
This code runs with Python 2.7. Required modules include NumPy and SciPy.

For trajectory detection the tracking algorithm and python module trackpy was used:
- trackpy 0.30 .  See https://github.com/soft-matter/trackpy

## Taking movies of swimming bacteria
Collect movies of swimming bacteria using phase contrast microscopy. Typical settings are a 10x objective, and image aquisition with a framerate of 20FPS. Observation is done with a properly diluted culture of cells ensuring the proper distinction of different cells while also allowing also for a good swimming statistics (typical OD600 value is 0.005). Dilution in preheated medium, observation within a temperature controlled enviornment. More details are provided in our work. Note the following when saving the aquired data:
- Script reads in acquired movie files, so save data as .avi file after aquisition.
- To make handling easier, all movies stored within one sub-folder are recognized by the script and analyzed sequentially. 
- Please store movies using a unique name (even if they are in different sub-folders), for example by having the date as part of the file name. All output files are stored within one central folder using the unique name as a prefix to label the specific movie.

## Analyze swimming behavior
The analysis is done in 3 consequetive steps:
- Detection of cells for each frame of the movie.
- Generation of trajectories. 
- Analysis of trajectories and output of swimming behavior (average swimming speed and additional statistics).
The details of trajectory generation are introduced in the Supplementary Text of our work. 

## Run script
Use "swimming.py" to run script. 

- Before running the first time, adjust the names of the folders raw movies are stored and where results should be stored, see beginning of the file. Same adjustments have to be made for module file gp_swimming.py as well. 
- Before running the first time, adjust the length calibration and the timing between two frames by adjusting gp_swimming.py. In particular adjust the 
- Consequetively run the different parts of the swimming.py script, as explained in the comments. Before running, adjust file names and detection and analysis parameters according to your file names and imageing settings.
- When running the script the first time, you can check for successfull detection by ploting movie-frames and detected cells within the same plot. See corresponding block in swimmg.py. 
- Make sure calibration settings in "gp_swimming.py" are adjusted to match the specifications of the microscope and objective used.

## Swimming analysis in soft agar
To analyze swimming cells in soft-agar images need to be aquired using a confocal microscope. We used a Leica Microsystems microscope and the corresponding aquisition software for this purpose. The code to handle the files, cell detection, trajectory generation and trajectory analysis is provided in the subfolder "tracking_softagar". For this purpose we historically used the MATLAB based tracking package 'u-track' from the Danuser lab: 
- utrack 2.0 - Copyright (C) 2007 Danuser Lab. See https://github.com/DanuserLab/u-track.

But the analysis can also be run using the Python module trackpy and only Python based code. See the subfolder and the provided README for additional details.

