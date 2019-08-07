# swimming_analysis
Analysis of movies to extract bacterial swimming behavior and statistics

Code for the swimming analysis as done in our work Cremer et al. 

## Taking movies of swimming bacteria

Collect movies of swimming bacteria using phase contrast microscopy. Typical settings are a 10x objective, and image aquisition with a rate of 20FPS. Observation is done with a properly diluted culture of cells allowing the proper distinction of different cells while allowing also for a good swimming statistics (Typical OD600 value 0.05). Dilution in preheated medium, observation within a temperature controlled enviornment. More details are provided in our work.

Note: Please store movies - even if they are in different folders - with a unique name, for example by having the date as part of the file-name. Movies stored within one folder are recognized by the script and can be analyzed sequentially. 

## Analyze swimming behavior

The analysis is done in 3 consequetive steps:
- Detection of cells for each frame of the movie.
- Generation of trajectories. 
- Analysis of trajectories and output of swimming behavior (average swimming speed and additional statistics).
The details of trajectory generation are introduced in the Supplementary Text of our work. 

## Run script

Use swimming.py to run script. 

- Before running the first time, adjust the names of the folders raw movies are stored and where results should be stored, see beginning of the file. Same adjustments have to be made for module file gp_swimming.py as well. 
- Before running the first time, adjust the length calibration and the timing between two frames by adjusting gp_swimming.py. In particular adjust the 
- Consequetively run the different parts of the swimming.py script, as explained in the comments. Before running, adjust file names and detection and analysis parameters according to your file names and imageing settings.
- When running the script the first time, you can check for successfull detection by ploting movie-frames and detected cells within the same plot. See corresponding block in swimmg.py. 

