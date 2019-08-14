###################################################################################
#Swimming analysis in agar
#Github: https://github.com/jonascremer/swimming_analysis
###################################################################################
#
#This Python code analyzes swimming behavior of motile bacteria. The code was used
#to analyze the movement characteristics of cells in liquid conditions. We used it
# for our manuscript:
#J.Cremer, T.Honda, Y.Tang, J.Wong-Ng, M.Vergassola, T.Hwa. 
#Chemotaxis as a navigation strategy to thrive in nutrient-replete environments
#
#August 2019, Jonas Cremer together with the other coauthors.
#
#Additional information provided in "README.md"



###################################################################################
#import of needed packages
import os
###################################################################################
from pylab import *
import colormaps
cmaps=colormaps.cmaps
import gp_swimming #import main module of swimming analysis


#####   
#: Folders where raw movies are located
#####
maindirc=os.getcwd() #output of file
datafolderrawvideo='/Volumes/xydrive'
datafolderrawvideo2='/Volumes/xydrive2'
    

#######
#the following define subfolders - they are generated if they don't exist
#######
data_detection=os.path.join(maindirc,"data_detection")
rawgrowthfolder=os.path.join(maindirc,"data_raw")
data_trajectories=os.path.join(datafolderrawvideo2,"data_trajectories")
data_analysisoutput=os.path.join(maindirc,"data_swimminganalysis")
folder_detectiontest=os.path.join(maindirc,"detection_test")
plot_controloutput=os.path.join(data_analysisoutput,"controloutput")
data_movieoutput=os.path.join(plot_controloutput,"movie_output")
#generate folders if they don't excist
if not os.path.exists(data_detection):
        os.makedirs(data_detection)
else:
        pass

if not os.path.exists(data_trajectories):
        os.makedirs(data_trajectories)
else:
        pass

if not os.path.exists(data_analysisoutput):
        os.makedirs(data_analysisoutput)
else:
        pass
    
if not os.path.exists(plot_controloutput):
        os.makedirs(plot_controloutput)
else:
        pass
if not os.path.exists(data_movieoutput):
        os.makedirs(data_movieoutput)
else:
        pass
        
#%% ########
#go through folder and generate short movies to check if detection settings are working
#run to check microscope settings, once settings are set this step can be fixed.
#the script generates a short movie to see the movement of detected cells.
###########

foldernamedata="2017.12.21"#give folder name with movies, script will go through all the movies within that folder
strcriteria="" #to select subset of movies within folder, emty string if eveyr file in folder should be used

#parameter settings detection
min_value_selectionbackground=17000
selectionfraction=0.15 #coarse detection
selectionfraction2=0.05 #fine detection, false positives possible.
gp_swimming.run_detection(foldernamedata,framemaxin=50,plot_detectionanalysis=True,num_avsteps=40,short=True,selectionfraction=selectionfraction, selectionfraction2=selectionfraction2,blur_parameter=20,min_value_selectionbackground=min_value_selectionbackground,strcriteria=strcriteria)  #if empty string main folder is searched

#%% ######
#detect cells 
##########
#this can also be done in parallel using a computation cluster
foldernamedata="2017.12.21" #give folder name with movies, script will go through all the movies within that folder
#output is stored in one common folder, so it is iportant that every movie has a unique name
strcriteria="" #to select subset of movies within folder, emty string if eveyr file in folder should be used

#parameter settings detection
min_value_selectionbackground=17000
selectionfraction=0.15#0.07 #0.07, 0.15(PDMS)
selectionfraction2=0.05
fluorescencedata=False #set if data is phase contrast or fluorescence...    
framemax=2400
shortmode=False
gp_swimming.run_detection(foldernamedata,framemaxin=framemax,plot_detectionanalysis=False,num_avsteps=40,short=shortmode,selectionfraction=selectionfraction, selectionfraction2=selectionfraction2,blur_parameter=20,min_value_selectionbackground=min_value_selectionbackground,strcriteria=strcriteria)  #if empty string main folder is searched
   
#%% ######
#go through detected runs and generate trajectories (and statistics)    
#######

#define criteria to select for which movies trajectories should be generated
crit1="2017.02.19" #only movie names with such a string will be analyzed.
crit2="" #define  a second criteria (like crit2="OD1" or leave empty).
excplude1=""
repeat=False #if repat is true, trajectory analysis is done even if file is already there, if False, only analysis of trajectories is done (if statisticalanalysis=True)
onlyselectiondisplay=False #if true, don't run it show only which files...
statisticalanalysis=True
gp_swimming.analysis(strcriteria1=crit1,strcriteria2=crit2,statisticalanalysis=statisticalanalysis,minlength_analysis=10,onlyselectiondisplay=onlyselectiondisplay,repeat=repeat,excplude1=excplude1)

#%% ######
#analyze trajectories and plot results
#######
namec="output" #this is the name string used for the output file
rit1="2017.02.19" #only movie names with such a string will be analyzed.
crit2="" #define  a second criteria (like crit2="OD1" or leave empty).
excplude1=""
namelist=[gp_swimming.get_folderentries(data_trajectories, strcriteria1=crit1,strcriteria2=crit2,excplude1=excplude1,ending=".pad")]


labellist=namelist
growthratelist=[]

#read in ....
filenameout="distributioncomparison"+namec

###the following options for plotting the distribution
vmaxhis=50
vhisbins=200 #200 is good for data with swimming cells
vsmode="fraction" #if growthlist is defined then "fraction" or "growth" is used

gp_swimming.plot_together(namelist,labellist,filenameout=filenameout,vmaxhis=vmaxhis,vhisbins=vhisbins,plot_trajectories=False,vsmode=vsmode)

