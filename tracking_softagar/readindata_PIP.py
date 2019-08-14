#readindata_PIP.py
###################################################################################
#Swimming analysis in agar
#Github: https://github.com/jonascremer/swimming_analysis
###################################################################################
#
#This collection of Python and Matlab code for single-cell detection, cell-tracking
#and the extraction of movement statistics (diffusion and drift). This code was used
#to analyze the movement characteristics of cells in expanding poputlations (within
# soft-agar), as published in our manuscript:
#J.Cremer, T.Honda, Y.Tang, J.Wong-Ng, M.Vergassola, T.Hwa. 
#Chemotaxis as a navigation strategy to thrive in nutrient-replete environments
#
#August 2019, Ying Tang and Jonas Cremer together with the other coauthors.
#
#additional information provided in "README_swimminganalysisinagar.md"

###################################################################################
#import needed packages
###################################################################################
import numpy as np
import scipy.io as sio
import os
import time
import colormaps
cmaps=colormaps.cmaps
import analyze_tools  #functions to read in meta files

###################################################################################
#Settings
###################################################################################
#set names of files to be analyzed
filename = "C:\Python\Experiment20180112HE486_WT_gly_UniformMixedCells600" #example
SquenceName='Sequence_002_Job 1_017' #prefix name of movies/images how it was stored by the microscope software.
movieoutput=False
listoutput=False 
thresholdv=15 #settings treshold analysis
opt_movie_fps_image=5 

#set names of folders, generate folders if needed
[rawdatafolder,digestdatafolder,movieoutputfolder,clusterfolder]=analyze_tools.get_folders(filename)

###################################################################################
#Start main program
###################################################################################
#go through image-data and determine number of tile-scans, z-steps, and t-steps, and color-channels
[maxT,maxS,maxC,maxM,jobnamelist,nameexpc]=analyze_tools.get_image_counts_movie(rawdatafolder)

###go through meta data first to determine min range...
[tstart,xmin,ymin,zmin]=analyze_tools.get_minvalues_movie(rawdatafolder,jobnamelist,nameexpc,maxS) #go through all timesteps and get smalles values

#S1=0f only one x, one z
T=0
S1=0
S2=S1+1 #0-300
T1=0#For this data, T is different z scan
T2=T1+5


##### select what to plot
generate_cluster=True# generate cluster when do_tracking is True
do_tracking=True
CalculateDiffusion=True
ratio_Cluster=True
ratio_fromtraj=False
RunTumble=False
Detection_matlab=False#start when first time, generate mat file for tracking: need to change dection_function mat file for different scripts
PlotMcherry=False
#If don't want to repeat drawing images of traj
NoRepeatDrawingImages=False
#save images to analyze tracking...
##################
 
if generate_cluster==True:


    for Sc in range(S1,S2): #go through different jobs (usually sequences)

       #for every S and Tcombination, one line of output...

      for Cc in range(0,maxC):#go through every color        

        counterimg=0

        for Tc in range(T1,T2):#maxT go through every time step 
       		[xpos,ypos,zpos,tlist,imgext]=analyze_tools.readin_moviemeta(rawdatafolder,jobnamelist[Sc],nameexpc,xmin=xmin,ymin=ymin,zmin=zmin,tstart=tstart)
            tcurrentfirst=tlist[maxM*Tc+0]
			clusterlist=[]#stores all cluster details
			avlist=[]#stores average intensity, number of clusters...
			MaxNumberCluster=0 #store largest number of cluster for different frames
            clusterarray = np.zeros((maxM,3), dtype=np.object)
            if Detection_matlab==True:
                analyze_tools.detection_matlab(maxM,Sc,Cc,Tc,maxT,maxC,maxM,jobnamelist,nameexpc,rawdatafolder,clusterfolder,thresholdv=thresholdv,Mc=1)
                time.sleep(120)
            #detect clusters and store them in a file...            
            NoCellDetected=analyze_tools.generate_clusterfile(maxM,Sc,Cc,Tc,maxT,maxC,maxM,jobnamelist,nameexpc,rawdatafolder,clusterfolder,thresholdv=thresholdv)
            print NoCellDetected
            #get trajectories buy running tracking algorithm in matlab
            if NoCellDetected==[0]:
                if do_tracking==True:
                    analyze_tools.do_tracking(Sc,Cc,Tc,clusterfolder,trackingradius=30)
                    time.sleep(120)
if CalculateDiffusion==True:
    for Sc in range(S1,S2): #go through different jobs (usually sequences)
      #for every S and Tcombination, one line of output...
      for Cc in range(0,maxC):#go through every color        
        counterimg=0
        for Tc in range(T1,T2):
            analyze_tools.Calculate_Diffusion(Sc,Cc,Tc,clusterfolder)
            time.sleep(90)
if RunTumble==True:
    for Sc in range(S1,S2): #go through different jobs (usually sequences)
       #for every S and Tcombination, one line of output...
      for Cc in range(0,maxC):#go through every color        
        counterimg=0
        for Tc in range(T1,T2):
            analyze_tools.RunandTumble(Sc,Cc,Tc,clusterfolder)
  
##########################################################
if ratio_Cluster==True:
 Ratio=np.zeros((S2-S1,T2-T1))
 RatioNon=np.zeros((S2-S1,T2-T1))
 Number=np.zeros((S2-S1,T2-T1))
 NumberNon=np.zeros((S2-S1,T2-T1))
 nonMotileDist=5 #set distance on judging nonMotile cells
 
 NumberCell=np.zeros((S2-S1,T2-T1))
 #RatioY=np.zeros((S2-S1,T2-T1))    
 #RatioNonY=np.zeros((S2-S1,T2-T1))
 for Sc in range(S1,S2): #go through different jobs (usually sequences)
  for Cc in range(0,maxC):#go through every color        
   for Tc in range(T1,T2):
        [NumberCell[Sc-S1,Tc-T1]]=analyze_tools.NumberCell_fromCluster(maxM,Sc,Cc,Tc,maxT,maxC,maxM,jobnamelist,nameexpc,rawdatafolder,clusterfolder,thresholdv=thresholdv,nonMotileDist=5)        
        if RatioNon[Sc-S1,Tc-T1]==-1:
            RatioNon[Sc-S1,Tc-T1]=0            
            Ratio[Sc-S1,Tc-T1]=0
            NumberNon[Sc-S1,Tc-T1]=0
            Number[Sc-S1,Tc-T1]=0
        else:
            Ratio[Sc-S1,Tc-T1]=1-RatioNon[Sc-S1,Tc-T1]
             
  filenameml=os.path.join(clusterfolder,"NumberCell"+".mat")
  sio.savemat(filenameml,mdict={'NumberCell': NumberCell})      #print "....."
##########################################################


            
            
    
