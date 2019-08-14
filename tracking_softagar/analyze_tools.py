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
#Import
###################################################################################
from datetime import datetime
import numpy as np
from scipy.ndimage import measurements
import scipy
import subprocess
import os
from PIL import Image
import scipy.io as sio
from pylab import *
import colormaps
cmaps=colormaps.cmaps



###################################################################################
#Preparation and Settings
###################################################################################
#
#to setup u-track, have u-track package in a subfolder. Copy all script files provided matlabfolder folder into main folder with scripts.
#when running first time, adjust rawdatafolder and mainfolder below such that data can be found.
#
#Before each run set name of Imagesequence to be analyzed
SquenceName='Sequence_002_Job 1_017'
tdelta_byhand=16

#generate needed folders...
def get_folders(filename):
    #folder for data
    rawdatafolder=os.path.join("~/movieanalysis/rawdata",filename)
    mainfolder=os.path.join("~/movieanalysis/",filename)
    try:
        os.makedirs(mainfolder)
    except:
        pass
    digestdatafolder=os.path.join(mainfolder,"digestdata")
    try:
        os.makedirs(digestdatafolder)
    except:
	pass
    movieoutputfolder=os.path.join(mainfolder,"movieoutput")
    try:
        os.makedirs(movieoutputfolder)
    except:
        pass
    clusterfolder=os.path.join(mainfolder,"clusterdata")
    try:
        os.makedirs(clusterfolder)
    except:
	pass
    return [rawdatafolder,digestdatafolder,movieoutputfolder,clusterfolder]
 

   
def generate_movie(clustermoviefolder,imgending="jpg",fpts=10,name_movieoutis=""):
            if name_movieoutis=="":
                name_movieoutis=clustermoviefolder+".mp4"
            
                process_name ="ffmpeg -y -r "+ str(fpts)+" -pattern_type glob -i '"+os.path.join(clustermoviefolder, "*."+imgending+"'")+" -c:v libx264 -pix_fmt yuv420p "+name_movieoutis
                               
                #process_name ="ffmpeg -framerate 2 -i"+".\movie_cluster_S1\0001%2d.jpg"+"-c:v libx264 -r 30 -pix_fmt yuv420p"+"out1.mp4"
                #outp=subprocess.check_output(process_name, shell=True)                    
                
                


#get full filename of images (moviemode)
def get_filenamecurimg_movie(rawdatafolder,Tc,maxT,Cc,maxC,Sc,jobnamelist,Mc,maxM,nameexpc):
            
            if maxT>100:
                if Tc==0:
                    Tcstr="000"
                elif Tc<10:
                    Tcstr="00"+str(Tc)
                elif Tc<100:
                    Tcstr="0"+str(Tc)
                else:
                    Tcstr=str(Tc)
            elif maxT>10:
                if Tc==0:
                    Tcstr="00"
                elif Tc<10:
                    Tcstr="0"+str(Tc)
                else:
                    Tcstr=str(Tc)        
            else:
                Tcstr=str(Tc)

            if Mc==0:
                Mcstr="00"
            elif Mc<10:
                Mcstr="0"+str(Mc)
            else:
                Mcstr=str(Mc)
            if maxM>100:
                if Mc==0:
                    Mcstr="000"
                elif Mc<10:
                    Mcstr="00"+str(Mc)
                elif Mc<100:
                    Mcstr="0"+str(Mc)
                else:
                    Mcstr=str(Mc)
            if maxM>1000:
                if Mc==0:
                    Mcstr="0000"
                elif Mc<10:
                    Mcstr="000"+str(Mc)
                elif Mc<100:
                    Mcstr="00"+str(Mc)
                elif Mc<1000:
                    Mcstr="0"+str(Mc)
                else:
                    Mcstr=str(Mc)
                
                
            #print maxM
                
            #for formath like blah 18_448_l49_t99_ch00.tif
             
            if rawdatafolder=='C:\Python\Experiment20161007':
                return jobnamelist[Sc]+"_l"+Tcstr+"_t"+Mcstr+".tif"
            else:
                #print 'name'+nameexpc
                #return jobnamelist[Sc]+"_z"+Tcstr+"_t"+Mcstr+"_ch0"+str(Cc)+".tif"
                return "l"+Tcstr+"_t"+Mcstr+"_ch0"+str(Cc)+".tif"

def readin_moviemeta(datafolder,jobname,nameexpc,xmin=0,ymin=0,zmin=0,tstart=0):
    ##############
    #read in 
    ##############


    
    ####
    ##read metadata
    #####
    #print "read in metatxt"
    position_xRI=[] #for tilescan
    position_yRI=[] #for tilescan
    timestamps=[]

    metadictRI={}
    datafolder=os.path.join(datafolder,'MetaData')
    metadataname=os.path.join(datafolder,SquenceName+".xml")
     
    with open(metadataname) as f:
        for lline in f:
#            
             if "StagePosX=" in lline:
                 xpos=float(lline.split('StagePosX="')[1].split('"')[0])
             if "StagePosY=" in lline:
                 ypos=float(lline.split('StagePosY="')[1].split('"')[0])
             if "ZPosition=" in lline:
                 zpos=float(lline.split('ZPosition="')[1].split('"')[0])
             if 'NumberOfElements="512" Origin="' in lline:
                 imgext=float(lline.split('NumberOfElements="512" Origin="')[1].split('" Length="')[1].split('" Unit=')[0])*1000000.#in mum
                
#            8.673617e-020" Length="1.923575e-004" Unit="m" BitInc="0" BytesInc="1"></DimensionDescription><DimensionDescription DimID="2" NumberOfElements="512" Origin="8.673617e-020" Length="1.923575e-004" Unit="m"
#                
             if "NumberOfTimeStamps=" in lline:
                 
                 timestamps=lline.split("NumberOfTimeStamps")[-1].split(">")[1].split(" ")[:-2]
                 #print timestamps
                 for il in range(0,len(timestamps)):
                     timestamps[il]=(int(timestamps[il],16)/10000.)/1000. #geturns time of image taken in seconds.
                     #print timestamps[il]
    xpos=0
    ypos=0
    zpos=0
    xpos=xpos-xmin
    ypos=ypos-ymin
    zpos=zpos-zmin
    timestamps.append(1000000.)
    timestamps=np.array(timestamps)
    timestamps=timestamps-tstart
    return [xpos,ypos,zpos,timestamps,imgext]
                
def get_image_movie(rawdatafolder,filenamcec,rgbformat=False):
                    filenamec=os.path.join(rawdatafolder,filenamcec)
                    try:
                        imgc=Image.open(filenamec)
                        imarr=np.array(imgc,dtype=np.float)
                        
                    except:
                        print "unable to open..."
                        print filenamec
                        #analyze images...
                    print "after open "
                    print datetime.datetime.now().time()
                    
                    if rgbformat:
                
                        imshowarr=np.zeros([imarr.shape[0],imarr.shape[1],3])
                        imshowarr[:,:,1]=imarr[:,:]
                        return imshowarr
                    else:
                        return imarr
                        
def get_cluster(Mc,Sc,Cc,Tc,maxT,maxC,maxM,jobnamelist,nameexpc,rawdatafolder,thresholdv=15):
                filenamcec=get_filenamecurimg_movie(rawdatafolder,Tc,maxT,Cc,maxC,Sc,jobnamelist,Mc,maxM,nameexpc)
                filenamcec=SquenceName+'_'+filenamcec
                filenamec=os.path.join(rawdatafolder,filenamcec)

                try:
                    imgc=Image.open(filenamec)
                    imarr=np.array(imgc,dtype=np.float)
                    #print imarr
                    imarr=imarr[:,:,1]
                    #print 'iiiiiiiiiiiiiiiiiiiiiiiiii'
                   # print imarr
                    
                except:
                    print "unable to open..."
                    print filenamec
                
                imarrts=np.copy(imarr)
                imarrts[imarrts < thresholdv] = 0
                imarrts[imarrts >= thresholdv] = 1
                imarr[:,:]=scipy.ndimage.filters.median_filter(imarr,3)[:,:]
                
                #get clusters on image
                sbindiagnoal=scipy.ndimage.morphology.generate_binary_structure(2,2)  ##This only generates a 3elements' 2d matrix!!
                labeled_array, numberclusters_curimg = measurements.label(imarr, structure=sbindiagnoal) #We can use imarr or imarrts!!
                centermass=np.zeros([numberclusters_curimg,6]) #x, y, mass, total size, extension, name
                Imgtreshhold=""
                for ilc in range(1,numberclusters_curimg):
                    #get center of mass for different clusters...
                    curpos=measurements.center_of_mass(imarr, labels=labeled_array,index=ilc)
                    centermass[ilc,0]=curpos[0]
                    centermass[ilc,1]=curpos[1]
                    centermass[ilc,2]=measurements.sum(imarr, labels=labeled_array, index=ilc)
                    centermass[ilc,3]=measurements.sum(imarrts, labels=labeled_array, index=ilc)
                    centermass[ilc,4]=np.sqrt(centermass[ilc,3])
                    
                    
                    
                    if np.nanmax(centermass[ilc,3])>40:
                        Imgtreshhold="*"
                
                
                disttreshold=10
                #go through clusters, give names and merge if too close...
                for ilc in range(0,numberclusters_curimg):
                    for ilc2 in range(centermass.shape[0]-1,ilc,-1):
                        distc=np.sqrt(np.power(centermass[ilc2,0]-centermass[ilc,0],2.)+np.power(centermass[ilc2,1]-centermass[ilc,1],2.))
                        if distc<disttreshold:
                            centermass[ilc,0]=(centermass[ilc,0]*centermass[ilc,2]+centermass[ilc2,0]*centermass[ilc2,2])/(centermass[ilc,2]+centermass[ilc2,2])
                            centermass[ilc,1]=(centermass[ilc,1]+centermass[ilc2,1])/2.
                            centermass[ilc,2]=(centermass[ilc,2]+centermass[ilc2,2])
                            centermass[ilc,3]=(centermass[ilc,3]+centermass[ilc2,3])
                            centermass[ilc,4]=np.sqrt(centermass[ilc,3])
                            
                            centermass = np.delete(centermass,(ilc2), axis=0)
                return centermass

def getintstr(Mc):                
                msspace=""
                if Mc<10:
                    msspace="00"
                elif Mc<100:
                    msspace="0"
                elif Mc<1000:
                    msspace=""
                return msspace+str(Mc)
                

def load_trackingdata(Sc,Cc,Tc,clusterfolder):
  filenameclusterstr="cluster_S"+str(Sc)+"_C"+str(Cc)+"_T"+str(Tc)
  xfilep=os.path.join(clusterfolder,filenameclusterstr+'_x.mat')
  yfilep=os.path.join(clusterfolder,filenameclusterstr+'_y.mat')
  #load trackinginfo....
  if os.path.isfile(xfilep)==False or os.path.isfile(yfilep)==False:
      print "file not found"
      return [False,0,0]
  else:
       
       TrajX=np.array(sio.loadmat(xfilep)['x'])
       TrajY=np.array(sio.loadmat(yfilep)['y']) 
       return [True,TrajX,TrajY]

def ratioOfmobile(Sc,Cc,Tc,clusterfolder):
    filenameclusterstr="cluster_S"+str(Sc)+"_C"+str(Cc)+"_T"+str(Tc)
    RfilepX=os.path.join(clusterfolder,filenameclusterstr+'RatioOfMobileX.mat')
    RfilepY=os.path.join(clusterfolder,filenameclusterstr+'RatioOfMobileY.mat')
    RfilepR=os.path.join(clusterfolder,filenameclusterstr+'RatioOfMobileR.mat')
    NfilepR=os.path.join(clusterfolder,filenameclusterstr+'NumberOfMobileR.mat')
    NNfilepR=os.path.join(clusterfolder,filenameclusterstr+'NumberOfNonMobileR.mat')
    if os.path.isfile(RfilepX)==False or os.path.isfile(RfilepY)==False or os.path.isfile(RfilepR)==False or os.path.isfile(NfilepR)==False:
      print "RatioOfMobile file not found"
      return [0,0,0,0,0,0,0,0]
      
    else:
      RatioX=np.array(sio.loadmat(RfilepX)['RatioMobileX'])
      RatioY=np.array(sio.loadmat(RfilepY)['RatioMobileY'])
      RatioR=np.array(sio.loadmat(RfilepR)['RatioMobileR'])
      RationonX=1-RatioX
      RationonY=1-RatioY
      RationonR=1-RatioR
      NumberR=np.array(sio.loadmat(NfilepR)['NumberMobileR'])
      NumberNonR=np.array(sio.loadmat(NNfilepR)['NumberNonMobileR'])
      return [RatioX,RationonX,RatioY,RationonY,RatioR,RationonR,NumberR,NumberNonR]




      
def Calculate_Diffusion(Sc,Cc,Tc,clusterfolder):
    filenameclusterstr="cluster_S"+str(Sc)+"_C"+str(Cc)+"_T"+str(Tc)
    xfilep=os.path.join(clusterfolder,filenameclusterstr+'_x.mat')
    yfilep=os.path.join(clusterfolder,filenameclusterstr+'_y.mat')
    if os.path.isfile(xfilep) and  os.path.isfile(yfilep):
        diffusionScript='DriftDiffusion' #DriftDiffusionBallFilter_20161112, DriftDiffusionGaussianFilter_20161112 
        process_name=matlabpath +" -nodisplay -nosplash -nodesktop -r \""+diffusionScript+"('"+xfilep+"','"+yfilep+"','"+clusterfolder+"','"+filenameclusterstr+"',"+str(Sc)+","+str(Cc)+","+str(Tc)+"); exit\""
        process_name="cd "+matlabtrackingscriptfolder+"; "+process_name
        outp=subprocess.call(process_name, shell=True)
        
        
def RunandTumble(Sc,Cc,Tc,clusterfolder):
    filenameclusterstr="cluster_S"+str(Sc)+"_C"+str(Cc)+"_T"+str(Tc)
    xfilep=os.path.join(clusterfolder,filenameclusterstr+'_x.mat')
    yfilep=os.path.join(clusterfolder,filenameclusterstr+'_y.mat')
    if os.path.isfile(xfilep) and  os.path.isfile(yfilep):
        diffusionScript='XYmatricesToRunTimes' #DriftDiffusionBallFilter_20161112, DriftDiffusionGaussianFilter_20161112 
        process_name=matlabpath +" -nodisplay -nosplash -nodesktop -r \""+diffusionScript+"('"+xfilep+"','"+yfilep+"','"+clusterfolder+"','"+filenameclusterstr+"',"+str(Sc)+","+str(Cc)+","+str(Tc)+"); exit\""
        process_name="cd "+matlabtrackingscriptfolder+"; "+process_name
        outp=subprocess.call(process_name, shell=True)        
        
        
        
        
def do_tracking(Sc,Cc,Tc,clusterfolder,trackingradius=30,reanalyze=False):
    filenameclusterstr="cluster_S"+str(Sc)+"_C"+str(Cc)+"_T"+str(Tc)
    xfilep=os.path.join(clusterfolder,filenameclusterstr+'_x.mat')
    yfilep=os.path.join(clusterfolder,filenameclusterstr+'_y.mat')
    if reanalyze==False:
        if os.path.isfile(xfilep) and  os.path.isfile(yfilep):
            trackingdetected=False
            print "tracking performed and trajectories detected: "+xfilep+"...."
        else:
            
            trackingdetected=True
        
        
    
    if reanalyze==True or trackingdetected:
       
        filenameclusterml=os.path.join(clusterfolder,filenameclusterstr+".mat")
        #check if file exists    
        if os.path.isfile(filenameclusterml)==False:
            print "Error: clusterfile not found: "+filenameclusterml
            print "File needs to be generated first before doing tracking."
            error
        
        radin1=trackingradius
        radin2=trackingradius
        process_name=matlabpath +" -nodisplay -nosplash -nodesktop -r \""+matlabtrackingscript+"('"+filenameclusterml+"','"+clusterfolder+"','"+filenameclusterstr+"',"+str(radin1)+","+str(radin2)+"); exit\""
              
        process_name="cd "+matlabtrackingscriptfolder+"; "+process_name
        
        outp=subprocess.call(process_name, shell=True)
        

    
def load_clusterfile(Sc,Cc,Tc,clusterfolder):
    filenameclusterstr="cluster_S"+str(Sc)+"_C"+str(Cc)+"_T"+str(Tc)
    filenameclusterpy=os.path.join(clusterfolder,filenameclusterstr+".npy")
    
    try:
        return np.load(filenameclusterpy)
    except:
        print "Error: clusterfile not found: "+filenameclusterpy
        print "File needs to be generated first."
        error
   
def generate_clusterfile(untilM,Sc,Cc,Tc,maxT,maxC,maxM,jobnamelist,nameexpc,rawdatafolder,clusterfolder,thresholdv=15,reanalyze=False):
    #check if file is there               
    filenameclusterstr="cluster_S"+str(Sc)+"_C"+str(Cc)+"_T"+str(Tc)
    filenameclusterml=os.path.join(clusterfolder,filenameclusterstr+".mat")
    filenameclusterpy=os.path.join(clusterfolder,filenameclusterstr+".npy")
    
    
    
    if reanalyze==False:
        if os.path.isfile(filenameclusterml) and  os.path.isfile(filenameclusterpy):
            clusterfilenotdetected=False
            print "clusterfile "+filenameclusterstr+" detected...."
        else:
            
            clusterfilenotdetected=True
        
        
    
    if reanalyze==True or clusterfilenotdetected:
            print "generating clusterfile "+filenameclusterstr+"...."
            
            clusterarray = np.zeros((maxM,3), dtype=np.object)
            NoCellDetected=1
            for Mc in range(0,maxM): #go through different frames within one movie
                
                centermass=get_cluster(Mc,Sc,Cc,Tc,maxT,maxC,maxM,jobnamelist,nameexpc,rawdatafolder,thresholdv=thresholdv)
                
                
                 
                if centermass.shape[0]>0:
                    NoCellDetected=0 
                    
                
                clusterarray[Mc,0] = np.transpose(np.vstack((centermass[0:centermass.shape[0],0]+0.0001,np.zeros(centermass.shape[0])))) #intensity or set uniformly as 255
                clusterarray[Mc,1] = np.transpose(np.vstack((centermass[0:centermass.shape[0],1]+0.0001,np.zeros(centermass.shape[0])))) #y
                clusterarray[Mc,2] = np.transpose(np.vstack((centermass[0:centermass.shape[0],2],np.zeros(centermass.shape[0])))) #x
            sio.savemat(filenameclusterml,mdict={'movieInfoCell': clusterarray})
            
            np.save(filenameclusterpy,clusterarray)
            #now we can start matlab script here....
            return[NoCellDetected]
            
            
            
def RatioOfNonMobile_fromCluster(untilM,Sc,Cc,Tc,maxT,maxC,maxM,jobnamelist,nameexpc,rawdatafolder,clusterfolder,thresholdv=15,reanalyze=False,nonMotileDist=5):
    clusterlist=[]
    avlist=[]
            
    filenameclusterstr="cluster_S"+str(Sc)+"_C"+str(Cc)+"_T"+str(Tc)
    filenameclustermDetected=os.path.join(clusterfolder,filenameclusterstr+"_x.mat")
           
    if os.path.isfile(filenameclustermDetected):
                    
            for Mc in range(0,maxM): #go through different frames within one movie
                centermass=get_cluster(Mc,Sc,Cc,Tc,maxT,maxC,maxM,jobnamelist,nameexpc,rawdatafolder,thresholdv=thresholdv)
            
                
                #give clusters names                
                if Mc==0:
                    for ilc in range(0,centermass.shape[0]):
                        centermass[ilc,5]=ilc+1
                else:
                    disttresholdtracking=30
                    for ilc in range(0,centermass.shape[0]):
                        alreadyfound=False
                        centermass[ilc,5]=0
                        for ilc2 in range(0,centermass_old.shape[0]):
                            distc=np.sqrt(np.power(centermass_old[ilc2,0]-centermass[ilc,0],2.)+np.power(centermass_old[ilc2,1]-centermass[ilc,1],2.))
                            if distc<disttresholdtracking:
                                if alreadyfound==False:
                                    centermass[ilc,5]=centermass_old[ilc2,5]   
                                    alreadyfound=True
                                    distbefore=distc
                                elif distc<distbefore:
                                    centermass[ilc,5]=centermass_old[ilc2,5]
                                    distbefore=distc
                      
                                
                clusterlist.append(centermass)
                centermass_old=centermass
                if Mc==0:
                    centermass_start=centermass
                
                if centermass.shape[0]==0:
                    nonCellDetected=1
                
                #go through cluster and determine if they were already observed at movie frame 0
            num_const=0
            for ilc in range(0,centermass.shape[0]):
                    namc=centermass[ilc,5]
                    if namc>0:
                        for ilc2 in range(0,centermass_start.shape[0]):
                            namc2=centermass_start[ilc2,5]
                            if namc2==namc:
                                distc=np.sqrt(np.power(centermass_start[ilc2,0]-centermass[ilc,0],2.)+np.power(centermass_start[ilc2,1]-centermass[ilc,1],2.))
                                if distc<nonMotileDist:
                                    num_const=num_const+1
            #go through every movie-frame and determine number of cells
            #use data: 
            #1: avlist[]#num clusters, and ave. intensit
            #2: #x, y, mass, total size, extension, name
            numcc=np.zeros([maxM])
            numsize=np.zeros([maxM])#max size for each image


            for Mc in range(0,maxM): 
                        
                        numcc[Mc]=clusterlist[Mc].shape[0]
                        
                        try:
                            numsize[Mc]=np.nanmax(clusterlist[Mc][:,2])
                        except:
                            numsize[Mc]=0
                        #avlist[Mc][0]
                        
            
            avnumclu=np.mean(numcc)
                       
            if  avnumclu>0:
                    fracc=num_const/avnumclu
            else:
                    fracc=0
               
            return[fracc,num_const,centermass.shape[0]-num_const]
    else:
         return [-1,-1,-1]
            
            
def NumberCell_fromCluster(untilM,Sc,Cc,Tc,maxT,maxC,maxM,jobnamelist,nameexpc,rawdatafolder,clusterfolder,thresholdv=15,reanalyze=False,nonMotileDist=5):

    filenameclusterstr="cluster_S"+str(Sc)+"_C"+str(Cc)+"_T"+str(Tc)
    filenameclustermDetected=os.path.join(clusterfolder,filenameclusterstr+".npy")
           
    if os.path.isfile(filenameclustermDetected):
            print "NumberCell_fromCluster "+filenameclusterstr+" is analyzing"
            #RatioMobileAnalysis=1
            
            NumberCell=0        
            for Mc in range(0,maxM): #go through different frames within one movie
                centermass=get_cluster(Mc,Sc,Cc,Tc,maxT,maxC,maxM,jobnamelist,nameexpc,rawdatafolder,thresholdv=thresholdv)

               
                if Mc==1:
                  NumberCell=centermass.shape[0]
                  print NumberCell
                  return [NumberCell]
   
    else:
        return [0]                       
            
            
            

def detection_matlab(untilM,Sc,Cc,Tc,maxT,maxC,maxM,jobnamelist,nameexpc,rawdatafolder,clusterfolder,thresholdv=15,reanalyze=False,Mc=1):
            filenameclusterstr="cluster_S"+str(Sc)+"_C"+str(Cc)+"_T"+str(Tc)
            filenameclusterml=os.path.join(clusterfolder,filenameclusterstr+".mat")  
            filenamcec=get_filenamecurimg_movie(rawdatafolder,Tc,maxT,Cc,maxC,Sc,jobnamelist,Mc,maxM,nameexpc)
            filenamcec=SquenceName+'_'+filenamcec
            filenamcec=filenamcec[0:25] #need to adjust for different data
            print filenamcec
            #ddd
            matlabdetectionscript="detection_function"    
            process_name=matlabpath +" -nodisplay -nosplash -nodesktop -r \""+matlabdetectionscript+"('"+filenameclusterml+"','"+filenamcec+"','"+rawdatafolder+"'); exit\""
              
           # process_name=matlabpath +" -nodisplay -nosplash -nodesktop -r \""+diffusionScript+"('"+xfilep+"','"+yfilep+"','"+clusterfolder+"','"+filenameclusterstr+"',"+str(Sc)+","+str(Cc)+","+str(Tc)+"); exit\""
            if machine=="jonas":
                process_name="cd "+matlabtrackingscriptfolder+"; "+process_name
            else:
                process_name="cd "+matlabtrackingscriptfolder+"& "+process_name
                
            outp=subprocess.call(process_name, shell=True)            
            
            
            
def get_minvalues_movie(rawdatafolder,jobnamelist,nameexpc,scmax):
    tstart=10000000000000000
    xmin=1000
    ymin=1000
    zmin=10000
    for Sc in range(0,scmax): #go through different jobs (usually sequences)
         [xpos,ypos,zpos,tlist,imgext]=readin_moviemeta(rawdatafolder,jobnamelist[Sc],nameexpc)
         if xpos<xmin:
             xmin=xpos
         if ypos<ymin:
             ymin=ypos
         if zpos<zmin:
             zmin=zpos
         if min(tlist)<tstart:
             tstart=min(tlist)
    return [tstart,xmin,ymin,zmin]

def get_image_counts_movie(rawdatafolder):
    #name of files
    #Sequence_Job 1_439_l00_t07_ch00.tif
    #Sequence_Job 18_448_l49_t99_ch00.tif
    #11st is sequence...take after blank
    #l: is run (so t)
    #t: is rpeat
    #ch: is channel
    listfiles=os.listdir(rawdatafolder)
    
    maxS=0 #tile scans
    maxT=0 #time scans
    maxC=0 #color channels
    maxZ=0 #z-scans
    maxM=0 #t counter for movies
    jobnamelist=[]
    nameexpclist=[]
    
    numimg=0
    for i in listfiles:
      if len(i.split("."))>1 and i.split(".")[-1]=="tif":
          numimg=numimg+1
          namestr=i.split(".t")[-2].split("_")
          
          
          if rawdatafolder=='C:\Python\Experiment20161007':
              #for formath like blah Sequence_S1_066_I01_t101.tif
              Cc=0 #color channel=0
              Mc=int(namestr[-1][1:]) #time counter for movies...e.g. 0-99 for 100 frame movie
              Tc=int(namestr[-2][1:]) #l is time counter of pattern loop
              jobcname=namestr[-5][:]+"_"+namestr[-4][:]
              jobccname=namestr[0][:]
          else:
          #for formath like blah 18_448_l49_t99_ch00.tif
              Cc=int(namestr[-1][2:]) #color channel
              Mc=int(namestr[-2][1:]) #time counter for movies...e.g. 0-99 for 100 frame movie
              Tc=int(namestr[-3][1:])#int(namestr[-3][12:]) now is z
              
             # print namestr[-3]
              jobcname=namestr[-4][:]
              jobccname=namestr[0][:]
              #print jobccname
          if jobcname in jobnamelist:
              pass
          else:
              jobnamelist.append(jobcname)
          if jobccname in nameexpclist:
              pass
          else:
              nameexpclist.append(jobccname)
          if maxT<Tc:
              maxT=Tc
          if maxC<Cc:
              maxC=Cc
          if maxM<Mc:
              maxM=Mc
          
    #
    #check length
    maxC=maxC+1
    #sort list
    
    
    jobnint=np.zeros([len(jobnamelist)])
    for il in range(0,len(jobnamelist)):
        if rawdatafolder=='C:\Python\Experiment20161007':
            jobnint[il]=int(jobnamelist[il].split("_")[-1].split("S")[-1])
            
        else:
            
            jobnint[il]=0
    sortind=np.argsort(jobnint)
    jobnamelistnew=[]
    for il in range(0,len(jobnamelist)):
        if rawdatafolder=='C:\Python\Experiment20161007':
            jobnamelistnew.append(jobnamelist[sortind[il]]+'_0'+str(66+il))
        else:
            jobnamelistnew.append(jobnamelist[sortind[il]])
        
    
    jobnamelist=jobnamelistnew    
    maxS=len(jobnamelist)
    maxT=maxT+1
    maxM=maxM+1
    
    
    if len(nameexpclist)>1:
        print "error: more than one sequence name found"
        error
    else:
        nameexpc=nameexpclist[0]
    print numimg
    if (numimg) != maxT*maxS*maxM*maxC:
             print "Wrong length"
             print maxT
             print maxS
             print maxM
             print maxC
                
    else:
    	print "Check successfull: correct number of images in folder"
    
    print "Read in of file with dimenions: M="+str(maxM)+" T="+str(maxT)+" S="+str(maxS)+" C="+str(maxC)
    
    return[maxT,maxS,maxC,maxM,jobnamelist,nameexpc]

