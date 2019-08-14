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

#import needed packages
import trackpy
import cv2
import numpy as np
import json
import subprocess
from scipy import ndimage
import os, sys
import scipy.io
import pandas as pd
import scipy
import csv

import matplotlib
import matplotlib.pyplot as plt
import colormaps
cmaps=colormaps.cmaps



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
data_trajectories=os.path.join(datafolderrawvideo2,"data_trajectories")
data_analysisoutput=os.path.join(maindirc,"data_swimminganalysis")
folder_detectiontest=os.path.join(maindirc,"detection_test")
plot_controloutput=os.path.join(data_analysisoutput,"controloutput")
data_movieoutput=os.path.join(plot_controloutput,"movie_output")
trjstat_output=os.path.join( plot_controloutput,"trjstat_output")


#make sure folders are generated if they don't exist...
try:
        
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
    if not os.path.exists(trjstat_output):
            os.makedirs(trjstat_output)
    else:
            pass  
except:
    print "CHECK PATH HARDDRIVE"


def tresholdanalysis(substracted,mintreshold=0):
    #go through different tresholds....take treshold with less than 1000 cells?
    tresholdlist=range(mintreshold,255)
    particlenum_treshold=np.zeros([255,2])
    particlenum_treshold[:]=np.nan
    
    tccount=-1
    determineconstanttreshold=True
    for treshold in tresholdlist:
        tccount=tccount+1
        #print "treshold test"
        #print treshold
        #image_to_be_labeled = ((difference_image > threshold) * 255).astype('uint8')  # not sure if it is necessary
        [haeh,image_to_be_labeled]=cv2.threshold(substracted,treshold,255,cv2.THRESH_BINARY)
        #print np.nanmin(image_to_be_labeled)
        #print np.nanmax(image_to_be_labeled)
        labelarray, particle_count = scipy.ndimage.measurements.label(image_to_be_labeled)
        #print particle_count
        particlenum_treshold[tccount,0]=treshold
        particlenum_treshold[tccount,1]=particle_count
        if particle_count==1:
            break
        
        tresholdconstant=treshold
        if tccount>1 and determineconstanttreshold:
            
            if particlenum_treshold[tccount,1]<500.:
                #print particlenum_treshold[tccount,1]
                #print particlenum_treshold[tccount-1,1]
                tresholdconstant=treshold
                determineconstanttreshold=False
                break
        
    return [particlenum_treshold,tresholdconstant]
          
    
def remove_doubles(particlelist,particlelist_weight,distancemegere=5):
    #go through particle list and merge close particles....
    
    particlelist_new=[[particlelist[0][0],particlelist[0][1],particlelist_weight[0]]]
    for ilc in range(1,len(particlelist)):
        
        xc=particlelist[ilc][0]
        yc=particlelist[ilc][1]
        curweight=particlelist_weight[ilc]
        
        nobreak=True
        indxsame=[]  #detected same particles...
        indxweight=[]
        for ilc2 in range(0,len(particlelist_new)):
            xc2=particlelist_new[ilc2][0]
            yc2=particlelist_new[ilc2][1]
            dist=np.sqrt((xc2-xc)*(xc2-xc)+(yc2-yc)*(yc2-yc))
            if dist<distancemegere:
                nobreak=False
                indxsame.append(ilc2)
                indxweight.append(particlelist_new[ilc2][2])
                
                
        
        
        if nobreak:
            particlelist_new.append([xc,yc,curweight])
        else:
            #merge all the indexes to one...
            indxheightest=indxweight.index(max(indxweight))
            
            #remove all lower weight entries
            for ilw in range(len(indxsame)-1,-1,-1):
                if ilw==indxheightest:
                    if max(indxweight)>=curweight:
                        pass
                    else:
                        del particlelist_new[indxsame[ilw]]
                        particlelist_new.append([xc,yc,curweight])
                else:
                    del particlelist_new[indxsame[ilw]]
                    #del indxweight[indxsame[ilw]]
                    
                
    
    #transform into array....
    particlelist=np.array(particlelist_new) #shape: num particles, 2(xy position)
    particlelist_weight=particlelist[:,2]
    return [particlelist,particlelist_weight]

def select_particles(particlelist,particlelist_weight,percentage=0.05,num_m=0,minvalue=0):
    if num_m==0:
        num_m=5
    num10=int(0.2*particlelist_weight.shape[0])
    if num10>num_m: #this is to make sure that analysis is not messed up in case many non-cellular particles are detected
        num10=num_m
    ind_highest=particlelist_weight.argsort()[-num10]
    avweight=particlelist_weight[ind_highest]#
    if minvalue>0 and avweight<minvalue:
        avweight=minvalue
    weightselection=avweight*percentage
    ind_selected=particlelist_weight>weightselection
    particlelist_selected=particlelist[ind_selected,:] 
    return [particlelist_selected,weightselection]

#go through differnt detected cell runs and detect trajectories
def analysis(strcriteria1="",strcriteria2="",statisticalanalysis=False,framemax=3000,minduration=3,minlength_analysis=50,onlyselectiondisplay=False,repeat=True,excplude1=""):

    filenmaelist=[]
    listdir=os.listdir(data_detection)
    
    for il in range(0,len(listdir)):
        curn=listdir[il][:-4]
        if listdir[il][-4:]==".npz":
            #print curn
            skip=False
                
            if strcriteria1=="":
                pass
            else:
                if strcriteria1 in curn:
                    pass
                else:
                    skip=True
            if strcriteria2=="":
                pass
            else:
                if strcriteria2 in curn:
                    pass
                else:
                    skip=True
            
            if excplude1=="":
                pass
            else:
                if excplude1 in curn:
                    skip=True
                    
            if "_background" in curn:
                skip=True
            if "_adaptive" in curn:
                skip=True
            if "_short" in curn:
                skip=True
            
            if skip==False:
                filenmaelist.append(curn)
                
    for ifile in filenmaelist:
        
        #print ifile
        dotra=False
        if 3>2:
        
            if 3>2:#==False:
                
                ifilec=os.path.basename(ifile)
                if ifilec[-4:]==".avi":
                    ifilec=ifilec[:-4]
                filenameoutc=os.path.join(data_trajectories,ifilec+".pad")
                if os.path.exists(filenameoutc)==False and onlyselectiondisplay==False:
                    
                    dotra=True
                elif not os.path.exists(filenameoutc): 
                    print "Trajectory file does not exist: "+ifilec
                    
                else:
                    print "Trajectory file does exist: "+ifilec
                    
           
                if repeat==True and onlyselectiondisplay==False:
                    dotra=True
                
                
            if dotra==True:
                print "Generating trajectory: "+ifile
                    
                generate_trajectories(ifile,framemaxin=framemax,search_range=20, memory=4)
            
            if statisticalanalysis:               
            
                analyze_trajectories(ifile,minduration=3)
                #trajectory_statistics(ifile,framemaxin=framemax,minlength_analysis=minlength_analysis)

        
def get_folderentries(folder, strcriteria1="",strcriteria2="",excplude1="",ending=".pad"):

    filenmaelist=[]
    listdir=os.listdir(folder)
    
    for il in range(0,len(listdir)):
        curn=listdir[il][:-4]
        if listdir[il][-4:]==ending:
            #print curn
            skip=False
                
            if strcriteria1=="":
                pass
            else:
                if strcriteria1 in curn:
                    pass
                else:
                    skip=True
            if strcriteria2=="":
                pass
            else:
                if strcriteria2 in curn:
                    pass
                else:
                    skip=True
            
            if excplude1=="":
                pass
            else:
                if excplude1 in curn:
                    skip=True
                    
            if "_background" in curn:
                skip=True
            if "_adaptive" in curn:
                skip=True
            if "_short" in curn:
                skip=True
            
            if skip==False:
                filenmaelist.append(curn)
                
    return filenmaelist

def detect_particles(imagein,imagein_timeav,blur_parameter=20,plot_controloutput=True,plot_outputname="",plot_annotatetext="",selectionfraction=0.05,selectionfraction2="",min_value_selectionbackground=0,fluorescencedata=False):
    if selectionfraction2=="":
        selectionfraction2=selectionfraction
    gray = cv2.cvtColor(imagein, cv2.COLOR_BGR2GRAY)
    
    gray_numpy = np.asarray(gray, dtype=float)
    
    blurred_grayscale_numpy = scipy.ndimage.filters.gaussian_filter(gray_numpy, blur_parameter)
    
    if fluorescencedata==False:
        rescaled = 255-gray_numpy#np.abs(original_grayscale - np.nanmin(original_grayscale))#(multiplier * blurred_grayscale));
        rescaled_blurred = 255-blurred_grayscale_numpy
    else:
        rescaled = gray_numpy#np.abs(original_grayscale - np.nanmin(original_grayscale))#(multiplier * blurred_grayscale));
        rescaled_blurred = blurred_grayscale_numpy
    
    
    #peaks removed
    #go through different sub areas and determin max values
    if fluorescencedata:
        arlist=np.hsplit(gray_numpy,32) #go through 80 subareas
    else:
        arlist=np.hsplit(gray_numpy,80)
    maxlist=[]
    for il in range(0,len(arlist)):
        maxlist.append(np.nanmax(arlist[il]))
    maxlist=np.array(maxlist)
    maxlist[::-1].sort()
    cutoffmax=maxlist[8] #take 8th heights array as max cufoff
    cutoffind = gray_numpy > cutoffmax
    gray_numpy_cutoff=np.copy(gray_numpy)
    if fluorescencedata==False: #do not use rescaling for fluorescence data
        gray_numpy_cutoff[cutoffind] = cutoffmax
    substracted=cv2.subtract(rescaled,rescaled_blurred)
    scalefactor=255./np.nanmax(substracted)
    scalefactor=255./65.
    substracted=substracted*scalefactor 
    [particlenum_treshold,tresholdconstant]=tresholdanalysis(substracted,mintreshold=10)

        
    #detect particles with determined treshold value
    [haeh,tresholdimage]=cv2.threshold(substracted,tresholdconstant,255,cv2.THRESH_BINARY)
    labelarray, particle_count = scipy.ndimage.measurements.label(tresholdimage)
    particlelist=ndimage.measurements.center_of_mass(rescaled,labels=labelarray,index=range(1,particle_count+1))
    particlelist_weight=ndimage.measurements.sum(rescaled,labels=labelarray,index=range(1,particle_count+1))
    particlelist_size=ndimage.measurements.sum(tresholdimage,labels=labelarray,index=range(1,particle_count+1))/255.
    [particlelist,particlelist_weight]=remove_doubles(particlelist,particlelist_weight)
    substracted2=cv2.subtract(rescaled,imagein_timeav)
    scalefactor2=255./np.nanmax(substracted2)
    scalefactor2=255./22.
    substracted2=substracted2*scalefactor2 
    [particlenum_treshold2,tresholdconstant2]=tresholdanalysis(substracted2,mintreshold=40)
    [haeh,tresholdimage2]=cv2.threshold(substracted2,tresholdconstant2,255,cv2.THRESH_BINARY)
    labelarray2, particle_count2 = scipy.ndimage.measurements.label(tresholdimage2)
    particlelist2=ndimage.measurements.center_of_mass(substracted2,labels=labelarray2,index=range(1,particle_count2+1))
    particlelist_weight2=ndimage.measurements.sum(substracted2,labels=labelarray2,index=range(1,particle_count2+1))
    [particlelist2,particlelist_weight2]=remove_doubles(particlelist2,particlelist_weight2)
    [particlelist_selected,weightselection]=select_particles(particlelist,particlelist_weight,selectionfraction) #percentage to sort out
    [particlelist_selected2,weightselection2]=select_particles(particlelist2,particlelist_weight2,selectionfraction2,num_m=1,minvalue=min_value_selectionbackground) #percentage to sort out
    strout1= "# "+str(particle_count)+"\n #_sel"+str(round(particlelist_selected.shape[0]))+"\n Wmax"+str(round(np.nanmax(particlelist_weight)))+"\n Wu: "+str(round(weightselection/selectionfraction))
    strout2= "# "+str(particle_count2)+"\n #_sel"+str(round(particlelist_selected2.shape[0]))+"\n Wmax"+str(round(np.nanmax(particlelist_weight2)))+"\n Wu: "+str(round(weightselection2/selectionfraction2))
    #combine
    particlelist_selected_combined=np.append(particlelist_selected,particlelist_selected2,axis=0)
    [particlelist_selected_combined,particlelist_selected_combined_weight]=remove_doubles(particlelist_selected_combined,particlelist_selected_combined[:,2])
    #plot 
    if plot_controloutput:
        imgsize= gray.shape
        
        fign=matplotlib.pyplot.figure(figsize=(15,15)) 
        #fign.set_canvas(matplotlib.pyplot.gcf().canvas)
        #fign.set_canvas(matplotlib.pyplot.gcf().canvas)
            
        if plot_outputname=="":
            ax1=fign.add_subplot(331)
            ax2=fign.add_subplot(332)
            ax3=fign.add_subplot(333)
            ax4=fign.add_subplot(334)
            ax5=fign.add_subplot(335)
            ax6=fign.add_subplot(336)
            ax7=fign.add_subplot(337)
            ax8=fign.add_subplot(338)
            ax9=fign.add_subplot(339)
            
            fign2=matplotlib.pylab.Figure(figsize=(30,30)) 
            axseparate=fign2.add_subplot(111)
            
            
            
            
            ax1.set_xticks([])   
            ax1.set_yticks([])   
            ax2.set_xticks([])   
            ax2.set_yticks([])  
            ax3.set_xticks([])   
            ax3.set_yticks([])   
            ax4.set_xticks([])   
            ax4.set_yticks([]) 
            ax5.set_xticks([])   
            ax5.set_yticks([]) 
            axseparate.set_xticks([])   
            axseparate.set_yticks([]) 
            #plot image
            ax1.set_ylim(0,imgsize[0])
            ax1.set_xlim(0,imgsize[1])
            ax2.set_ylim(0,imgsize[0])
            ax2.set_xlim(0,imgsize[1])
            ax3.set_ylim(0,imgsize[0])
            ax3.set_xlim(0,imgsize[1])
            ax4.set_ylim(0,imgsize[0])
            ax4.set_xlim(0,imgsize[1])
            ax5.set_ylim(0,imgsize[0])
            ax5.set_xlim(0,imgsize[1])
            axseparate.set_ylim(0,imgsize[0])
            axseparate.set_xlim(0,imgsize[1])
        else:
            #gs = matplotlib.gridspec.GridSpec(4, 3)
            
            gs = matplotlib.gridspec.GridSpec(4, 3,#numer rows, number  colums
                                   width_ratios=[1,1,1],
                                   height_ratios=[1,1,0.5,0.5]
                                   )               
            gs.update(wspace=0.4, hspace=0.4) #wspace left/right, hspace top/bottom
    
            
            
            
            ax1mout = fign.add_subplot(gs[0:2, :])
            ax2mout = fign.add_subplot(gs[2, 0])
            ax3mout = fign.add_subplot(gs[2,1])
            ax4mout = fign.add_subplot(gs[2,2])
            ax2bmout = fign.add_subplot(gs[3, 0])
            ax3bmout = fign.add_subplot(gs[3,1])
            ax4bmout = fign.add_subplot(gs[3,2])
                        
            ax1mout.set_xticks([])   
            ax1mout.set_yticks([])   
            ax1mout.set_ylim(0,imgsize[0])
            ax1mout.set_xlim(0,imgsize[1])
            
            
            plot_annotatetext=plot_annotatetext+", cells: "+str(particlelist_selected.shape[0])
            ax1mout.annotate(plot_annotatetext, xy=(0.,1.02),xycoords='axes fraction', color='k',fontsize=32)
    
            ax1mout.annotate(strout1, xy=(1.01,0.8),xycoords='axes fraction', color='k',fontsize=20)
            ax1mout.annotate(strout2, xy=(1.01,0.5),xycoords='axes fraction', color='k',fontsize=20)

    
            ax2mout.set_xlabel("intensity")
            ax3mout.set_xlabel("treshold")
            ax4mout.set_xlabel("cluster weight")
            
            ax2bmout.set_xlabel("intensity")
            ax3bmout.set_xlabel("treshold")
            ax4bmout.set_xlabel("cluster weight")
            
            ax2mout.set_ylabel("abundance")
            ax3mout.set_ylabel("# particles")
            ax4mout.set_ylabel("abundance")
            
            ax2bmout.set_ylabel("abundance")
            ax3bmout.set_ylabel("# particles")
            ax4bmout.set_ylabel("abundance")
            
            ax2mout.annotate("adapted contrast method (all cells)", xy=(0.,1.05),xycoords='axes fraction', color='k',fontsize=20)
            ax2bmout.annotate("average backround (only moving cells)", xy=(0.,1.05),xycoords='axes fraction', color='k',fontsize=20)
    
        
        sizemarker=500
        if plot_outputname=="":
        
        
            ax1.imshow(gray)
            ax2.imshow(gray,cmap='gray',vmin=0,vmax=255)
            ax3.imshow(blurred_grayscale,cmap='gray',vmin=0,vmax=255)
            axseparate.imshow(substracted)
            #ax4.imshow(blurred_grayscale)
        
        
        
            ax8.hist(gray.flatten(),bins=255,range=(0,256))#np.nanmax(gray)+1))
            ax9.hist(substracted.flatten(),bins=256,range=(0,256))#np.nanmax(gray)+1))
        
            
            #plot detected particles into 
            ax5.imshow(gray)
            for ip in range(0,particlelist_selected.shape[0]):
                ax5.scatter(particlelist_selected[ip,1],particlelist_selected[ip,0],s=sizemarker,edgecolors='w', facecolors='none')
                ax4.scatter(particlelist_selected[ip,1],particlelist_selected[ip,0],s=sizemarker,edgecolors='w', facecolors='none')
                axseparate.scatter(particlelist_selected[ip,1],particlelist_selected[ip,0],s=sizemarker,edgecolors='w', facecolors='none')
            #facecolors='none', edgecolors='r'
            ax7.plot(particlenum_treshold[:,0],particlenum_treshold[:,1])
            ax7.set_xlim(0,np.nanmax(particlenum_treshold[:,0]))
            ax7.set_ylim(0,1000)
            ax7.axvline(tresholdconstant,color='k',ls='--')
            ax6.hist(particlelist_weight,bins=100,range=(0,np.nanmax(particlelist_weight)))#np.nanmax(gray)+1))
            ax6.axvline(weightselection,color='k',ls='--')
        else:
            
        
            ax1mout.imshow(gray_numpy_cutoff)
            ax2mout.hist(gray.flatten(),bins=255,range=(0,256),color='red')#np.nanmax(gray)+1))
            ax2mout.hist(substracted.flatten(),bins=256,range=(0,256),color='blue')#np.nanmax(gray)+1))
            #ax2bmout.hist(gray2.flatten(),bins=255,range=(0,256),color='red')#np.nanmax(gray)+1))
            ax2bmout.hist(substracted2.flatten(),bins=256,range=(0,256),color='blue')#np.nanmax(gray)+1))
                        
            for ip in range(0,particlelist_selected.shape[0]):
                ax1mout.scatter(particlelist_selected[ip,1],particlelist_selected[ip,0],s=sizemarker,edgecolors='w', facecolors='none')
            
            for ip in range(0,particlelist_selected2.shape[0]):
                ax1mout.scatter(particlelist_selected2[ip,1],particlelist_selected2[ip,0],s=sizemarker,edgecolors='b', facecolors='none',linestyle=':')
            
            if 3>4: #plot details 
                for ip in range(0,particlelist.shape[0]):
                    ax1mout.scatter(particlelist[ip,1],particlelist[ip,0],s=sizemarker*0.2,edgecolors='r', facecolors='none')
                    labeltext=str(particlelist_weight[ip])+", "+str(particlelist_size[ip])
                    print labeltext
                    ax1mout.text(particlelist[ip,1],particlelist[ip,0],labeltext)
                
                
            #facecolors='none', edgecolors='r'
            ax3mout.plot(particlenum_treshold[:,0],particlenum_treshold[:,1])
            ax3mout.set_xlim(0,np.nanmax(particlenum_treshold[:,0]))
            ax3mout.set_ylim(0,1000)
            
            ax3mout.axvline(tresholdconstant,color='k',ls='--')
            ax4mout.hist(particlelist_weight,bins=300,range=(0,10*weightselection))#np.nanmax(particlelist_weight)))#np.nanmax(gray)+1))
            ax4mout.axvline(weightselection,color='k',ls='--')
            
            
            
            ax3bmout.plot(particlenum_treshold2[:,0],particlenum_treshold2[:,1])
            ax3bmout.set_xlim(0,np.nanmax(particlenum_treshold2[:,0]))
            ax3bmout.set_ylim(0,1000)
            
            ax3bmout.axvline(tresholdconstant2,color='k',ls='--')
            ax4bmout.hist(particlelist_weight2,bins=300,range=(0,10*weightselection2))#np.nanmax(particlelist_weight2)))#np.nanmax(gray)+1))
            ax4bmout.axvline(weightselection2,color='k',ls='--')
            
        if plot_outputname=="":
            fign.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,hspace=.1,wspace=.1)#,wspace=5.0, hspace=.1)
            figname="plot_detectedcells.png"
            #fign.show()
            #fign.draw()
            fign.savefig(os.path.join(folder_detectiontest,figname))
            
            fign2.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,hspace=.1,wspace=.1)#,wspace=5.0, hspace=.1)
            figname="plot_detectedcells_zoom.pdf"
            fign2.savefig(os.path.join(folder_detectiontest,figname))
        
            cv2.imwrite(os.path.join(folder_detectiontest,"testimage.jpg"),gray)
            cv2.imwrite(os.path.join(folder_detectiontest,"testimage_substracted.jpg"),substracted)
        else:
            fign.subplots_adjust(left=0.1, bottom=0.07, right=0.95, top=0.95,hspace=.1,wspace=.1)#,wspace=5.0, hspace=.1)
            figname=plot_outputname
            #fign.show()
            #fign.canvas.draw_idle()
            fign.savefig(figname)
            print plot_annotatetext.split(",")[0]
            cv2.imwrite(os.path.join(folder_detectiontest,plot_annotatetext.split(",")[0].replace(": ","")+".jpg"),gray)
            
        matplotlib.pylab.close()
         
    return [particlelist_selected_combined,particlelist_selected,particlelist_selected2]

def analyze_movie(filenamein,framemaxin=2000,plot_detectionanalysis=True,num_avsteps=20, blur_parameter=20, selectionfraction=0.05, selectionfraction2="",short=False,min_value_selectionbackground=0,fulloutput=True,fluorescencedata=False):
    if selectionfraction2=="":
        selectionfraction2=selectionfraction
    if short==True:
        adds="_short"
    else:
        adds=""
    movieoutputdir2=os.path.join(plot_controloutput,filenamein+"_detection"+adds)
    if not os.path.exists(movieoutputdir2):
        os.makedirs(movieoutputdir2)
    else:
        pass
    filename=os.path.join(datafolderrawvideo2,filenamein)
    #check if filename exists...
    if os.path.exists(filename)==False:
        filename=os.path.join(datafolderrawvideo1,filenamein)
        if os.path.exists(filename)==False:
            print "File not found: "+filename
            errorfilenotfound
    
    filenameinm=os.path.basename(filenamein[:-4])
    filenameout=os.path.join(data_detection,filenameinm+adds+".npz")
    filenameout_adaptedcontrast=os.path.join(data_detection,filenameinm+adds+"_adaptive.npz")
    filenameout_background=os.path.join(data_detection,filenameinm+adds+"_background.npz")
    
    filenameoutmat=filenameout[:-4]+adds+".mat"
    filenameout_adaptedcontrastmat=filenameout_adaptedcontrast[:-4]+adds+".mat"
    filenameout_backgroundmat=filenameout_background[:-4]+adds+".mat"

    particlesdetect_list={}
    particlesdetect_adaptedcontrast_list={}
    particlesdetect_background_list={}
    
    print "load avi file: "+filename
    vc = cv2.VideoCapture(filename)
    framecounter=0
    length_movie = int(vc.get(cv2.CAP_PROP_FRAME_COUNT))
    if length_movie<framemaxin:
        framemaxin=length_movie
    
    ####load movie first time to calculate average
    imglist=[]
    if vc.isOpened():
        rval , frame = vc.read()
        imglist.append(frame)
    else:
        rval = False
        
    if os.path.exists(filename)==False:
        print "File not found: "+filename
        here78

    while rval:
        if framecounter>num_avsteps:
            break
        #print framecounter
        rval, frame = vc.read()
        imglist.append(frame)
        framecounter = framecounter + 1    
    
    #get image average
    if fluorescencedata:
        imgc_av=np.zeros([512,512,num_avsteps])
    else:
        imgc_av=np.zeros([1024,1280,num_avsteps])
    for ila in range(0,num_avsteps):
        imgc=imglist[ila].copy()
        
        gray = cv2.cvtColor(imgc, cv2.COLOR_BGR2GRAY)
        if fluorescencedata==False:
            imgc_av[:,:,ila] = 255.-np.asarray(gray, dtype=float)
        else:
            imgc_av[:,:,ila] = np.asarray(gray, dtype=float)
    imgc_avv=np.average(imgc_av,axis=2)

    #end movie....close again....
    #release memmory for avi file...
    vc.release()
    
    
    #load movie second time to detect cells
    vc = cv2.VideoCapture(filename)
    
    #make sure frame max is not too long...
    framemaxinvideo = int(vc.get(cv2.CAP_PROP_FRAME_COUNT))
    if framemaxinvideo-1<=framemaxin:
        framemaxin=framemaxinvideo-2
    framecounter=0
    
    if vc.isOpened():
        rval , frame = vc.read()
        imgc=frame
    else:
        rval = False
    
    while rval:
        if framecounter>framemaxin:
            break
        print(framecounter)
        rval, frame = vc.read()
        imgc=frame
        framecounter = framecounter + 1
    
        #generate name (for movie)
        framecounterstr=str(framecounter)
        if framecounter<10:
            mvstr="000"+framecounterstr
        elif framecounter<100:
            mvstr="00"+framecounterstr
        elif framecounter<1000:
            mvstr="0"+framecounterstr
        else:
            mvstr=""+framecounterstr
        figname=mvstr+".jpg"
        
        plot_annotatetext="frame: "+str(framecounter)+", time: "+str(framecounter/20.)+"s"
        
        if plot_detectionanalysis==True:
            plot_outputname_detection=os.path.join(movieoutputdir2,figname)
        else:
            plot_outputname_detection=""
        [particlesdetect,particlesdetect_adaptedcontrast,particlesdetect_background]=detect_particles(imgc,imgc_avv,plot_controloutput=plot_detectionanalysis,plot_outputname=plot_outputname_detection,plot_annotatetext=plot_annotatetext,selectionfraction=selectionfraction,blur_parameter=blur_parameter,selectionfraction2=selectionfraction2,min_value_selectionbackground=min_value_selectionbackground,fluorescencedata=fluorescencedata)
        #output particlesdetect is [numper particles, 3 (x,y,weight)]
        particlesdetect_list["frame"+str(framecounter)]=particlesdetect
        particlesdetect_adaptedcontrast_list["frame"+str(framecounter)]=particlesdetect_adaptedcontrast
        particlesdetect_background_list["frame"+str(framecounter)]=particlesdetect_background
        
        
        
    #release memmory for avi file...
    vc.release()
    
    #save detected cells
    
    diroutc=os.path.dirname(filenameout)
    if not os.path.exists(diroutc):
        os.makedirs(diroutc)
    else:
        pass
    
    np.savez(filenameout,**particlesdetect_list)
    if fulloutput==True:
        np.savez(filenameout_adaptedcontrast,**particlesdetect_adaptedcontrast_list)
        np.savez(filenameout_background,**particlesdetect_background_list)
    
    
    
    if fulloutput==True:
        convert_detectedcellsintomat(particlesdetect_list,filenameoutmat)
        convert_detectedcellsintomat(particlesdetect_adaptedcontrast_list,filenameout_adaptedcontrastmat)
        convert_detectedcellsintomat(particlesdetect_background_list,filenameout_backgroundmat)
    
    if plot_detectionanalysis:
        filenameoutdetmovie=filenameout[:-4].split('/')[-1]
    
        name_movieoutis=os.path.join(data_movieoutput,"movie_detection_"+filenameoutdetmovie+".mp4")
        opt_movie_fps=20
        process_name ="ffmpeg -y -r "+ str(opt_movie_fps)+" -pattern_type glob -i '"+os.path.join(movieoutputdir2, "*.jpg'")+" -c:v libx264 -pix_fmt yuv420p "+name_movieoutis
        subprocess.call(process_name, shell=True)
        print ("detection movie saved to: "+name_movieoutis)
    print("data detected cells saved to "+filenameout)

#if matlab file is wanted    
def convert_detectedcellsintomat(listin,filenout):
    numt=len(listin)
    FrameStack = np.empty((numt-1,3,), dtype=np.object)
    for i in range(numt-1):
        #print "frame"+str(i+1)
        x=listin["frame"+str(i+1)].copy()[:,0:2]+0.01
        x[:,1]=0
        #print str(np.nanmax(x))+" "+str(np.nanmin(x[:,0]))
        y=listin["frame"+str(i+1)].copy()[:,0:2]
        y[:,0]=listin["frame"+str(i+1)][:,1]+0.01
        y[:,1]=0
        weight=listin["frame"+str(i+1)].copy()[:,0:2]
        weight[:,1]=0
        weight[:,0]=listin["frame"+str(i+1)][:,2]
        FrameStack[i,0] = x
        FrameStack[i,1] = y
        FrameStack[i,2] = weight
        
    scipy.io.savemat(filenout, {"movieInfoCell":FrameStack})

#merge trajectories
def mergetrajectory(trj_xdata,trj_ydata,mergdistance=20.):
    num_traj=trj_xdata.shape[0] #trajectories and duration
    deleterows=[]
    print "num traj before merging"
    print num_traj
    #merge trajectories
    
    for il in range(num_traj-1,-1,-1):
        #find where trajectory begins....
        valuesc= np.isfinite(trj_ydata[il,:])
        indt=np.argwhere(valuesc==True)[0][0]
        
        if indt<3:
            break
        
        xc1=trj_xdata[il,indt]
        yc1=trj_ydata[il,indt]
        #go through other trajectories and check if they end here  
        for il2 in range(il-1,-1,-1):
            
            #check if il2 is nan
            if np.isnan(trj_xdata[il2,indt]):
               
                #check if timepoint indt-1 is number
                if np.isfinite(trj_xdata[il2,indt-1]):
                    xc2=trj_xdata[il2,indt-1]
                    yc2=trj_ydata[il2,indt-1]
                    
                    distance=np.sqrt(np.power(xc2-xc1,2.)+np.power(yc2-yc1,2.))
                    if(distance<mergdistance):
                        trj_xdata[il,:indt-1]=trj_xdata[il2,:indt-1]
                        trj_ydata[il,:indt-1]=trj_ydata[il2,:indt-1]
                        "merging2 "+str(il2)
                        trj_xdata[il,:]
                        deleterows.append(il2)
        
                        break
                #check if timepoint indt-1 is number
                elif np.isfinite(trj_xdata[il2,indt-2]):
                    
                    xc2=trj_xdata[il2,indt-2]
                    yc2=trj_ydata[il2,indt-2]
                    distance=np.sqrt(np.power(xc2-xc1,2.)+np.power(yc2-yc1,2.))
                    if(distance<mergdistance):
                        trj_xdata[il,:indt-2]=trj_xdata[il2,:indt-2]
                        trj_ydata[il,:indt-2]=trj_ydata[il2,:indt-2]
                        "merging2 "+str(il2)
                        deleterows.append(il2)
                        break
                #check if timepoint indt-1 is number
                elif np.isfinite(trj_xdata[il2,indt-3]):
                    
                    xc2=trj_xdata[il2,indt-3]
                    yc2=trj_ydata[il2,indt-3]
                    distance=np.sqrt(np.power(xc2-xc1,2.)+np.power(yc2-yc1,2.))
                    if(distance<mergdistance):
                        trj_xdata[il,:indt-3]=trj_xdata[il2,:indt-3]
                        trj_ydata[il,:indt-3]=trj_ydata[il2,:indt-3]
                        "merging3 "+str(il2)
                        deleterows.append(il2)
                        break
                    
            #check if il2-2 is number
    #delete 
    trj_xdata=np.delete(trj_xdata,deleterows,axis=0)
    trj_ydata=np.delete(trj_ydata,deleterows,axis=0)
    # x = numpy.delete(x,(2), axis=1)
    return [trj_xdata,trj_ydata]

def generate_movie(filenamein,plot_detectedcells=True,plot_trajectories=False,framemaxin=2000,trajectory_color='trajectories'):

    if filenamein[-4:]==".avi":
        filenamein=filenamein[:-4]
    filenameinbase=os.path.basename(filenamein)
           
    #load avi file using imagej
    movieoutputdir=os.path.join(plot_controloutput,filenamein)
    if not os.path.exists(movieoutputdir):
        os.makedirs(movieoutputdir)
    else:
        pass
    
    movieoutputdir2=os.path.join(plot_controloutput,filenamein+"_detection")
    if not os.path.exists(movieoutputdir2):
        os.makedirs(movieoutputdir2)
    else:
        pass

    try:
        filename=os.path.join(datafolderrawvideo2,filenamein+".avi")
    except:
                filename=os.path.join(datafolderrawvideo,filenamein+".avi")
                print os.path.join(datafolderrawvideo,filenamein+".avi")
                haeh
    
    print "load avi file: "+filename
    vc = cv2.VideoCapture(filename)
    
    length = int(vc.get(cv2.CAP_PROP_FRAME_COUNT))
    print length
    framecounter=1
    
    
    if vc.isOpened():
        rval , frame = vc.read()
        imgc=frame
    else:
        rval = False
        
    if os.path.exists(filename)==False:
        print "File not found: "+filename
        here78

    
    #load detected trajectories
    if plot_trajectories:
             
            filenametrj=os.path.join(data_trajectories,filenameinbase+".npy")
            trj_data=np.load(filenametrj)
            num_traj=trj_data.shape[0]
            
            
            num_trjovertime=np.zeros([trj_data.shape[1]])
            
            trajectorylength=np.zeros([num_traj])
            for it in range(0,num_traj):
                trajectorylength[it]=np.count_nonzero(~np.isnan(trj_data[it,:,0]))
          
            for it in range(0,trj_data.shape[1]):
                num_trjovertime[it]=np.count_nonzero(~np.isnan(trj_data[:,it,0]))
                
            #0 index
            #1 x
            #2 y
            #3 time
            #4 velocity
            #5 theta
            #6 tubmling event....
            #7 theta
            #8 tubmling event....
                
                
            #get ranges
            vmin=np.nanmin(trj_data[:,:,4])
            vmax=np.nanmax(trj_data[:,:,4])
            dvmin=np.nanmin(trj_data[:,:,5])
            dvmax=np.nanmax(trj_data[:,:,5])
            dthetamin=np.nanmin(trj_data[:,:,7])
            dthetamax=np.nanmax(trj_data[:,:,7])
            
            
           
    #load detected particles...
    if plot_detectedcells:
            
            filenamedetection=os.path.join(data_detection,filenamein+".npz")
            filenamedetection_adaptive=os.path.join(data_detection,filenamein+"_adaptive.npz")
            filenamedetection_background=os.path.join(data_detection,filenamein+"_background.npz")
            
            
            detected=np.load(filenamedetection)
            detected_adaptive=np.load(filenamedetection_adaptive)
            detected_background=np.load(filenamedetection_background)
            


    while rval:
        if framecounter>framemaxin:
            break
        print framecounter
        rval, frame = vc.read()
        imgc=frame
        framecounter = framecounter + 1
    
        gray = cv2.cvtColor(imgc, cv2.COLOR_BGR2GRAY)
        #go through every frame....
        framecounterstr=str(framecounter)
        if framecounter<10:
            mvstr="000"+framecounterstr
        elif framecounter<100:
            mvstr="00"+framecounterstr
        elif framecounter<1000:
            mvstr="0"+framecounterstr
        else:
            mvstr=""+framecounterstr
        figname=mvstr+".jpg"
        
        fign=figure(figsize=(20,12)) 
        
        gs = matplotlib.gridspec.GridSpec(2, 3,#numer rows, number  colums
                                   width_ratios=[1,1,0.1],
                                   height_ratios=[1,0.5]
                                   )               
        gs.update(wspace=0.4, hspace=0.4) #wspace left/right, hspace top/bottom
    
        ax1=fign.add_subplot(gs[0, 0])
        ax2=fign.add_subplot(gs[0, 1])
        ax3=fign.add_subplot(gs[1,0])
        
        if trajectory_color in ["speed","acceleration","dtheta"]:
            ax2c=fign.add_subplot(gs[0,2])
            cb=matplotlib.colorbar.ColorbarBase(ax2c, cmap=cm.jet)
            cb.set_ticks([0,vmax/2.,vmax])
            
        ax1.set_xticks([])   
        ax1.set_yticks([])   
        ax2.set_xticks([])   
        ax2.set_yticks([])  
        imgsize= imgc.shape
        
        #plot image
        ax1.set_ylim(0,imgsize[0])
        ax1.set_xlim(0,imgsize[1])
        ax2.set_ylim(0,imgsize[0])
        ax2.set_xlim(0,imgsize[1])

        ax1.imshow(gray)
        
        if plot_detectedcells:
            
            plot_annotatetext="frame: "+str(framecounter)+", time: "+str(framecounter/20.)+"s"
            try:
                xc=detected["frame"+str(framecounter-1)][:,0]
                yc=detected["frame"+str(framecounter-1)][:,1]
            except:
                "cell detection not found"
                xc[:]=np.nan
                yc[:]=np.nan
                
            
            
            detected=np.load(filenamedetection)
            detected_adaptive=np.load(filenamedetection_adaptive)
            detected_background=np.load(filenamedetection_background)
                            
                
            try:
                xc_background=detected_background["frame"+str(framecounter-1)][:,0]
                yc_background=detected_background["frame"+str(framecounter-1)][:,1]
            except:
                "cell detection not found"
                xc_background[:]=np.nan
                yc_background[:]=np.nan

            try:
                xc_adaptive=detected["frame"+str(framecounter-1)][:,0]
                yc_adaptive=detected["frame"+str(framecounter-1)][:,1]
            except:
                "cell detection not found"
                xc_adaptive[:]=np.nan
                yc_adaptive[:]=np.nan

            ax1.scatter(yc_adaptive,xc_adaptive,s=100,edgecolors='w', facecolors='none')

            ax1.scatter(yc_background,xc_background,s=100,edgecolors='b', facecolors='none',linestyle=":")
        timestr=str(round(framecounter/20.,2))+" seconds"
        ax2.annotate(timestr, xy=(0.1,1.1),xycoords='axes fraction', color='k',fontsize=16)

        if plot_trajectories:
            
            numbertimeframestrajectory=1000
            minind=framecounter-numbertimeframestrajectory  
            
            if minind<0:
                minind=0
            colorlist=['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']
            num_c=len(colorlist)
                        
            #plot tumbling events
                        
            for it in range(0,num_traj):
                if framecounter>2:
                    for frr in range(minind,framecounter):
                        if trj_data[it,frr,8]==1 and np.isnan(trj_data[it,frr-1,8]):
                            ax2.scatter(trj_data[it,frr,2],trj_data[it,frr,1],facecolors='none',edgecolors='k',s=100,alpha=0.5)
               
                        
            if trajectory_color=='trajectories':
                
                for it in range(0,num_traj):
                    colorc=colorlist[it%num_c]
                    ax2.plot(trj_data[it,minind:framecounter,2],trj_data[it,minind:framecounter,1],color=colorc)
                
                    
            
            elif trajectory_color=='speed':
                for it in range(0,num_traj):
                    colorc=colorlist[it%num_c]
                    colorv=(trj_data[it,minind:framecounter,4]+0)/(0+vmax)
                    colorc=cm.jet(colorv)
                    print vmin
                    print vmax
                    #here78
                    ax2.scatter(trj_data[it,minind:framecounter,2],trj_data[it,minind:framecounter,1],color=colorc,edgecolors='none',s=10)
                           
            ax3.hist(trajectorylength,bins=100,range=(0,1000))            
            
        #save
        fign.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,hspace=.1,wspace=.1)#,wspace=5.0, hspace=.1)
        fign.savefig(os.path.join(movieoutputdir,figname))
        close()
    #release memmory for avi file...
    vc.release()
    close()
    
    name_movieoutis=os.path.join(data_movieoutput,"movie_"+filenameinbase+".mp4")
    opt_movie_fps=10
    print("...generating movie file....")
    process_name ="ffmpeg -y -r "+ str(opt_movie_fps)+" -pattern_type glob -i '"+os.path.join(movieoutputdir, "*.jpg'")+" -c:v libx264 -pix_fmt yuv420p "+name_movieoutis
    subprocess.call(process_name, shell=True)
    print("Movie savet to: "+name_movieoutis)
    
def weighted_sdt(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return np.sqrt(variance)
    
    
    
    

def plot_together(listin,labellistin,odlist="",growthratelist=[],timelist=[],filenameout="distributioncomparison",minlength_analysis=0,minduration=3,vmaxhis=50,vhisbins=20,plot_trajectories=False,vsmode="growth"):#,framemaxin=2000,minlength_analysis=100):
    
    num_runs=len(listin)
    num_t=[]
    resultslist=[]
    for ir in range(0,num_runs):
        print listin[ir]
        num_t.append(len(listin[ir]))
        resultslist.append(np.zeros([num_t[-1],22]))
        resultslist[-1][:]=np.nan
        #0: time
        #1: od
        #2: frac_motility
        #3: mean velocity
        #4: mean velocity mobile cells
        #5: mean run time.... 
        #6: mean run time mobile cells
    
    num_tmax=max(num_t)
    
    #load od information form file....
    if odlist=="" or timelist=="": #alternative, set index
        for ir in range(0,num_runs):
            for it in range(0,num_t[ir]):
                resultslist[ir][it,0]=it
                resultslist[ir][it,1]=it
    else:
        for ir in range(0,num_runs):
            for it in range(0,num_t[ir]):
                resultslist[ir][it,0]=timelist[ir][it]/60.
                resultslist[ir][it,1]=odlist[ir][it]
    fign=plt.figure(figsize=(num_runs*8,num_tmax*7)) 
            
    gs = matplotlib.gridspec.GridSpec(num_tmax*2, num_runs,#numer rows, number  colums
                                       width_ratios=[1]*num_runs,
                                       height_ratios=[1,0.5]*num_tmax
                                       )               
    gs.update(wspace=1.0, hspace=0.4) #wspace left/right, hspace top/bottom
          
    listcvv=[]
    #go through different runs 
    for ir in range(0,num_runs):
        for iTT in range(0,num_t[ir]):
            filenamec=listin[ir][iTT]
            if filenamec[-4:]==".avi":
                filenamec=filenamec[:-4]
            filenameinbase=os.path.basename(filenamec)
    

            #load trajectory information.....

            filenametrj=os.path.join(data_trajectories,filenameinbase+".npy")
            filenametrjpad=os.path.join(data_trajectories,filenameinbase+".pad")
            
            if not os.path.exists(filenametrj):
               if os.path.exists(filenametrjpad):
                   print filenametrj
                   print "trajectories not analyzed. Analyze file first."
                   #analyze_trajectories(filenametrjpad,minduration=minduration)

               else:
                   print "!!!Trajectory file not found: "+filenameinbase
               #here78
            else:
                print "file found: "+filenameinbase
                trj_data=np.load(filenametrj)
                
                
                
                #### analysis trajectories...
                num_traj=trj_data.shape[0]
                num_trjovertime=np.zeros([trj_data.shape[1]])
                trajectorylength=np.zeros([num_traj])
                for it in range(0,num_traj):
                    trajectorylength[it]=np.count_nonzero(~np.isnan(trj_data[it,:,0]))
                for it in range(0,trj_data.shape[1]):
                    num_trjovertime[it]=np.count_nonzero(~np.isnan(trj_data[:,it,0]))
                    
                #selected trajectories.....
                indselect=[]
                indnonselect=[]
                for it in range(0,num_traj):
                    if trajectorylength  [it]>minlength_analysis:
                        indselect.append(it)
                    else:
                        indnonselect.append(it)
                
                #get basic properties of trajectories...
                num_trajS=len(indselect)
                TDspeed=np.zeros([num_trajS])
                #TDrunningspeed=np.zeros([num_trajS])
                TDlength=np.zeros([num_trajS])
                
                #add more variables....
                dt=0.05
                its=-1
                for it in indselect:
                    its=its+1
                    lengthc=np.nansum(trj_data[it,:,4])
                    timec=np.count_nonzero(~np.isnan(trj_data[it,:,0]))*dt
                    TDspeed[its]=lengthc/timec
                    TDlength[its]=timec
 
 
                ####get run-time
                itc=-1
                firstfast=True
                for it in indselect:
                    itc=itc+1
                    tumblind=np.where(trj_data[it,:,9]==1)
                    #ediff is only taking time between differences...                        
                    if itc==0:
                        tumblingdiff=np.ediff1d(tumblind)*dt #assumes camera at 20hz
                    else:
                        tumblingdiff=np.append(tumblingdiff,np.ediff1d(tumblind)*dt)
                    if TDspeed[itc]>10:
                        if firstfast:
                            tumblingdifffast=np.ediff1d(tumblind)*dt
                            firstfast=False
                        else:
                            tumblingdifffast=np.append(tumblingdifffast,np.ediff1d(tumblind)*dt)
                            
                            
                tumblingdiff=tumblingdiff[tumblingdiff>0.05]#select only non-nearest neighbours
                tumblingdifffast=tumblingdifffast[tumblingdifffast>0.05]                
                tubltav=np.average(tumblingdiff)
                tubltavfast=np.average(tumblingdifffast)
                
                    
                totaltime=trj_data[0,:,0].shape[0]/20.
                #0 index
                #1 x
                #2 y
                #3 time
                #4 velocity
                #5 theta
                #6 tubmling event....
                #7 theta
                #8 tubmling event....
                #9 tumbling....
                    
                    
                #get ranges
                vmin=np.nanmin(trj_data[:,:,4])
                vmax=np.nanmax(trj_data[:,:,4])
                print vmin
                print vmax
                dvmin=np.nanmin(trj_data[:,:,5])
                dvmax=np.nanmax(trj_data[:,:,5])
                print dvmin
                print dvmax
                dthetamin=np.nanmin(trj_data[:,:,7])
                dthetamax=np.nanmax(trj_data[:,:,7])
                print dthetamin
                print dthetamax
                        
                weight=TDlength/np.sum(TDlength)
                #so now we have trajectories and lenght of those. Do statistics of those

                indxfast=TDspeed>10.
                weightfast=weight[indxfast]/np.sum(weight[indxfast])
                meanspeed=np.average(TDspeed,weights=weight)
                
                avnumtrj=num_trajS*np.average(TDlength,weights=weight)/totaltime #average number of traj. each timepoint
                stdspeed=weighted_sdt(TDspeed,weight)
                meanspeedfast=np.average(TDspeed[indxfast],weights=weightfast)
                stdspeedfast=weighted_sdt(TDspeed[indxfast],weightfast)
                fractionMobile=np.sum(weight[indxfast])
                print "-------------------"
                print "Trajectory analysis for: "+filenameinbase
                print "-------------------"
                print "Number of trajectories: "+str(num_trajS)+" [min frames "+str(minlength_analysis)+"; total traj.:"+str(num_traj)+"]"
                print "Mean trajectory length [s]:"+str(np.average(TDlength,weights=weight))
                print "Mean speed [mu m/s]: "+str(meanspeed)+" +- "+str(stdspeed)
                print "Mean speed (fast trj., >10) [mu m/s]: "+str(meanspeedfast)+" +- "+str(stdspeedfast)
                print "-------------------"
                
                #for output with plot
                txlist=[]
                textan=filenameinbase
                txlist.append(textan)        
                textan="#Trj: "+str(num_trajS)+" [min frames "+str(minlength_analysis)+"; total traj.:"+str(num_traj)+"]" 
                txlist.append(textan)
                textan="mean length:"+str(np.average(TDlength,weights=weight))+"s"
                txlist.append(textan)
                textan="mean speed: "+str(round(meanspeed,2))+" +- "+str(round(stdspeed,2))+"$\\mu m/s$"
                txlist.append(textan)   
                textan="mean speed (>10): "+str(round(meanspeedfast,2))+" +- "+str(round(stdspeedfast,2))+"$\\mu m/s$"
                txlist.append(textan)    
                textan="fraction motile (>10): "+str(round(fractionMobile*100,0))+"%"
                txlist.append(textan) 
                
                
                
                #update 
                resultslist[ir][iTT,2]=fractionMobile*100.
                resultslist[ir][iTT,3]=meanspeed
                resultslist[ir][iTT,4]=meanspeedfast
                resultslist[ir][iTT,5]=tubltav
                resultslist[ir][iTT,6]=tubltavfast
                
                ###put 
    
                #WEIGHTED HISTOGRAM
                ax1=fign.add_subplot(gs[iTT*2,ir])
                if vmaxhis==0:
                    vmaxuse=np.nanmax(TDspeed)
                else:
                    vmaxuse=vmaxhis
                
                ax1.set_xlim(0,vmaxuse)                    
                ax1.hist(TDspeed,weights=TDlength,normed=True,bins=vhisbins,range=(0,vmaxuse))#np.max(TDspeed)*1.05))
                ax1.set_xlabel("speed [mu m/s]")
                ax1.axvline(meanspeed,alpha=0.5,color='k',lw=5)
                ax1.axvline(meanspeedfast,ls='--',alpha=0.5,color='k',lw=5)
                
                if iTT==0:
                    ax1.annotate(labellistin[ir], xy=(0.3,1.2),xycoords='axes fraction', color='k',fontsize=25)
                if odlist=="" or timelist=="":
                    pass
                else:
                    ax1.annotate("OD="+str(odlist[ir][iTT])+", time="+str(round(timelist[ir][iTT]/60.,2))+"h", xy=(0.,1.05),xycoords='axes fraction', color='k',fontsize=20)                
                
                ###give info on what was detected
                fsu=15
                
                ax1.plot(distributionfit[0],distributionfit[1],ls='-',color='orange')           
                
                ax1.annotate("$N_T=$"+str(int(round(avnumtrj))), xy=(1.02,0.9),xycoords='axes fraction', color='k',fontsize=fsu)
                ax1.annotate("$\\tau_T=$"+str(round(np.average(TDlength,weights=weight),1))+"$s$", xy=(1.02,0.8),xycoords='axes fraction', color='k',fontsize=fsu)
                ax1.annotate("$\\langle v\\rangle=$"+str(round(meanspeed,1))+"$\\pm$"+str(round(stdspeed,1))+"$um/s$", xy=(1.02,0.7),xycoords='axes fraction', color='k',fontsize=fsu)
                ax1.annotate("$\\langle v\\rangle_{f}=$"+str(round(meanspeedfast,1))+"$\\pm$"+str(round(stdspeedfast,1))+"$um/s$", xy=(1.02,0.6),xycoords='axes fraction', color='k',fontsize=fsu)
                ax1.annotate("$\\alpha_m=$"+str(int(round(fractionMobile*100,0)))+"$\\%$",xy=(1.02,0.5),xycoords='axes fraction', color='k',fontsize=fsu)
                ax1.annotate("$\\tau_r=$"+str(int(round(tubltav,2)))+"$s$",xy=(1.02,0.4),xycoords='axes fraction', color='k',fontsize=fsu)
                ax1.annotate("$\\tau_{r,f}=$"+str(int(round(tubltavfast,2)))+"$s$",xy=(1.02,0.3),xycoords='axes fraction', color='k',fontsize=fsu)
                     
                
                #plot distriution running time...
                ax2=fign.add_subplot(gs[iTT*2+1,ir])
                numtbins=20
                xatr=np.linspace(0,1.2,100)
                                
                ax2.hist(tumblingdiff,range=(0,1.2),bins=numtbins,color='r',alpha=0.5)
                ax2.axvline(tubltav,ls='-',alpha=0.5,color='r')
                scalec=tumblingdiff.shape[0]/numtbins
                ax2.plot(xatr,(scalec/tubltav)*np.exp(-(1./tubltav)*xatr),ls='-',alpha=0.5)
                ax2.set_xlabel("time between tumbling $[s]$")

                ax2.hist(tumblingdifffast,range=(0,1.2),bins=numtbins,color='b',alpha=0.5)
                ax2.axvline(tubltavfast,ls='-',alpha=0.5,color='b')
                scalec=tumblingdifffast.shape[0]/numtbins
                ax2.plot(xatr,(scalec/tubltavfast)*np.exp(-(1./tubltavfast)*xatr),ls='-',alpha=0.5)
                ax2.set_xlabel("time between tumbling $[s]$")

                
                if plot_trajectories: #plot trajecories and indicate different properties...
                    print "plot trajectories...."
                    its=-1
                    minlengthplot=3
                    numplotedtrj=0
                    colorlist=['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']
                    #generate plot....

                    trajectorymode=["trajectories","speed"]       
                    numT=len(trajectorymode)
                    figtrj=plt.figure(figsize=(len(trajectorymode)*5,7)) 
                    axltr=[]
                    axltr2=[]
                                   
                    gsT = matplotlib.gridspec.GridSpec(2, len(trajectorymode)*2,#numer rows, number  colums
                                       width_ratios=[1,0.1]*len(trajectorymode),
                                       height_ratios=[1,0.5]
                                       )               
                    gsT.update(wspace=0.3, hspace=0.4) #wspace left/right, hspace top/bottom
    
    
                    for irj in range(0,len(trajectorymode)):
                        axltr.append(figtrj.add_subplot(gsT[0,2*irj]))
                        axltr2.append(figtrj.add_subplot(gsT[1,2*irj]))
                    num_c=len(colorlist)    
                    numbers_trajectoriestoplot=30
                    firsttrj=True
                    #use first traectories:
                    #indselectTR=indselect
                    #use sorted trajectories
                    lengthTR=np.nansum(trj_data[indselect,:,0],axis=1)
                    cursel= np.argsort(lengthTR[:,])[::-1]
                    indselectTR = [ indselect[i] for i in cursel.tolist()]
                    #indselectTR=indselect[cursel.tolist()[:]]
                    #indselectTR=indselect[]
                    
                    
                    for it in indselectTR:
                        
                        its=its+1
                        lengthc=np.nansum(trj_data[it,:,4])
                        indx=~np.isnan(trj_data[it,:,0])
                        timec=np.count_nonzero(indx)*dt
                        if numplotedtrj>numbers_trajectoriestoplot:
                            break
                        if timec>minlengthplot: #selection criteria minimal trajectory length
                             
                            
                            #get trajectory and plot
                            numplotedtrj=numplotedtrj+1
                             
                            iT=-1
                            for trajectory_color in   trajectorymode:
                                iT=iT+1
                                
                                minind=0
                                framecounter=-1
                                #plot trajectory properties
                                
                                if trajectory_color=='trajectories':
                                
                                        colorc=colorlist[it%num_c]
                                        #minind:framecounter
                                        axltr[iT].scatter(trj_data[it,minind:framecounter,2],trj_data[it,minind:framecounter,1],color=colorc,s=10)
                                    
                                    
                                        tumblind=np.where(trj_data[it,minind:framecounter,9]==1)
                                        tumblindrot=np.where(trj_data[it,minind:framecounter,10]==1)
                                        tumblindspeed=np.where(trj_data[it,minind:framecounter,11]==1)
                                        
                                        
                                        
                                        #axltr[iT].scatter(trj_data[it,tumblind,2],trj_data[it,tumblind,1],color='k',s=10,alpha=0.2)
                                        axltr[iT].scatter(trj_data[it,tumblindrot,2],trj_data[it,tumblindrot,1],s=20,alpha=0.5,facecolors='none', edgecolors='r')
                                        axltr[iT].scatter(trj_data[it,tumblindspeed,2],trj_data[it,tumblindspeed,1],s=40,alpha=0.5,facecolors='none', edgecolors='b')
                                                                                                                        
                                        #ax2.plot(trj_xdataA[it,minind:framecounter],trj_ydataA[it,minind:framecounter],color=colorc)
                                        #ax2.plot(trj_xdataB[it,minind:framecounter],trj_ydataB[it,minind:framecounter],color=colorc)
                                        
                                        
                                        
                                                                                
                                        #trj_xdata
                                    
                            
                                elif trajectory_color=='speed':
                                        vmaxtr=40
                                        #colorc=colorlist[it%num_c]
                                        colorv=(trj_data[it,minind:framecounter,5]+0)/(0+vmaxtr)
                                        colorc=matplotlib.cm.jet(colorv)
                                        #here78
                                        axltr[iT].scatter(trj_data[it,minind:framecounter,2],trj_data[it,minind:framecounter,1],color=colorc,edgecolors='none',s=10)
                    
                                        
                                        
                                else:
                                    print trajectoryplotmodnotknown
                                    
                                if firsttrj:
                                    firsttrj=False
                    
                    #numpy.polyfit(x, numpy.log(y), 1, w=numpy.sqrt(y))
                    numtbins=20
                    axltr2[0].hist(tumblingdiff,range=(0,1.2),bins=numtbins)
                    axltr2[0].axvline(tubltav,ls='-',alpha=0.5)
                    xatr=np.linspace(0,1.2,100)
                    scalec=tumblingdiff.shape[0]/numtbins
                    axltr2[0].plot(xatr,(scalec/tubltav)*np.exp(-(1./tubltav)*xatr),ls='-',alpha=0.5)
                    
                    axltr2[0].set_xlabel("time between tumbling $[s]$")

                    print np.nanmin(trj_data[:,:,5])  
                    print np.nanmax(trj_data[:,:,5])
                    vlc=trj_data[:,:,5]
                    vlc = vlc[~np.isnan(vlc)]
                    axltr2[1].hist(vlc,range=(0,np.nanmax(vlc)),bins=20)
                    axltr2[1].set_xlabel("speed $[\mu m/s]$")
                                        #axltr2[iT].set_ylabel("")
                    figtrj.subplots_adjust(left=0.1/float(numT), bottom=0.25/float(numT), right=1.-0.3/float(numT), top=1.-0.2/float(numT),hspace=.1,wspace=.1)#,wspace=5.0, hspace=.1)
                    folderouttr=os.path.join(trjstat_output,filenameout)
                    if not os.path.exists(folderouttr):
                        os.makedirs(folderouttr)
                    figname=os.path.join(folderouttr,filenameout+"_run"+str(ir)+"_t"+str(iTT)+".pdf")
                    figtrj.savefig(figname)
            
            
    fign.subplots_adjust(left=0.1/float(num_runs), bottom=0.25/float(num_tmax), right=1.-0.4/float(num_runs), top=1.-0.2/float(num_tmax),hspace=.1,wspace=.1)#,wspace=5.0, hspace=.1)
    figname=os.path.join(trjstat_output,filenameout+".pdf")
    fign.savefig(figname)
        
    print("Distribution saved to: "+trjstat_output)

    ### get average values
    resultslistav=[]
    resultsav=np.zeros([num_runs,7])
    resultsav[:]=np.nan
    for ir in range(0,num_runs):
        resultslistav.append([])
        for icc in range(0,7):
            resultslistav[-1].append(np.nanmean(resultslist[ir][:,icc]))
            resultsav[ir,icc]=np.nanmean(resultslist[ir][:,icc])

    ####
    #save to csv file
    ####
    listcvv=[]
    listcvv.append([])
    for ir in range(0,num_runs):
            
                listcvv[-1].append(labellistin[ir])
                listcvv[-1].append("")
                listcvv[-1].append("")
                listcvv[-1].append("")
                listcvv[-1].append("")
                listcvv[-1].append("")
                listcvv[-1].append("")
    listcvv.append([])
    for ir in range(0,num_runs):
            
                listcvv[-1].append("time [h]")
                listcvv[-1].append("density [OD]")
                listcvv[-1].append("motile fraction")
                listcvv[-1].append("mean speed [mum/s]")
                listcvv[-1].append("mean speed fast [mum/s]")
                listcvv[-1].append("run time [s]")
                listcvv[-1].append("run time fast [s]")
                listcvv[-1].append("3 distr fit - frac ns")
                listcvv[-1].append("3 distr fit - mean ns [mum/s]")
                listcvv[-1].append("3 distr fit - width ns [mum/s]")
                listcvv[-1].append("3 distr fit - frac sw")
                listcvv[-1].append("3 distr fit - mean sw [mum/s]")
                listcvv[-1].append("3 distr fit - width sw [mum/s]")
                listcvv[-1].append("3 distr fit - frac sl")
                listcvv[-1].append("3 distr fit - mean sl [mum/s]")
                listcvv[-1].append("3 distr fit - width sl [mum/s]")
                listcvv[-1].append("2 distr fit - frac ns")
                listcvv[-1].append("2 distr fit - mean ns [mum/s]")
                listcvv[-1].append("2 distr fit - width ns [mum/s]")
                listcvv[-1].append("2 distr fit - frac s")
                listcvv[-1].append("2 distr fit - mean s [mum/s]")
                listcvv[-1].append("2 distr fit - width s [mum/s]")
              
                    
            
    for iTT in range(0,max(num_t)):
        listcvv.append([])
        for ir in range(0,num_runs):
                    try:
                        listcvv[-1].append(resultslist[ir][iTT,0])
                        listcvv[-1].append(resultslist[ir][iTT,1])
                        listcvv[-1].append(resultslist[ir][iTT,2])
                        listcvv[-1].append(resultslist[ir][iTT,3])
                        listcvv[-1].append(resultslist[ir][iTT,4])
                        listcvv[-1].append(resultslist[ir][iTT,5])
                        listcvv[-1].append(resultslist[ir][iTT,6])
                        listcvv[-1].append(resultslist[ir][iTT,7])
                        listcvv[-1].append(resultslist[ir][iTT,8])
                        listcvv[-1].append(resultslist[ir][iTT,9])
                        listcvv[-1].append(resultslist[ir][iTT,10])
                        listcvv[-1].append(resultslist[ir][iTT,11])
                        listcvv[-1].append(resultslist[ir][iTT,12])
                        listcvv[-1].append(resultslist[ir][iTT,13]) 
                        listcvv[-1].append(resultslist[ir][iTT,14])
                        listcvv[-1].append(resultslist[ir][iTT,15])
                        listcvv[-1].append(resultslist[ir][iTT,16])
                        listcvv[-1].append(resultslist[ir][iTT,17])
                        listcvv[-1].append(resultslist[ir][iTT,18])
                        listcvv[-1].append(resultslist[ir][iTT,19])
                        listcvv[-1].append(resultslist[ir][iTT,20])
                        listcvv[-1].append(resultslist[ir][iTT,21])
                        listcvv[-1].append(resultslist[ir][iTT,22]) 
                    except:
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
                        listcvv[-1].append("")
    

                    
    #####add lines as list format
    entries=["time","density","motilefrac","vmean","vmeanfast","runtime","runtimefast","dis frac 1","dis mean 1","dis width 1","dis frac 2","dis mean 2","dis width 2","dis frac 3","dis mean 3","dis width 3","dis frac 1","dis mean 1","dis width 1","dis frac 2","dis mean 2","dis width 2"]
    for ie in range(0,len(entries    )):
        listcvv.append([])
        for ir in range(0,num_runs):
            print ir
            print ie
            listc=','.join(map(str,resultslist[ir][:,ie]))
            #listc=str(resultslist[ir][:,0])
            listcvv[-1].append(entries[ie]+"=["+listc+"]")
            for icc in range(0,6):
                listcvv[-1].append("")
                      
                      
    #### save array as txt output for plotting
    txtlineout=[]
    for ir in range(0,num_runs):
        strout=labellistin[ir]
        txtlineout.append(strout)
        print txtlineout[-1]
        strout="np.array(["
        for iTT in range(0,max(num_t)):
            #### add output as txt file
    #order: time, od, swimming speed, run time, swimming speed fast, fraction motile
            try:
                curadd=str(resultslist[ir][iTT,[0,1,3,5,4,2]].tolist())
                strout=strout+curadd+","
            except:
                pass
        strout=strout[:-1]
        strout=strout+"])"
        
        txtlineout.append(strout)
        print txtlineout[-1]
        #swimmdatagly3mexpont=np.array([[0.05,	0.038,	18.9269601,	0.56002865],[0.533333333,	0.059,	14.64381348,	0.40611354],[1,	0.086,	15.23267813	,0.45327696],[1.25,	0.106,	16.45361601	,0.44051095],[1.633333333,	0.138,	15.7036132,	0.51632653],[1.966666667,	0.177,	15.64225541	,0.49982394],[2.383333333,	0.243,	18.13289145,	0.67414773],[2.683333333,	0.3,	18.99446723	,0.58994709],[3.333333333,	0.316,	11.54530894	,0.46246418],[4.133333333,	0.315,	11.42053891,	0.70454545],[5.15,	0.315,	7.58828032,	0.65690476],[5.966666667,	0.315,	6.79386399,	0.54921569]])
    nptxtout=os.path.join(trjstat_output,filenameout+"_arrayout.txt")
        
    nptxtoutfile = open(nptxtout, 'w')
    for item in txtlineout:
        print>>nptxtoutfile, item
        
    
   #output array: time, od, swimming speed, run time, swimming speed fast, ration motile
                    
                    
                    
    
    listcvvshort=[] 
    listcvvshort.append([])
    
    listcvvshort[-1].append("run")
    #listcvv[-1].append("time [h]")
    #listcvv[-1].append("density [OD]")
    listcvvshort[-1].append("motile fraction")
    listcvvshort[-1].append("mean speed [mum/s]")
    listcvvshort[-1].append("mean speed fast [mum/s]")
    listcvvshort[-1].append("run time [s]")
    listcvvshort[-1].append("run time fast [s]")
    for ir in range(0,num_runs):
        listcvvshort.append([])
        listcvvshort[-1].append(labellistin[ir])
        
        listcvvshort[-1].append(resultslistav[ir][2])
        listcvvshort[-1].append(resultslistav[ir][3])
        listcvvshort[-1].append(resultslistav[ir][4])
        listcvvshort[-1].append(resultslistav[ir][5])
        listcvvshort[-1].append(resultslistav[ir][6])
                        
    
    
                        
                
    csvname=os.path.join(trjstat_output,filenameout+"_full.csv")
    with open(csvname, "wb") as f:
        writer = csv.writer(f)
        writer.writerows(listcvv)
    print("full output written to:")
    print(csvname)
    
    csvname=os.path.join(trjstat_output,filenameout+".csv")
    with open(csvname, "wb") as f:
        writer = csv.writer(f)
        writer.writerows(listcvvshort)
    print("short output written to:")
    print(csvname)
                

    
def smooth(a,WSZ):
    out0 = np.convolve(a,np.ones(WSZ,dtype=int),'valid')/WSZ    
    r = np.arange(1,WSZ-1,2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return np.concatenate((  start , out0, stop  ))

def analyze_trajectories(filenamein,minduration=3):
    
    filenamein=os.path.basename(filenamein)
    if filenamein[-4:]==".avi":
        filenamein=filenamein[:-4]
    if filenamein[-4:]==".pad":
        filenamein=filenamein[:-4]
    
    
    filenametrj=os.path.join(data_trajectories,filenamein+".pad")
    filenameout2=filenametrjx[:-6]+".npy"
    #load pandas....to numpy file...to have it in same format as the matlab trajectories...
    trjpd = pd.read_pickle(filenametrj) 
    num_trj=trjpd['particle'].nunique() 
    num_frame=trjpd['frame'].nunique()
    
    trj_xdata=np.zeros([num_trj,num_frame])
    trj_xdata[:]=np.nan
    trj_ydata=np.zeros([num_trj,num_frame])
    trj_ydata[:]=np.nan
    for ip in range(0,trjpd['particle'].max()+1):
        dcc=trjpd.loc[trjpd['particle'] == ip]
        
        for index, row in dcc.iterrows():
            #print row["frame"]
            frc=int(row["frame"]-1)
            trj_xdata[ip,frc]=row["x"]
            trj_ydata[ip,frc]=row["y"]
    
    Traj=[]
    #0: time
            #1: x
            #2: y
            #3: real time
            #4: velocity
            #5: dv
            #6: theta
            #7: dtheta
            #8:  tumbling evet 
    #Analysis of trajectories...
    #Traj has length of number of trajectories (sturcutred array).....
    #every entry has a matrix....time, xy before it was (-,-,time, x, y,-)
    
    
    #setting parameters of trajectory analysis...
    
    ####
    #length calibration    
    ####    
    muparpix = 0.5123 #calibration when camer on the right
    dt=0.05# %to calibrate time, frame rate 20Hz
    a=2 
    b=2
    for j in range(0,trj_xdata.shape[0]):
        indx=np.isfinite(trj_xdata[j,:]) #1:length(Matrix_x(:,1))
        curt=np.count_nonzero(~np.isnan(trj_xdata[j,:]))

        curnp=np.zeros([curt,12])
        curnp[:,1]=trj_xdata[j,indx]
        curnp[:,2]=trj_ydata[j,indx]
        curnp[:,0]=np.arange(0,curt,1)*dt
        curnp[:,3]=np.array(np.nonzero(indx[:]))[0,:] #time
        curnp[:,4]=np.nan #dr
        curnp[:,5]=np.nan #velocity
        curnp[:,6]=np.nan #dv
        curnp[:,7]=np.nan #theta event....
        curnp[:,8]=np.nan #dtheta
        curnp[:,9]=np.nan #tubmling event....
        curnp[:,10]=np.nan #tubmling event rotation....
        curnp[:,11]=np.nan #tubmling event speed....
        
        if curt<minduration:
            pass
        else:
            Traj.append(curnp)
    v_tot=[]
    dr_dtot=[]
    runsTot=[]
    
    #add these numbers to traje list
    Nrtotal_trjlist=[]# %length of trajectory (with both swimming and tumbling)
    Nrrun_trjlist=[]# %length of trajectories, with only running
    Nrtumbling_trjlist=[]#%length of trajectories with only tumbling
    NewTraj=[]
    
    
    for i in range(0, len(Traj)):#go through every trajectory....
        
        if Traj[i].shape[0]>10: #only consder trajectories > 0.5 seconds...
    
            #smoth trajectories...(averaging window)
            xNew=smooth(Traj[i][:,1]*muparpix,3)##m(:,3)*muparpix,3); %muparpix is calibrating length...
            yNew=smooth(Traj[i][:,2]*muparpix,3)#smooth(Traj(i).m(:,4)*muparpix,3);
            f=Traj[i][:,0]*dt #f=Traj(i).m(:,2)*dt; %dt is calibrating time
            num_x=xNew.shape[0]
            num_y=yNew.shape[0]
            #calculate dx and dy for whole trajectory....
            dx=xNew[a:num_x]-xNew[:num_x-a] #calculate dx over several timesteps, parameter a)
            dy=yNew[a:num_y]-yNew[:num_x-a]
    
            inda=np.array(range(int(a/2.),num_x-int(a/2.)))
            #calculate change in angle...
            theta=np.zeros([dx.shape[0]]) #prepare vector for angle variation....
            for j in range(0,theta.shape[0]):
                if dx[j]<0:
                    
                    theta[j]=np.pi+np.arctan(dy[j]/dx[j]);
                elif dx[j]==0:
                    if dy[j]==0:
                        theta[j]=np.pi/2;
                    else:
                        theta[j]=3*np.pi/2;
                elif dx[j]>0:
                    if dy[j]>0:
                        theta[j]=np.arctan(dy[j]/dx[j]);
                    elif dy[j]<0:
                        theta[j]=2*np.pi+np.arctan(dy[j]/dx[j]);
                    else:
                        theta[j]=0;
            indb=np.array(inda[int(b/2.):-int(b/2.)])
            
            #adjust theta if larger than pi
            dtheta=theta[b:theta.shape[0]]-theta[0:theta.shape[0]-b];
            idx1=np.where( dtheta > np.pi )
            idx2=np.where( dtheta < -np.pi )
            dtheta[idx1]=dtheta[idx1]-np.pi
            dtheta[idx2]=dtheta[idx2]+np.pi
    
            dr=np.sqrt(dx*dx+dy*dy)# %calculate dr
            v=dr/(dt*a)# %calculate v....this is taking several timesteps into account (parameter a)
            dv=v[b+0:v.shape[0]]-v[0:v.shape[0]-b]
            vtest=dv/v[0:dv.shape[0]]#; %normalized change in velocity...
            #find tumbling events, defined by:
            #1. change in angle (dtheta larger than traghold)
            #2. change in speed, relative change smaller than 0.5)
    
            #condition 1 (if change is larger than a certain angle)
            tumbleAngle=indb[np.where(np.absolute(dtheta)>=1.2)]#; %jeromes other skript had here a tvalue defined 
            
            #condition 2 (if relative change is smaller)
            tumbleSpeed=indb[np.where(vtest<-.5)]#; %jeromes other skript had treshold defined...
    
            tumble=np.unique(np.append(tumbleAngle,tumbleSpeed))#logical_or(tumbleAngle,tumbleSpeed)
    
            Traj[i][inda,4]=dr/float(a)
            Traj[i][inda,5]=v
            Traj[i][indb,6]=dv
            Traj[i][inda,7]=theta
            Traj[i][indb,8]=dtheta
            Traj[i][tumble,9]=1
            Traj[i][tumbleAngle,10]=1
            Traj[i][tumbleSpeed,11]=1
            #0: time
            #1: x
            #2: y
            #3: real time
            #4: dr            
            #5: velocity
            #6: dv
            #7: theta
            #8: dtheta
            #9:  tumbling evet 
            
            if np.std(dtheta)<1.2 :
            #take only trajectories which are not changing angle too much.
                NewTraj.append(Traj[i])
                
                #difference between tumble events...
                dtumble=np.diff(f[tumble])# %time difference between tumbling events (in dt units)
                Nrtotal=xNew.shape[0]-a #get length of trajectory with both, swimming and tumbling
                Nr=xNew.shape[0]-a-tumble.shape[0]#%get length of trajectory with no tumblin
                Nrt=tumble.shape[0]#;  %number tumbling events....
                Nrtotal_trjlist.append(Nrtotal)
                Nrrun_trjlist.append(Nr)
                Nrtumbling_trjlist.append(Nrt)
    
    
                if tumble.shape[0]>0:
                # %in case there are tumbling events to analyze...
                    runsTot.append(dtumble)# %add times of tumble...calibration of time
                    v_tot.append(v)#[v_tot; v];
                    dr_dtot.append(dr)#     dr_tot=[dr_tot; dr];
                    
    
    #take trajectories which have been analyzed ...
    
    #calculate average values....\
    Nrrun_trjlist=np.array(Nrrun_trjlist)
    Nrtumbling_trjlist=np.array(Nrtumbling_trjlist)
    v_tot=np.array(v_tot)
    
    #save trajectory analysis....each trajectory as single array..... 
    outarray=np.zeros([len(Traj),trj_xdata.shape[1],Traj[0].shape[1]])
    outarray[:]=np.nan
    
    for i in range(0, len(Traj)):
        outarray[i,Traj[i][:,3].astype(int)[:],:]=Traj[i]
    np.save(filenameout2,outarray)
    

def get_avifiles(foldernamedatain):
    
    try: #try hdd first...
        if len(foldernamedatain)>0:
            foldernamedatainfullpath=os.path.join(datafolderrawvideo2,foldernamedatain)
        else:
            foldernamedatainfullpath=datafolderrawvideo2
        listdir=os.listdir(foldernamedatainfullpath)
    except:
         if len(foldernamedatain)>0:
            foldernamedatainfullpath=os.path.join(datafolderrawvideo,foldernamedatain)
         else:
            foldernamedatainfullpath=datafolderrawvideo
         listdir=os.listdir(foldernamedatainfullpath)
    for il in range(0,len(listdir)):
        if listdir[il][-4:]==".avi":
            print '#filenamein="'+os.path.join(foldernamedatain,listdir[il])+'"'
    namepure=[]
    for il in range(0,len(listdir)):
            if listdir[il][-4:]==".avi":
                namepure.append(listdir[il][:-4])
    return namepure
     
       
def generate_trajectories(filenamein,framemaxin=500, search_range=10, memory=5):
    filenamein=os.path.basename(filenamein)
    if filenamein[-4:]==".avi":
        filenamein=filenamein[:-4]
    filenamedetection=os.path.join(data_detection,filenamein+".npz")
    try: 
        detected=np.load(filenamedetection)
        
        numc=len(detected.keys())
        if framemaxin>numc:
            framemaxin=numc
            print "Less analyzed frames than framemax. Used frames: "+str(framemaxin)
        
        c1=detected["frame"+str(1)][:,:3]
        c1[:,2]=1
        frame=c1[:,2].astype(int)
        mass=(c1[:,2]*0+1).astype(int)
        c1f=np.ones(c1.shape[0]).astype(int)
        c1m=np.ones(c1.shape[0]).astype(int)
        
        for il in range(2,framemaxin):
            c2=detected["frame"+str(il)][:,:3]
            c2[:,2]=il
            #x=c2[:,0]
            #y=c2[:,1]
            frame=c2[:,2].astype(int)
            mass=(c2[:,2]*0+1).astype(int)
            #mass=c1[:,2].astype(int)
            
            c1=np.append(c1,c2,axis=0)
            c1f=np.append(c1f,frame,axis=0)
            c1m=np.append(c1m,mass,axis=0)
        
        numdet=c1.shape[0]/float(framemaxin)
        print "number cells detected: "+str(numdet)
        
        if numdet>10000:
            print "too many detected cells"
        else:
            
            dataframe = pd.DataFrame({'x':c1[:,0], 'y':c1[:,1], 'frame':c1f, 'mass':c1m, 'size':c1m, 'ecc': c1m, 'raw_mass': c1m, 'ep':c1m})
            dataframeout=trackpy.link_df(dataframe, search_range=search_range, memory=memory)#, neighbor_strategy='KDTree', link_strategy='auto', predictor=None, adaptive_stop=None, adaptive_step=0.95, copy_features=False, diagnostics=False, pos_columns=None, t_column=None, hash_size=None, box_size=None, verify_integrity=True, retain_index=False)
            num_trj=dataframeout['particle'].nunique()    
            print("trajectories found: "+str(num_trj))
            filenameout=os.path.join(data_trajectories,filenamein+".pad")
            dataframeout.to_pickle(filenameout)  # where to save it, usually as a .pkl
            print "trajectories save to: "+filenameout
            
    except:
        print "!!!!!could not open "+filenamedetection

#go through whole oflder, if only one file....start analyze_movie directly
def run_detection(foldernamedatain,framemaxin=2000,plot_detectionanalysis=True,num_avsteps=10,short=False,selectionfraction=0.1, selectionfraction2=0.15,blur_parameter=20,min_value_selectionbackground=0,strcriteria="",fluorescencedata=False):
    
    if len(foldernamedatain)>0:
        foldernamedatainfullpath=os.path.join(datafolderrawvideo2,foldernamedatain)
    else:
        foldernamedatainfullpath=datafolderrawvideo

    filenmaelist=[]
    #go through filenmaes.....
    listdir=os.listdir(foldernamedatainfullpath)
    for il in range(0,len(listdir)):
        if listdir[il][-4:]==".avi":
            if strcriteria=="":
                filenmaelist.append(listdir[il])
            else:
                if strcriteria in listdir[il]:
                    filenmaelist.append(listdir[il])
    for il in range(0,len(filenmaelist)):
        filename=os.path.join(foldernamedatain,filenmaelist[il])
        analyze_movie(filename,framemaxin=framemaxin,plot_detectionanalysis=plot_detectionanalysis,num_avsteps=num_avsteps,short=short,selectionfraction=selectionfraction, selectionfraction2=selectionfraction2,blur_parameter=blur_parameter,min_value_selectionbackground=min_value_selectionbackground,fluorescencedata=fluorescencedata)

    
    
