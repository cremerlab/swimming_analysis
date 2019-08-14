% This file is used to run the u-track package. A similar file is provided by the u-track package. We adjusted it to integrate it into our
% file reading and cell detection routine.
%
% u-track is available: https://github.com/DanuserLab/u-track
% 
% 
% Licence of the u-track package asfollows: u-track is free software: 
% you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

function track_function(filenamein,dirnameout,filenameout,maxSearchRadiusin,maxSearchRadiusingap)
disp(filenamein)
load(filenamein)
%disp(strcat('Running tracking algorithm in matlab...',filenamein))

movieInfo = cell2struct(movieInfoCell, {'xCoord','yCoord','amp'}, 2);
save movieInfo.mat
load movieInfo.mat
%filenameout='outm'
%maxSearchRadiusin=30
%maxSearchRadiusingap=30
%dirnameout='~/Desktop'; %directory where to save input and output
dirnemeoutc=dirnameout;
filenameoutc='tracksTest4.mat'; %name of file where input and output are saved

strcat(filenameout,'x.m')

filenameoutMobileX=fullfile(dirnameout,strcat(filenameout,'RatioOfMobileX.mat'));
filenameoutMobileY=fullfile(dirnameout,strcat(filenameout,'RatioOfMobileY.mat'));
filenameoutMobileR=fullfile(dirnameout,strcat(filenameout,'RatioOfMobileR.mat'));
filenameoutNumberMobileR=fullfile(dirnameout,strcat(filenameout,'NumberOfMobileR.mat'));
filenameoutNumberNonMobileR=fullfile(dirnameout,strcat(filenameout,'NumberOfNonMobileR.mat'));
filenameoutcx=fullfile(dirnameout,strcat(filenameout,'_x.mat'));
filenameoutcy=fullfile(dirnameout,strcat(filenameout,'_y.mat'));

filenameoutcxmat=fullfile(dirnameout,strcat(filenameout,'_x.mat'));
filenameoutcymat=fullfile(dirnameout,strcat(filenameout,'_y.mat'));

filenameoutcfig1=fullfile(dirnameout,strcat(filenameout,'_Distribution.fig'));
filenameoutcfig2=fullfile(dirnameout,strcat(filenameout,'_SingleTraj.fig'));
filenameoutcfig3=fullfile(dirnameout,strcat(filenameout,'_Scatter.fig'));
filenameoutcfig4=fullfile(dirnameout,strcat(filenameout,'_ScatterDistribution.fig'));
filenameoutcfig5=fullfile(dirnameout,strcat(filenameout,'_TrajDuration.fig'));

gapCloseParam.timeWindow = 3; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
gapCloseParam.mergeSplit = 1; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
gapCloseParam.minTrackLen = 1; %minimum length of track segments from linking to be used in gap closing: reduce isolated gap points

%optional input:
gapCloseParam.diagnostics = 1; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.

%% cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';

%parameters

parameters.linearMotion = 0; %use linear motion Kalman filter.
parameters.minSearchRadius = 0; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
parameters.maxSearchRadius = maxSearchRadiusin; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.

parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).

parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
% parameters.kalmanInitParam.searchRadiusFirstIteration = 10; %Kalman filter initialization parameters.

%optional input
parameters.diagnostics = []; %if you want to plot the histogram of linking distances up to certain frames, indicate their numbers; 0 or empty otherwise. Does not work for the first or last frame of a movie.

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing

%function name
costMatrices(2).funcName = 'costMatRandomDirectedSwitchingMotionCloseGaps';

%parameters

%needed all the time
parameters.linearMotion = 0; %use linear motion Kalman filter.

parameters.minSearchRadius = 0; %minimum allowed search radius.
parameters.maxSearchRadius = maxSearchRadiusingap; %maximum allowed search radius.
parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.

parameters.brownScaling = [0 0.01]; %power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
% parameters.timeReachConfB = 3; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).
parameters.timeReachConfB = gapCloseParam.timeWindow; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).

parameters.ampRatioLimit = [0.7 4]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.

parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.

parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).

parameters.linStdMult = 1*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.

parameters.linScaling = [0.25 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
% parameters.timeReachConfL = 4; %similar to timeReachConfB, but for the linear part of the motion.
parameters.timeReachConfL = gapCloseParam.timeWindow; %similar to timeReachConfB, but for the linear part of the motion.

parameters.maxAngleVV = 30; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.

%optional; if not input, 1 will be used (i.e. no penalty)
parameters.gapPenalty = 1.5; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).

%optional; to calculate MS search radius
%if not input, MS search radius will be the same as gap closing search radius
parameters.resLimit = []; %resolution limit, which is generally equal to 3 * point spread function sigma.

costMatrices(2).parameters = parameters;
clear parameters

%% Kalman filter function names

kalmanFunctions.reserveMem  = 'kalmanResMemLM';
kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

%% additional input

%saveResults
saveResults.dir = dirnemeoutc

saveResults.filename = filenameoutc
% saveResults = 0; %don't save results

%verbose state
verbose = 1;

%problem dimension
probDim = 2;

%% tracking function call

[tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,...
    costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);

% for i = 1 : 12
%     movieInfoTmp((i-1)*1200+1:i*1200) = movieInfo((i-1)*1200+1:i*1200);
%     saveResults.filename = ['tracks1All_' sprintf('%02i',i) '.mat'];
%     [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfoTmp,...
%         costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
%     clear movieInfoTmp
% end
if isempty(tracksFinal)
    exit
end
%% ~~~ the end ~~~

[trackedFeatureInfo, trackedFeatureIndx, trackStartRow, numSegments] = convStruct2MatIgnoreMS(tracksFinal);

row=size(trackedFeatureInfo,1);
col=size(trackedFeatureInfo,2);

x=NaN(row,col/8);
y=NaN(row,col/8);

for i=1:row
for j=1:col/8
    if ~isnan(trackedFeatureInfo(i,1+8*(j-1)))
        y(i,j)=trackedFeatureInfo(i,1+8*(j-1));
    end
    if ~isnan(trackedFeatureInfo(i,2+8*(j-1)))
        x(i,j)=trackedFeatureInfo(i,2+8*(j-1));
    end
end
end
% save x.mat
% save y.mat
% save('x.mat',filenameoutcx);
% save('y.mat',filenameoutcy);

save(filenameoutcx,'x');
save(filenameoutcy,'y');
save(filenameoutcxmat,'x');
save(filenameoutcymat,'y');



%% Distribution of mean displacement and mean square for a movie


PositionX=zeros(size(x,1),1);
PositionY=zeros(size(y,1),1);

LengthX=zeros(size(x,1),1);
LengthY=zeros(size(y,1),1);
DiffusionX=zeros(size(x,1),1);
DiffusionY=zeros(size(y,1),1);
Xn=zeros(size(x,1),1);
Yn=zeros(size(y,1),1);
X0=zeros(size(x,1),1);
Y0=zeros(size(y,1),1);

for i=1:size(x,1)
    for j=1:size(x,2)
        if ~isnan(x(i,j)) 
            X0(i,1)=x(i,j);
            break
        end
    end
    for j=1:size(x,2)
        if  ~isnan(y(i,j))
            Y0(i,1)=y(i,j);
            break
        end
    end
end

for j=1:size(x,2)
    for i=1:size(x,1)
        if ~isnan(x(i,j))  && ~isnan(y(i,j)) %&& ~isnan(x(i,j+1)) && ~isnan(y(i,j+1))
        %B=plot(y(i,j),x(i,j),'b.');hold on;
        %axis([0 512 0 512])
        %hold on;
        Xn(i,1)=Xn(i,1)+1;Yn(i,1)=Yn(i,1)+1;
        PositionX(i,1)=PositionX(i,1)+(x(i,j)-X0(i,1));
        PositionY(i,1)=PositionY(i,1)+(y(i,j)-Y0(i,1));
        DiffusionX(i,1)=DiffusionX(i,1)+(x(i,j)-X0(i,1))^2;
        DiffusionY(i,1)=DiffusionY(i,1)+(y(i,j)-Y0(i,1))^2;
        
        
        LengthX(i,1)=LengthX(i,1)+sqrt((x(i,j)-X0(i,1))^2);
        end
    end
    %drawnow
    %mov(j) = getframe(gcf);
end
PositionX=PositionX./Xn;
PositionY=PositionY./Yn;
DiffusionX=DiffusionX./Xn;
DiffusionY=DiffusionY./Yn;

% figure 
% x1=sort(PositionX);y1=sort(PositionY);x1(x1==0)=[];y1(y1==0)=[];
% x2=sort(DiffusionX);y2=sort(DiffusionY);x2(x2==0)=[];y2(y2==0)=[];
% [a1,b1]=hist(x1,50); a1=a1/length(x1); subplot(2,2,1);bar(b1,a1);title('\langle(x-x_0)\rangle');xlabel('pixels');
% [a2,b2]=hist(x2,50); a2=a2/length(x2); subplot(2,2,2);bar(b2,a2);title('\langle(x-x_0)^2\rangle');xlabel('pixels^2');
% [c1,d1]=hist(y1,50); c1=c1/length(y1); subplot(2,2,3);bar(d1,c1);title('\langle(y-y_0)\rangle');xlabel('pixels');
% [c2,d2]=hist(y2,50); c2=c2/length(y2); subplot(2,2,4);bar(d2,c2);title('\langle(y-y_0)^2\rangle');xlabel('pixels^2');

% savefig(filenameoutcfig1);







%% Statistics for single trajectory
SingleTraj=0; %save image 1 or not 0
I=1; %denote how many trajectories wanted to look into

XnT=zeros(size(x,1),1);
YnT=zeros(size(y,1),1);

for i=1:I %size(x,1)  
    for j=1:size(x,2)
        if ~isnan(x(i,j))  && ~isnan(y(i,j)) %&& ~isnan(x(i,j+1)) && ~isnan(y(i,j+1))
        XnT(i,1)=XnT(i,1)+1;YnT(i,1)=YnT(i,1)+1;
        end
    end 
    DriftXT=zeros(1,XnT(i,1));
    DriftYT=zeros(1,YnT(i,1));
    PositionXT=zeros(1,XnT(i,1));
    PositionYT=zeros(1,YnT(i,1));
    %PositionXT=zeros(1,XnT(i,1));
    DiffusionXT=zeros(1,XnT(i,1));
    DiffusionYT=zeros(1,YnT(i,1));
    XnTline=zeros(1,XnT(i,1));
    YnTline=zeros(1,YnT(i,1));
    XnT(i,1)=0;
    YnT(i,1)=0;
    count=0;
    for j=1:size(x,2)
        XnT(i,1)=XnT(i,1)+1;YnT(i,1)=YnT(i,1)+1;
        if ~isnan(x(i,j))  && ~isnan(y(i,j)) %&& ~isnan(x(i,j+1)) && ~isnan(y(i,j+1))
        count=count+1;%XnT(i,1)=XnT(i,1)+1;YnT(i,1)=YnT(i,1)+1;
        XnTline(1,count)=XnT(i,1);YnTline(1,count)=YnT(i,1);
        PositionXT(1,count)=x(i,XnT(i,1));
        PositionYT(1,count)=y(i,XnT(i,1));
        
        if count==1
        DriftXT(1,count)=(x(i,XnT(i,1))-X0(i,1));
        DriftYT(1,count)=(y(i,YnT(i,1))-Y0(i,1));
        XPrevious=x(i,XnT(i,1));YPrevious=y(i,YnT(i,1));
        DiffusionXT(1,count)=(x(i,XnT(i,1))-X0(i,1))^2;
        DiffusionYT(1,count)=(y(i,YnT(i,1))-Y0(i,1))^2;
        end
        if count~=1
        DriftXT(1,count)=(x(i,XnT(i,1))-XPrevious);
        DriftYT(1,count)=(y(i,YnT(i,1))-YPrevious);
        XPrevious=x(i,XnT(i,1));YPrevious=y(i,YnT(i,1));
        DiffusionXT(1,count)=DiffusionXT(1,count-1)+(x(i,XnT(i,1))-X0(i,1))^2;
        DiffusionYT(1,count)=DiffusionYT(1,count-1)+(y(i,YnT(i,1))-Y0(i,1))^2;
        end
        end
    end
%     figure
% %subplot(3,2,1);plot(XnTline(1,:),PositionXT(1,:));title('x-t');xlabel('t');hold on;
% %subplot(3,2,2);plot(YnTline(1,:),PositionYT(1,:));title('y-t');xlabel('t');hold on;
% subplot(2,2,1);plot(XnTline(1,:),DriftXT(1,:));title('\Deltax-t');xlabel('t');hold on;
% subplot(2,2,2);plot(YnTline(1,:),DriftYT(1,:));title('\Deltay-t');xlabel('t');hold on;
% subplot(2,2,3);plot(XnTline(1,:),DiffusionXT(1,:));title('\Sigma_n^N(x_n-x_0)^2-t');xlabel('t');hold on;
% subplot(2,2,4);plot(YnTline(1,:),DiffusionYT(1,:));title('\Sigma_n^N(y_n-y_0)^2-t');xlabel('t');hold on;
end

% if SingleTraj==1
% savefig(filenameoutcfig2);
% end
%% Scatter
XN=zeros(size(x,1),1);
YN=zeros(size(y,1),1);
TrajDuration=zeros(size(x,1),1);
%YN=zeros(size(y,1),1);
R=zeros(size(x,1),1);
Rmax=0;
for i=1:size(x,1)

    for j=1:size(x,2)
        if ~isnan(x(i,j))  && ~isnan(y(i,j)) %&& ~isnan(x(i,j+1)) && ~isnan(y(i,j+1))
        XN(i,1)=x(i,j)-X0(i,1);
        YN(i,1)=y(i,j)-Y0(i,1);
        TrajDuration(i,1)=TrajDuration(i,1)+1;
        if Rmax<XN(i,1) 
            Rmax=XN(i,1);
        end
        if Rmax<YN(i,1)
            Rmax=YN(i,1);
        end
        end
    end
    R(i,1)=sqrt(XN(i,1)^2+YN(i,1)^2);
    if Rmax<R(i,1) 
            Rmax=R(i,1);
    end
end
DriftTotalX=sum(XN)/length(XN);
DriftTotalY=sum(YN)/length(YN);
%DriftTotalR=sum(R)/length(R);
% figure
% plot(XN,YN,'.');hold on;
% r=Rmax; 
% n=10;
% a=0;
% b=0;
% t=linspace(-pi,pi); 
% Rx=sin(t)'*linspace(0,r,n+1)+a;
% Ry=cos(t)'*linspace(0,r,n+1)+b;
% plot(Rx,Ry);title('Scatter Plot');xlabel('x in pixels');ylabel('y in pixels');text(Rmax-Rmax/4,Rmax-Rmax/4,['\Deltax=',num2str(DriftTotalX),'\newline\Deltay=',num2str(DriftTotalY),'\newlineTraj#=',num2str(size(x,1))])%,'\newline\DeltaR=',num2str(DriftTotalR)])
% axis equal
%savefig(filenameoutcfig3);


% figure
% RR=sort(R);[r1,R1]=hist(RR,25);r1=r1/length(R); bar(R1,r1);title('Distribution of R');xlabel('R in pixels');
% % savefig(filenameoutcfig4);

% figure
% TrajDurationS=sort(TrajDuration);[traj1,Traj1]=hist(TrajDurationS,100);traj1=traj1/length(TrajDuration); bar(Traj1,traj1);title('Distribution of TrajDuration');xlabel('frames');
% %savefig(filenameoutcfig5);

%% Ratio of Mobile and non-Mobile
I=4;%Set the threshold of non-Mobile: for hundreds of traj, Guanssian fit is suitble for deleting inner 4 displancements.
XN=sqrt(XN.^2);YN=sqrt(YN.^2);
XS=sort(XN);[xs1,XS1]=hist(XS,25);xs1=xs1/length(XS);
YS=sort(YN);[ys1,YS1]=hist(YS,25);ys1=ys1/length(YS);
RS=sort(R);[Rs1,RS1]=hist(RS,25);Rs1=Rs1/length(RS);
RatioMobile=zeros(3,1);%NumberMobile=zeros(1,1);
for i=1:I  % Some problem on definition of Mobile
    RatioMobile(1,1)=RatioMobile(1,1)+xs1(i);
    RatioMobile(2,1)=RatioMobile(2,1)+ys1(i);
    RatioMobile(3,1)=RatioMobile(3,1)+Rs1(i);
end
%RatioNon(1,1)=1-RatioMobile(1,1);
%RatioNon(2,1)=1-RatioMobile(2,1);

RatioMobileX=RatioMobile(1,1);
RatioMobileY=RatioMobile(2,1);
RatioMobileR=RatioMobile(3,1);
NumberMobileR=RatioMobile(3,1)*length(RS);
NumberNonMobileR=length(RS)-NumberMobileR;
% save(filenameoutMobile,'RatioMobile');
% save(filenameoutMobileX,'RatioMobileX');
% save(filenameoutMobileY,'RatioMobileY');
% save(filenameoutMobileR,'RatioMobileR');
% save(filenameoutNumberMobileR,'NumberMobileR');
% save(filenameoutNumberNonMobileR,'NumberNonMobileR');
end

