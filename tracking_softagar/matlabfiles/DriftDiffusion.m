function DriftDiffusion(filex,filey,dirnameout,filefolder,Sc,Cc,Tc)

%disp(filex)
load(filex)
%disp(filey)
load(filey)
PlotDetail=0;%To plot the distribution of statistics from single cells
PlotDriftDiffusionFit=1;%To plot the linear fit between mean square displacement and time.
ThrowNonMotile=1;%Delete non-motile cells to get correlation: set creteria below


PhysicalLength=192.36;%length size of images in the unit of mum
pixel=512;%pixel number
tLength=0.068; %time step between images in the unit of second
num=20;% binsize for the histogram
TimeWindow=10;%choose traj last onger than TimeWindow frames only
Howlongtolook=10;%Number of consecutive images for making the linear plot between mean square displacement and time.
CoObservationLength=50;%Number of consecutive images for calculating the auto-correlation function

xAll=x;% positions of cells in the x-direction
yAll=y;% positions of cells in the y-direction


xjudge=x;
yjudge=y;
xjudge(isnan(xjudge)==1)=[];
yjudge(isnan(yjudge)==1)=[];
if size(xjudge,2)==1
    return;
end   
if size(yjudge,2)==1
    return;
end
diffusionout=[dirnameout,'\Diffusion',filefolder]; 
dirnameout=[dirnameout,'\',filefolder];  
if exist(dirnameout)==0
    mkdir(dirnameout);
end

%%
%Pick traj with length longer than TimeWindow of images
NumWindow=1;
DiffusionTotal=zeros(NumWindow,NumWindow*10);
DriftXTotal=zeros(NumWindow,NumWindow*10);
DriftYTotal=zeros(NumWindow,NumWindow*10);
TrajLengthX=[];
TrajLengthY=[];
for window=1:NumWindow
changeWindow=1;
if changeWindow==1
NumLength=window*TimeWindow;
xnew=0;ynew=0;
for i=1:size(x,1)
    %xtemp=x(i,:);
    TempLength=0;TempPosition=0;
    for j=1:size(x,2)
        if ~isnan(x(i,j))
            if j>1 && ~isnan(x(i,j-1))
            TempLength=TempLength+1;
            else
            TempLength=0;TempPosition=j;    
            end
            if TempLength>NumLength-1 && ((j<size(x,2) && isnan(x(i,j+1))) || j==size(x,2))
                xsave=NaN(1,size(x,2));
                xsave(1,TempPosition:TempPosition+TempLength)=x(i,TempPosition:TempPosition+TempLength);
               
                if xnew==0
                xnew=xsave(1,:);
                TrajLengthX=xsave(1,TempPosition:TempPosition+9);
                else
                xnew=[xnew',xsave(1,:)']';
                TrajLengthX=[TrajLengthX',xsave(1,TempPosition:TempPosition+9)']';
                end
            end
        end
    end
end
x=xnew;
for i=1:size(y,1)
    %xtemp=x(i,:);
    TempLength=0;TempPosition=0;
    for j=1:size(y,2)
        if ~isnan(y(i,j))
            if j>1 && ~isnan(y(i,j-1))
            TempLength=TempLength+1;
            else
            TempLength=0;TempPosition=j;    
            end
            if TempLength>NumLength-1 && ((j<size(y,2) && isnan(y(i,j+1))) || j==size(y,2))
                ysave=NaN(1,size(y,2));
                ysave(1,TempPosition:TempPosition+TempLength)=y(i,TempPosition:TempPosition+TempLength);
                if ynew==0
                ynew=ysave(1,:);
                TrajLengthY=ysave(1,TempPosition:TempPosition+9);
                else
                ynew=[ynew',ysave(1,:)']';
                TrajLengthY=[TrajLengthY',ysave(1,TempPosition:TempPosition+9)']';
                end
            end
        end
    end
end
y=ynew;
end
if CoObservationLength>size(x,2)-1
    CoObservationLength=size(x,2)-1;
end

%%
%use ball filter to remove the non-motile cells
if ThrowNonMotile==1 %Ball filter
TrajLength=zeros(size(x,1),1);
LengthThreshold=1*pixel/PhysicalLength;%mum The larger, the remaining path is less 
howlong=10;%The larger, the remaining path is more. set stay in a region for continuous 10 times


Xstart=zeros(size(x,1),1);
Ystart=zeros(size(y,1),1);
for i=1:size(x,1)
    for j=1:size(x,2)
        if ~isnan(x(i,j)) && ~isnan(y(i,j))
            Xstart(i,1)=x(i,j);Ystart(i,1)=y(i,j);
            break
        end
    end
end
count=0;xnew=0;ynew=0;xThrow=0;yThrow=0;
for i=1:size(x,1)
    times=0;continuous=-1;out=1;countStayIn=0;%displacementX=0;displacementY=0;
    for j=1:size(x,2)-1 %displacment
        displacement=Inf(1);
        if ~isnan(x(i,j)) && ~isnan(y(i,j)) &&~isnan(x(i,j+1)) && ~isnan(y(i,j+1)) 
            displacement=sqrt((x(i,j+1)-x(i,j))^2+(y(i,j+1)-y(i,j))^2);
             displacementX=sqrt((x(i,j+1)-x(i,j))^2);
             displacementY=sqrt((y(i,j+1)-y(i,j))^2);
            times=times+1;
        end
        %if (displacement<LengthThreshold || (displacementY<LengthThreshold && displacementX<LengthThreshold))  && displacement~=0
        if displacement<LengthThreshold && displacement~=Inf(1)
                if continuous==-1 
                    continuous=times;
                    countStayIn=1;
                    if countStayIn==howlong
                          out=0;break
                    end
                else
                    if continuous==times-1
                         continuous=continuous+1;
                         countStayIn=countStayIn+1;        
                         if countStayIn==howlong
                            out=0;break
                         end
                    else
                        continuous=times;countStayIn=1;
                        if countStayIn==howlong
                          out=0;break
                        end
                    end
                end
        end 

    end
    if out==1
        if xnew==0
                xnew=x(i,:);
                ynew=y(i,:);
            else
                xnew=[xnew',x(i,:)']';
                ynew=[ynew',y(i,:)']';
        end        
    else
        count=count+1; 
        if xThrow==0
                xThrow=x(i,:);
                yThrow=y(i,:);
            else
                xThrow=[xThrow',x(i,:)']';
                yThrow=[yThrow',y(i,:)']';
        end    
    end
end
fractionDrop=count/size(x,1);
x2=xnew;y2=ynew;

%% 
%To align the trajectories starting from different frames
TrajLengthX2=[];
TrajLengthY2=[];
if changeWindow==1 
NumLength=window*TimeWindow;
xnew=0;ynew=0;
for i=1:size(x2,1)
    %xtemp=x(i,:);
    TempLength=0;TempPosition=0;
    for j=1:size(x2,2)
        if ~isnan(x2(i,j))
            if j>1 && ~isnan(x2(i,j-1))
            TempLength=TempLength+1;
            else
            TempLength=0;TempPosition=j;    
            end
            if TempLength>NumLength-1 && ((j<size(x2,2) && isnan(x2(i,j+1))) || j==size(x2,2))
                xsave=NaN(1,size(x2,2));
                xsave(1,TempPosition:TempPosition+TempLength)=x2(i,TempPosition:TempPosition+TempLength);
               
                if xnew==0
                xnew=xsave(1,:);
                TrajLengthX2=xsave(1,TempPosition:TempPosition+9);
                else
                xnew=[xnew',xsave(1,:)']';
                TrajLengthX2=[TrajLengthX2',xsave(1,TempPosition:TempPosition+9)']';
                end
            end
        end
    end
end
x2=xnew;
for i=1:size(y2,1)
    %xtemp=x(i,:);
    TempLength=0;TempPosition=0;
    for j=1:size(y2,2)
        if ~isnan(y2(i,j))
            if j>1 && ~isnan(y2(i,j-1))
            TempLength=TempLength+1;
            else
            TempLength=0;TempPosition=j;    
            end
            if TempLength>NumLength-1 && ((j<size(y2,2) && isnan(y2(i,j+1))) || j==size(y2,2))
                ysave=NaN(1,size(y2,2));
                ysave(1,TempPosition:TempPosition+TempLength)=y2(i,TempPosition:TempPosition+TempLength);
                if ynew==0
                ynew=ysave(1,:);
                TrajLengthY2=ysave(1,TempPosition:TempPosition+9);
                else
                ynew=[ynew',ysave(1,:)']';
                TrajLengthY2=[TrajLengthY2',ysave(1,TempPosition:TempPosition+9)']';
                end
            end
        end
    end
end
y2=ynew;


TrajLengthX3=[];
TrajLengthY3=[];
xThrownew=0;yThrownew=0;
for i=1:size(xThrow,1)
    %xtemp=x(i,:);
    TempLength=0;TempPosition=0;
    for j=1:size(xThrow,2)
        if ~isnan(xThrow(i,j))
            if j>1 && ~isnan(xThrow(i,j-1))
            TempLength=TempLength+1;
            else
            TempLength=0;TempPosition=j;    
            end
            if TempLength>NumLength-1 && ((j<size(xThrow,2) && isnan(xThrow(i,j+1))) || j==size(xThrow,2))
                xsave=NaN(1,size(xThrow,2));
                xsave(1,TempPosition:TempPosition+TempLength)=xThrow(i,TempPosition:TempPosition+TempLength);
               
                if xThrownew==0
                xThrownew=xsave(1,:);
                TrajLengthX3=xsave(1,TempPosition:TempPosition+9);
                else
                xThrownew=[xThrownew',xsave(1,:)']';
                TrajLengthX3=[TrajLengthX3',xsave(1,TempPosition:TempPosition+9)']';
                end
            end
        end
    end
end
xThrow=xThrownew;
for i=1:size(yThrow,1)
    %xtemp=x(i,:);
    TempLength=0;TempPosition=0;
    for j=1:size(yThrow,2)
        if ~isnan(yThrow(i,j))
            if j>1 && ~isnan(yThrow(i,j-1))
            TempLength=TempLength+1;
            else
            TempLength=0;TempPosition=j;    
            end
            if TempLength>NumLength-1 && ((j<size(yThrow,2) && isnan(yThrow(i,j+1))) || j==size(yThrow,2))
                ysave=NaN(1,size(yThrow,2));
                ysave(1,TempPosition:TempPosition+TempLength)=yThrow(i,TempPosition:TempPosition+TempLength);
                if yThrownew==0
                yThrownew=ysave(1,:);
                TrajLengthY3=ysave(1,TempPosition:TempPosition+9);
                else
                yThrownew=[yThrownew',ysave(1,:)']';
                TrajLengthY3=[TrajLengthY3',ysave(1,TempPosition:TempPosition+9)']';
                end
            end
        end
    end
end
yThrow=yThrownew;
end



%%
%Ball filter plots

TrajDuration=zeros(size(x2,1),1);
for i=1:size(x2,1)
    for j=1:size(x2,2)
        if ~isnan(x2(i,j))
            TrajDuration(i,1)=TrajDuration(i,1)+1;
        end
    end
end
x1=sort(TrajDuration(:,1));
[a1,b1]=hist(x1,num);
figure
bar(b1,a1);
DataTrajLength=[b1',a1'];%hhhhhhhhhhhhhhhhh
title('Distribution of traj length');xlabel('frames');ylabel('Ratio');
h=text(max(b1)*0.2,max(a1)*0.8,['#>',num2str(NumLength),'frames =',num2str(size(x2,1))]);
set(h,'FontSize',14);
set(0,'DefaultFigureVisible','off');
set(gca,'FontSize',18);%figurename=[dirnameout,'\',filefolder,'ScatterFinalPosition.pdf'];
figurename=[dirnameout,'\',filefolder,'Ball_TrajDuration.jpg'];
saveas(gcf,figurename)


%%
%Plot mean displacement versus time

TrajLength2=zeros(size(TrajLengthX2,1),10);
for j=2:10
for i=1:size(TrajLengthX2,1)
    if j==1
        TrajLength2(i,j)=sqrt((TrajLengthX2(i,j)-TrajLengthX2(i,j-1)).^2+(TrajLengthY2(i,j)-TrajLengthY2(i,j-1)).^2)/pixel*PhysicalLength;
    else
        
        TrajLength2(i,j)=TrajLength2(i,j-1)+sqrt((TrajLengthX2(i,j)-TrajLengthX2(i,j-1)).^2+(TrajLengthY2(i,j)-TrajLengthY2(i,j-1)).^2)/pixel*PhysicalLength;
    end
end
end

%%display(size(TrajLengthX2,1));%%display(size(TrajLengthY2,1));
xplot=[0:9]*tLength;
err=zeros(1,10);
TrajLengthAverage2=zeros(1,10);
c=jet(5);
figure
    for i=1:10
        TrajLengthAverage2(1,i)=sum(TrajLength2(:,i),1)/size(TrajLength2,1);
        err(1,i)= sqrt(sum((TrajLength2(:,i)-TrajLengthAverage2(1,i)).^2)/size(TrajLength2,1));
    end
    e=errorbar(xplot,TrajLengthAverage2(1,:),err(1,:),'o','MarkerSize',5,...
    'MarkerEdgeColor',c(1,:),'MarkerFaceColor',c(1,:));hold on;
  

yplot=TrajLengthAverage2(1,:);
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
plot(xplot,yfit,'b-','LineWidth',1);
 axis([0,10*tLength, 0, 30])
 h=legend(['Mean velocity= ',num2str(roundn(P(1),-2)),'\mum/s'],'Location','Best');
set(h,'FontSize',14);
set(gca,'FontSize',18);
xlabel('\delta t (s)');
    ylabel('Mean trajectory length (\mu m)');
    title('Mean trajectory length')
    figurename=[dirnameout,'\',filefolder,'Ball_MeanTrajLength.jpg'];
    saveas(gcf,figurename)

VelocityMatBall=roundn(P(1),-2);




TrajLength3=zeros(size(TrajLengthX3,1),10);
for j=2:10
for i=1:size(TrajLengthX3,1)
    if j==1
        TrajLength3(i,j)=sqrt((TrajLengthX3(i,j)-TrajLengthX3(i,j-1)).^2+(TrajLengthY3(i,j)-TrajLengthY3(i,j-1)).^2)/pixel*PhysicalLength;
    else
        
        TrajLength3(i,j)=TrajLength3(i,j-1)+sqrt((TrajLengthX3(i,j)-TrajLengthX3(i,j-1)).^2+(TrajLengthY3(i,j)-TrajLengthY3(i,j-1)).^2)/pixel*PhysicalLength;
    end
end
end
xplot=[0:9]*tLength;
err=zeros(1,10);
TrajLengthAverage3=zeros(1,10);
c=jet(5);
figure
    for i=1:10
        TrajLengthAverage3(1,i)=sum(TrajLength3(:,i),1)/size(TrajLength3,1);
        err(1,i)= sqrt(sum((TrajLength3(:,i)-TrajLengthAverage3(1,i)).^2)/size(TrajLength3,1));
    end
    e=errorbar(xplot,TrajLengthAverage3(1,:),err(1,:),'o','MarkerSize',5,...
    'MarkerEdgeColor',c(1,:),'MarkerFaceColor',c(1,:));hold on;
yplot=TrajLengthAverage3(1,:);
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
plot(xplot,yfit,'b-','LineWidth',1);
 axis([0,10*tLength, 0, 30])
 h=legend(['Mean velocity= ',num2str(roundn(P(1),-2)),'\mum/s'],'Location','Best');
set(h,'FontSize',14);
set(gca,'FontSize',18);
xlabel('\delta t (s)');
    ylabel('Mean trajectory length (\mu m)');
    title('Mean trajectory length')
    figurename=[dirnameout,'\',filefolder,'Remain_MeanTrajLength.jpg'];
    saveas(gcf,figurename)

VelocityMatRemain=roundn(P(1),-2);




if CoObservationLength>size(x2,2)-1
    CoObservationLength=size(x2,2)-1;
end


%%
%Auto-Correlation function with deleting non-motile cells

CorrelationDisplacementX=NaN(size(x2,1),size(x2,2)-1);
CorrelationDisplacementY=NaN(size(x2,1),size(x2,2)-1);
CorrelationDisplacementR=NaN(size(x2,1),size(x2,2)-1);
Correlation=zeros(1,size(x2,2)-1);
CorrelationDivided=zeros(1,size(x2,2)-1);
numberCorrelationR=zeros(1,size(x2,2)-1);

if size(x2,2)==1
    deltaRCorrelation=0;
end

for j=1:size(x2,2)-1
   for i=1:size(x2,1)
        if ~isnan(x2(i,j))  && ~isnan(y2(i,j)) && ~isnan(x2(i,j+1)) && ~isnan(y2(i,j+1))
        CorrelationDisplacementX(i,j)=x2(i,j+1)-x2(i,j);
        CorrelationDisplacementY(i,j)=y2(i,j+1)-y2(i,j);
        CorrelationDisplacementR(i,j)=sqrt(CorrelationDisplacementX(i,j)^2+CorrelationDisplacementY(i,j)^2);
        end
   end
   DeltaRCorrelation=CorrelationDisplacementR(:,j);
    DeltaRCorrelation(isnan(DeltaRCorrelation)==1)=[];
    if j==1
        deltaRCorrelation=DeltaRCorrelation';
    else
    deltaRCorrelation=[deltaRCorrelation,DeltaRCorrelation'];
    end
end

for k=1:size(x2,2)-2
for j=1:size(x2,2)-k-1
   for i=1:size(x2,1)
        if ~isnan(CorrelationDisplacementX(i,j))  && ~isnan(CorrelationDisplacementY(i,j)) && ~isnan(CorrelationDisplacementX(i,j+k)) && ~isnan(CorrelationDisplacementY(i,j+k))
        Correlation(1,k)=Correlation(1,k)+CorrelationDisplacementX(i,j+k)*CorrelationDisplacementX(i,j)+CorrelationDisplacementY(i,j+k)*CorrelationDisplacementY(i,j);        
        numberCorrelationR(1,k)=numberCorrelationR(1,k)+1;
        end
   end
end
Correlation(1,k)=Correlation(1,k)/pixel*PhysicalLength/pixel*PhysicalLength/((size(x2,2)-k-1)*tLength);

end
figure
plot((1:CoObservationLength)*tLength,Correlation(1,1:CoObservationLength),'ro','MarkerSize',4,'MarkerFaceColor','r'); hold on;
set(0,'DefaultFigureVisible','off')
set(gca,'FontSize',18);
%axis([0 size(x2,2)*tLength min(Correlation(1,:)) max(Correlation(1,:))])
xlabel('\tau (s)');
ylabel('\Sigma_{t=0}^{T-\tau}\langle\Deltar(t+\tau)\cdot\Deltar(t)\rangle/(T-\tau) (\mum^2/s)');
title('Auto-correlation')
figurename=[dirnameout,'\',filefolder,'Ball_Correlation.jpg'];
saveas(gcf,figurename)

DataCorrelation=[(1:size(Correlation,2))'*tLength,Correlation(1,1:size(Correlation,2))'];%hhhhhhhhhhhhhhhhh



for BallNum=1:2
if BallNum==2 %To show both ball filter and the remaining
x2=xThrow;y2=yThrow;filefolderBall=[filefolder,'remain_'];
else
filefolderBall=filefolder;    
end

MaxSteplength=Howlongtolook; %how long to look at
NumberOfTraj=zeros(MaxSteplength,1);
DisplacementRb=zeros(MaxSteplength,num);
DisplacementRa=zeros(MaxSteplength,num);
DisplacementXb=zeros(MaxSteplength,num);
DisplacementXa=zeros(MaxSteplength,num);
DisplacementYb=zeros(MaxSteplength,num);
DisplacementYa=zeros(MaxSteplength,num);
DisplacementXbNew=zeros(MaxSteplength,num);
DisplacementXaNew=zeros(MaxSteplength,num);
DisplacementYbNew=zeros(MaxSteplength,num);
DisplacementYaNew=zeros(MaxSteplength,num);
if BallNum==1
DisplacementXbNew2=zeros(MaxSteplength,num);
DisplacementXaNew2=zeros(MaxSteplength,num);
DisplacementYbNew2=zeros(MaxSteplength,num);
DisplacementYaNew2=zeros(MaxSteplength,num);
end
DisplacementRbNew=zeros(MaxSteplength,num);
DisplacementRaNew=zeros(MaxSteplength,num);
Angleb=zeros(MaxSteplength,num);
Anglea=zeros(MaxSteplength,num);

for  Steplength=1:MaxSteplength
    MotileDisThreshold=0;%0.3*pixel/PhysicalLength+0.02*Steplength*pixel/PhysicalLength;%cut non-motile
DisplacementX=NaN(size(x2,1),size(x2,2)-1);
DisplacementY=NaN(size(x2,1),size(x2,2)-1);
DisplacementR=NaN(size(x2,1),size(x2,2)-1);

if size(x2,2)==1
    deltaX=0;deltaR=0;deltaY=0;
    savematV=0;
else
deltaR=[];deltaX=[];deltaY=[];
for i=1:size(x2,1)
    for j=1:size(x2,2)-Steplength
        if ~isnan(x2(i,j))  && ~isnan(y2(i,j)) && ~isnan(x2(i,j+Steplength)) && ~isnan(y2(i,j+Steplength))
        DisplacementX(i,j)=x2(i,j+Steplength)-x2(i,j);
        DisplacementY(i,j)=y2(i,j+Steplength)-y2(i,j);
        end
    end
   for j=1:size(x2,2)-Steplength
        
        if ~isnan(x2(i,j))  && ~isnan(y2(i,j)) && ~isnan(x2(i,j+Steplength)) && ~isnan(y2(i,j+Steplength))
        deltaX=[deltaX,x2(i,j+Steplength)-x2(i,j)];
        deltaY=[deltaY,y2(i,j+Steplength)-y2(i,j)];
        deltaR=[deltaR,sqrt((x2(i,j+Steplength)-x2(i,j))^2+(y2(i,j+Steplength)-y2(i,j))^2)];
        break;
        end
   end

end
end

    deltaX=deltaX';
    x1=sort((deltaX(:,1))/pixel*PhysicalLength);
    [a1,b1]=hist(x1,num);a1=a1/length(x1);
    DisplacementXb(Steplength,:)=b1;
    DisplacementXa(Steplength,:)=a1;
    NumberOfTraj(Steplength,1)=size(deltaX,1);
    
    deltaXNew=deltaX;
    deltaXNew(abs(deltaXNew)<MotileDisThreshold)=[];
    x1=sort((deltaXNew(:,1))/pixel*PhysicalLength);
    [a1,b1]=hist(x1,num);a1=a1/length(x1); 
    DisplacementXbNew(Steplength,:)=b1;
    DisplacementXaNew(Steplength,:)=a1;
    DisplacementXaNew(Steplength,:)=DisplacementXaNew(Steplength,:)/sum(DisplacementXaNew(Steplength,:));
    if BallNum==1
    DisplacementXbNew2(Steplength,:)=b1;
    DisplacementXaNew2(Steplength,:)=a1;
    DisplacementXaNew2(Steplength,:)=DisplacementXaNew2(Steplength,:)/sum(DisplacementXaNew2(Steplength,:));
    end
    
      deltaY=deltaY';
    x1=sort((deltaY(:,1))/pixel*PhysicalLength);
    [a1,b1]=hist(x1,num);a1=a1/length(x1);
   
    DisplacementYb(Steplength,:)=b1;
    DisplacementYa(Steplength,:)=a1;
    deltaYNew=deltaY;
    deltaYNew(abs(deltaYNew)<MotileDisThreshold)=[];
    x1=sort((deltaYNew(:,1))/pixel*PhysicalLength);
    [a1,b1]=hist(x1,num);a1=a1/length(x1);   
    DisplacementYbNew(Steplength,:)=b1;
    DisplacementYaNew(Steplength,:)=a1;
    DisplacementYaNew(Steplength,:)=DisplacementYaNew(Steplength,:)/sum(DisplacementYaNew(Steplength,:));
     if BallNum==1
        DisplacementYbNew2(Steplength,:)=b1;
    DisplacementYaNew2(Steplength,:)=a1; 
    DisplacementYaNew2(Steplength,:)=DisplacementYaNew2(Steplength,:)/sum(DisplacementYaNew2(Steplength,:));
   
     end
    
     
        deltaR=deltaR';
    x1=sort((deltaR(:,1))/pixel*PhysicalLength); 
    [a1,b1]=hist(x1,num);a1=a1/length(x1);
    DisplacementRb(Steplength,:)=b1;
    DisplacementRa(Steplength,:)=a1;
 DisplacementRbNew(Steplength,:)=b1;
    DisplacementRaNew(Steplength,:)=a1;
    DisplacementRaNew(Steplength,:)=DisplacementRaNew(Steplength,:)/sum(DisplacementRaNew(Steplength,:));
 
 
if Steplength==1 && BallNum==1
    CorrelationDisplacementX=NaN(size(x2,1),size(x2,2)-1);
CorrelationDisplacementY=NaN(size(x2,1),size(x2,2)-1);
CorrelationDisplacementR=NaN(size(x2,1),size(x2,2)-1);
Correlation=zeros(1,size(x2,2)-1);
CorrelationDivided=zeros(1,size(x2,2)-1);
numberCorrelationR=zeros(1,size(x2,2)-1);

if size(x2,2)==1
    deltaRCorrelation=0;
else

for j=1:size(x2,2)-1
   for i=1:size(x2,1)
        if ~isnan(x2(i,j))  && ~isnan(y2(i,j)) && ~isnan(x2(i,j+1)) && ~isnan(y2(i,j+1))
        CorrelationDisplacementX(i,j)=x2(i,j+1)-x2(i,j);
        CorrelationDisplacementY(i,j)=y2(i,j+1)-y2(i,j);
        CorrelationDisplacementR(i,j)=sqrt(CorrelationDisplacementX(i,j)^2+CorrelationDisplacementY(i,j)^2);
        end
   end
   DeltaRCorrelation=CorrelationDisplacementR(:,j);
    DeltaRCorrelation(isnan(DeltaRCorrelation)==1)=[];
    if j==1
        deltaRCorrelation=DeltaRCorrelation';
    else
    deltaRCorrelation=[deltaRCorrelation,DeltaRCorrelation'];
    end
end

for k=1:size(x2,2)-2
for j=1:size(x2,2)-k-1
   for i=1:size(x2,1)
        if ~isnan(CorrelationDisplacementX(i,j))  && ~isnan(CorrelationDisplacementY(i,j)) && ~isnan(CorrelationDisplacementX(i,j+k)) && ~isnan(CorrelationDisplacementY(i,j+k))
        Correlation(1,k)=Correlation(1,k)+CorrelationDisplacementX(i,j+k)*CorrelationDisplacementX(i,j)+CorrelationDisplacementY(i,j+k)*CorrelationDisplacementY(i,j);        
        numberCorrelationR(1,k)=numberCorrelationR(1,k)+1;
        end
   end
end
end

 AnglesNew=NaN(size(x2,1),size(x2,2)-2);
for j=1:size(x2,2)-Steplength-1
   for i=1:size(x2,1)
        if ~isnan(CorrelationDisplacementX(i,j))  && ~isnan(CorrelationDisplacementY(i,j)) && ~isnan(CorrelationDisplacementX(i,j+Steplength)) && ~isnan(CorrelationDisplacementY(i,j+Steplength))
            A=[CorrelationDisplacementX(i,j),CorrelationDisplacementY(i,j)];
            B=[CorrelationDisplacementX(i,j+Steplength),CorrelationDisplacementY(i,j+Steplength)];
            AnglesNew(i,j)=acos(dot(A,B)/(norm(A)*norm(B)))*180/pi;
        end
   end
    Delta=AnglesNew(:,j);
    Delta(isnan(Delta)==1)=[];
    if j==1
        deltaAngle=Delta';
    else
    deltaAngle=[deltaAngle,Delta'];
    end
end
    deltaAngle=deltaAngle';
    if size(deltaAngle,2)>0 && size(deltaAngle,1)>0       
        x1=sort((deltaAngle(:,1)));
        [a1,b1]=hist(x1,num);a1=a1/length(x1);
    else
        a1=0;b1=0;
    end
    [a1,b1]=hist(x1,num);a1=a1/length(x1);
    Angleb(Steplength,:)=b1;
    Anglea(Steplength,:)=a1;
    
    figure
    bar(b1,a1);
    a=get(gca);
    xxx=a.XLim;%
yyy=a.YLim;%
k=[0.2 0.9];%
x0=xxx(1)+k(1)*(xxx(2)-xxx(1));%
y0=yyy(1)+k(2)*(yyy(2)-yyy(1));%
set(gca,'FontSize',18);
    text(x0,y0,['#displacements=',num2str(size(deltaAngle,1)),'\newline\deltat=',num2str(Steplength*tLength),'s'],'FontSize',18);
    xlabel('\Delta\theta=\theta(t+\deltat)-\theta(t)');
    ylabel('Ratio');
    title('Angular distribution')
    
    figurename=[dirnameout,'\',filefolder,'Step',num2str(Steplength),'Ball_AngleDistribution.jpg'];
    saveas(gcf,figurename)
    end
end
end

%%
%Plot X-direction's mean displacement versus time and mean square displacement versus time
MeanDisplacementX=sum(DisplacementXb.*DisplacementXa,2)';
DiffusionCoeffcientX=sum(DisplacementXb.^2.*DisplacementXa,2)'-MeanDisplacementX.^2;

MeanDisplacementXNew=sum(DisplacementXbNew.*DisplacementXaNew,2)';
DiffusionCoeffcientXNew=sum(DisplacementXbNew.^2.*DisplacementXaNew,2)'-MeanDisplacementXNew.^2;
DriftXTotal(window,1:MaxSteplength)=MeanDisplacementXNew(1,1:MaxSteplength);

if PlotDriftDiffusionFit==1 && window==1
figure
xplot=(0:MaxSteplength)*tLength;
yplot=[0,MeanDisplacementX(1,:)];
p1=plot(xplot,yplot,'ro','MarkerSize',10,'MarkerFaceColor','r'); hold on;
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
P1=P(1);P2=P(2);
p2=plot(xplot,yfit,'r-','LineWidth',1);
yplot=[0,MeanDisplacementXNew(1,:)];
p3=plot(xplot,yplot,'b>','MarkerSize',10,'MarkerFaceColor','b'); hold on;
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
p4=plot(xplot,yfit,'b-','LineWidth',1);
set(0,'DefaultFigureVisible','off')
set(gca,'FontSize',18);
a=get(gca);
xxx=a.XLim;%
yyy=a.YLim;%
k=[0.5 0.9];%
x0=xxx(1)+k(1)*(xxx(2)-xxx(1));%
y0=yyy(1)+k(2)*(yyy(2)-yyy(1));%

set(gca,'FontSize',18);
%axis([0,MaxSteplength*tLength, min(min(MeanDisplacementXNew(1,:)),min(MeanDisplacementX(1,:))), max(max(MeanDisplacementXNew(1,:)),max(MeanDisplacementX(1,:)))+0.5])
axis([0,MaxSteplength*tLength, -8, 8])
xlabel('\deltat (s)');
ylabel('\langlex(t+\deltat)-x(t)\rangle (\mum)');
title('First moment of X-displacement')
 h=legend('For all trajectories',['Linear fit: x=',num2str(roundn(P1,-2)),'*\deltat+(',num2str(roundn(P2,-2)),')'],'For motile',['Linear fit: x=',num2str(roundn(P(1),-2)),'*\deltat+(',num2str(roundn(P(2),-2)),')'],'Location','Best');
set(h,'FontSize',14);
figurename=[dirnameout,'\',filefolderBall,'Ball_DisplacementX.jpg'];
saveas(gcf,figurename)


if BallNum==1
Ball_DriftXmat=roundn(P(1),-2)/2;
else
Ball_DriftXmat_Remain=roundn(P(1),-2)/2;
end



figure
xplot=(0:MaxSteplength)*tLength;
yplot=[0,DiffusionCoeffcientX(1,:)];
p1=plot(xplot,yplot,'ro','MarkerSize',10,'MarkerFaceColor','r'); hold on;
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
P1=P(1);P2=P(2);
p2=plot(xplot,yfit,'r-','LineWidth',1);
yplot=[0,DiffusionCoeffcientXNew(1,:)];
p3=plot(xplot,yplot,'b>','MarkerSize',10,'MarkerFaceColor','b'); hold on;
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
p4=plot(xplot,yfit,'b-','LineWidth',1);
set(0,'DefaultFigureVisible','off')
set(gca,'FontSize',18);
a=get(gca);
xxx=a.XLim;%
yyy=a.YLim;%
k=[0.5 0.9];%
x0=xxx(1)+k(1)*(xxx(2)-xxx(1));%
y0=yyy(1)+k(2)*(yyy(2)-yyy(1));%
set(gca,'FontSize',18);
%axis([0,MaxSteplength*tLength, -8, 8])
xlabel('\deltat (s)');
ylabel('Varianc (\mum^2)');
title('Varianc of X-displacement')
 h=legend('For all trajectories',['Linear fit: V_x=',num2str(roundn(P1,-2)),'*\deltat+(',num2str(roundn(P2,-2)),')'],'For motile',['Linear fit: V_x=',num2str(roundn(P(1),-2)),'*\deltat+(',num2str(roundn(P(2),-2)),')'],'Location','Best');
set(h,'FontSize',14);
figurename=[dirnameout,'\',filefolderBall,'Ball_DiffusionX.jpg'];
%saveas(gcf,figurename)
end

%%
%Plot Y-direction's mean displacement versus time and mean square displacement versus time

MeanDisplacementY=sum(DisplacementYb.*DisplacementYa,2)';
DiffusionCoeffcientY=sum(DisplacementYb.^2.*DisplacementYa,2)'-MeanDisplacementY.^2;

MeanDisplacementYNew=sum(DisplacementYbNew.*DisplacementYaNew,2)';
DiffusionCoeffcientYNew=sum(DisplacementYbNew.^2.*DisplacementYaNew,2)'-MeanDisplacementYNew.^2;
DriftYTotal(window,1:MaxSteplength)=MeanDisplacementYNew(1,1:MaxSteplength);

if PlotDriftDiffusionFit==1 && window==1
    figure
xplot=(0:MaxSteplength)*tLength;
yplot=[0,MeanDisplacementY(1,:)];
p1=plot(xplot,yplot,'ro','MarkerSize',10,'MarkerFaceColor','r'); hold on;
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
P1=P(1);P2=P(2);
p2=plot(xplot,yfit,'r-','LineWidth',1);
yplot=[0,MeanDisplacementYNew(1,:)];
p3=plot(xplot,yplot,'b>','MarkerSize',10,'MarkerFaceColor','b'); hold on;
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
p4=plot(xplot,yfit,'b-','LineWidth',1);
set(0,'DefaultFigureVisible','off')
set(gca,'FontSize',18);
a=get(gca);
    xxx=a.XLim;%
yyy=a.YLim;%
k=[0.5 0.9];%
x0=xxx(1)+k(1)*(xxx(2)-xxx(1));%
y0=yyy(1)+k(2)*(yyy(2)-yyy(1));%
set(gca,'FontSize',18);
axis([0,MaxSteplength*tLength, -8, 8])
xlabel('\deltat (s)');
ylabel('\langley(t+\deltat)-y(t)\rangle (\mum)');
title('First moment of Y-displacement')
 h=legend('For all trajectories',['Linear fit: y=',num2str(roundn(P1,-2)),'*\deltat+(',num2str(roundn(P2,-2)),')',],'For motile',['Linear fit: y=',num2str(roundn(P(1),-2)),'*\deltat+(',num2str(roundn(P(2),-2)),')'],'Location','Best');
set(h,'FontSize',14);figurename=[dirnameout,'\',filefolder,'DisplacementY.pdf'];
figurename=[dirnameout,'\',filefolderBall,'Ball_DisplacementY.jpg'];
saveas(gcf,figurename)


if BallNum==1
Ball_DriftYmat=roundn(P(1),-2)/2;
else
Ball_DriftYmat_Remain=roundn(P(1),-2)/2;
end



figure
xplot=(0:MaxSteplength)*tLength;
yplot=[0,DiffusionCoeffcientY(1,:)];
p1=plot(xplot,yplot,'ro','MarkerSize',10,'MarkerFaceColor','r'); hold on;
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
P1=P(1);P2=P(2);
p2=plot(xplot,yfit,'r-','LineWidth',1);
yplot=[0,DiffusionCoeffcientYNew(1,:)];
p3=plot(xplot,yplot,'b>','MarkerSize',10,'MarkerFaceColor','b'); hold on;
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
p4=plot(xplot,yfit,'b-','LineWidth',1);
set(0,'DefaultFigureVisible','off')
set(gca,'FontSize',18);
a=get(gca);
xxx=a.XLim;%
yyy=a.YLim;%
k=[0.5 0.9];%
x0=xxx(1)+k(1)*(xxx(2)-xxx(1));%
y0=yyy(1)+k(2)*(yyy(2)-yyy(1));%
set(gca,'FontSize',18);
%axis([0,MaxSteplength*tLength, -8, 8])
xlabel('\deltat (s)');
ylabel('Varianc (\mum^2)');
title('Varianc of Y-displacement')
 h=legend('For all trajectories',['Linear fit: V_y=',num2str(roundn(P1,-2)),'*\deltat+(',num2str(roundn(P2,-2)),')'],'For motile',['Linear fit: V_y=',num2str(roundn(P(1),-2)),'*\deltat+(',num2str(roundn(P(2),-2)),')'],'Location','Best');
set(h,'FontSize',14);
figurename=[dirnameout,'\',filefolderBall,'Ball_DiffusionY.jpg'];
%saveas(gcf,figurename)
end

%MeanSquare and Diffusion
MeanDisplacement=sqrt(MeanDisplacementXNew.^2+MeanDisplacementYNew.^2);
MeanDisplacementOld=sqrt(MeanDisplacementX.^2+MeanDisplacementY.^2);
MeanSquare=sum(DisplacementRb.^2.*DisplacementRa,2)';
DiffusionCoeffcient=sum(DisplacementRb.^2.*DisplacementRa,2)'-MeanDisplacementOld.^2;
MeanSquareNew=sum(DisplacementRb.^2.*DisplacementRaNew,2)';
DiffusionCoeffcientNew=sum(DisplacementRb.^2.*DisplacementRaNew,2)'-MeanDisplacement.^2;
DiffusionTotal(window,1:MaxSteplength)=DiffusionCoeffcientNew(1,1:MaxSteplength);

if PlotDriftDiffusionFit==1 && window==1
 
%Diffusion-----------------------
figure
xplot=(0:MaxSteplength)*tLength;
yplot=[0,DiffusionCoeffcient(1,:)];
p1=plot(xplot,yplot,'ro','MarkerSize',10,'MarkerFaceColor','r'); hold on;
set(gca,'FontSize',18);
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
P1=P(1);P2=P(2);
p2=plot(xplot,yfit,'r-','LineWidth',1);
yplot=[0,DiffusionCoeffcientNew(1,:)];
display(yplot);
if BallNum==1
DiffusionRx=xplot;
DiffusionRy=[0,DiffusionCoeffcientNew(1,:)];
end

p3=plot(xplot,yplot,'b>','MarkerSize',10,'MarkerFaceColor','b');
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
p4=plot(xplot,yfit,'b-','LineWidth',1);
%text(max(xplot)*0.05,max(yplot)*0.6,['all trajectories: V_r=',num2str(roundn(P1,-2)),'*\deltat+(',num2str(roundn(P2,-2)),'\newlinemotile: V_r=',num2str(roundn(P(1),-2)),'*\deltat+(',num2str(roundn(P(2),-2)),')'],'FontSize',18);
xlabel('\deltat (s)');
axis([min(xplot), max(xplot), min(yplot), max(yplot)+300]);
ylabel('Variance of displacement')
set(gca,'FontSize',18);
 h=legend('For all trajectories',['Linear fit: V_r=',num2str(roundn(P1,-2)),'*\deltat+(',num2str(roundn(P2,-2)),')'],'For trajecories of motile',['Linear fit: V_r=',num2str(roundn(P(1),-2)),'*\deltat+(',num2str(roundn(P(2),-2)),')'],'Location','best');
set(h,'FontSize',14);figurename=[dirnameout,'\',filefolder,'Diffusion.pdf'];
%saveas(gcf,figurename)
figurename=[dirnameout,'\',filefolderBall,'Ball_Diffusion.jpg'];
saveas(gcf,figurename)

if BallNum==1
Ball_DiffusionRmat=roundn(P(1),-2)/2;
else
Ball_DiffusionRmat_Remain=roundn(P(1),-2)/2;
end


end

%% 
%Angular-distributed displacement
if size(x2,2)>1
for Steplength=[Howlongtolook]
NumOfAngle=60;
NumOfTraj=zeros(NumOfAngle,1);
StepLengthInAngle=cell(NumOfAngle,1);StepLengthInAnglePlot=zeros(NumOfAngle,1);
AngleDistri=[-180:360/NumOfAngle:180-360/NumOfAngle];
%AngleDistriPlot=[-180:360/NumOfAngle:180-360/NumOfAngle];
for i=1:size(y2,1)
    for j=1:size(y2,2)-Steplength
        if ~isnan(y2(i,j)) && ~isnan(x2(i,j)) &&  ~isnan(y2(i,j+Steplength)) && ~isnan(x2(i,j+Steplength))
            AngleTemp=rad2deg(atan2((y2(i,j+Steplength)-y2(i,j)),(x2(i,j+Steplength)-x2(i,j))));
            indexTemp=find(abs(AngleDistri-AngleTemp)==min(abs(AngleDistri-AngleTemp)));
            %%display(AngleTemp);%display(indexTemp);
            if size(indexTemp,2)>1
                indexTemp=indexTemp(1,1);
            end
            if AngleTemp-AngleDistri(1,indexTemp)<0
                indexTemp=indexTemp-1;
            end
            %NumOfTraj(indexTemp,1)=NumOfTraj(indexTemp,1)+1;
            StepLengthInAngle{indexTemp,1}=[StepLengthInAngle{indexTemp,1},sqrt((y2(i,j+Steplength)-y2(i,j))^2+(x2(i,j+Steplength)-x2(i,j))^2)/pixel*PhysicalLength];
        end
    end
end
for i=1:NumOfAngle
    NumOfTraj(i,1)=size(StepLengthInAngle{i,1},2);
    %display(StepLengthInAngle);%display(NumOfTraj);
    if NumOfTraj(i,1)==0
        StepLengthInAnglePlot(i,1)=0;
    else
    StepLengthInAnglePlot(i,1)=sum(StepLengthInAngle{i,1},2)/NumOfTraj(i,1);
    end
end


AngleDistri=[AngleDistri/360*2*pi, pi];
NumOfTrajNew=[NumOfTraj',NumOfTraj(1,1)];
StepLengthInAngleNew=[StepLengthInAnglePlot',StepLengthInAnglePlot(1,1)];
StepLengthInAngleSumNew=[StepLengthInAnglePlot'.*NumOfTraj',StepLengthInAnglePlot(1,1).*NumOfTraj(1,1)];

figure
polarplot(AngleDistri,NumOfTrajNew,'LineWidth',1.5);hold on;
ax = gca;ax.LineWidth = 1.5;ax.ThetaTick = [0 45 90 135 180 225 270 315];  ax.RAxisLocation = 90;%ax.RAxisLocation='manual';
title('# of displacement');
set(gca,'FontSize',15);
figurename=[dirnameout,'\',filefolderBall,'Ball_Step',num2str(Steplength),'NumAngluar.jpeg'];
hold off; print(gcf, '-djpeg', '-r600',figurename)


figure
polarplot(AngleDistri,StepLengthInAngleNew,'LineWidth',1.5);hold on;
ax = gca;ax.LineWidth = 1.5;ax.ThetaTick = [0 45 90 135 180 225 270 315];ax.RAxisLocation = 90;
title('Average displacement (\mum)');
set(gca,'FontSize',15);
figurename=[dirnameout,'\',filefolderBall,'Ball_Step',num2str(Steplength),'LengthAngluar.jpeg'];
hold off; print(gcf, '-djpeg', '-r600',figurename)


figure
polarplot(AngleDistri,StepLengthInAngleSumNew,'LineWidth',1.5);hold on;
ax = gca;ax.LineWidth = 1.5;ax.ThetaTick = [0 45 90 135 180 225 270 315];ax.RAxisLocation = 90;
title('Sum of displacement (\mum)');
set(gca,'FontSize',15);
figurename=[dirnameout,'\',filefolderBall,'Ball_Step',num2str(Steplength),'TotalLengthAngluar.jpeg'];
hold off; print(gcf, '-djpeg', '-r600',figurename)

end

if BallNum==1
QuadrantXNeg=[[1:1:ceil(NumOfAngle/8)-1],[ceil(NumOfAngle/8*7):1:NumOfAngle]];
QuadrantXPos=ceil(NumOfAngle/8*3):1:ceil(NumOfAngle/8*5)-1;
QuadrantYNeg=ceil(NumOfAngle/8):1:ceil(NumOfAngle/8*3)-1;
QuadrantYPos=ceil(NumOfAngle/8*5):1:ceil(NumOfAngle/8*7)-1;
    
YPositivePro=sum(NumOfTrajNew(1,QuadrantYPos),2)/(sum(NumOfTrajNew(1,QuadrantYPos),2)+sum(NumOfTrajNew(1,QuadrantYNeg),2));
XPositivePro=sum(NumOfTrajNew(1,QuadrantXPos),2)/(sum(NumOfTrajNew(1,QuadrantXPos),2)+sum(NumOfTrajNew(1,QuadrantXNeg),2));
XPositiveLength=sum(StepLengthInAngleNew(1,QuadrantXPos),2)/size(QuadrantXPos,2);
XNegativeLength=sum(StepLengthInAngleNew(1,QuadrantXNeg),2)/size(QuadrantXNeg,2);
YPositiveLength=sum(StepLengthInAngleNew(1,QuadrantYPos),2)/size(QuadrantYPos,2);
YNegativeLength=sum(StepLengthInAngleNew(1,QuadrantYNeg),2)/size(QuadrantYNeg,2);
else
QuadrantXNeg=[[1:1:ceil(NumOfAngle/8)-1],[ceil(NumOfAngle/8*7):1:NumOfAngle]];
QuadrantXPos=ceil(NumOfAngle/8*3):1:ceil(NumOfAngle/8*5)-1;
QuadrantYNeg=ceil(NumOfAngle/8):1:ceil(NumOfAngle/8*3)-1;
QuadrantYPos=ceil(NumOfAngle/8*5):1:ceil(NumOfAngle/8*7)-1;

YPositiveProTotal=sum(NumOfTrajNew(1,QuadrantYPos),2)/(sum(NumOfTrajNew(1,QuadrantYPos),2)+sum(NumOfTrajNew(1,QuadrantYNeg),2));
XPositiveProTotal=sum(NumOfTrajNew(1,QuadrantXPos),2)/(sum(NumOfTrajNew(1,QuadrantXPos),2)+sum(NumOfTrajNew(1,QuadrantXNeg),2));
XPositiveLengthTotal=sum(StepLengthInAngleNew(1,QuadrantXPos),2)/size(QuadrantXPos,2);
XNegativeLengthTotal=sum(StepLengthInAngleNew(1,QuadrantXNeg),2)/size(QuadrantXNeg,2);
YPositiveLengthTotal=sum(StepLengthInAngleNew(1,QuadrantYPos),2)/size(QuadrantYPos,2);
YNegativeLengthTotal=sum(StepLengthInAngleNew(1,QuadrantYNeg),2)/size(QuadrantYNeg,2);
end


else
 QuadrantXNeg=0;
QuadrantXPos=0;
QuadrantYNeg=0;
QuadrantYPos=0;  
YPositivePro=0;
XPositivePro=0;
YPositiveProTotal=0;
XPositiveProTotal=0;XPositiveLength=0;XNegativeLength=0;YPositiveLength=0;YNegativeLength=0;XPositiveLengthTotal=0;
XNegativeLengthTotal=0;YPositiveLengthTotal=0;YNegativeLengthTotal=0;
end

end

end


%%
%Statistics after using the Gaussian Filter to remove the non-motile cells 


Analysis=1;
%Plot traj duration distribution
TrajDuration=zeros(size(x,1),1);
for i=1:size(x,1)
    for j=1:size(x,2)
        if ~isnan(x(i,j))
            TrajDuration(i,1)=TrajDuration(i,1)+1;
        end
    end
end

x1=sort(TrajDuration(:,1));
[a1,b1]=hist(x1,num);

figure

bar(b1,a1);

DataTrajLength=[b1',a1'];%

%DataTrajLength=[{'TrajLength_frames',b1}',{'TrajLength_Ratio',a1}'];%hhhhhhhhhhhhhhhhh
title('Distribution of traj length');xlabel('frames');ylabel('Ratio');
h=text(max(b1)*0.2,max(a1)*0.8,['#>',num2str(NumLength),'frames =',num2str(size(x,1))]);
set(h,'FontSize',14);
set(0,'DefaultFigureVisible','off');
set(gca,'FontSize',18);%figurename=[dirnameout,'\',filefolder,'ScatterFinalPosition.pdf'];
%saveas(gcf,figurename)   
figurename=[dirnameout,'\',filefolder,'TrajDuration.jpg'];
saveas(gcf,figurename)


%%
%plot traj length distribution
TrajLength=zeros(size(TrajLengthX,1),10);
for j=2:10
for i=1:size(TrajLengthX,1)
    if j==1
        TrajLength(i,j)=sqrt((TrajLengthX(i,j)-TrajLengthX(i,j-1)).^2+(TrajLengthY(i,j)-TrajLengthY(i,j-1)).^2)/pixel*PhysicalLength;
    else
        
        TrajLength(i,j)=TrajLength(i,j-1)+sqrt((TrajLengthX(i,j)-TrajLengthX(i,j-1)).^2+(TrajLengthY(i,j)-TrajLengthY(i,j-1)).^2)/pixel*PhysicalLength;
    end
end
end
xplot=[0:9]*tLength;
err=zeros(1,10);
TrajLengthAverage=zeros(1,10);
c=jet(5);
figure
    for i=1:10
        TrajLengthAverage(1,i)=sum(TrajLength(:,i),1)/size(TrajLength,1);
        err(1,i)= sqrt(sum((TrajLength(:,i)-TrajLengthAverage(1,i)).^2)/size(TrajLength,1));
    end
    e=errorbar(xplot,TrajLengthAverage(1,:),err(1,:),'o','MarkerSize',5,...
    'MarkerEdgeColor',c(1,:),'MarkerFaceColor',c(1,:));hold on;
  

yplot=TrajLengthAverage(1,:);
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
plot(xplot,yfit,'b-','LineWidth',1);
 axis([0,10*tLength, 0, 50])
 h=legend(['Mean velocity= ',num2str(roundn(P(1),-2)),'\mum/s'],'Location','Best');
set(h,'FontSize',14);
set(gca,'FontSize',18);
xlabel('\delta t (s)');
    ylabel('Mean trajectory length (\mu m)');
    title('Mean trajectory length')
    figurename=[dirnameout,'\',filefolder,'MeanTrajLength.jpg'];
    saveas(gcf,figurename)
    figurename=[dirnameout,'\',filefolder,'MeanTrajLength.fig'];
    %saveas(gcf,figurename)
VelocityMat=roundn(P(1),-2);

%%
MaxSteplength=Howlongtolook; %how long to look at
if Analysis==1
NumberOfTraj=zeros(MaxSteplength,1);
DisplacementRb=zeros(MaxSteplength,num);
DisplacementRa=zeros(MaxSteplength,num);
DisplacementXb=zeros(MaxSteplength,num);
DisplacementXa=zeros(MaxSteplength,num);
DisplacementYb=zeros(MaxSteplength,num);
DisplacementYa=zeros(MaxSteplength,num);
DisplacementXbNew=zeros(MaxSteplength,num);
DisplacementXaNew=zeros(MaxSteplength,num);
DisplacementYbNew=zeros(MaxSteplength,num);
DisplacementYaNew=zeros(MaxSteplength,num);
DisplacementRbNew=zeros(MaxSteplength,num);
DisplacementRaNew=zeros(MaxSteplength,num);
Angleb=zeros(MaxSteplength,num);
Anglea=zeros(MaxSteplength,num);

%%
%Auto-Correlation function
CorrelationDisplacementX=NaN(size(x,1),size(x,2)-1);
CorrelationDisplacementY=NaN(size(x,1),size(x,2)-1);
Correlation=zeros(1,size(x,2)-1);
CorrelationDivided=zeros(1,size(x,2)-1);
numberCorrelationR=zeros(1,size(x,2)-1);


for j=1:size(x,2)-1
   for i=1:size(x,1)
        if ~isnan(x(i,j))  && ~isnan(y(i,j)) && ~isnan(x(i,j+1)) && ~isnan(y(i,j+1))
        CorrelationDisplacementX(i,j)=x(i,j+1)-x(i,j);
        CorrelationDisplacementY(i,j)=y(i,j+1)-y(i,j);
        end
   end
end

for k=1:size(x,2)-2
for j=1:size(x,2)-k-1
   for i=1:size(x,1)
        if ~isnan(CorrelationDisplacementX(i,j))  && ~isnan(CorrelationDisplacementY(i,j)) && ~isnan(CorrelationDisplacementX(i,j+k)) && ~isnan(CorrelationDisplacementY(i,j+k))
        Correlation(1,k)=Correlation(1,k)+CorrelationDisplacementX(i,j+k)*CorrelationDisplacementX(i,j)+CorrelationDisplacementY(i,j+k)*CorrelationDisplacementY(i,j);        
        numberCorrelationR(1,k)=numberCorrelationR(1,k)+1;
        end
   end
end
Correlation(1,k)=Correlation(1,k)/pixel*PhysicalLength/pixel*PhysicalLength/((size(x,2)-k-1)*tLength);

end
figure
plot((1:CoObservationLength)*tLength,Correlation(1,1:CoObservationLength),'ro','MarkerSize',4,'MarkerFaceColor','r'); hold on;
set(0,'DefaultFigureVisible','off')
set(gca,'FontSize',18);
%axis([0 size(x,2)*tLength min(Correlation(1,:)) max(Correlation(1,:))])
xlabel('\tau (s)');
ylabel('\Sigma_{t=0}^{T-\tau}\langle\Deltar(t+\tau)\cdot\Deltar(t)\rangle/(T-\tau) (\mum^2/s)');
title('Auto-correlation')
figurename=[dirnameout,'\',filefolder,'Correlation.jpg'];
saveas(gcf,figurename)

%%

for  Steplength=1:MaxSteplength
    MotileDisThreshold=0;%0.3*pixel/PhysicalLength+0.02*Steplength*pixel/PhysicalLength;%cut non-motile
DisplacementX=NaN(size(x,1),size(x,2)-1);
DisplacementY=NaN(size(x,1),size(x,2)-1);
DisplacementR=NaN(size(x,1),size(x,2)-1);

if size(x,2)==1
    deltaX=0;deltaR=0;deltaY=0;
    savematV=0;
else
deltaR=[];deltaX=[];deltaY=[];
for i=1:size(x,1)
    for j=1:size(x,2)-Steplength
        if ~isnan(x(i,j))  && ~isnan(y(i,j)) && ~isnan(x(i,j+Steplength)) && ~isnan(y(i,j+Steplength))
        DisplacementX(i,j)=x(i,j+Steplength)-x(i,j);
        DisplacementY(i,j)=y(i,j+Steplength)-y(i,j);
        end
    end
   for j=1:size(x,2)-Steplength
        
        if ~isnan(x(i,j))  && ~isnan(y(i,j)) && ~isnan(x(i,j+Steplength)) && ~isnan(y(i,j+Steplength))
        deltaX=[deltaX,x(i,j+Steplength)-x(i,j)];
        deltaY=[deltaY,y(i,j+Steplength)-y(i,j)];
        deltaR=[deltaR,sqrt((x(i,j+Steplength)-x(i,j))^2+(y(i,j+Steplength)-y(i,j))^2)];
        break;
        end
   end

end
end 
    deltaX=deltaX';
    x1=sort((deltaX(:,1))/pixel*PhysicalLength);
    [a1,b1]=hist(x1,num);a1=a1/length(x1);
    DisplacementXb(Steplength,:)=b1;
    DisplacementXa(Steplength,:)=a1;
    NumberOfTraj(Steplength,1)=size(deltaX,1);
    
    deltaXNew=deltaX;
    deltaXNew(abs(deltaXNew)<MotileDisThreshold)=[];
    x1=sort((deltaXNew(:,1))/pixel*PhysicalLength);
    [a1,b1]=hist(x1,num);a1=a1/length(x1);
    if PlotDetail==1  
    figure
    set(0,'DefaultFigureVisible','off')
    set(gca,'FontSize',18);
    plot(b1,a1,'ro','MarkerSize',10,'MarkerFaceColor','r');
    
    
    
    [m1, ind1]=max(a1);
    aa=a1;bb=b1;
    aa(ind1)=0;
    [m2, ind2]=max(aa);
    aa([ind1,ind2])=[];bb([ind1,ind2])=[];
    %%display(size(aa));
    if size(aa,2)>1 || size(aa,1)>1
    opts=fitoptions('Method','NonlinearLeastSquares');
        opts.Lower=[-100000000000000 -1000000000000000];
    f = fit(bb.',aa.','gauss1',opts);
    va1=f(b1(ind1));va2=f(b1(ind2));
    fraction=a1(ind1)-va1+a1(ind2)-va2;
    a1(ind1)=va1;a1(ind2)=va2;
    hold on;plot(b1,f(b1),'b-','LineWidth',1);
    else
        fraction=0;
        
    end

    fractionX=fraction;
  
    DisplacementXbNew(Steplength,:)=b1;
    DisplacementXaNew(Steplength,:)=a1;
    DisplacementXaNew(Steplength,:)=DisplacementXaNew(Steplength,:)/sum(DisplacementXaNew(Steplength,:));
   hold on;
    
    axis([b1(1,1),b1(1,size(b1,2)), 0.0001, m1+2])
    
    a=get(gca);
    xxx=a.XLim;%
    yyy=a.YLim;%
    k=[0.01 0.5];%
    x0=xxx(1)+k(1)*(xxx(2)-xxx(1));%
    y0=yyy(1)+k(2)*(yyy(2)-yyy(1));%
    set(gca,'FontSize',18);
    text(x0,y0,['#disp=',num2str(size(deltaX,1)),', \deltat=',num2str(Steplength*tLength),'s','\newlineFraction=',num2str(roundn(fraction,-2))],'FontSize',16);   
    xlabel('\Deltax=x(t+\deltat)-x(t) (\mum)');
    ylabel('Ratio');
    
    title('X-displacement')
    if Steplength==1

     set(gca,'YScale','log');
 
    figurename=[dirnameout,'\',filefolder,'Step',num2str(Steplength),'XDistributionLog.jpg'];
    saveas(gcf,figurename)
    end

    end
%     
    deltaY=deltaY';
    x1=sort((deltaY(:,1))/pixel*PhysicalLength);
    [a1,b1]=hist(x1,num);a1=a1/length(x1);
   
    DisplacementYb(Steplength,:)=b1;
    DisplacementYa(Steplength,:)=a1;
    deltaYNew=deltaY;
    deltaYNew(abs(deltaYNew)<MotileDisThreshold)=[];
    x1=sort((deltaYNew(:,1))/pixel*PhysicalLength);
    [a1,b1]=hist(x1,num);a1=a1/length(x1);
    if PlotDetail==1
    figure
    set(0,'DefaultFigureVisible','off')
    set(gca,'FontSize',18);
    plot(b1,a1,'ro','MarkerSize',10,'MarkerFaceColor','r');
    
    [m1, ind1]=max(a1);
    aa=a1;bb=b1;
    aa(ind1)=0;
    [m2, ind2]=max(aa);
    aa([ind1,ind2])=[];bb([ind1,ind2])=[];
    %%display(aa);%%display(bb);
    if size(aa,2)>1 || size(aa,1)>1
        opts=fitoptions('Method','NonlinearLeastSquares');
        opts.Lower=[-100000000000000 -1000000000000000];
    f = fit(bb.',aa.','gauss1',opts);
    %%display(f);
    va1=f(b1(ind1));va2=f(b1(ind2));
    fraction=a1(ind1)-va1+a1(ind2)-va2;
    a1(ind1)=va1;a1(ind2)=va2;
    hold on;plot(b1,f(b1),'b-','LineWidth',1);
    else
        fraction=0;
        
    end
    

    fractionY=fraction;
   
    
    DisplacementYbNew(Steplength,:)=b1;
    
    DisplacementYaNew(Steplength,:)=a1;
    DisplacementYaNew(Steplength,:)=DisplacementYaNew(Steplength,:)/sum(DisplacementYaNew(Steplength,:));
    
    hold on;
    axis([b1(1,1),b1(1,size(b1,2)), 0.0001, m1+2])
    
    a=get(gca);
    xxx=a.XLim;%
yyy=a.YLim;%
k=[0.01 0.5];%
x0=xxx(1)+k(1)*(xxx(2)-xxx(1));%
y0=yyy(1)+k(2)*(yyy(2)-yyy(1));%
set(gca,'FontSize',18);
    text(x0,y0,['#disp=',num2str(size(deltaY,1)),', \deltat=',num2str(Steplength*tLength),'s','\newlineFraction=',num2str(roundn(fraction,-2))],'FontSize',16);
    xlabel('\Deltay=y(t+\deltat)-y(t) (\mum)');
    ylabel('Ratio');
    title('Y-displacement')
    if Steplength==1 

    set(gca,'YScale','log');

    figurename=[dirnameout,'\',filefolder,'Step',num2str(Steplength),'YDistributionLog.jpg'];
    saveas(gcf,figurename)
    end

    end
   
    deltaR=deltaR';
    x1=sort((deltaR(:,1))/pixel*PhysicalLength); 
    [a1,b1]=hist(x1,num);a1=a1/length(x1);
    DisplacementRb(Steplength,:)=b1;
    DisplacementRa(Steplength,:)=a1;
    if PlotDetail==1
    figure
    set(0,'DefaultFigureVisible','off')
    set(gca,'FontSize',18);
    plot(b1,a1,'ro','MarkerSize',10,'MarkerFaceColor','r');
    
    
    if Steplength==1
    DataDisplacementR=[b1',a1'];%
    end
   
    [m1, ind1]=max(a1);
    aa=a1;bb=b1;
    aa(ind1)=0;
    [m2, ind2]=max(aa);
    aa([ind1,ind2])=[];bb([ind1,ind2])=[];
    if size(aa,2)>1 || size(aa,1)>1
    opts=fitoptions('Method','NonlinearLeastSquares');
        opts.Lower=[-100000000000000 -1000000000000000];
    f = fit(bb.',aa.','gauss1',opts);
    va1=f(b1(ind1));va2=f(b1(ind2));
    fraction=a1(ind1)-va1+a1(ind2)-va2;
    a1(ind1)=va1;a1(ind2)=va2;
    hold on;plot(b1,f(b1),'b-','LineWidth',1);
    else
        fraction=0;
        
    end

    fractionR=fraction;

    DisplacementRbNew(Steplength,:)=b1;
    DisplacementRaNew(Steplength,:)=a1;
    DisplacementRaNew(Steplength,:)=DisplacementRaNew(Steplength,:)/sum(DisplacementRaNew(Steplength,:));
    
    
    hold on;
    
    axis([b1(1,1),b1(1,size(b1,2)), 0.0001, m1+2])
    
    set(gca,'YScale','log');
    a=get(gca);
    xxx=a.XLim;%
yyy=a.YLim;%
k=[0.5 0.3];%
x0=xxx(1)+k(1)*(xxx(2)-xxx(1));%
y0=yyy(1)+k(2)*(yyy(2)-yyy(1));%
set(gca,'FontSize',18);
   
    text(x0,y0,['\deltat=',num2str(Steplength*tLength),'s','\newlineFraction=',num2str(roundn(fraction,-2))],'FontSize',18);
    xlabel('\Deltar=r(t+\deltat)-r(t) (\mum)');
    ylabel('Ratio');
    title('Displacement distribution')
    if Steplength==1

    figurename=[dirnameout,'\',filefolder,'Step',num2str(Steplength),'RDistribution.jpg'];
    saveas(gcf,figurename)
    end

    end

    if Steplength==1% || Steplength==MaxSteplength
if size(x,2)<3
    deltaAngle=0;
else
Angles=NaN(size(x,1),size(x,2)-2);
end

for j=1:size(x,2)-Steplength-1
   for i=1:size(x,1)
        if ~isnan(DisplacementX(i,j))  && ~isnan(DisplacementY(i,j)) && ~isnan(DisplacementX(i,j+Steplength)) && ~isnan(DisplacementY(i,j+Steplength))
            A=[DisplacementX(i,j),DisplacementY(i,j)];
            B=[DisplacementX(i,j+Steplength),DisplacementY(i,j+Steplength)];
            Angles(i,j)=acos(dot(A,B)/(norm(A)*norm(B)))*180/pi;
        end
   end
    Delta=Angles(:,j);
    Delta(isnan(Delta)==1)=[];
    if j==1
        deltaAngle=Delta';
    else
    deltaAngle=[deltaAngle,Delta'];
    end
end
    deltaAngle=deltaAngle';
   % %%display(deltaAngle);
    if size(deltaAngle,2)>0 && size(deltaAngle,1)>0       
    x1=sort((deltaAngle(:,1)));
    [a1,b1]=hist(x1,num);a1=a1/length(x1);
    else
        a1=0;b1=0;
    end
    Angleb(Steplength,:)=b1;
    Anglea(Steplength,:)=a1;
    if PlotDetail==1
    figure
    bar(b1,a1);
      a=get(gca);
    xxx=a.XLim;%
yyy=a.YLim;%
k=[0.2 0.9];%
x0=xxx(1)+k(1)*(xxx(2)-xxx(1));%
y0=yyy(1)+k(2)*(yyy(2)-yyy(1));%
set(gca,'FontSize',18);
    text(x0,y0,['#displacements=',num2str(size(deltaAngle,1)),'\newline\deltat=',num2str(Steplength*tLength),'s'],'FontSize',18);
    xlabel('\Delta\theta=\theta(t+\deltat)-\theta(t)');
    ylabel('Ratio');
    title('Angular distribution')
    
    figurename=[dirnameout,'\',filefolder,'Step',num2str(Steplength),'AngleDistribution.jpg'];
    saveas(gcf,figurename)
    end

    end
end
%%
%Plot X-direction's mean displacement versus time and mean square displacement versus time
MeanDisplacementX=sum(DisplacementXb.*DisplacementXa,2)';
DiffusionCoeffcientX=sum(DisplacementXb.^2.*DisplacementXa,2)'-MeanDisplacementX.^2;

MeanDisplacementXNew=sum(DisplacementXbNew.*DisplacementXaNew,2)';
DiffusionCoeffcientXNew=sum(DisplacementXbNew.^2.*DisplacementXaNew,2)'-MeanDisplacementXNew.^2;
DriftXTotal(window,1:MaxSteplength)=MeanDisplacementXNew(1,1:MaxSteplength);

if PlotDriftDiffusionFit==1 && window==1
figure
xplot=(0:MaxSteplength)*tLength;
yplot=[0,MeanDisplacementX(1,:)];
p1=plot(xplot,yplot,'ro','MarkerSize',10,'MarkerFaceColor','r'); hold on;
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
P1=P(1);P2=P(2);
p2=plot(xplot,yfit,'r-','LineWidth',1);
yplot=[0,MeanDisplacementXNew(1,:)];
p3=plot(xplot,yplot,'b>','MarkerSize',10,'MarkerFaceColor','b'); hold on;

DataDriftX=[xplot',yplot'];%hhhhhhhhhhhhhhhhh

P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
p4=plot(xplot,yfit,'b-','LineWidth',1);
set(0,'DefaultFigureVisible','off')
set(gca,'FontSize',18);
a=get(gca);
xxx=a.XLim;%
yyy=a.YLim;%
k=[0.5 0.9];%
x0=xxx(1)+k(1)*(xxx(2)-xxx(1));%
y0=yyy(1)+k(2)*(yyy(2)-yyy(1));%

set(gca,'FontSize',18);
%axis([0,MaxSteplength*tLength, min(min(MeanDisplacementXNew(1,:)),min(MeanDisplacementX(1,:))), max(max(MeanDisplacementXNew(1,:)),max(MeanDisplacementX(1,:)))+0.5])
axis([0,MaxSteplength*tLength, -8, 8])
xlabel('\deltat (s)');
ylabel('\langlex(t+\deltat)-x(t)\rangle (\mum)');
title('First moment of X-displacement')
 h=legend('For all trajectories',['Linear fit: x=',num2str(roundn(P1,-2)),'*\deltat+(',num2str(roundn(P2,-2)),')'],'For motile',['Linear fit: x=',num2str(roundn(P(1),-2)),'*\deltat+(',num2str(roundn(P(2),-2)),')'],'Location','Best');
set(h,'FontSize',14);
figurename=[dirnameout,'\',filefolder,'DisplacementX.jpg'];
saveas(gcf,figurename)
DriftXmat=roundn(P(1),-2)/2;

figure
xplot=(0:MaxSteplength)*tLength;
yplot=[0,DiffusionCoeffcientX(1,:)];
p1=plot(xplot,yplot,'ro','MarkerSize',10,'MarkerFaceColor','r'); hold on;
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
P1=P(1);P2=P(2);
p2=plot(xplot,yfit,'r-','LineWidth',1);
yplot=[0,DiffusionCoeffcientXNew(1,:)];
p3=plot(xplot,yplot,'b>','MarkerSize',10,'MarkerFaceColor','b'); hold on;

DataDriftY=[xplot',yplot'];%hhhhhhhhhhhhhhhhh

P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
p4=plot(xplot,yfit,'b-','LineWidth',1);
set(0,'DefaultFigureVisible','off')
set(gca,'FontSize',18);
a=get(gca);
xxx=a.XLim;%
yyy=a.YLim;%
k=[0.5 0.9];%
x0=xxx(1)+k(1)*(xxx(2)-xxx(1));%
y0=yyy(1)+k(2)*(yyy(2)-yyy(1));%
set(gca,'FontSize',18);
%axis([0,MaxSteplength*tLength, -8, 8])
xlabel('\deltat (s)');
ylabel('Varianc (\mum^2)');
title('Varianc of X-displacement')
 h=legend('For all trajectories',['Linear fit: V_x=',num2str(roundn(P1,-2)),'*\deltat+(',num2str(roundn(P2,-2)),')'],'For motile',['Linear fit: V_x=',num2str(roundn(P(1),-2)),'*\deltat+(',num2str(roundn(P(2),-2)),')'],'Location','Best');
set(h,'FontSize',14);
figurename=[dirnameout,'\',filefolder,'DiffusionX.jpg'];
%saveas(gcf,figurename)


DiffusionXmat=roundn(P(1),-2)/2;
DiffusionXmatAll=roundn(P1,-2)/2;

end

%%
%Plot X-direction's mean displacement versus time and mean square displacement versus time
MeanDisplacementY=sum(DisplacementYb.*DisplacementYa,2)';
DiffusionCoeffcientY=sum(DisplacementYb.^2.*DisplacementYa,2)'-MeanDisplacementY.^2;

MeanDisplacementYNew=sum(DisplacementYbNew.*DisplacementYaNew,2)';
DiffusionCoeffcientYNew=sum(DisplacementYbNew.^2.*DisplacementYaNew,2)'-MeanDisplacementYNew.^2;
DriftYTotal(window,1:MaxSteplength)=MeanDisplacementYNew(1,1:MaxSteplength);

if PlotDriftDiffusionFit==1 && window==1
    figure
xplot=(0:MaxSteplength)*tLength;
yplot=[0,MeanDisplacementY(1,:)];
p1=plot(xplot,yplot,'ro','MarkerSize',10,'MarkerFaceColor','r'); hold on;
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
P1=P(1);P2=P(2);
p2=plot(xplot,yfit,'r-','LineWidth',1);
% set(p1,'MarkerSize',10);set(p3,'MarkerSize',10);
% set(p2,'LineWidth',2);set(p4,'LineWidth',2);
yplot=[0,MeanDisplacementYNew(1,:)];
p3=plot(xplot,yplot,'b>','MarkerSize',10,'MarkerFaceColor','b'); hold on;
%yplot=[0,DiffusionCoeffcientNew(1,:)];
%p2=plot(xplot,yplot,'b>','MarkerSize',10,'MarkerFaceColor','b');
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
p4=plot(xplot,yfit,'b-','LineWidth',1);
set(0,'DefaultFigureVisible','off')
set(gca,'FontSize',18);
a=get(gca);
    xxx=a.XLim;%
yyy=a.YLim;%
k=[0.5 0.9];%
x0=xxx(1)+k(1)*(xxx(2)-xxx(1));%
y0=yyy(1)+k(2)*(yyy(2)-yyy(1));%
%text(x0,y0,['For motile: y=',num2str(roundn(P1,-2)),'*\deltat+(',num2str(roundn(P2,-2)),')','\newlineFor all trajectories: y=',num2str(roundn(P(1),-2)),'*\deltat+(',num2str(roundn(P(2),-2)),')'],'FontSize',18);
set(gca,'FontSize',18);
%axis([0,MaxSteplength*tLength, min(min(MeanDisplacementYNew(1,:)),min(MeanDisplacementY(1,:))), max(max(MeanDisplacementYNew(1,:)),max(MeanDisplacementY(1,:)))+0.5])
axis([0,MaxSteplength*tLength, -8, 8])
xlabel('\deltat (s)');
ylabel('\langley(t+\deltat)-y(t)\rangle (\mum)');
title('First moment of Y-displacement')
 h=legend('For all trajectories',['Linear fit: y=',num2str(roundn(P1,-2)),'*\deltat+(',num2str(roundn(P2,-2)),')',],'For motile',['Linear fit: y=',num2str(roundn(P(1),-2)),'*\deltat+(',num2str(roundn(P(2),-2)),')'],'Location','Best');
set(h,'FontSize',14);figurename=[dirnameout,'\',filefolder,'DisplacementY.pdf'];
%saveas(gcf,figurename)
figurename=[dirnameout,'\',filefolder,'DisplacementY.jpg'];
saveas(gcf,figurename)


DriftYmat=roundn(P(1),-2)/2;


figure
xplot=(0:MaxSteplength)*tLength;
yplot=[0,DiffusionCoeffcientY(1,:)];
p1=plot(xplot,yplot,'ro','MarkerSize',10,'MarkerFaceColor','r'); hold on;
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
P1=P(1);P2=P(2);
p2=plot(xplot,yfit,'r-','LineWidth',1);
yplot=[0,DiffusionCoeffcientYNew(1,:)];
p3=plot(xplot,yplot,'b>','MarkerSize',10,'MarkerFaceColor','b'); hold on;
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
p4=plot(xplot,yfit,'b-','LineWidth',1);
set(0,'DefaultFigureVisible','off')
set(gca,'FontSize',18);
a=get(gca);
xxx=a.XLim;%
yyy=a.YLim;%
k=[0.5 0.9];%
x0=xxx(1)+k(1)*(xxx(2)-xxx(1));%
y0=yyy(1)+k(2)*(yyy(2)-yyy(1));%
set(gca,'FontSize',18);
%axis([0,MaxSteplength*tLength, -8, 8])
xlabel('\deltat (s)');
ylabel('Varianc (\mum^2)');
title('Varianc of Y-displacement')
 h=legend('For all trajectories',['Linear fit: V_y=',num2str(roundn(P1,-2)),'*\deltat+(',num2str(roundn(P2,-2)),')'],'For motile',['Linear fit: V_y=',num2str(roundn(P(1),-2)),'*\deltat+(',num2str(roundn(P(2),-2)),')'],'Location','Best');
set(h,'FontSize',14);
figurename=[dirnameout,'\',filefolder,'DiffusionY.jpg'];
%saveas(gcf,figurename)



DiffusionYmat=roundn(P(1),-2)/2;
DiffusionYmatAll=roundn(P1,-2)/2;
end




%MeanSquare and Diffusion
MeanDisplacement=sqrt(MeanDisplacementXNew.^2+MeanDisplacementYNew.^2);
MeanDisplacementOld=sqrt(MeanDisplacementX.^2+MeanDisplacementY.^2);
MeanSquare=sum(DisplacementRb.^2.*DisplacementRa,2)';
DiffusionCoeffcient=sum(DisplacementRb.^2.*DisplacementRa,2)'-MeanDisplacementOld.^2;
%DiffusionCoeffcient=MeanSquare+MeanDisplacement.^2;


%MeanDisplacementNew=sum(DisplacementRb.*DisplacementRaNew,2)';
MeanSquareNew=sum(DisplacementRb.^2.*DisplacementRaNew,2)';
DiffusionCoeffcientNew=sum(DisplacementRb.^2.*DisplacementRaNew,2)'-MeanDisplacement.^2;
DiffusionTotal(window,1:MaxSteplength)=DiffusionCoeffcientNew(1,1:MaxSteplength);

if PlotDriftDiffusionFit==1 && window==1
 
%Diffusion-----------------------
figure
xplot=(0:MaxSteplength)*tLength;
yplot=[0,DiffusionCoeffcient(1,:)];
p1=plot(xplot,yplot,'ro','MarkerSize',10,'MarkerFaceColor','r'); hold on;
set(gca,'FontSize',18);
P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
P1=P(1);P2=P(2);
p2=plot(xplot,yfit,'r-','LineWidth',1);
yplot=[0,DiffusionCoeffcientNew(1,:)];
p3=plot(xplot,yplot,'b>','MarkerSize',10,'MarkerFaceColor','b');

DataDiffusionR=[xplot',yplot'];%hhhhhhhhhhhhhhhhh



P = polyfit(xplot,yplot,1);
yfit = P(1)*xplot+P(2);
p4=plot(xplot,yfit,'b-','LineWidth',1);
%text(max(xplot)*0.05,max(yplot)*0.6,['all trajectories: V_r=',num2str(roundn(P1,-2)),'*\deltat+(',num2str(roundn(P2,-2)),'\newlinemotile: V_r=',num2str(roundn(P(1),-2)),'*\deltat+(',num2str(roundn(P(2),-2)),')'],'FontSize',18);
xlabel('\deltat (s)');
%ylabel('V_r=\langle[r(t+\deltat)-r(t)]^2-\langler(t+\deltat)-r(t)\rangle^2\rangle (\mum^2)');
logy=set(gca,'YScale','log');
logy=set(gca,'XScale','log');
axis([min(xplot), max(xplot), min(yplot), max(yplot)+300]);
ylabel('Variance of displacement')
set(gca,'FontSize',18);
 h=legend('For all trajectories',['Linear fit: V_r=',num2str(roundn(P1,-2)),'*\deltat+(',num2str(roundn(P2,-2)),')'],'For trajecories of motile',['Linear fit: V_r=',num2str(roundn(P(1),-2)),'*\deltat+(',num2str(roundn(P(2),-2)),')'],'Location','best');
set(h,'FontSize',14);figurename=[dirnameout,'\',filefolder,'Diffusion.pdf'];
%saveas(gcf,figurename)
figurename=[dirnameout,'\',filefolder,'Diffusion.jpg'];
saveas(gcf,figurename)
% figurename=[dirnameout,'\',filefolder,'Diffusion.fig'];
% saveas(gcf,figurename)

DiffusionRmat=roundn(P(1),-2)/2;
DiffusionRmatAll=roundn(P1,-2)/2;


end

%% 
%Angular-distributed displacement: Use non-filter traj (both the above and here need to use all 10-step for enough statistics)
if size(x,2)>1
for Steplength=[Howlongtolook]
NumOfAngle=60;
NumOfTraj=zeros(NumOfAngle,1);
StepLengthInAngle=cell(NumOfAngle,1);StepLengthInAnglePlot=zeros(NumOfAngle,1);
AngleDistri=[-180:360/NumOfAngle:180-360/NumOfAngle];
%AngleDistriPlot=[-180:360/NumOfAngle:180-360/NumOfAngle];
for i=1:size(y,1)
    for j=1:size(y,2)-Steplength
        if ~isnan(y(i,j)) && ~isnan(x(i,j)) &&  ~isnan(y(i,j+Steplength)) && ~isnan(x(i,j+Steplength))
            AngleTemp=rad2deg(atan2((y(i,j+Steplength)-y(i,j)),(x(i,j+Steplength)-x(i,j))));
            indexTemp=find(abs(AngleDistri-AngleTemp)==min(abs(AngleDistri-AngleTemp)));
            %%display(AngleTemp);%display(indexTemp);
            if size(indexTemp,2)>1
                indexTemp=indexTemp(1,1);
            end
            if AngleTemp-AngleDistri(1,indexTemp)<0
                indexTemp=indexTemp-1;
            end
            %NumOfTraj(indexTemp,1)=NumOfTraj(indexTemp,1)+1;
            StepLengthInAngle{indexTemp,1}=[StepLengthInAngle{indexTemp,1},sqrt((y(i,j+Steplength)-y(i,j))^2+(x(i,j+Steplength)-x(i,j))^2)/pixel*PhysicalLength];
        end
    end
end
for i=1:NumOfAngle
    NumOfTraj(i,1)=size(StepLengthInAngle{i,1},2);
    if NumOfTraj(i,1)==0
        StepLengthInAnglePlot(i,1)=0;
    else
    StepLengthInAnglePlot(i,1)=sum(StepLengthInAngle{i,1},2)/NumOfTraj(i,1);
    end
end



AngleDistri=[AngleDistri/360*2*pi, pi];
NumOfTrajNew=[NumOfTraj',NumOfTraj(1,1)];
StepLengthInAngleNew=[StepLengthInAnglePlot',StepLengthInAnglePlot(1,1)];
StepLengthInAngleSumNew=[StepLengthInAnglePlot'.*NumOfTraj',StepLengthInAnglePlot(1,1).*NumOfTraj(1,1)];

figure
polarplot(AngleDistri,NumOfTrajNew,'LineWidth',1.5);hold on;
ax = gca;ax.LineWidth = 1.5;ax.ThetaTick = [0 45 90 135 180 225 270 315];  ax.RAxisLocation = 90;%ax.RAxisLocation='manual';
title('# of displacement');
set(gca,'FontSize',15);
figurename=[dirnameout,'\',filefolder,'Step',num2str(Steplength),'NumAngluar.jpeg'];
hold off; print(gcf, '-djpeg', '-r600',figurename)


figure
polarplot(AngleDistri,StepLengthInAngleNew,'LineWidth',1.5);hold on;
ax = gca;ax.LineWidth = 1.5;ax.ThetaTick = [0 45 90 135 180 225 270 315];ax.RAxisLocation = 90;
title('Average displacement (\mum)');
set(gca,'FontSize',15);
figurename=[dirnameout,'\',filefolder,'Step',num2str(Steplength),'LengthAngluar.jpeg'];
hold off; print(gcf, '-djpeg', '-r600',figurename)


figure
polarplot(AngleDistri,StepLengthInAngleSumNew,'LineWidth',1.5);hold on;
ax = gca;ax.LineWidth = 1.5;ax.ThetaTick = [0 45 90 135 180 225 270 315];ax.RAxisLocation = 90;
title('Sum of displacement (\mum)');
set(gca,'FontSize',15);
figurename=[dirnameout,'\',filefolder,'Step',num2str(Steplength),'TotalLengthAngluar.jpeg'];
hold off; print(gcf, '-djpeg', '-r600',figurename)


QuadrantXNeg=[[1:1:ceil(NumOfAngle/8)-1],[ceil(NumOfAngle/8*7):1:NumOfAngle]];
QuadrantXPos=ceil(NumOfAngle/8*3):1:ceil(NumOfAngle/8*5)-1;
QuadrantYNeg=ceil(NumOfAngle/8):1:ceil(NumOfAngle/8*3)-1;
QuadrantYPos=ceil(NumOfAngle/8*5):1:ceil(NumOfAngle/8*7)-1;
    
YPositiveProOrigin=sum(NumOfTrajNew(1,QuadrantYPos),2)/(sum(NumOfTrajNew(1,QuadrantYPos),2)+sum(NumOfTrajNew(1,QuadrantYNeg),2));
XPositiveProOrigin=sum(NumOfTrajNew(1,QuadrantXPos),2)/(sum(NumOfTrajNew(1,QuadrantXPos),2)+sum(NumOfTrajNew(1,QuadrantXNeg),2));
XPositiveLengthOrigin=sum(StepLengthInAngleNew(1,QuadrantXPos),2)/size(QuadrantXPos,2);
XNegativeLengthOrigin=sum(StepLengthInAngleNew(1,QuadrantXNeg),2)/size(QuadrantXNeg,2);
YPositiveLengthOrigin=sum(StepLengthInAngleNew(1,QuadrantYPos),2)/size(QuadrantYPos,2);
YNegativeLengthOrigin=sum(StepLengthInAngleNew(1,QuadrantYNeg),2)/size(QuadrantYNeg,2);


end
end

%%
%Save calculated statistics to mat files
% if savematV==1
 filenamenew=[dirnameout,'AllTrajectory.mat'];   
 save(filenamenew,'xAll','yAll');
 filenamenew=[dirnameout,'LongTrajectory.mat'];   
 save(filenamenew,'x','y');

 filenamenew=[diffusionout,'DisplacmentWOFilter.mat'];   
 save(filenamenew,'DisplacementXb','DisplacementXa','DisplacementYb','DisplacementYa');
  filenamenew=[diffusionout,'DisplacmentWithFilter.mat'];   
 save(filenamenew,'DisplacementXbNew2','DisplacementXaNew2','DisplacementYbNew2','DisplacementYaNew2','DiffusionRx','DiffusionRy');
 
    
filename=[diffusionout,'.mat'];
save(filename,'DriftXmat','DriftYmat','DiffusionXmat','DiffusionXmatAll','DiffusionYmat','DiffusionYmatAll',...
    'DiffusionRmat','DiffusionRmatAll','fractionX','fractionY','fractionR','fractionDrop','VelocityMat',...
    'VelocityMatBall','VelocityMatRemain','XPositiveProTotal','YPositiveProTotal',...
    'XPositivePro','YPositivePro','XPositiveLength','XNegativeLength','YPositiveLength','YNegativeLength',...
    'XPositiveProOrigin','YPositiveProOrigin','XPositiveLengthOrigin','XNegativeLengthOrigin','YPositiveLengthOrigin','YNegativeLengthOrigin',...
    'XPositiveLengthTotal','XNegativeLengthTotal','YPositiveLengthTotal','YNegativeLengthTotal',...
    'Ball_DriftXmat','Ball_DriftXmat_Remain','Ball_DriftYmat','Ball_DriftYmat_Remain','Ball_DiffusionRmat',...
    'Ball_DiffusionRmat_Remain');
end

close all;
end

clear
end