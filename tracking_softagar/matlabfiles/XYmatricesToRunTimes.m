function XYmatricesToRunTimes(filex,filey,dirnameout,filefolder,Sc,Cc,Tc)
folderRes=filefolder;
filenameinNew=[dirnameout,'\clusterdata'];
filename=dirnameout;

disp(filex)
load(filex)
disp(filey)
load(filey)
PhysicalLength=192.36;%mum
pixel=512;
tLength=0.068; %0.13 for old data %0.068;%second
changeWindowWithoutBall=0;
changeWindow=0;
Ballfilter=0;
%folder='C:\Python\20170307_Tomo\2017-03-07\WT_glu_OD011_1-Res';
% folder='/Volumes/DIVE_FAT/2016-11-30/Trajectories/';
folder=folderRes;
TimeWindow=10;
a=2;% parameter to differentiate run and tumble
b=2;

%%
%Pick traj with length---------
if changeWindowWithoutBall==1

NumLength=1*TimeWindow;
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


%%
if Ballfilter==1 %Ball filter
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
x=xnew;y=ynew;

end




%%
%cd(folder)

% filenameRadical='Cluster_WT';
% list=dir;
% 
% %display(list);
% 
% mask = ones(1,length(list));
% for i=1:length(list)
%     if isempty(findstr(filenameRadical,list(i).name))
%         mask(i)=0;
%     end
% end
% 
% list2=list(find(mask));
% 
% %display(list2(1).name);display(list2(2).name);display(list2(3).name);
% clear('list');



    
 
%% Explorating Trajerctories
dt=tLength;
% thresh_t=1/dt;
% thresh_gyr=10;
% thresh_runfrac_lb=.7;
% thresh_runfrac_ub=1;
% fracToPlot=.3;
% thresh_tumble_angle=60;
% thresh_tumble_speed=-.7;


 
showtotalTraj=1;
showOnebyOne=0;
% for nn=1:length(list)
for nn=1:1
    
    

    % for nn=length(list)-3:length(list)
    % for nn=1:1
   
    
     Matrix_x=x;
    Matrix_y=y;
    for j=1:length(Matrix_x(:,1))
        mask=~isnan(Matrix_x(j,:));
        index=find(mask);
        temp=zeros(length(index),6);
        temp(1:length(index),2)=index;
        temp(1:length(index),3)=Matrix_x(j,index);
        temp(1:length(index),4)=Matrix_y(j,index);
        
        Traj(j).m =temp;
        
    end
    
    
    muparpix = PhysicalLength/pixel; %camera 1322M
    
    
    MedianRuntime_matrix=zeros(6,5,9,5,5);
    MedianSpeed_matrix=zeros(6,5,9,5,5);
    NonMotileFrac_matrix=zeros(6,5,9,5,5);
    
    iteration=0;
    
    for idx_t=1:3
        thresh_t=(idx_t-1)/dt;
        for idx_gyr=1:3
            thresh_gyr=2+2*(idx_gyr-1);% in mum unit
            for idx_runfrac_lb=1:9
                thresh_runfrac_lb=idx_runfrac_lb*.1;
                for idx_tumble_angle=1:5
                    thresh_tumble_angle=30+idx_tumble_angle*10;
                    for idx_tumble_speed=1:5
                        iteration=iteration+1;
                        
                        thresh_tumble_speed=-1+idx_tumble_speed*.1;
                        %     thresh_t=1/dt;
                        % thresh_gyr=10;
                        % thresh_runfrac_lb=.7;
                        thresh_runfrac_ub=1;
                        fracToPlot=.1;
                        % thresh_tumble_angle=60;
                        % thresh_tumble_speed=-.7;
                        
                        
                        
                       % muparpix = 0.4682963379; %camera 1322M
                        
                        dr_tot=[];
                        a=2;
                        b=2;
                        dt=tLength; %for confocal images
                        v_tot=[];
                        runsTot=[];
                        Vmeantot=[];
                        dthetatot=[];
                        Angles_Ying=[];
                        Angles_Ying_tot=[];
                        TrajLength=[];
                        V_hist=[];
                        Nr_tot=[];
                        Nrt_tot=[];
                        NewTraj=[];
                        RunTimeTot=[];
                        sdTot=[];
                        flag=0;
                        nm=0;
                        nmbase=0;
                        for i=1:length(Traj)
                           % xNew=smooth(Traj(i).m*muparpix,3);
                            %yNew=smooth(Traj(i).m*muparpix,4);
                            xNew=smooth(Traj(i).m(:,3)*muparpix,3);%whether this smooth function should for (x,y)???
                            yNew=smooth(Traj(i).m(:,4)*muparpix,3);
                            f=Traj(i).m*dt;
                            dx=xNew(2:length(xNew))-xNew(1:length(xNew)-1);%diaplacement vector
                            dy=yNew(2:length(yNew))-yNew(1:length(yNew)-1);
                            theta=zeros(1,length(dx));
                            for j=1:length(theta) %angle in [0,2pi]
                                if dx(j)<0
                                    theta(j)=pi+atan(dy(j)/dx(j));
                                end
                                if dx(j)==0
                                    if dy(j)==0
                                        theta(j)=pi/2;
                                    else
                                        theta(j)=3*pi/2;
                                    end
                                end
                                if dx(j)>0
                                    if dy(j)>0
                                        theta(j)=atan(dy(j)/dx(j));
                                    elseif dy(j)<0
                                        theta(j)=2*pi+atan(dy(j)/dx(j));
                                    else
                                        theta(j)=0;
                                    end
                                end
                                
                                
                            end
                            
                            dtheta=theta(b+1:length(theta))-theta(1:length(theta)-b);%angle changes between vectors with displacement b
                            idx1=find(dtheta>pi); %absolute angle difference [0,\pi]
                            idx2=find(dtheta<-pi);
                            dtheta(idx1)=dtheta(idx1)-pi;
                            dtheta(idx2)=dtheta(idx2)+pi;
                            dtheta_c=pi-abs(dtheta);
                            dtheta=min(abs(dtheta),pi-abs(dtheta));% use arc-cos is better???
                            
                            
                            
                            dr=sqrt(dx.*dx+dy.*dy);
                            v=dr/(dt);
                            dv=v(b+1:length(v))-v(1:length(v)-b);%velocity changes between vectors with displacement b
                            vtest=dv(:)./v(1:length(dv));%relative velocity changes
                            
                            tumbleAngle=find(abs(dtheta)>=thresh_tumble_angle/180*pi);
                            tumbleSpeed=find(vtest<thresh_tumble_speed);
                            tumble=[tumbleAngle(:); tumbleSpeed(:)];
                            tumble=unique(sort(tumble));
                            
                            run_idx1=find(abs(dtheta)<=thresh_tumble_angle/180*pi);
                            run_idx2=find(vtest>=thresh_tumble_speed);
                            run_idx=intersect(run_idx1,run_idx2);%satisfy both two conditions
                           % showOnebyOne=0;
                            
                            if showOnebyOne==1
                                figure(3)
                                subplot(1,2,1)
                                
                                plot(xNew,yNew);
                                hold on
                                plot(xNew(1),yNew(1),'or','MarkerSize',10);
                                
                                plot(xNew(tumbleAngle+b/2),yNew(tumbleAngle+b/2),'or');
                                plot(xNew(tumble+b/2),yNew(tumble+b/2),'xr');
                                hold off
                                theta_bins=[-180:10:180];
                                n_theta=hist(dtheta*180/pi,theta_bins);
                                subplot(1,2,2)
                                plot(theta);
                                plot(theta_bins,n_theta/sum(n_theta));
                                pause
                                figurename2=[filenameinNew,'Do1.jpeg'];
                                print(gcf, '-djpeg', '-r800',figurename2)
                            end
                            
                            dtumbletemp=diff(tumble);
                            Nr=length(xNew)-a-length(tumble);
                            Nrt=length(find(dtumbletemp>1));
                            idxtumb=find(dtumbletemp>1);
                            
                            Nr_tot=[Nr_tot Nr];
                            Nrt_tot=[Nrt_tot Nrt];
                            RunTime=dt*sum(Nr)/sum(Nrt);
                            x1=diff(xNew);
                            y1=diff(yNew);
                            sd=x1'*x1+y1'*y1;
                            
                            %         if ~isempty(idxtumb)
                            %                 dtumble=dtumbletemp(idxtumb);
                            %                 runsTot=[runsTot ; dtumble];
                            %                 v_tot=[v_tot; v];
                            %         end
                            gyration=sqrt((max(xNew)-min(xNew)).^2+(max(yNew)-min(yNew)).^2);%want to box filter
                            
                            
                            if(length(xNew)>thresh_t)
                                nmbase=nmbase+1;
                                if (gyration<thresh_gyr)
                                    nm=nm+1;
                                else
                                    if (Nr/length(xNew)>thresh_runfrac_lb)&&(Nr/length(xNew)<thresh_runfrac_ub) % filtering trajectories of certain run ratio
                                        flag=flag+1;
                                        dtumble=dtumbletemp(idxtumb);
                                        runsTot=[runsTot ; dtumble];%index of run in xNew, what's this used for???
                                        v_tot=[v_tot; mean(v(run_idx+b/2))];%why use +b/2: index moved
                                        RunTimeTot=[RunTimeTot RunTime];
                                        sdTot=[sdTot sd];
                                        TrajLength=[TrajLength length(xNew)];
                                        %         else
                                        if showtotalTraj==1 && iteration==1
                                            
                                            if rand<fracToPlot % This plots a fraction of the trajectories, all in the same graph
                                                
                                                figure(2)
                                                hold on
                                                plot(xNew,yNew,'-b');
                                                %                 plot(xNew(tumbleSpeed+b/2),yNew(tumbleSpeed+b/2),'or');
                                                %                 plot(xNew(tumbleAngle+b/2),yNew(tumbleAngle+b/2),'og');
                                                plot(xNew(tumble+b/2),yNew(tumble+b/2),'or','MarkerSize',0.8);
                                                hold off
                                               
                                                figurename2=[filenameinNew,'Tumble.jpeg'];
                                                print(gcf, '-djpeg', '-r800',figurename2)
                                                figurename2=[filenameinNew,'Tumble.fig'];
                                                saveas(gcf,figurename2)
                                            end
                                            
                                        end
                                        
                                        
                                        
                                    end
                                end
                            end
                        end
                      
    
    %     x_bins=[0:.1:20];
    %     n=hist(TrajLength*dt,x_bins);
    %     n2=hist(RunTimeTot,x_bins);
    %
    %     figure(1)
    %     subplot(2,2,1)
    %     plot(x_bins,n)
    %     subplot(2,2,2)
    %     plot(x_bins,n2)
    %     subplot(2,2,3)
    %     plot(RunTimeTot,TrajLength*dt,'.')
    %     subplot(2,2,4)
    %     plot(RunTimeTot,sdTot,'.')
    %
    % %       pause
    %
    %     x_bins=[0:.1:5];
    %     n=hist(runsTot*dt,x_bins);
    %     s=n;
    %     for kk=1:length(s)
    %         s(kk)=sum(n(1:kk));
    %     end
    %     x_tauToFit=x_bins;
    %     tauToFit=1-s/sum(n);
    %     options1 = optimset('MaxFunEvals',2000, 'Display','final','TolX',1e-6, 'TolFun',1e-6);
    %     options1 = optimset('Display','final');
    %     x_init=[.001; 1];
    %     x_lb=[.001,0.01];
    %     x_ub=[10, 10];
    %     paraminit=[.15; 1];
    %     FunctionToFit= @(param,x_data) param(1)*exp(-x_data/param(2));
    %     paramfound=fminsearch(FunctionToFit,paraminit,x_tauToFit(10:end),tauToFit(10:end),x_lb,x_ub,options1);
    %     FunctionToMinimize= @(T) max((tauToFit(3:end)-T(1)*exp(-x_tauToFit(3:end)/T(2))).^2);
    %    FunctionToMinimize= @(T) max((n(3:end)-T(1)*exp(-x_tauToFit(3:end)/T(2))).^2./(T(1)*exp(-x_tauToFit(3:end)/T(2))).^2);
    
    %     (tauToFit(3:end)-T(1)*exp(-x_tauToFit(3:end)/T(2))).^2
    %     (tauToFit(3:end)-T(1)*exp(-x_tauToFit(3:end)/T(2))).^2./(T(1)*exp(-x_tauToFit(3:end)/T(2))).^2
    %     (T(1)*exp(-x_tauToFit(3:end)/T(2))).^2
    %     [x0,fval,exitFlag,output]=fminsearch(FunctionToMinimize,x_init,options1);
    %     paramfound=fminsearch(FunctionToFit,paraminit);
    %     paramfound=lsqcurvefit(FunctionToFit,paraminit,x_tauToFit(3:end),tauToFit(3:end),x_lb,x_ub,options1);
    
    
    %     figure(4)
    %     hold on
    %     semilogy(x_tauToFit,tauToFit,'k')
    % %     semilogy(x_tauToFit(10:end),paramfound(1)*exp(-x_tauToFit(10:end)/paramfound(2)),'g')
    %     semilogy(x_tauToFit(3:end),x0(1)*exp(-x_tauToFit(3:end)/x0(2)),'g')
    %
    %     hold off
    %     figure(5)
    
    
    disp([filename ' - MedianRunTime = ' num2str(median(RunTimeTot)) ' on ' num2str(flag) ' trajectories'])
    disp([filename ' - MedianSpeed = ' num2str(median(v_tot)) ' on ' num2str(flag) ' trajectories'])
    disp([filename ' - NonMotileFraction = ' num2str(nm/nmbase) ' on ' num2str(nmbase) ' trajectories'])
    disp([''])
    
    
    
    MedianRuntime_matrix(idx_t,idx_gyr,idx_runfrac_lb,idx_tumble_angle,idx_tumble_speed)=median(RunTimeTot);
%Tarj>frames; box filter; runfraction; angular change; speed decrease ratio  
    MedianSpeed_matrix(idx_t,idx_gyr,idx_runfrac_lb,idx_tumble_angle,idx_tumble_speed)=(median(v_tot));
    NonMotileFrac_matrix(idx_t,idx_gyr,idx_runfrac_lb,idx_tumble_angle,idx_tumble_speed)=nm/nmbase;
    
       
                    end
                    
                    
                    
                end
            end
        end
        
    end     
        
   
    
    %     Mean_RuntimeToPlot(nn)=(median(RunTimeTot));
    %     MeanSpeedToPlot(nn)=mean(v_tot);
    %     disp([filename ' - Fit Mean = ' num2str(x0(2))])
    %     disp([filename ' - Pooled Mean = ' num2str(dt*sum(Nr_tot)/sum(Nrt_tot))])
     figurename2=[filenameinNew,'Data.mat'];
    save(figurename2,'MedianRuntime_matrix','MedianSpeed_matrix','NonMotileFrac_matrix');
end


% for nn=1:1
% % for nn=1:1
%     %filename=list(nn).name;
%     Matrix_x=xNew;
%     Matrix_y=yNew;
%     for j=1:length(Matrix_x(:,1))
%         mask=~isnan(Matrix_x(j,:));
%         index=find(mask);
%         temp=zeros(length(index),6);
%         temp(1:length(index),2)=index;
%         temp(1:length(index),3)=Matrix_x(j,index);
%         temp(1:length(index),4)=Matrix_y(j,index);
%         
%         Traj(j).m =temp;
%         
%     end
%     
%     
%     muparpix = PhysicalLength/pixel; %camera 1322M
%     
%     dr_tot=[];
%     dt=tLength; %for confocal images
%     v_tot=[];
%     runsTot=[];
%     Vmeantot=[];
%     dthetatot=[];
%     Angles_Ying=[];
%     Angles_Ying_tot=[];
%     TrajLength=[];
%     V_hist=[];
%     Nr_tot=[];
%     Nrt_tot=[];
%     NewTraj=[];
%     RunTimeTot=[];
%     sdTot=[];
% flag=0;
%     for i=1:length(Traj)
% %         xNew=smooth(Traj(i).m*muparpix,3);
% %         yNew=smooth(Traj(i).m*muparpix,3);
% %         f=Traj(i).m*dt;
%         xNew=smooth(Traj(i).m(:,3)*muparpix,3);
%         yNew=smooth(Traj(i).m(:,4)*muparpix,3);
%         f=Traj(i).m(:,2)*dt;
%         dx=xNew(2:length(xNew))-xNew(1:length(xNew)-1);
%         dy=yNew(2:length(yNew))-yNew(1:length(yNew)-1);
%         theta=zeros(1,length(dx));
%         for j=1:length(theta)
%             if dx(j)<0
%                 theta(j)=pi+atan(dy(j)/dx(j));
%             end
%             if dx(j)==0
%                 if dy(j)==0
%                     theta(j)=pi/2;
%                 else
%                     theta(j)=3*pi/2;
%                 end
%             end
%             if dx(j)>0
%                 if dy(j)>0
%                     theta(j)=atan(dy(j)/dx(j));
%                 elseif dy(j)<0
%                     theta(j)=2*pi+atan(dy(j)/dx(j));
%                 else
%                     theta(j)=0;
%                 end
%             end
%             
%             
%         end
%         
%         dtheta=theta(b+1:length(theta))-theta(1:length(theta)-b);
%         idx1=find(dtheta>pi);
%         idx2=find(dtheta<-pi);
%         dtheta(idx1)=dtheta(idx1)-pi;
%         dtheta(idx2)=dtheta(idx2)+pi;
%         dtheta_c=pi-abs(dtheta);
%         dtheta=min(abs(dtheta),pi-abs(dtheta));
%         
%         
%         
%         dr=sqrt(dx.*dx+dy.*dy);
%         v=dr/(dt*a);
%         dv=v(b+1:length(v))-v(1:length(v)-b);
%         vtest=dv(:)./v(1:length(dv));
%         
%         tumbleAngle=find(abs(dtheta)>=60/180*pi);
%         tumbleSpeed=find(vtest<-.7);
%         tumble=[tumbleAngle(:); tumbleSpeed(:)];
% 
%         tumble=unique(sort(tumble));
% 
%         
%         showOnebyOne=0;
%         if showOnebyOne==1
%             figure(3)
%             subplot(1,2,1)
%             
%             plot(xNew,yNew);
%             hold on
%             plot(xNew(1),yNew(1),'or','MarkerSize',10);
%             
%             plot(xNew(tumbleAngle+b/2),yNew(tumbleAngle+b/2),'or');
%             plot(xNew(tumble+b/2),yNew(tumble+b/2),'xr');
%             hold off
%             theta_bins=[-180:10:180];
%             n_theta=hist(dtheta*180/pi,theta_bins);
%             subplot(1,2,2)
%             plot(theta);
%             plot(theta_bins,n_theta/sum(n_theta));
%             pause
%             figurename2=[filenameinNew,'Do.jpeg'];
%              print(gcf, '-djpeg', '-r800',figurename2)
%         end
%         
%         dtumbletemp=diff(tumble);
%         Nr=length(xNew)-a-length(tumble);
%         Nrt=length(find(dtumbletemp>1));
%         idxtumb=find(dtumbletemp>1);
%         
%         Nr_tot=[Nr_tot Nr];
%             Nrt_tot=[Nrt_tot Nrt];
%         RunTime=dt*sum(Nr)/sum(Nrt);
%         x1=diff(xNew);
%         y1=diff(yNew);
%         sd=x1'*x1+y1'*y1;
%         
%         if ~isempty(idxtumb)
%                 dtumble=dtumbletemp(idxtumb);
%                 runsTot=[runsTot ; dtumble];
%                 v_tot=[v_tot; v];
%         end
%         
%         if (length(xNew)>1/dt)&&(Nr/length(xNew)>.8) % filtering trajectories
%             flag=flag+1;
%         RunTimeTot=[RunTimeTot RunTime];   
%         sdTot=[sdTot sd];
%         TrajLength=[TrajLength length(xNew)];
%             
%             if rand<.1 % This plots a fraction of the trajectories, all in the same graph
%                 
%                 figure(2)
%                 hold on
%                 plot(xNew,yNew,'-b');
%                 plot(xNew(tumbleSpeed+b/2),yNew(tumbleSpeed+b/2),'or');
%                 plot(xNew(tumbleAngle+b/2),yNew(tumbleAngle+b/2),'og');
% %                 plot(xNew(tumble+b/2),yNew(tumble+b/2),'ok');
%                 hold off
%                 figurename2=[filenameinNew,'Traj.jpeg'];
%                 print(gcf, '-djpeg', '-r800',figurename2)
%                 
%             end
%             
%         
%         
%         end
%         
%     end
%     
%     x_bins=[0:.1:20];
%     n=hist(TrajLength*dt,x_bins);
%     n2=hist(RunTimeTot,x_bins);
%     
%     figure(1)
%     subplot(2,2,1)
%     plot(x_bins,n)
%     subplot(2,2,2)
%     plot(x_bins,n2)
%     subplot(2,2,3)
%     plot(RunTimeTot,TrajLength*dt,'.')
%     subplot(2,2,4)
%     plot(RunTimeTot,sdTot,'.')
%     figurename2=[filenameinNew,'Subplot.jpeg'];
%     print(gcf, '-djpeg', '-r800',figurename2)
% %       pause
%     
%     x_bins=[0:.1:5];
%     n=hist(runsTot*dt,x_bins);
%     s=n;
%     for kk=1:length(s)
%         s(kk)=sum(n(1:kk));
%     end
%     x_tauToFit=x_bins;
%     tauToFit=1-s/sum(n);
%     options1 = optimset('MaxFunEvals',2000, 'Display','final','TolX',1e-6, 'TolFun',1e-6);
%     x_init=[1; 1];
%     x_lb=[.1,0.01];
%     x_ub=[10, 10];
%     paraminit=[1,.5];    
%     FunctionToFit= @(param,x_data) param(1)*exp(-x_data/param(2));
%     paramfound=lsqcurvefit(FunctionToFit,paraminit,x_tauToFit(2:end),tauToFit(2:end),x_lb,x_ub,options1);
%      
%     FunctionToMinimize= @(T) mean((tauToFit(3:end-3)-T(1)*exp(-x_tauToFit(3:end-3)/T(2))).^2);
%     [x0,fval,exitFlag,output]=fminsearch(FunctionToMinimize,x_init,options1);
% 
%     
%     figure(4)
%     hold on
%     semilogy(x_tauToFit,tauToFit,'k')
%     plot(x_tauToFit,x0(1)*exp(-x_tauToFit/x0(2)))
%     plot(x_tauToFit,paramfound(1)*exp(-x_tauToFit/paramfound(2)),'g')
%     hold off
%     MeanSpeed=mean(v_tot);
%     
%      set(gca,'YScale','log');
%      axis([0 5 0.001 2])
%     legend(['MedianRunTime = ' num2str(median(RunTimeTot)) ' on ' num2str(flag) ' trajectories','\newlineMeanSpead= ',num2str(MeanSpeed)],['Fit Mean = ' num2str(paramfound(2))],['Pooled Mean = ' num2str(dt*sum(Nr_tot)/sum(Nrt_tot))],'Location','Best')
% %     disp([filename ' - MedianRunTime = ' num2str(median(RunTimeTot)) ' on ' num2str(flag) ' trajectories'])
% %     disp([filename ' - Fit Mean = ' num2str(paramfound(2))])
% %     disp([filename ' - Pooled Mean = ' num2str(dt*sum(Nr_tot)/sum(Nrt_tot))])
%     figurename2=[filenameinNew,'Runtime.jpeg'];
%     print(gcf, '-djpeg', '-r800',figurename2)
% end
close all;
end