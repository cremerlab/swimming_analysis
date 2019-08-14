%% AVI to particles

alpha=.00001; % quantile fraction to normalize intensity
PSF=fspecial('gaussian',9,.3); % ps for gaussian blur
tau=100;
    tau_inv=1/tau;
% Find all the folders to analyze

folder1='/Volumes/DIVE_FAT/20170307_Tomo/2017-03-07/';
cd(folder1);
listdir=dir;
mask = ones(1,length(listdir));
for i=1:length(listdir)
    if listdir(i).isdir||isempty(findstr('OD',listdir(i).name))||isempty(findstr('.avi',listdir(i).name))
        mask(i)=0;
    end
end
listdir2=listdir(find(mask));
clear('listdir');

% c=cell(3000,3)
%
for j=1length(listdir2)
    

%     disp(j)
    text=listdir2(j).name;
    text2=text(1:length(text)-4);

%     tic
    folder=[folder1 listdir2(j).name '/'];
    folderRes=[folder1 text2 '-Res/'];
    mkdir(folderRes);
    cd(folder1);
    mov=VideoReader([folder1 listdir2(j).name]);
    i=1;
    Imean=double(read(mov,i));  
    initialmeanlength=20
    for i=2:20
    Imean=Imean+double(read(mov,i));
    end
    Imean=Imean/initialmeanlength;
    
    Nimages=uint16(mov.Duration*mov.FrameRate);
    movieInfoCell=cell(Nimages,3);
    
    for i=1:Nimages
%     for i=max([initial_index(j) 2]):uint16(mov.Duration*mov.FrameRate)

%     for i=2:length(list2)
        tic
%         disp(i)
        I=double(read(mov,i));
        
        
        I1=imabsdiff(I,Imean);
%         imshoçw(I1)
%         break
        I2=(imfilter((I1),PSF,'conv'));
        minI=quantile(double(I2(:)),alpha);
        maxI=quantile(double(I2(:)),1-alpha);
        I3=((double(I2)-minI)/(maxI-minI));
        
        thresh=mean(I3(:))+3*std(I3(:));
        BW=im2bw(I3,thresh);
        stats=regionprops(BW,'centroid','area');
        p_centroid=vec2mat([stats.Centroid],2);
        p_area=vec2mat([stats.Area],1);
        big_index=find(p_area>10);
        p_export=p_centroid(big_index,:);
%         length(p(big_index))
        Imean=I*tau_inv+Imean*(1-tau_inv);
%         imshow((BW))
%         hold on
%         for k=1:length(length(p_export))
%             plot(p_export(:,1),p_export(:,2),'+');
%         end
%         hold off
%         pause(.001)        
        cd(folderRes)
        filenameexport='Frame_';
        fid=fopen([filenameexport num2str(i-1)],'w+');
        fprintf(fid, ['frame ' num2str(i-1)]);
        fclose(fid);
        if ~isempty(p_export)
            p2=zeros(length(p_export(:,1)),6);
            p2(1:length(p2),2:3)=p_export;
            dlmwrite([filenameexport num2str(i-1)],p2,'-append','delimiter',' ','roffset',1);
        end
        zerovector=zeros(length(p2(:,2)),1);
        movieInfoCell{i,1}=[p2(:,2) zerovector];
        movieInfoCell{i,2}=[p2(:,3) zerovector];
        movieInfoCell{i,3}=[p2(:,3) zerovector];
        
        disp(['Film : ' text2 '  image : ' num2str(i) '  time elapsed : ' num2str(toc)])
    end
    save(['Cluster_' text2 ],'movieInfoCell');
    
end

filenamein=[folderRes,'Cluster_' text2 ];
dirnameout=filenamein;filenameout=filenamein;
maxSearchRadiusin=30;maxSearchRadiusingap=30;
track_function_liquid(filenamein,dirnameout,filenameout,maxSearchRadiusin,maxSearchRadiusingap);


