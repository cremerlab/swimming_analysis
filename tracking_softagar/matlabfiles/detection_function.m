function detection_function(filenamein,dirnameout,rawdatafolder)
%% Confocal images
%folder='/Volumes/DIVE_FAT/2016-12-02/Experiment20161201_compare1_Agar_1/';
cd(rawdatafolder);
filenameRadical=dirnameout;


% Make a list of the images to open
list=dir;
mask = ones(1,length(list));
for i=1:length(list)
    if isempty(findstr(filenameRadical,list(i).name))
        mask(i)=0;
    end
end
list2=list(find(mask));

PSF=fspecial('gaussian',9,.3); %PSF used to smooth the images
ptot=[];
l=[];

%Open the images and extract the position of the bacteria
movieInfoCell=cell(length(list2),3);
display(length(list2));
display(list2);
for i=1:length(list2)
    I=imread(list2(i).name); %Open image i
    I1=imfilter(double(I),PSF,'conv'); % smooth image using psf
    
    I2=imadjust(I1);
    I2=imclearborder(I2);
    
    BW =imregionalmax(I2,26); % Transform into a binary using regional max
    BW2=bwmorph(BW,'spur'); % Next lines are to clean the binary from hot pixels
    BW2=bwmorph(BW2,'hbreak');
    BW2=bwmorph(BW2,'spur');
    BW2=bwmorph(BW2,'clean');
    BW2=bwmorph(BW2,'open',1);
    BW2=bwpropfilt(BW2,'area',[10 1000]); % remove small objects, adjust if size of objects changes a lot
    stats=regionprops(BW2,'centroid'); %get the centroid of the objects
    
    
     p=vec2mat([stats.Centroid],2); % p is the position matrix for the current image
    
        movieInfoCell{i,1}(:,1)=p(:,2);%x,y inter-changes
        movieInfoCell{i,1}(:,2)=0;
        movieInfoCell{i,2}(:,1)=p(:,1);
        movieInfoCell{i,2}(:,2)=0;
        movieInfoCell{i,3}(:,1)=255*ones(size(p,1),1);
        movieInfoCell{i,3}(:,2)=0;

end



save(filenamein,'movieInfoCell');


end