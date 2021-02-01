clear all; clc;
tic
%% User Inputs %%

Folder = 'C:\Users\dplad\Desktop\megdata\male gels 20x, cacl2\'; %Tells the script where to find your images (must leave the trailing slash).

BotCrop = 16; %Number of pixels to remove from the bottom of the image to avoid the weird camera error.
ShowImages = 1; %Do you want to display images of the analyzed images? (Yes=1, No=0)
WriteImages = 1; %Do you want to save PDFs of the analyzed images? (Yes=1, No=0)
MinSize = 5; %Minimum number of pixels for VK spot to be detected.
SplitSize = 60000; %Area threshold (pixels) where big and small are split.
EccentricityMAX = 0.97; % (0-1) How non-circular can VK spots be? This helps filter out out-of-focus cell edges. The closer to 1, the more linear objects (cell edges) will in included.
StdDevThresh = 3; %This determines the threshold for identifying VK spots. The higher the number, the more stringent the threshold.

%Result_Cell variable contains the following information for each image in the data set: 
%(1) Filename, (2) Number of VK objects,(3) Average area of VK objects, (4) VK area threshold, (5) VK Integrated Density, (6) Areas of VK all VK objects.

%Result_SplitDemographic variable contains: 
%(1) Integrated area of small VK objects (below SplitSize threshold), (2) Integrated area of large VK objects (above SplitSize threshold), and (3) Integrated area of ALL VK objects.

%Result_VKAREAS_ALL variable contains a concatenated list of all VK objects detected in the data set (for histograms, for example).

%% Everything else (don't touch) %%
varNames = {'Filename','NumObj','AvgArea','AreaMinUsed','IntDen'};
varTypes = {'string','double','double','double','double'};
inDir = append(Folder,'*.jpg'); %Don't touch this.
srcFiles = dir(inDir);
Result_Table = table('Size',[length(srcFiles) length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
Area_Distro = cell(length(srcFiles),1);

for g = 1:length(srcFiles)
    rgbprog = (g/length(srcFiles)*100);
    rgbprog2 = sprintf('Calculating mean RGB levels for data set, image %d of %d. %0.2f%c complete.',g,length(srcFiles),rgbprog,'%');
    if g == 1, disp(rgbprog2);
    else, clc; disp(rgbprog2); end
    filename = strcat(Folder,srcFiles(g).name);
    I = imread(filename); I2 = I(1:end-BotCrop,:,:);
    [ResY ResX ResZ] = size(I,[1 2 3]);
    I_Red = I2(:,:,1); I_Red_inv = imcomplement(I_Red);
    I_Green = I2(:,:,2); I_Green_inv = imcomplement(I_Green);
    I_Blue = I2(:,:,3); I_Blue_inv = imcomplement(I_Blue);
    I_AVG(g,1) = mean(I_Red_inv(:));
    I_AVG(g,2) = mean(I_Green_inv(:));
    I_AVG(g,3) = mean(I_Blue_inv(:));
end

I_THRESH(1,1:3) = [mean(I_AVG(:,1)) mean(I_AVG(:,2)) mean(I_AVG(:,3))];
I_overSTD(1,1:3) = [std(I_AVG(:,1))*StdDevThresh std(I_AVG(:,2))*StdDevThresh std(I_AVG(:,3))*StdDevThresh];
I_THRESH2(1,1:3) = [I_THRESH(1,1)+I_overSTD(1,1) I_THRESH(1,2)+I_overSTD(1,2) I_THRESH(1,3)+I_overSTD(1,3)];
if ShowImages ==1, figure(); else end
for f = 1:3 %length(srcFiles)
    clf;
    progress = (f/length(srcFiles)*100);
    progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,length(srcFiles),progress,'%');
    if f==1, disp(progress2);
    else, clc; disp(progress2); end
    filename = strcat(Folder,srcFiles(f).name);
    
    cd(Folder);
    
    I = imread(filename); I2 = I(1:end-BotCrop,:,:);
    
    [ResY ResX ResZ] = size(I,[1 2 3]);
    I_Red = I2(:,:,1); I_Red_inv = imcomplement(I_Red);
    I_Green = I2(:,:,2); I_Green_inv = imcomplement(I_Green);
    I_Blue = I2(:,:,3); I_Blue_inv = imcomplement(I_Blue);
        
    I_BW = zeros([ResY-BotCrop ResX]);
    for y = 1:ResY-BotCrop,
        for x = 1:ResX,
            if I_Red_inv(y,x)>I_THRESH2(1,1) & I_Green_inv(y,x)>I_THRESH2(1,2) & I_Blue_inv(y,x)>I_THRESH2(1,3)
                I_BW(y,x) = 1;
            else I_BW(y,x) = 0;
            end
        end
    end
    
    I_bwfilt = bwareaopen(I_BW,MinSize);
    I_bwfilt2 = bwpropfilt(I_bwfilt,'Eccentricity',[0 EccentricityMAX]);
    I_cc = bwconncomp(I_bwfilt2,4);
    I_regionprops = regionprops(I_cc,'Area','Circularity','Solidity');
    
    
    Areas = [I_regionprops.Area]; Areas = Areas'; RowDist = length(Areas);
    Area_Distro(f,1) = mat2cell(Areas,RowDist);
    
    Result_Table.Filename(f) = filename;
    Result_Table.NumObj(f) = I_cc.NumObjects;
    Result_Table.AvgArea(f) = mean(Areas(:));
    Result_Table.AreaMinUsed(f) = MinSize;
    
    Gray1 = im2uint8(rgb2gray(I2));
    Gray2 = imcomplement(Gray1);
    IntDen1 = regionprops(I_cc,Gray2,'PixelValues');
    
    for id = 1:length(IntDen1)
        IntDen2(id,1) = sum(IntDen1(id,1).PixelValues);
    end
    
    Result_Table.IntDen(f) = sum(IntDen2(:));
    
    %%% Visualization of segmentation %%%
if ShowImages == 1 && WriteImages == 0
    B = bwboundaries(I_BW,8);
    BB = bwboundaries(I_bwfilt2,8);
    imshow(I2);
    hold on
    for k = 1:size(B,1)
        b = B{k};
        plot(b(:,2),b(:,1),'g','LineWidth',1);
    end
    for kk = 1:size(BB,1)
        bb = BB{kk};
        plot(bb(:,2),bb(:,1),'r','LineWidth',1);
    end
    hold off
elseif ShowImages == 1 && WriteImages == 1    
    B = bwboundaries(I_BW,8);
    BB = bwboundaries(I_bwfilt2,8);
    imshow(I2);
    hold on
    for k = 1:size(B,1)
        b = B{k};
        plot(b(:,2),b(:,1),'g','LineWidth',1);
    end
    for kk = 1:size(BB,1)
        bb = BB{kk};
        plot(bb(:,2),bb(:,1),'r','LineWidth',1);
    end
    hold off
    
    figfilename = strcat(srcFiles(f).name(1:end-4),' Analysis.pdf');
    Folder_Outputs = append(Folder,'VKFociSegmentation\');
    if f==1,mkdir VKFociSegmentation; 
        Folder_Outputs = append(Folder,'VKFociSegmentation\');
        cd(Folder_Outputs);
    else cd(Folder_Outputs);
    end
    ax = gca;
    exportgraphics(ax,figfilename,'ContentType','image','Resolution','400');
elseif ShowImages == 0 && WriteImages == 1
    figfilename = strcat(srcFiles(f).name(1:end-4),' Analysis.pdf');
    Folder_Outputs = append(Folder,'VKFociSegmentation\');
    if f==1,mkdir VKFociSegmentation; 
        Folder_Outputs = append(Folder,'VKFociSegmentation\');
        cd(Folder_Outputs);
    else cd(Folder_Outputs);
    end
    ax = gca;
    exportgraphics(ax,figfilename,'ContentType','image','Resolution','400');
else
end

end

Result_Cell = table2cell(Result_Table);
Result_Cell(1:length(srcFiles),6) = Area_Distro;
%Total_NumObj = sum(Result_Table.NumObj(:));
% Area_Distro_Long = zeros([Total_NumObj 1]);

for d = 1:length(srcFiles)
if d == 1,
    Area_Distro_Long(1:Result_Table.NumObj(d),1) = cell2mat(Area_Distro(d,1));
else Area_Distro_Long(end+1:end+Result_Table.NumObj(d),1) = cell2mat(Area_Distro(d,1));
end
end

Area_Distro_BigFilter = Area_Distro_Long>SplitSize;
Area_Distro_SmallFilter = Area_Distro_Long<SplitSize;

Area_Distro_Big = Area_Distro_Long(Area_Distro_BigFilter);   
Area_Distro_Small = Area_Distro_Long(Area_Distro_SmallFilter);
Result_SplitDemographics(1,1) = sum(Area_Distro_Small(:));
Result_SplitDemographics(1,2) = sum(Area_Distro_Big(:));
Result_SplitDemographics(1,3) = sum(Result_SplitDemographics(1,1:2));

close all


for r = 1:length(Result_Cell)
    if r == 1, 
        Result_VKAREAS_ALL(1:length(Result_Cell{r,5}),1) = Result_Cell{r,5}(:);
    else Result_VKAREAS_ALL(end+1:end+length(Result_Cell{r,5}),1) = Result_Cell{r,5}(:);
    end
end

matfilename = strcat(Folder,' Combined Analysis Workspace.mat');
save(matfilename);

toc
