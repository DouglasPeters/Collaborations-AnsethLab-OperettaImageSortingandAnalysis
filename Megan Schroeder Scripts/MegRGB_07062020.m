clear all; clc;
tic
%cd(userpath);

inDir = 'G:\For Doug analysis_20201203\Humans\Male\rep3\*.jpg';
inFolder = 'G:\For Doug analysis_20201203\Humans\Male\rep3\';

BotCrop = 16; %Number of pixels to remove from the bottom of the image to avoid the weird camera error.
WriteImages = 1;
MinSize = 5; %Minimum number of pixels for VK spot to be detected.
SplitSize = 20000; %Area threshold where big and small are split.
EccentricityMAX = 0.99; %How non-circular can VK spots be? This helps filter out out-of-focus cell edges.
StdDevThresh = 3;


varNames = {'Filename','NumObj','AvgArea','AreaMinUsed'};
varTypes = {'string','double','double','double'};

srcFiles = dir(inDir);

Result_Table = table('Size',[length(srcFiles) length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
Area_Distro = cell(length(srcFiles),1);
for g = 1:length(srcFiles)
    disp('Calculating mean RGB levels for data set.')
    filename = strcat(inFolder,srcFiles(g).name);
    I = imread(filename);
    [ResY ResX ResZ] = size(I,[1 2 3]);
    I_Red = I(1:end-BotCrop,:,1); I_Red_inv = imcomplement(I_Red);
    I_Green = I(1:end-BotCrop,:,2); I_Green_inv = imcomplement(I_Green);
    I_Blue = I(1:end-BotCrop,:,3); I_Blue_inv = imcomplement(I_Blue);
    I_AVG(g,1) = mean(I_Red_inv(:));
    I_AVG(g,2) = mean(I_Green_inv(:));
    I_AVG(g,3) = mean(I_Blue_inv(:));
end
    
I_THRESH(1,1:3) = [mean(I_AVG(:,1)) mean(I_AVG(:,2)) mean(I_AVG(:,3))];
I_2xSTD(1,1:3) = [std(I_AVG(:,1))*StdDevThresh std(I_AVG(:,2))*StdDevThresh std(I_AVG(:,3))*StdDevThresh];
I_THRESH2(1,1:3) = [I_THRESH(1,1)+I_2xSTD(1,1) I_THRESH(1,2)+I_2xSTD(1,2) I_THRESH(1,3)+I_2xSTD(1,3)];


for f = 1:length(srcFiles)
    progress = (f/length(srcFiles)*100);
    progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,length(srcFiles),progress,'%');
    disp(progress2)
    filename = strcat(inFolder,srcFiles(f).name);
    
    cd(inFolder);
    
    I = imread(filename);
    
    [ResY ResX ResZ] = size(I,[1 2 3]);
    I_Red = I(1:end-BotCrop,:,1); I_Red_inv = imcomplement(I_Red);
    I_Green = I(1:end-BotCrop,:,2); I_Green_inv = imcomplement(I_Green);
    I_Blue = I(1:end-BotCrop,:,3); I_Blue_inv = imcomplement(I_Blue);
        
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
    
    %%% Visualization of segmentation %%%
if WriteImages == 1
    for slash = 1:length(strfind(filename,'\'))
        if slash ==1, filename2 = extractAfter(filename,'\');
        else filename2 = extractAfter(filename2,'\');
        end
    end
    filename2 = filename2(1:end-4);
    filename3 = append(filename2,'.pdf');
    Folder_Outputs = append(inFolder,'VKFociSegmentation\');
    
    if f==1,mkdir VKFociSegmentation; 
        Folder_Outputs = append(inFolder,'VKFociSegmentation\');
        cd(Folder_Outputs);
    else cd(Folder_Outputs);
    end
    
    B = bwboundaries(I_BW,8);
    BB = bwboundaries(I_bwfilt2,8);
    imshow(I(1:end-15,:,:));
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
    ax = gca;
    exportgraphics(ax,filename3,'ContentType','image','Resolution','400');
else 
end    
end

Result_Cell = table2cell(Result_Table);
Result_Cell(1:length(srcFiles),5) = Area_Distro;
Total_NumObj = sum(Result_Table.NumObj(:));
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
        VK_AREAS(1:length(Result_Cell{r,5}),1) = Result_Cell{r,5}(:);
    else VK_AREAS(end+1:end+length(Result_Cell{r,5}),1) = Result_Cell{r,5}(:);
    end
end

cd(Folder_Outputs);
save('MatLabWorkspace.mat');
