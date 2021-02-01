clear all; clc;
tic
cd(userpath);

%% USER INPUTS %%

Folder = 'C:\Users\dplad\Desktop\test\'; %This tells the script where to find your images (must keep the trailing slash).

Channels = 2; %How many channels are in the image? (Inputs can be 2, 3, or 4). This script assumes Channel 1 is DAPI and will use it to ID nuclei.
Ch2_Analysis = 1; %Do you want to analyze fluorescence in Channel 2? (Yes=1, No=0)
Ch3_Analysis = 0; %Do you want to analyze fluorescence in Channel 3? (Yes=1, No=0)
Ch4_Analysis = 0; %Do you want to analyze fluorescence in Channel 4? (Yes=1, No=0)

Radius_Erode = 1; %Adjust the erosion radius to trim more signal from edges of nuclei (prior to dilation). This can help clean up small unwanted pixels.
Radius_Dilate = 5; %Adjust the dilation radius that is used to identify "NEAR" (i.e., cellular) pixels around nuclei. This can make your NEAR:FAR ratio more or less stringent.

ShowImages = 0; %Do you want to render images of the dilation analysis? (Yes=1, No=0)
Ch_Show = 2; %Which channel do you want to render? (Inputs can be 2, 3, or 4).

WriteImages = 0; %Do you want to save images of the dilation analysis as PDFs? (Yes=1, No=0)






%% Everything else (don't touch) %%
cd(Folder);
srcFiles = dir(Folder);

if Ch2_Analysis == 1, Ch2_NEARFAR_MEANS = table(); else end
if Ch3_Analysis == 1, Ch3_NEARFAR_MEANS = table(); else end
if Ch4_Analysis == 1, Ch4_NEARFAR_MEANS = table(); else end

for f = 1:length(srcFiles)
cd(Folder);    
    for dots = 1:length(srcFiles)
        isimage(dots,1) = double(srcFiles(dots).isdir);
    end
    ff = sum(isimage);
if srcFiles(f).isdir == 0
    progress = ((f-ff)/(length(srcFiles)-ff)*100);
    progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f-ff,(length(srcFiles)-ff),progress,'%');
    clc; disp(progress2);

    filename = srcFiles(f).name;
    cd(Folder);
    I = bfopen(filename);
    Res = length(I{1,1}{1,1});
    Reshape_lin = Res*Res;
    ReshapeV = [Reshape_lin 1];
    Slices = (length(I{1,1})/Channels);
    Blank = zeros(Res,Res,Slices);
    Planes = Channels*Slices;

    Ch1 = uint16(Blank);
    Ch2 = uint16(Blank);
    if Channels>2, Ch3 = uint16(Blank); else end
    if Channels>3, Ch4 = uint16(Blank); else end

    for i = 1:Slices
        Ch1_planes(i,1) = 1+(Channels*i-Channels);
        Ch2_planes(i,1) = 2+(Channels*i-Channels);
        if Channels>2, Ch3_planes(i,1) = 3+(Channels*i-Channels); else end
        if Channels>3, Ch4_planes(i,1) = 4+(Channels*i-Channels); else end
    end

    for m = 1:Slices
        Ch1(:,:,m) = I{1,1}{Ch1_planes(m,1),1};
            Ch1_lin(1:Reshape_lin,m) = reshape(Ch1(:,:,m),ReshapeV);
        Ch2(:,:,m) = I{1,1}{Ch2_planes(m,1),1};
            Ch2_lin(1:Reshape_lin,m) = reshape(Ch2(:,:,m),ReshapeV);
        if Channels>2, Ch3(:,:,m) = I{1,1}{Ch3_planes(m,1),1};
            Ch3_lin(1:Reshape_lin,m) = reshape(Ch3(:,:,m),ReshapeV); else end
        if Channels>3, Ch4(:,:,m) = I{1,1}{Ch4_planes(m,1),1};
            Ch4_lin(1:Reshape_lin,m) = reshape(Ch4(:,:,m),ReshapeV); else end
    end

    Ch1_IntDen(f,1) = sum(Ch1_lin(:)); Ch1_MeanInt(f,1) = mean(Ch1_lin(:));
    Ch2_IntDen(f,1) = sum(Ch2_lin(:)); Ch2_MeanInt(f,1) = mean(Ch2_lin(:));
    if Channels>2, Ch3_IntDen(f,1) = sum(Ch3_lin(:)); Ch3_MeanInt(f,1) = mean(Ch3_lin(:)); else end
    if Channels>3, Ch4_IntDen(f,1) = sum(Ch4_lin(:)); Ch4_MeanInt(f,1) = mean(Ch4_lin(:)); else end

    for sumx = 1:Res
        for sumy = 1:Res
        Ch1_med(sumx,sumy) = max(Ch1(sumx,sumy,:));
        Ch2_med(sumx,sumy) = median(Ch2(sumx,sumy,:));
        if Channels>2, Ch3_med(sumx,sumy) = median(Ch3(sumx,sumy,:)); else end
        if Channels>3, Ch4_med(sumx,sumy) = median(Ch4(sumx,sumy,:)); else end
        end
    end

    Ch1_medd = uint16(Ch1_med);
    Ch2_medd = uint16(Ch2_med);
    if Channels>2, Ch3_medd = uint16(Ch3_med); else end
    if Channels>3, Ch4_medd = uint16(Ch3_med); else end

    erodestrel = strel('disk',Radius_Erode);
    dilatestrel = strel('disk',Radius_Dilate);
    Ch1_BW1 = imbinarize(Ch1_medd);
    Ch1_BW2 = imerode(Ch1_BW1,erodestrel);
    Ch1_BW3 = imdilate(Ch1_BW2,dilatestrel);
   
if Ch2_Analysis == 1,
    Ch2_NEAR = Ch2_medd(Ch1_BW3); Ch2_NEAR_MEAN(f-ff,1) = mean(Ch2_NEAR);
    Ch2_FAR = Ch2_medd(~Ch1_BW3); Ch2_FAR_MEAN(f-ff,1) = mean(Ch2_FAR);
    Ch2_NEARFAR_MEANS.NEAR(f-ff) = Ch2_NEAR_MEAN(f-ff,1);
    Ch2_NEARFAR_MEANS.FAR(f-ff) = Ch2_FAR_MEAN(f-ff,1);
    Ch2_NEARFAR_MEANS.NFRatio(f-ff) = Ch2_NEARFAR_MEANS.NEAR(f-ff)/Ch2_NEARFAR_MEANS.FAR(f-ff);
else end

if Ch3_Analysis == 1,
    Ch3_NEAR = Ch3_medd(Ch1_BW3); Ch3_NEAR_MEAN(f-ff,1) = mean(Ch3_NEAR);
    Ch3_FAR = Ch3_medd(~Ch1_BW3); Ch3_FAR_MEAN(f-ff,1) = mean(Ch3_FAR);
    Ch3_NEARFAR_MEANS.NEAR(f-ff) = Ch3_NEAR_MEAN(f-ff,1);
    Ch3_NEARFAR_MEANS.FAR(f-ff) = Ch3_FAR_MEAN(f-ff,1);
    Ch3_NEARFAR_MEANS.NFRatio(f-ff) = Ch3_NEARFAR_MEANS.NEAR(f-ff)/Ch3_NEARFAR_MEANS.FAR(f-ff);
else end

if Ch4_Analysis == 1,
    Ch4_NEAR = Ch4_medd(Ch1_BW3); Ch4_NEAR_MEAN(f-ff,1) = mean(Ch4_NEAR);
    Ch4_FAR = Ch4_medd(~Ch1_BW3); Ch4_FAR_MEAN(f-ff,1) = mean(Ch4_FAR);
    Ch4_NEARFAR_MEANS.NEAR(f-ff) = Ch4_NEAR_MEAN(f-ff,1);
    Ch4_NEARFAR_MEANS.FAR(f-ff) = Ch4_FAR_MEAN(f-ff,1);
    Ch4_NEARFAR_MEANS.NFRatio(f-ff) = Ch4_NEARFAR_MEANS.NEAR(f-ff)/Ch4_NEARFAR_MEANS.FAR(f-ff);
else end

%%% Visualization of segmentation %%%
if ShowImages == 1 && WriteImages == 0
    B = bwboundaries(Ch1_BW3,8);
    if Ch_Show ==2, imshow(Ch2_medd.*5);
    elseif Ch_Show ==3, imshow(Ch3_medd.*5);
    elseif Ch_Show ==4, imshow(Ch4_medd.*5);
    else imshow(Ch1_medd.*5); end
    hold on
    for k = 1:size(B,1)
        b = B{k};
        plot(b(:,2),b(:,1),'r','LineWidth',1);
    end
    hold off
elseif ShowImages == 1 && WriteImages == 1    
    B = bwboundaries(Ch1_BW3,8);
    if Ch_Show ==2, imshow(Ch2_medd.*5);
    elseif Ch_Show ==3, imshow(Ch3_medd.*5);
    elseif Ch_Show ==4, imshow(Ch4_medd.*5);
    else imshow(Ch1_medd.*5); end
    hold on
    for k = 1:size(B,1)
        b = B{k};
        plot(b(:,2),b(:,1),'r','LineWidth',1);
    end
    hold off
    
    figfilename = strcat(srcFiles(f).name(1:end-4),' DilationAnalysis.pdf');
    Folder_Outputs = append(Folder,'DilationAnalysis\');
    if f-ff==1,mkdir DilationAnalysis; 
        Folder_Outputs = append(Folder,'DilationAnalysis\');
        cd(Folder_Outputs);
    else cd(Folder_Outputs);
    end
    ax = gca;
    exportgraphics(ax,figfilename,'ContentType','image','Resolution','400');
elseif ShowImages == 0 && WriteImages == 1
    figfilename = strcat(srcFiles(f).name(1:end-4),' DilationAnalysis.pdf');
    Folder_Outputs = append(Folder,'DilationAnalysis\');
    if f-ff==1,mkdir DilationAnalysis; 
        Folder_Outputs = append(Folder,'DilationAnalysis\');
        cd(Folder_Outputs);
    else cd(Folder_Outputs);
    end
    ax = gca;
    exportgraphics(ax,figfilename,'ContentType','image','Resolution','400');
elseif ShowImages == 0 && WriteImages == 0 && f-ff==1,
    mkdir DilationAnalysis;
    Folder_Outputs = append(Folder,'DilationAnalysis\');
    cd(Folder_Outputs);
else cd(Folder_Outputs);
end 

Workspace_Image = strcat(srcFiles(f).name(1:end-4),' Analysis Workspace.mat');
save(Workspace_Image);

else
end

end

close all;
if Ch2_Analysis == 1, save('Ch2_NearFarRatios.mat','Ch2_NEARFAR_MEANS'); else end
if Ch3_Analysis == 1, save('Ch3_NearFarRatios.mat','Ch3_NEARFAR_MEANS'); else end
if Ch4_Analysis == 1, save('Ch4_NearFarRatios.mat','Ch4_NEARFAR_MEANS'); else end

toc