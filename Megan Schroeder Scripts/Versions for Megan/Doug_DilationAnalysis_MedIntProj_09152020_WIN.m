clear all; clc;
tic
cd(userpath);

%% USER INPUTS %%

Folder = 'C:\Users\dplad\Desktop\megdata\Ch1 nuclear to _cyto_ ratio (ch1 red RUNX2, ch2 Vimentin, ch3 DAPI)\'; %This tells the script where to find your images (must keep the trailing slash).

Channels = 3; %How many channels are in the image? (Inputs can be 2, 3, or 4). This script assumes Channel 1 is DAPI and will use it to ID nuclei.

Ch_DAPI = 3; %Which channel is being used for nuclear segmentation (i.e., DAPI)? (1-4)
Ch1_Analysis = 1; %Do you want to analyze fluorescence in Channel 1? (Yes=1, No=0)
Ch2_Analysis = 0; %Do you want to analyze fluorescence in Channel 2? (Yes=1, No=0)
Ch3_Analysis = 0; %Do you want to analyze fluorescence in Channel 3? (Yes=1, No=0)
Ch4_Analysis = 0; %Do you want to analyze fluorescence in Channel 4? (Yes=1, No=0)

Radius_Erode = 1; %Adjust the erosion radius to trim more signal from edges of nuclei (prior to dilation). This can help clean up small unwanted pixels.
Radius_Dilate = 5; %Adjust the dilation radius that is used to identify "NEAR" (i.e., cellular) pixels around nuclei. This can make your NEAR:FAR ratio more or less stringent.

ShowImages = 0; %Do you want to render images of the dilation analysis? (Yes=1, No=0)
Ch_Show = 3; %Which channel do you want to render? (Inputs can be 2, 3, or 4).
WriteImages = 0; %Do you want to save images of the dilation analysis as PDFs? (Yes=1, No=0)


%% Everything else (don't touch) %%
cd(Folder);
srcFiles = dir(Folder);

if Ch1_Analysis == 1, RESULTS_Ch1 = table(); else end
if Ch2_Analysis == 1, RESULTS_Ch2 = table(); else end
if Ch3_Analysis == 1, RESULTS_Ch3 = table(); else end
if Ch4_Analysis == 1, RESULTS_Ch4 = table(); else end

for dots = 1:length(srcFiles)
    if srcFiles(dots).name(1) == '.' || srcFiles(dots).isdir == 1
        isimage(dots,1) = 1;
    else isimage(dots,1) = 0; 
    end
end
ff = sum(isimage);

for f = 1+ff:length(srcFiles)
cd(Folder);    
    
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

    Ch1_IntDen(f-ff,1) = sum(Ch1_lin(:)); Ch1_MeanInt(f-ff,1) = mean(Ch1_lin(:));
    Ch2_IntDen(f-ff,1) = sum(Ch2_lin(:)); Ch2_MeanInt(f-ff,1) = mean(Ch2_lin(:));
    if Channels>2, Ch3_IntDen(f-ff,1) = sum(Ch3_lin(:)); Ch3_MeanInt(f-ff,1) = mean(Ch3_lin(:)); else end
    if Channels>3, Ch4_IntDen(f-ff,1) = sum(Ch4_lin(:)); Ch4_MeanInt(f-ff,1) = mean(Ch4_lin(:)); else end

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
    
    if Ch_DAPI == 1, DAPI_medd = Ch1_medd;
    elseif Ch_DAPI == 2, DAPI_medd = Ch2_medd;
    elseif Ch_DAPI == 3, DAPI_medd = Ch3_medd;
    elseif Ch_DAPI == 4, DAPI_medd = Ch4_medd;
    else
    end
    DAPI_BW1 = imbinarize(DAPI_medd);
    DAPI_BW2 = imerode(DAPI_BW1,erodestrel);
    DAPI_BW3 = imdilate(DAPI_BW2,dilatestrel);
    DAPI_Donut = logical(DAPI_BW3-DAPI_BW2);

if Ch1_Analysis == 1,
    Ch1_CELLULAR = Ch1_medd(DAPI_BW3); Ch1_CELLULAR_MEAN(f-ff,1) = mean(Ch1_CELLULAR);
    Ch1_ECM = Ch1_medd(~DAPI_BW3); Ch1_ECM_MEAN(f-ff,1) = mean(Ch1_ECM);
    Ch1_NUC = Ch1_medd(DAPI_BW2); Ch1_NUC_MEAN(f-ff,1) = mean(Ch1_NUC);
    Ch1_PERINUC = Ch1_medd(DAPI_Donut); Ch1_PERINUC_MEAN(f-ff,1) = mean(Ch1_PERINUC);
    RESULTS_Ch1.CELLULAR(f-ff) = Ch1_CELLULAR_MEAN(f-ff,1);
    RESULTS_Ch1.ECM(f-ff) = Ch1_ECM_MEAN(f-ff,1);
    RESULTS_Ch1.CERatio(f-ff) = RESULTS_Ch1.CELLULAR(f-ff)/RESULTS_Ch1.ECM(f-ff);
    RESULTS_Ch1.IntDen(f-ff) = sum(Ch1_medd,'all');
    RESULTS_Ch1.NUC(f-ff) = Ch1_NUC_MEAN(f-ff,1);
    RESULTS_Ch1.PERINUC(f-ff) = Ch1_PERINUC_MEAN(f-ff,1);
    RESULTS_Ch1.NPRatio(f-ff) = RESULTS_Ch1.NUC(f-ff)/RESULTS_Ch1.PERINUC(f-ff);
else end    
    
if Ch2_Analysis == 1,
    Ch2_CELLULAR = Ch2_medd(DAPI_BW3); Ch2_CELLULAR_MEAN(f-ff,1) = mean(Ch2_CELLULAR);
    Ch2_ECM = Ch2_medd(~DAPI_BW3); Ch2_ECM_MEAN(f-ff,1) = mean(Ch2_ECM);
    Ch2_NUC = Ch2_medd(DAPI_BW2); Ch2_NUC_MEAN(f-ff,1) = mean(Ch2_NUC);
    Ch2_PERINUC = Ch2_medd(DAPI_Donut); Ch2_PERINUC_MEAN(f-ff,1) = mean(Ch2_PERINUC);
    RESULTS_Ch2.CELLULAR(f-ff) = Ch2_CELLULAR_MEAN(f-ff,1);
    RESULTS_Ch2.ECM(f-ff) = Ch2_ECM_MEAN(f-ff,1);
    RESULTS_Ch2.CERatio(f-ff) = RESULTS_Ch2.CELLULAR(f-ff)/RESULTS_Ch2.ECM(f-ff);
    RESULTS_Ch2.IntDen(f-ff) = sum(Ch2_medd,'all');
    RESULTS_Ch2.NUC(f-ff) = Ch2_NUC_MEAN(f-ff,1);
    RESULTS_Ch2.PERINUC(f-ff) = Ch2_PERINUC_MEAN(f-ff,1);
    RESULTS_Ch2.NPRatio(f-ff) = RESULTS_Ch2.NUC(f-ff)/RESULTS_Ch2.PERINUC(f-ff);
else end

if Ch3_Analysis == 1,
    Ch3_CELLULAR = Ch3_medd(DAPI_BW3); Ch3_CELLULAR_MEAN(f-ff,1) = mean(Ch3_CELLULAR);
    Ch3_ECM = Ch3_medd(~DAPI_BW3); Ch3_ECM_MEAN(f-ff,1) = mean(Ch3_ECM);
    Ch3_NUC = Ch3_medd(DAPI_BW2); Ch3_NUC_MEAN(f-ff,1) = mean(Ch3_NUC);
    Ch3_PERINUC = Ch3_medd(DAPI_Donut); Ch3_PERINUC_MEAN(f-ff,1) = mean(Ch3_PERINUC);
    RESULTS_Ch3.CELLULAR(f-ff) = Ch3_CELLULAR_MEAN(f-ff,1);
    RESULTS_Ch3.ECM(f-ff) = Ch3_ECM_MEAN(f-ff,1);
    RESULTS_Ch3.CERatio(f-ff) = RESULTS_Ch3.CELLULAR(f-ff)/RESULTS_Ch3.ECM(f-ff);
    RESULTS_Ch3.IntDen(f-ff) = sum(Ch3_medd,'all');
    RESULTS_Ch3.NUC(f-ff) = Ch3_NUC_MEAN(f-ff,1);
    RESULTS_Ch3.PERINUC(f-ff) = Ch3_PERINUC_MEAN(f-ff,1);
    RESULTS_Ch3.NPRatio(f-ff) = RESULTS_Ch3.NUC(f-ff)/RESULTS_Ch3.PERINUC(f-ff);
else end

if Ch4_Analysis == 1,
    Ch4_CELLULAR = Ch4_medd(DAPI_BW3); Ch4_CELLULAR_MEAN(f-ff,1) = mean(Ch4_CELLULAR);
    Ch4_ECM = Ch4_medd(~DAPI_BW3); Ch4_ECM_MEAN(f-ff,1) = mean(Ch4_ECM);
    Ch4_NUC = Ch4_medd(DAPI_BW2); Ch4_NUC_MEAN(f-ff,1) = mean(Ch4_NUC);
    Ch4_PERINUC = Ch4_medd(DAPI_Donut); Ch4_PERINUC_MEAN(f-ff,1) = mean(Ch4_PERINUC);
    RESULTS_Ch4.CELLULAR(f-ff) = Ch4_CELLULAR_MEAN(f-ff,1);
    RESULTS_Ch4.ECM(f-ff) = Ch4_ECM_MEAN(f-ff,1);
    RESULTS_Ch4.CERatio(f-ff) = RESULTS_Ch4.CELLULAR(f-ff)/RESULTS_Ch4.ECM(f-ff);
    RESULTS_Ch4.IntDen(f-ff) = sum(Ch4_medd,'all');
    RESULTS_Ch4.NUC(f-ff) = Ch4_NUC_MEAN(f-ff,1);
    RESULTS_Ch4.PERINUC(f-ff) = Ch4_PERINUC_MEAN(f-ff,1);
    RESULTS_Ch4.NPRatio(f-ff) = RESULTS_Ch4.NUC(f-ff)/RESULTS_Ch4.PERINUC(f-ff);
else end

%%% Visualization of segmentation %%%
if ShowImages == 1 && WriteImages == 0
    B = bwboundaries(DAPI_BW3,8);
    if Ch_Show ==1, imshow(Ch1_medd.*5);
    elseif Ch_Show ==2, imshow(Ch2_medd.*5);
    elseif Ch_Show ==3, imshow(Ch3_medd.*5);
    elseif Ch_Show ==4, imshow(Ch4_medd.*5);
    else
        imshow(Ch1_medd.*5); 
    end
    hold on
    for k = 1:size(B,1)
        b = B{k};
        plot(b(:,2),b(:,1),'r','LineWidth',1);
    end
    hold off
    drawnow;
elseif ShowImages == 1 && WriteImages == 1    
    B = bwboundaries(DAPI_BW3,8);
    if Ch_Show ==1, imshow(Ch1_medd.*5);
    elseif Ch_Show ==2, imshow(Ch2_medd.*5);
    elseif Ch_Show ==3, imshow(Ch3_medd.*5);
    elseif Ch_Show ==4, imshow(Ch4_medd.*5);
    else
        imshow(Ch1_medd.*5); 
    end
    hold on
    for k = 1:size(B,1)
        b = B{k};
        plot(b(:,2),b(:,1),'r','LineWidth',1);
    end
    hold off
    drawnow;
    figfilename = strcat(srcFiles(f).name(1:end-4),' DilationAnalysis.pdf');
    Folder_Outputs = append(Folder,'DilationAnalysis\');
    if f-ff==1,mkdir DilationAnalysis; 
        Folder_Outputs = append(Folder,'DilationAnalysis\');
        cd(Folder_Outputs);
    else
        cd(Folder_Outputs);
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

if f == length(srcFiles)
    Workspace_Image = strcat(srcFiles(f).name(1:end-4),' Analysis Workspace.mat');
    save(Workspace_Image);
else
end

end

close all;
cd(Folder_Outputs);
if Ch1_Analysis == 1, save('Ch1 Results.mat','RESULTS_Ch1'); else end
if Ch2_Analysis == 1, save('Ch2 Results.mat','RESULTS_Ch2'); else end
if Ch3_Analysis == 1, save('Ch3 Results.mat','RESULTS_Ch3'); else end
if Ch4_Analysis == 1, save('Ch4 Results.mat','RESULTS_Ch4'); else end

toc