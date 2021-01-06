clear; clc; tic; cd(userpath);
%% Editables %%

Folder = 'F:\MatLab_Exports\2020-11-11_MaleVICs_Soft_Y-27632\2020-11-11_MaleVICs_Soft_Y27632__2020-11-11T08_18_16-Measurement1\Images\ImageStacks\'; 
%Where are the images located? IMPORTANT: Use format 'FILEPATH\'. The apostrophes and ending slash are important.

FigShow = 1; %Do you want to show the Figure with segmentation overlay? (1=yes, 0=no)
FigSave = 0; %Do you want to save the Figure that is generated during analysis? (1=yes, 0=no)

Channels = 3; %How many fluorescent channels are in the image?

CH_DAPI = 1; %Which channel corresponds to DAPI signal?
CH_CellMask = 3; %Which channel corresponds to CellMask (CM) signal?
CM_SizeFilt = 1000;

TF_Analysis = 0; %Will do TF translocalization analysis (1=yes, 0=no)
CH_TF = 3; %Which channel corresponds to transcription factor (TF) signal?

aSMA_Analysis = 1; %Will do gradient-based aSMA analysis (1=yes, 0=no)
CH_aSMA = 2; %Which channel corresponds to aSMA signal?

SplitClusterResults = 1; %Will split results for each well into clusters (>MinClusterSize) and non-clusters (<MinClusterSize). (1=yes, 0=no)
MinClusterSize = 8; %This dictates how "clusters" are identified for result sorting. 
%%The number refers to how many nuclei must be detected for a "group" to be called a "cluster.

%% Analysis Pre-Reqs and Metadata %%

cd(Folder);
srcFiles = dir('*.tiff');

START = 1;
FINISH = length(srcFiles);

for u = 1:length(srcFiles)
    srcFiles_FOV(u,1) = str2num(srcFiles(u).name(8:9));
    srcFiles_WELL(u,1) = string(srcFiles(u).name(1:6));
end
Fields = max(srcFiles_FOV);
Wells = numel(srcFiles)/Fields;

if FigShow == 1, figure,
else end

%% Analysis %%
for f = START:FINISH
    time(f,1).ElapsedSeconds = toc;
    try
clc
filename = strcat(Folder,srcFiles(f).name);
progress = (((FINISH-START+1)-(FINISH-f))/FINISH)*100;
progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,FINISH,progress,'%');
disp(progress2);
if progress < 10,
    disp('Estimated time remaining will display after 10% of images are analyzed...');
else
    time(f).AverageSecondsPerLoop = time(f).ElapsedSeconds/((FINISH-START+1)-(FINISH-f));
    time(f).EstimatedTotalSeconds = time(f).AverageSecondsPerLoop*(FINISH-START+1);
    time(f).EstimatedSecondsLeft = time(f).EstimatedTotalSeconds-time(f).ElapsedSeconds;
    time(f).EstimatedMinutesLeft = time(f).EstimatedSecondsLeft/60;
    time(f).EstimatedMinutesElapsed = time(f).ElapsedSeconds/60;
    estimate = sprintf('Run time: %0.2f minutes.',time(f).EstimatedMinutesElapsed);
    estimate2 = sprintf('Estimated time remaining: %0.2f minutes.',time(f).EstimatedMinutesLeft);
    disp(estimate);
    disp(estimate2);
end
    
Results(f).FileName = srcFiles(f).name;
Results(f).Well = srcFiles_WELL{f};
Results(f).FieldofView = srcFiles_FOV(f);

I = bfopen(filename);

ResY = size(I{1,1}{1,1},1);
ResX = size(I{1,1}{1,1},2);
ZPlanes = (length(I{1,1})/Channels);
Blank3D = zeros(ResY,ResX,ZPlanes);
Slices = Channels*ZPlanes; 

%% Parsing Focal Plane and Creating In-Focus Image %%

disp('Generating in-focus images...');
[Ch1,Ch2,Ch3,Ch4,Focus_Ch1,Focus_Ch2,Focus_Ch3,Focus_Ch4,ZFocus] = InFocusImage(I,ZPlanes,ResY,ResX,Channels,f);

Results(f).ZFocus = ZFocus(f,1);

if Results(f).ZFocus == 1 || Results(f).ZFocus == ZPlanes,
warning('No in-focus plane found. Continuing to next image.');
pause(2);
    if f == START,
        cd(Folder); mkdir('Analysis'); cd Analysis;
    else end;
continue
else
end

%% Cell Mask Segmentation %%

disp('Segmenting Cell Mask signal...');
if CH_CellMask == 1, [CM_IndCells,CMseg_cc,CMseg_props,CMseg_CentroidXY,CM_Watershed_BW2,CM_Watershed_Perim] = CellMaskSegmentation(Focus_Ch1,CM_SizeFilt,ResY,ResX);
elseif CH_CellMask == 2, [CM_IndCells,CMseg_cc,CMseg_props,CMseg_CentroidXY,CM_Watershed_BW2,CM_Watershed_Perim] = CellMaskSegmentation(Focus_Ch2,CM_SizeFilt,ResY,ResX);
elseif CH_CellMask == 3, [CM_IndCells,CMseg_cc,CMseg_props,CMseg_CentroidXY,CM_Watershed_BW2,CM_Watershed_Perim] = CellMaskSegmentation(Focus_Ch3,CM_SizeFilt,ResY,ResX);
elseif CH_CellMask == 4, [CM_IndCells,CMseg_cc,CMseg_props,CMseg_CentroidXY,CM_Watershed_BW2,CM_Watershed_Perim] = CellMaskSegmentation(Focus_Ch4,CM_SizeFilt,ResY,ResX); 
else end

%% Nuclear Segmentation %%

disp('Segmenting DAPI signal...');
if CH_DAPI == 1, [DAPI_75Percentile,DAPI_Watershed_BW2,DAPI_Watershed_Perim,DAPIseg_cc,DAPIseg_props] = NucleiSegmentation(Focus_Ch1,ResY,CM_Watershed_BW2);
elseif CH_DAPI == 2, [DAPI_75Percentile,DAPI_Watershed_BW2,DAPI_Watershed_Perim,DAPIseg_cc,DAPIseg_props] = NucleiSegmentation(Focus_Ch2,ResY,CM_Watershed_BW2);
elseif CH_DAPI == 3, [DAPI_75Percentile,DAPI_Watershed_BW2,DAPI_Watershed_Perim,DAPIseg_cc,DAPIseg_props] = NucleiSegmentation(Focus_Ch3,ResY,CM_Watershed_BW2);
elseif CH_DAPI == 4, [DAPI_75Percentile,DAPI_Watershed_BW2,DAPI_Watershed_Perim,DAPIseg_cc,DAPIseg_props] = NucleiSegmentation(Focus_Ch4,ResY,CM_Watershed_BW2); 
else end

%% Cellular Analysis %%

disp('Quantifying image properties for each cell/nucleus...');
if CH_DAPI == 1, DAPI = Focus_Ch1;
elseif CH_DAPI == 2, DAPI = Focus_Ch2;
elseif CH_DAPI == 3, DAPI = Focus_Ch3;
elseif CH_DAPI == 4, DAPI = Focus_Ch4;
else end

if CH_CellMask == 1, CellMask = Focus_Ch1;
elseif CH_CellMask == 2, CellMask = Focus_Ch2;
elseif CH_CellMask == 3, CellMask = Focus_Ch3;
elseif CH_CellMask == 4, CellMask = Focus_Ch4;
else end

clearvars NearestNuc;
[Results_CellAnalysis,NearestNucDistanceFiltered] = CellularAnalysis(DAPI_75Percentile,DAPI_Watershed_BW2,Channels,CMseg_props,CM_IndCells,DAPI,Focus_Ch2,Focus_Ch3,CellMask);

clearvars NucleiNumber AdjustedNucleiNumber;
for q = 1:size(Results_CellAnalysis,2)
    NucleiNumber(q,1) = Results_CellAnalysis(q).NucleiNumber; 
    AdjustedNucleiNumber(q,1) = Results_CellAnalysis(q).AdjustedNucleiNumber;
end

%% TF Nuclear Translocalization Analysis %%

clearvars NucCytoRatio_sorted
if TF_Analysis == 1,
    disp('Quantifying TF localization in each nucleus...');
    if CH_TF == 1, [TFAnalysis,NucCytoRatio_sorted] = NucTranslocation(Focus_Ch1,Results_CellAnalysis,Blank3D);
    elseif CH_TF == 2, [TFAnalysis,NucCytoRatio_sorted] = NucTranslocation(Focus_Ch2,Results_CellAnalysis,Blank3D);
    elseif CH_TF == 3, [TFAnalysis,NucCytoRatio_sorted] = NucTranslocation(Focus_Ch3,Results_CellAnalysis,Blank3D);
    elseif CH_TF == 4, [TFAnalysis,NucCytoRatio_sorted] = NucTranslocation(Focus_Ch4,Results_CellAnalysis,Blank3D);
    else
    end
else
end

%% aSMA Activation Analysis (work in progress) %%

clearvars Grad_Stats
if aSMA_Analysis == 1,
    disp('Quantifying aSMA gradient intensities for each cell...');
    if CH_aSMA == 1, aSMA = Focus_Ch1;
    elseif CH_aSMA == 2, aSMA = Focus_Ch2;
    elseif CH_aSMA == 3, aSMA = Focus_Ch3;
    elseif CH_aSMA == 4, aSMA = Focus_Ch4;
    else
    end
    
    [Grad_Stats,aSMAGradImageTotal] = aSMAActivation(aSMA,Results_CellAnalysis);
else
end


%% Figure %%

if FigShow > 0,
    disp('Generating Figure...');
totalfilteredcellperims = false(ResY,ResX);
for v = 1:size(Results_CellAnalysis,2)
    if v == 1, totalfilteredcellperims = Results_CellAnalysis(v).LogPerim;
    else totalfilteredcellperims = or(totalfilteredcellperims,Results_CellAnalysis(v).LogPerim);
    end
end

figimage1 = imoverlay(DAPI.*100,imdilate(totalfilteredcellperims,strel('disk',1)),'r'); %hold on;
figimage1 = imoverlay(figimage1,imdilate(DAPI_Watershed_Perim,strel('disk',1)),'b');

figimage2 = imoverlay(CellMask.*25,imdilate(totalfilteredcellperims,strel('disk',1)),'r'); %hold on;
figimage2 = imoverlay(figimage2,imdilate(DAPI_Watershed_Perim,strel('disk',1)),'b');

if aSMA_Analysis ==1, 
    figimage3 = aSMA.*100;
    figimage3 = imfuse(figimage3,aSMAGradImageTotal.*100,'ColorChannels',[2 1 2]);
else
    figimage3 = Blank3D(:,:,1);
end

figimage4 = Blank3D(:,:,1:3);

C = [figimage1 figimage2; figimage3 figimage4];
imshow(C);
title('Segmentation and Analysis Summary');

drawnow; hold off;

else
end

%% Results %%

clearvars CellularMeans CellularSums CMAreas aSMAGradientMean;
disp('Collating results...')

for r = 1:size(Results_CellAnalysis,2)
    CellularMeans(r,1) = Results_CellAnalysis(r).MeanCh1;
    CellularMeans(r,2) = Results_CellAnalysis(r).MeanCh2;
    if Channels>2, CellularMeans(r,3) = Results_CellAnalysis(r).MeanCh3; else end;
    if Channels>3, CellularMeans(r,4) = Results_CellAnalysis(r).MeanCh4; else end;
    CellularSums(r,1) = Results_CellAnalysis(r).SumCh1;
    CellularSums(r,2) = Results_CellAnalysis(r).SumCh2;
    if Channels >2, CellularSums(r,3) = Results_CellAnalysis(r).SumCh3; else end;
    if Channels >3, CellularSums(r,4) = Results_CellAnalysis(r).SumCh4; else end;
    CMAreas(r,1) = Results_CellAnalysis(r).CellArea;    
       
    if aSMA_Analysis == 1, aSMAGradientMean(r,1) = Grad_Stats(r).MeanGradientValue;
    else
    end
end

if Results(f).ZFocus == 1 || Results(f).ZFocus == ZPlanes,
    Results(f).NucNearestNeighborDistance = "";
    Results(f).TotalNuclei = 0;
    Results(f).NucleiPerGroup = 0;
    Results(f).AdjustedTotalNuclei = 0;
    Results(f).AdjustedNucGroupSizes = "";
    Results(f).MeanIntensities = ""; %Each column is a channel.
    Results(f).SumIntensities = ""; %Each column is a channel.
    Results(f).CMNumber = 0;
    Results(f).CMAreas = "";
    Results(f).CMTotalArea = 0;
    if TF_Analysis==1, Results(f).TFNucCytoRatios = 0;
    else end
    if aSMA_Analysis>0, Results(f).aSMAGradientMeans = 0;
    else end
else
    Results(f).NucNearestNeighborDistance = NearestNucDistanceFiltered;
    Results(f).TotalNuclei = sum(NucleiNumber,'all');
    Results(f).NucleiPerGroup = NucleiNumber;
    Results(f).AdjustedTotalNuclei = sum(AdjustedNucleiNumber,'all');
    Results(f).AdjustedNucGroupSizes = AdjustedNucleiNumber;
    Results(f).MeanIntensities = CellularMeans; %Each column is a channel.
    Results(f).SumIntensities = CellularSums; %Each column is a channel.
    Results(f).CMNumber = size(Results_CellAnalysis,2);
    Results(f).CMAreas = CMAreas;
    Results(f).CMTotalArea = sum(CMAreas);
    if TF_Analysis>0, Results(f).TFNucCytoRatios = NucCytoRatio_sorted;
    else end
    if aSMA_Analysis>0, Results(f).aSMAGradientMeans = aSMAGradientMean;
    else end
end

if FigSave == 1
    disp('Saving Figure...');
    if f == 1,
        cd(Folder); mkdir('Analysis'); cd Analysis;
    else cd(Folder); cd Analysis;
    end
    ax = gca;
    Fig_Name = append(filename(end-17:end-9),' Segmentation.tif');
    exportgraphics(ax,Fig_Name,'ContentType','image','Resolution','400');
elseif f == START,
    cd(Folder); mkdir('Analysis'); cd Analysis;
elseif f > START
    cd(Folder); cd Analysis;
end
    catch
         warning('An error occurred during analysis. Saving Results and skipping to next image.'); pause(2);
         Results = Results(START:f);
         cd(Folder); cd Analysis;
         save('AnalysisResults.mat','Results', '-v7.3');
         if FigShow == 1,
             imshow(imadjust(Focus_Ch1));
         else end
    end
end %End of Analysis Loop

cd(Folder); cd Analysis;
save('AnalysisResults.mat','Results', '-v7.3');

    %% Result Sorting By Well%%

disp('Finished analyzing images, now sorting results by well...');
[Results_Sorted,WellNumber] = WellSort(Channels,Results,ZPlanes,TF_Analysis,aSMA_Analysis);
cd(Folder); cd Analysis;
save('AnalysisResultsSortedByWell.mat','Results_Sorted','-v7.3');
    
    %% Sorting Results by Cluster Status %%

if SplitClusterResults == 1,
    disp('Sorting results by cluster status...');
    [Results_Sorted_InClusters,Results_Sorted_OutClusters] = ClusterSort(Channels,Results_Sorted,WellNumber,MinClusterSize,TF_Analysis,aSMA_Analysis);

cd(Folder); cd Analysis;
    save('AnalysisResultsSortedByWellandWithinClusters.mat','Results_Sorted_InClusters','-v7.3');
    save('AnalysisResultsSortedByWellandOutsideClusters.mat','Results_Sorted_OutClusters','-v7.3');
else end

toc