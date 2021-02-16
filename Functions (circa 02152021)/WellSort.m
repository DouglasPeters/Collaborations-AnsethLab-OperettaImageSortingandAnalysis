function [Results_Sorted,WellNumber] = WellSort(Channels,Results,ZPlanes,TF_Analysis,aSMA_Analysis)

for s = 1:size(Results,2)
    NoCellFilter(s,1) = ZPlanes>Results(s).ZFocus && Results(s).ZFocus>1;
    if Results(s).TotalNuclei > 0
        NoCellFilter(s,1) = 1;
    else
        NoCellFilter(s,1) = 0;
    end
end

Results_Filtered = Results(NoCellFilter);

for i = 1:size(Results_Filtered,2)
    FieldsofView(i,1) = string(Results_Filtered(i).FieldofView);
    WellsMeasured(i,1) = string(Results_Filtered(i).Well);
end

WellsUnique = unique(WellsMeasured);
WellNumber = length(WellsUnique);

for w = 1:WellNumber
    Results_Sorted(w).Well = WellsUnique(w,1);
    NumFields(w,1) = sum(WellsMeasured == WellsUnique(w,1));
    Results_Sorted(w).NumFields = NumFields(w,1);
    for f = (sum(NumFields(1:(w-1)))+1):sum(NumFields(1:w,1))
        if f==(sum(NumFields(1:(w-1)))+1)
            Results_Sorted(w).TotalNumberofCMObjects = Results_Filtered(f).CMNumber;
            Results_Sorted(w).TotalNumberofNuclei = Results_Filtered(f).TotalNuclei;
            Results_Sorted(w).TotalAdjustedNumberofNuclei = Results_Filtered(f).AdjustedTotalNuclei;
            Results_Sorted(w).CMNumberObjects(Results_Filtered(f).FieldofView,1) = Results_Filtered(f).CMNumber;
            Results_Sorted(w).CMTotalArea(Results_Filtered(f).FieldofView,1) = Results_Filtered(f).CMTotalArea;
            Results_Sorted(w).CMAreas(1:size(Results_Filtered(f).CMAreas,1),1) = Results_Filtered(f).CMAreas;
            Results_Sorted(w).NearestNucDistance(1:size(Results_Filtered(f).NucNearestNeighborDistance,1),1:2) = Results_Filtered(f).NucNearestNeighborDistance(:,1:2);
            Results_Sorted(w).TotalNuclei(Results_Filtered(f).FieldofView,1) = Results_Filtered(f).TotalNuclei;
            Results_Sorted(w).NucleiPerGroup(1:numel(Results_Filtered(f).NucleiPerGroup),1) = Results_Filtered(f).NucleiPerGroup;
            Results_Sorted(w).AdjustedTotalNuclei(Results_Filtered(f).FieldofView,1) = Results_Filtered(f).AdjustedTotalNuclei;
            Results_Sorted(w).AdjustedNucGroupSizes(1:numel(Results_Filtered(f).AdjustedNucGroupSizes),1) = Results_Filtered(f).AdjustedNucGroupSizes;
            Results_Sorted(w).MeanCh1(1:size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,1);
            Results_Sorted(w).MeanCh2(1:size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,2);
            if Channels >2, Results_Sorted(w).MeanCh3(1:size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,3); else end
            if Channels >3, Results_Sorted(w).MeanCh4(1:size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,4); else end
            Results_Sorted(w).SumCh1(1:size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,1);
            Results_Sorted(w).SumCh2(1:size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,2);
            if Channels >2, Results_Sorted(w).SumCh3(1:size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,3);else end
            if Channels >3, Results_Sorted(w).SumCh4(1:size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,4); else end
            if TF_Analysis == 1, Results_Sorted(w).TFNucCytoRatios = Results_Filtered(f).TFNucCytoRatios; else end
            if aSMA_Analysis == 1
                Results_Sorted(w).aSMAGradientMeans = Results_Filtered(f).aSMAGradientMeans; 
%                 for c = 1:Results_Filtered(f).CMNumber
%                     if c == 1
%                         Results_Sorted(w).aSMAGradientMeansNucAdjusted(1:Results_Sorted(w).AdjustedNucGroupSizes(c),1) = Results_Sorted(w).aSMAGradientMeans(c);
%                     else
%                         Results_Sorted(w).aSMAGradientMeansNucAdjusted(end+1:end+Results_Sorted(w).AdjustedNucGroupSizes(c),1) = Results_Sorted(w).aSMAGradientMeans(c);
%                     end
%                 end
            else
            end
            
        else
            Results_Sorted(w).TotalNumberofCMObjects = Results_Filtered(f).CMNumber + Results_Sorted(w).TotalNumberofCMObjects;
            Results_Sorted(w).TotalNumberofNuclei = Results_Filtered(f).TotalNuclei + Results_Sorted(w).TotalNumberofNuclei;
            Results_Sorted(w).TotalAdjustedNumberofNuclei = Results_Filtered(f).AdjustedTotalNuclei + Results_Sorted(w).TotalAdjustedNumberofNuclei;
            Results_Sorted(w).CMNumberObjects(Results_Filtered(f).FieldofView,1) = Results_Filtered(f).CMNumber;
            Results_Sorted(w).CMTotalArea(Results_Filtered(f).FieldofView,1) = Results_Filtered(f).CMTotalArea;
            Results_Sorted(w).CMAreas(end+1:end+size(Results_Filtered(f).CMAreas,1),1) = Results_Filtered(f).CMAreas;
            Results_Sorted(w).NearestNucDistance(end+1:end+size(Results_Filtered(f).NucNearestNeighborDistance,1),1:2) = Results_Filtered(f).NucNearestNeighborDistance(:,1:2);
            Results_Sorted(w).TotalNuclei(Results_Filtered(f).FieldofView,1) = Results_Filtered(f).TotalNuclei;
            Results_Sorted(w).NucleiPerGroup(end+1:end+numel(Results_Filtered(f).NucleiPerGroup),1) = Results_Filtered(f).NucleiPerGroup;
            Results_Sorted(w).AdjustedTotalNuclei(Results_Filtered(f).FieldofView,1) = Results_Filtered(f).AdjustedTotalNuclei;
            Results_Sorted(w).AdjustedNucGroupSizes(end+1:end+numel(Results_Filtered(f).AdjustedNucGroupSizes),1) = Results_Filtered(f).AdjustedNucGroupSizes;
            Results_Sorted(w).MeanCh1(end+1:end+size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,1);
            Results_Sorted(w).MeanCh2(end+1:end+size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,2);
            if Channels >2, Results_Sorted(w).MeanCh3(end+1:end+size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,3); else end
            if Channels >3, Results_Sorted(w).MeanCh4(end+1:end+size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,4); else end
            Results_Sorted(w).SumCh1(end+1:end+size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,1);
            Results_Sorted(w).SumCh2(end+1:end+size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,2);
            if Channels >2, Results_Sorted(w).SumCh3(end+1:end+size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,3); else end
            if Channels >3, Results_Sorted(w).SumCh4(end+1:end+size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,4); else end
            if TF_Analysis == 1, Results_Sorted(w).TFNucCytoRatios(end+1:end+size(Results_Filtered(f).TFNucCytoRatios,1),1:2) = Results_Filtered(f).TFNucCytoRatios; else end
            if aSMA_Analysis == 1
                Results_Sorted(w).aSMAGradientMeans(end+1:end+size(Results_Filtered(f).aSMAGradientMeans,1),1) = Results_Filtered(f).aSMAGradientMeans;
%                 for c = 1:Results_Filtered(f).CMNumber
%                     Results_Sorted(w).aSMAGradientMeansNucAdjusted(end+1:end+Results_Sorted(w).AdjustedNucGroupSizes(c),1) = Results_Sorted(w).aSMAGradientMeans(c);
%                 end
            else
            end
        end
    end
end

if aSMA_Analysis == 1
    for w = 1:WellNumber
        for f = (sum(NumFields(1:(w-1)))+1):sum(NumFields(1:w,1))
            for c = 1:Results_Filtered(f).CMNumber
                if c==1 && f == (sum(NumFields(1:(w-1)))+1)
                    Results_Sorted(w).aSMAGradientMeansNucAdjusted(1:Results_Filtered(f).AdjustedNucGroupSizes(c),1) = Results_Filtered(f).aSMAGradientMeans(c);
                else
                    Results_Sorted(w).aSMAGradientMeansNucAdjusted(end+1:end+Results_Filtered(f).AdjustedNucGroupSizes(c),1) = Results_Filtered(f).aSMAGradientMeans(c);
                end
            end
        end
    end
else
end
    
    clearvars ZeroMask;
    for e = 1:size(Results_Sorted(w).CMNumberObjects,1)
    ZeroMask(e,1) = Results_Sorted(w).CMNumberObjects(e,1) > 0;
    end
    Results_Sorted(w).CMNumberObjects = Results_Sorted(w).CMNumberObjects(ZeroMask);
    Results_Sorted(w).CMTotalArea = Results_Sorted(w).CMTotalArea(ZeroMask);
    Results_Sorted(w).TotalNuclei = Results_Sorted(w).TotalNuclei(ZeroMask);
    Results_Sorted(w).AdjustedTotalNuclei = Results_Sorted(w).AdjustedTotalNuclei(ZeroMask);
end