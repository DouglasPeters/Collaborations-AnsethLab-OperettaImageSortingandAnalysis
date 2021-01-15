function [Results_Sorted_InClusters,Results_Sorted_OutClusters] = ClusterSort(Channels,Results_Sorted,WellNumber,MinClusterSize,TF_Analysis,aSMA_Analysis)

for w = 1:WellNumber    
    clearvars Results_InClusters_Filter Results_OutClusters_Filter NNClusterINFilt NNClusterOUTFilt;
    Results_InClusters_Filter = Results_Sorted(w).AdjustedNucGroupSizes >= MinClusterSize;
    Results_InClusters_Find = find(Results_InClusters_Filter);
    Results_OutClusters_Filter = Results_Sorted(w).AdjustedNucGroupSizes < MinClusterSize;
    Results_OutClusters_Find = find(Results_OutClusters_Filter);
        
    Results_Sorted_InClusters(w).Well = Results_Sorted(w).Well;
    if sum(Results_InClusters_Filter)>0,
        for i = 1:sum(Results_InClusters_Filter)
            Results_Sorted_InClusters(w).CMAreas(i,1) = Results_Sorted(w).CMAreas(Results_InClusters_Find(i,1),1);
            Results_Sorted_InClusters(w).NumberCMObjects = "";
            Results_Sorted_InClusters(w).PercentofAllCMObjects = "";
            Results_Sorted_InClusters(w).TotalCMArea = "";
            Results_Sorted_InClusters(w).PercentofAllCMArea = "";
            Results_Sorted_InClusters(w).NucleiPerGroup(i,1) = Results_Sorted(w).NucleiPerGroup(Results_InClusters_Find(i,1),1);
            Results_Sorted_InClusters(w).AdjustedNucGroupSizes(i,1) = Results_Sorted(w).AdjustedNucGroupSizes(Results_InClusters_Find(i,1),1);
            Results_Sorted_InClusters(w).TotalNuclei = "";
            Results_Sorted_InClusters(w).PercentofAllNuclei = "";
            Results_Sorted_InClusters(w).NearestNucDistanceFiltered = "";
            Results_Sorted_InClusters(w).MeanCh1(i,1) = Results_Sorted(w).MeanCh1(Results_InClusters_Find(i,1),1);
            Results_Sorted_InClusters(w).MeanCh2(i,1) = Results_Sorted(w).MeanCh2(Results_InClusters_Find(i,1),1);
            if Channels >2, Results_Sorted_InClusters(w).MeanCh3(i,1) = Results_Sorted(w).MeanCh3(Results_InClusters_Find(i,1),1); else end;
            if Channels >3, Results_Sorted_InClusters(w).MeanCh4(i,1) = Results_Sorted(w).MeanCh4(Results_InClusters_Find(i,1),1); else end;
            Results_Sorted_InClusters(w).SumCh1(i,1) = Results_Sorted(w).SumCh1(Results_InClusters_Find(i,1),1);
            Results_Sorted_InClusters(w).SumCh2(i,1) = Results_Sorted(w).SumCh2(Results_InClusters_Find(i,1),1);
            if Channels >2, Results_Sorted_InClusters(w).SumCh3(i,1) = Results_Sorted(w).SumCh3(Results_InClusters_Find(i,1),1); else end;
            if Channels >3, Results_Sorted_InClusters(w).SumCh4(i,1) = Results_Sorted(w).SumCh4(Results_InClusters_Find(i,1),1); else end;
            if TF_Analysis == 1, Results_Sorted_InClusters(w).TFNucCytoRatios = Results_Sorted(w).TFNucCytoRatios(Results_Sorted(w).NearestNucDistance(:,2)>=MinClusterSize); else end;
            if aSMA_Analysis == 1, Results_Sorted_InClusters(w).aSMAGradientMeans(i,1) = Results_Sorted(w).aSMAGradientMeans(Results_InClusters_Find(i,1),1); else end;
            
            Results_Sorted_InClusters(w).NumberCMObjects = numel(Results_Sorted_InClusters(w).CMAreas);
            Results_Sorted_InClusters(w).PercentofAllCMObjects = Results_Sorted_InClusters(w).NumberCMObjects/Results_Sorted(w).TotalNumberofCMObjects*100;
            Results_Sorted_InClusters(w).TotalCMArea = sum(Results_Sorted_InClusters(w).CMAreas);
            Results_Sorted_InClusters(w).PercentofAllCMArea = Results_Sorted_InClusters(w).TotalCMArea/sum(Results_Sorted(w).CMTotalArea)*100;
            Results_Sorted_InClusters(w).TotalNuclei = sum(Results_Sorted(w).NucleiPerGroup(Results_InClusters_Find));
            Results_Sorted_InClusters(w).PercentofAllNuclei = Results_Sorted_InClusters(w).TotalNuclei/sum(Results_Sorted(w).TotalNuclei)*100;
            Results_Sorted_InClusters(w).NearestNucDistanceFiltered = Results_Sorted(w).NearestNucDistance(Results_Sorted(w).NearestNucDistance(:,2)>=MinClusterSize);
        end
    else end
        
    Results_Sorted_OutClusters(w).Well = Results_Sorted(w).Well;
    if sum(Results_OutClusters_Filter)>0
        for j = 1:sum(Results_OutClusters_Filter)
            Results_Sorted_OutClusters(w).CMAreas(j,1) = Results_Sorted(w).CMAreas(Results_OutClusters_Find(j,1),1);
            Results_Sorted_OutClusters(w).NumberCMObjects = "";
            Results_Sorted_OutClusters(w).PercentofAllCMObjects = "";
            Results_Sorted_OutClusters(w).TotalCMArea = "";
            Results_Sorted_OutClusters(w).PercentofAllCMArea = "";
            Results_Sorted_OutClusters(w).NucleiPerGroup(j,1) = Results_Sorted(w).NucleiPerGroup(Results_OutClusters_Find(j,1),1);
            Results_Sorted_OutClusters(w).AdjustedNucGroupSizes(j,1) = Results_Sorted(w).AdjustedNucGroupSizes(Results_OutClusters_Find(j,1),1);
            Results_Sorted_OutClusters(w).TotalNuclei = "";
            Results_Sorted_OutClusters(w).PercentofAllNuclei = "";
            Results_Sorted_OutClusters(w).NearestNucDistanceFiltered = "";
            Results_Sorted_OutClusters(w).MeanCh1(j,1) = Results_Sorted(w).MeanCh1(Results_OutClusters_Find(j,1),1);
            Results_Sorted_OutClusters(w).MeanCh2(j,1) = Results_Sorted(w).MeanCh2(Results_OutClusters_Find(j,1),1);
            if Channels >2, Results_Sorted_OutClusters(w).MeanCh3(j,1) = Results_Sorted(w).MeanCh3(Results_OutClusters_Find(j,1),1); else end;
            if Channels >3, Results_Sorted_OutClusters(w).MeanCh4(j,1) = Results_Sorted(w).MeanCh4(Results_OutClusters_Find(j,1),1); else end;
            Results_Sorted_OutClusters(w).SumCh1(j,1) = Results_Sorted(w).SumCh1(Results_OutClusters_Find(j,1),1);
            Results_Sorted_OutClusters(w).SumCh2(j,1) = Results_Sorted(w).SumCh2(Results_OutClusters_Find(j,1),1);
            if Channels >2, Results_Sorted_OutClusters(w).SumCh3(j,1) = Results_Sorted(w).SumCh3(Results_OutClusters_Find(j,1),1); else end;
            if Channels >3, Results_Sorted_OutClusters(w).SumCh4(j,1) = Results_Sorted(w).SumCh4(Results_OutClusters_Find(j,1),1); else end;
            if TF_Analysis == 1, Results_Sorted_OutClusters(w).TFNucCytoRatios = Results_Sorted(w).TFNucCytoRatios(Results_Sorted(w).NearestNucDistance(:,2)<MinClusterSize); else end;
            if aSMA_Analysis == 1, Results_Sorted_OutClusters(w).aSMAGradientMeans(j,1) = Results_Sorted(w).aSMAGradientMeans(Results_OutClusters_Find(j,1),1); else end;
        
            Results_Sorted_OutClusters(w).NumberCMObjects = numel(Results_Sorted_OutClusters(w).CMAreas);
            Results_Sorted_OutClusters(w).PercentofAllCMObjects = Results_Sorted_OutClusters(w).NumberCMObjects/Results_Sorted(w).TotalNumberofCMObjects*100;
            Results_Sorted_OutClusters(w).TotalCMArea = sum(Results_Sorted_OutClusters(w).CMAreas);
            Results_Sorted_OutClusters(w).PercentofAllCMArea = Results_Sorted_OutClusters(w).TotalCMArea/sum(Results_Sorted(w).CMTotalArea)*100;
            Results_Sorted_OutClusters(w).TotalNuclei = sum(Results_Sorted(w).NucleiPerGroup(Results_OutClusters_Find));
            Results_Sorted_OutClusters(w).PercentofAllNuclei = Results_Sorted_OutClusters(w).TotalNuclei/sum(Results_Sorted(w).TotalNuclei)*100;
            Results_Sorted_OutClusters(w).NearestNucDistanceFiltered = Results_Sorted(w).NearestNucDistance(Results_Sorted(w).NearestNucDistance(:,2)<MinClusterSize);
        end
    else
    end
end