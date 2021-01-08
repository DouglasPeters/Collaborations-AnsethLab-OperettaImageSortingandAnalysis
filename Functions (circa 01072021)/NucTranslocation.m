function [TFAnalysis,NucCytoRatio_sorted] = NucTranslocation(TF,Results_CellAnalysis,Blank3D)

for r = 1:size(Results_CellAnalysis,2)
    BlankALLNUC = false(size(Blank3D(:,:,1)));
    for n = 1:size(Results_CellAnalysis(r).NucleiCC.PixelIdxList,2)
        BlankNUC = false(size(Blank3D(:,:,1)));
        BlankNUC(Results_CellAnalysis(r).NucleiCC.PixelIdxList{n}) = true;
        TFAnalysis(r).NUCmean(n,1) = mean(TF(BlankNUC),'all');
        
        BlankALLNUC(Results_CellAnalysis(r).NucleiCC.PixelIdxList{n}) = true;
    end
    
    BlankCYTO = Results_CellAnalysis(r).LogImage;
    CytoMASK = logical(BlankCYTO - BlankALLNUC);
    TFAnalysis(r).CYTOmean = mean(TF(CytoMASK),'all');
    
    for m = 1:size(TFAnalysis(r).NUCmean,1)
        TFAnalysis(r).NucCytoRatio(m,1) = TFAnalysis(r).NUCmean(m,1)/TFAnalysis(r).CYTOmean;
    end
    
    if r == 1
        NucCytoRatio_sorted(1:size(TFAnalysis(r).NUCmean,1),1) = TFAnalysis(r).NucCytoRatio(:);
        NucCytoRatio_sorted(1,2) = Results_CellAnalysis(r).NucleiNumber;
    else
        NucCytoRatio_sorted(end+1:end+size(TFAnalysis(r).NUCmean,1),1) = TFAnalysis(r).NucCytoRatio(:);
        NucCytoRatio_sorted(r,2) = Results_CellAnalysis(r).NucleiNumber; 
    end
end