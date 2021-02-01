function [Grad_Stats,aSMAGradImageTotal] = aSMAActivation(aSMA,Results_CellAnalysis)

    aSMAGradImageTotal = uint16(Results_CellAnalysis(1).LogImage);
for r = 1:size(Results_CellAnalysis,2)
    BlankCYTO = logical(Results_CellAnalysis(r).LogImage);
    aSMA_gradient = imgradient(aSMA);
    aSMA_masked = aSMA_gradient(BlankCYTO);
    maskforimage = uint16(Results_CellAnalysis(r).LogImage);
    Grad_Stats(r).GradientArray = aSMA_masked;
    Grad_Stats(r).MaskedImage = aSMA.*maskforimage;
    Grad_Stats(r).MaskedGradientImage = uint16(aSMA_gradient).*maskforimage;
    Grad_Stats(r).MeanGradientValue = mean(Grad_Stats(r).GradientArray,'all');
    Grad_Stats(r).MaxGradientValue = max(Grad_Stats(r).GradientArray);
    aSMAGradImageTotal = aSMAGradImageTotal + Grad_Stats(r).MaskedGradientImage;
end
    
