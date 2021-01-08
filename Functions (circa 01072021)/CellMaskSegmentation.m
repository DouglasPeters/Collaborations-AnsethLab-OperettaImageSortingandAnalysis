

function [CM_IndCells,CMseg_cc,CMseg_props,CMseg_CentroidXY,CM_Watershed_BW2,CM_Watershed_Perim] = CellMaskSegmentation(Ch_CM,CM_SizeFilt,ResY,ResX)

% I_CM_bg = imopen(Ch_CM,strel('disk',50));
% Ch_CM2 = Ch_CM - I_CM_bg;
% I_CM_adapthisteq = imlocalbrighten(Ch_CM);
% I_CM_thresh = multithresh(I_CM_adapthisteq,3);
I_CM_thresh = mean(imopen(Ch_CM,strel('disk',10)),'all');
%I_CM_BW = I_CM_adapthisteq>I_CM_thresh(1);
I_CM_BW = Ch_CM>I_CM_thresh;
I_CM_BW2 = bwareaopen(imdilate(imerode(I_CM_BW,strel('disk',3)),strel('disk',2)),CM_SizeFilt,4);

CM_DistTransf1 = bwdist(~I_CM_BW2);
CM_DistTransf1 = imcomplement(CM_DistTransf1);
CM_DistTransfMask = imextendedmin(CM_DistTransf1,10);
CM_DistTransf2 = imimposemin(CM_DistTransf1,CM_DistTransfMask);
CM_Watershed = watershed(CM_DistTransf2);
CM_Watershed_BW = I_CM_BW2;
CM_Watershed_BW(CM_Watershed == 0) = 0;
CM_Watershed_BW2 = bwareaopen(CM_Watershed_BW,200,4);
CM_Watershed_Perim = bwperim(CM_Watershed_BW2);

CMseg_cc = bwconncomp(CM_Watershed_BW2,4);
CMseg_props = regionprops(CMseg_cc,Ch_CM,'basic','Centroid');

for h = 1:CMseg_cc.NumObjects
    CM_IndCells(h).LogImage = zeros(ResY,ResX);
    CM_IndCells(h).LogImage(CMseg_cc.PixelIdxList{h}) = 1;
    CM_IndCells(h).LogFilledHoles = logical(imfill(CM_IndCells(h).LogImage,'holes'));
    CM_IndCells(h).LogPerim = bwperim(CM_IndCells(h).LogFilledHoles);
    [row,col] = find(CM_IndCells(h).LogPerim);
    CM_IndCells(h).PerimRow = row;
    CM_IndCells(h).PerimCol = col;
    CM_IndCells(h).alphaShape = alphaShape(CM_IndCells(h).PerimCol,CM_IndCells(h).PerimRow);
%     CM_IndCells(h).inShape = inShape(CM_IndCells(h).alphaShape,48,741);
end

for g = 1:length(CMseg_props)
    CMseg_CentroidXY(g,1:2) = [CMseg_props(g).Centroid(1) CMseg_props(g).Centroid(2)];
end