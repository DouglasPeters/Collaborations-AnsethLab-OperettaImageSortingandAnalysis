% For von kossa images
% 
% clear %clear workspace
% clc %clear command window

filenames = dir('G:\For Doug analysis_20201203\Female\CaCl2\rep 1\*.jpg');

Resultlist = {'Pixel reduction factor','PixRedFactor',' ',' ','','',;...
    'Sobel threshold','SobelThresh',' ',' ','','',;...
    'Filename','Cell Number','Vonkossa','IntegratedDensity','factin (mean int)','aSMA/factin'}; %headerlines in results list for output

%%
hl = 3;
for q = 1:numel(filenames)
    img = imread(filenames(q).name);
%     aSMA = img(:,:,2);
%     factin = img(:,:,1);
%     dapi = img(:,:,3);
%     img = imread(filenames(q).name);
%     [r,b,g] = imsplit(imread(filenames(q).name));
%         ri = uint8(255) - r;
%         bi = uint8(255) - b;
%         gi = uint8(255) - g;
%     img = ri - gi;

      
    Resultlist{q+hl,1} = filenames(q,1).name;
    filenames(q,1).name
    %% find nuclui
    img1 = rgb2gray(img);
    img2 = im2uint8(img1);
    img3 = img2 > 750;
    img4 = uint8(255) - im2uint8(img3);
  %  bw = im2bw(dapi, graythresh(dapi));
    %img = rgb2gray(img);
    %bw = im2bw(img, graythresh(img));
   % bw = img > mean(img);
    %bw1 = imerode(bw,strel('disk',1));
    %bw2 = imclose(bw,strel('disk',5));
    %bw3 = imdilate(bw2, strel('disk',5));
    %inverseGrayImage = uint8(255) - im2uint8(bw3);
    %bw4 = im2bw(inverseGrayImage);
    %stats = regionprops(bw3,'Area','Centroid');
   % stats2 = regionprops(bw3,bi,'Area','MeanIntensity');
    %area = [stats.Area];
   % cells = area(area>50);
    
    %cellNumber = length(cells);
   % Resultlist{q+hl,2} = cellNumber;
  %  imshowpair(img,bw3);
%     subplot(2,2,2), imshow(EdU)
    %% Black around
  %  Int = [stats2.MeanIntensity];
    
    Black = mean2(img4);
    Int = sum(img4(:));
    Resultlist{q+hl,3} = Black;
    Resultlist{q+hl,4} = Int;

    %% EdU
%     bw2 = EdU > median(EdU,2)*3;
%     bw3 = imopen(bw2, strel('disk',2));
%     stats = regionprops(bw3,'Area');
%     pos = find([stats.Area] < 200);
%     EDUpositive = length(pos);
%     subplot(2,2,2), imshow(bw2)
%     subplot(2,2,3), imshow(imadjust(EdU))
%     subplot(2,2,4), imshow(bw3)
%     imshow(imadjust(EdU))
    %% Percent EdU positive

%     percentEDU = (EDUpositive/cellNumber)*100;
%     Resultlist{q+hl,5} = percentEDU;
%     Resultlist{q+hl,6} = EDUpositive;

    
end

