clear all; clc;
tic
cd(userpath);

inDir = 'C:\Users\dplad\Desktop\megdapi\*.czi';
Folder = 'C:\Users\dplad\Desktop\megdapi\';
srcFiles = dir(inDir);

BWAreaMin = 100;

for f = 1:length(srcFiles)
    progress = (f/length(srcFiles)*100);
    progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,length(srcFiles),progress,'%');
    disp(progress2)
    filename = strcat(Folder,srcFiles(f).name);
    filetype = filename(end-3:end);
    cd(Folder);
   
    if filetype == '.lsm',
        Channels = 1;
        I = bfopen(filename);
    elseif filetype == '.czi',
        Channels = 2;
        I = bfopen(filename);
    else
    end
              
    Res = length(I{1,1}{1,1});
    Slices = (length(I{1,1})/Channels);
    Blank = zeros(Res,Res,Slices);
    Planes = Channels*Slices;
    Ch1 = uint16(Blank);    

    for i = 1:Slices
        Ch1_planes(i,1) = 1+(Channels*i-Channels);
        Ch1(:,:,i) = I{1,1}{Ch1_planes(i,1),1};
    end
    
    for x = 1:Res
        for y = 1:Res
            Ch1_MIP(x,y) = max(Ch1(x,y,:));
        end
    end
        
    BW1 = imbinarize(Ch1_MIP);
    BW2 = imdilate(imerode(BW1,strel('disk',2)),strel('disk',2));
    BW3 = bwareaopen(BW2,BWAreaMin);
    %imshow(BW3);
    cc = bwconncomp(BW3,4);
    Nuclei_ImageCount(f,1) = cc.NumObjects;
end
toc