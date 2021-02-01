clear all; clc;
tic
%% User Inputs %%

Folder = 'C:\Users\dplad\Desktop\megdata\male gels 20x, cacl2\';  %Tells the script where to find your images (must keep the trailing slash).

BotCrop = 16; %Number of pixels to remove from the bottom of the image to avoid the weird camera error.
Threshold = 100; %8-bit threshold value (0-255) to segment VK signal from background. Arbitrary...
ShowImages = 0; %Do you want to show images of the VK segmentation? (Yes=1, No=0)

%% Don't touch this %%
inDir = append(Folder,'*.jpg'); %Don't touch this.
srcFiles = dir(inDir);
varNames = {'Filename','VKArea','VKMean','VKIntDen'};
varTypes = {'string','double','double','double'};
Result_Table = table('Size',[length(srcFiles) length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);

for q = 1:2 %length(srcFiles)
    filename = strcat(Folder,srcFiles(q).name);
    I = imread(filename);
    I2 = I(1:end-BotCrop,:,:);
      
    Result_Table.Filename(q) = filename;
   
    I3 = rgb2gray(I2);
    I4 = imcomplement(im2uint8(I3));
    I5 = I4 > Threshold;
    I6 = I4(I5);

    Result_Table.VKArea(q) = sum(I5,'all');
    Result_Table.VKMean(q) = mean(I6,'all');
    Result_Table.VKIntDen(q) = sum(I6,'all');  

if ShowImages == 1
    B = bwboundaries(I5,8);
    imshow(I3); 
    hold on
    for k = 1:size(B,1)
        b = B{k};
        plot(b(:,2),b(:,1),'r','LineWidth',1);
    end
    hold off
    drawnow;
else
end
end

matfilename = strcat(Folder,' VK Analysis Results.mat');
save(matfilename,'Result_Table');

toc

