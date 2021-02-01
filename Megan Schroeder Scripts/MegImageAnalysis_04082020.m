clear all; clc;
tic
cd(userpath);

inDir = 'C:\Users\dplad\Desktop\meg test\*.czi';
Folder = 'C:\Users\dplad\Desktop\meg test\';
Channels = 2;
%BinWidth = 500;

erodestrel = strel('disk',1);
dilatestrel = strel('disk',5);

srcFiles = dir(inDir);
%BinIdx = 1:BinWidth:65535;

%Ch1_DATA = zeros([length(BinIdx) length(srcFiles)]);
%Ch2_DATA = zeros([length(BinIdx) length(srcFiles)]);
%if Channels>2, Ch3_DATA = zeros([length(BinIdx) length(srcFiles)]); else end
%if Channels>3, Ch4_DATA = zeros([length(BinIdx) length(srcFiles)]); else end

for f = 1:1 %length(srcFiles)
    progress = (f/length(srcFiles)*100);
    progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,length(srcFiles),progress,'%');
    disp(progress2)
    filename = strcat(Folder,srcFiles(f).name);
    
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

Ch1_IntDen(f,1) = sum(Ch1_lin(:)); Ch1_MeanInt(f,1) = mean(Ch1_lin(:));
Ch2_IntDen(f,1) = sum(Ch2_lin(:)); Ch2_MeanInt(f,1) = mean(Ch2_lin(:));
if Channels>2, Ch3_IntDen(f,1) = sum(Ch3_lin(:)); Ch3_MeanInt(f,1) = mean(Ch3_lin(:)); else end
if Channels>3, Ch4_IntDen(f,1) = sum(Ch4_lin(:)); Ch4_MeanInt(f,1) = mean(Ch4_lin(:)); else end

for sumx = 1:Res
    for sumy = 1:Res
    Ch1_sum(sumx,sumy) = sum(Ch1(sumx,sumy,:));
    Ch2_sum(sumx,sumy) = sum(Ch2(sumx,sumy,:));
    if Channels>2, Ch3_sum(sumx,sumy) = sum(Ch3(sumx,sumy,:)); else end
    if Channels>3, Ch4_sum(sumx,sumy) = sum(Ch4(sumx,sumy,:)); else end
    end
end

Ch1_sumd = uint16(Ch1_sum);
Ch2_sumd = uint16(Ch2_sum);
if Channels>2, Ch3_sumd = uint16(Ch3_sum); else end
if Channels>3, Ch4_sumd = uint16(Ch3_sum); else end

%%%%%%%%%%%%%%%%%%%
Ch1_BW1 = imbinarize(Ch1_sumd);
Ch1_BW2 = imerode(Ch1_BW1,erodestrel);
Ch1_BW3 = imdilate(Ch1_BW2,dilatestrel);

Ch1_NEAR = Ch1_sumd(Ch1_BW3); Ch1_NEAR_MEAN(f,1) = mean(Ch1_NEAR);
Ch1_FAR = Ch1_sumd(~Ch1_BW3); Ch1_FAR_MEAN(f,1) = mean(Ch1_FAR);
Ch2_NEAR = Ch2_sumd(Ch1_BW3); Ch2_NEAR_MEAN(f,1) = mean(Ch2_NEAR);
Ch2_FAR = Ch2_sumd(~Ch1_BW3); Ch2_FAR_MEAN(f,1) = mean(Ch2_FAR);


%Ch1_Histo = histcounts(Ch1_lin,BinIdx,'Normalization','probability');
 %   Ch1_DATA(1:length(BinIdx),1) = reshape(BinIdx,[length(BinIdx) 1]);
  %  Ch1_DATA(1:length(Ch1_Histo),f+1) = reshape(Ch1_Histo,[length(Ch1_Histo) 1]);
%Ch2_Histo = histcounts(Ch2_lin,BinIdx,'Normalization','probability');
 %   Ch2_DATA(1:length(BinIdx),1) = reshape(BinIdx,[length(BinIdx) 1]);
  %  Ch2_DATA(1:length(Ch2_Histo),f+1) = reshape(Ch2_Histo,[length(Ch2_Histo) 1]);
%if Channels>2, Ch3_Histo = histcounts(Ch3_lin,BinIdx,'Normalization','probability'); 
 %   Ch3_DATA(1:length(BinIdx),1) = reshape(BinIdx,[length(BinIdx) 1]);
  %  Ch3_DATA(1:length(Ch3_Histo),f+1) = reshape(Ch3_Histo,[length(Ch3_Histo) 1]);
   % else end
%if Channels>3, Ch4_Histo = histcounts(Ch4_lin,BinIdx,'Normalization','probability'); 
 %   Ch4_DATA(1:length(BinIdx),1) = reshape(BinIdx,[length(BinIdx) 1]);
  %  Ch4_DATA(1:length(Ch4_Histo),f+1) = reshape(Ch4_Histo,[length(Ch4_Histo) 1]);
   % else end

end

Ch2_NEARFAR_MEANS(:,1) = Ch2_NEAR_MEAN(:,1);
Ch2_NEARFAR_MEANS(:,2) = Ch2_FAR_MEAN(:,1);
Ch2_NEARFAR_MEANS(:,3) = Ch2_NEARFAR_MEANS(:,1)./Ch2_NEARFAR_MEANS(:,2);

