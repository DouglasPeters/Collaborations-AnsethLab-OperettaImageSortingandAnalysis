clear all; clc;
tic
cd(userpath);

Folder = 'F:\MatLab_Exports\2020-11-10_FemaleVICs_Soft_Y-27632\2020-11-10_FemaleVICs_Soft_Y27632__2020-11-10T10_40_54-Measurement1\Images\';
%Where are the images located? This filepath should end with a folder name, not a file name.
%IMPORTANT: Use format 'FILEPATH/'. The apostrophes and ending slash are important.

cd(Folder);
srcFiles = dir('*.tiff');
FINISH = length(srcFiles);

for f = 1:FINISH
    if f == 1, disp('Collecting image information...'); else end;
    [pathstr,name,ext] = fileparts(srcFiles(f).name);
    
    location = strfind(name,'-');
    FileInfo(f,:).PlateLocation = name(1:12);
    FileInfo(f,:).Row = str2double(name(2:3));
    FileInfo(f,:).Column = str2double(name(5:6));
    FileInfo(f,:).Field = str2double(name(8:9));
    FileInfo(f,:).RCF = name(1:9);
    FileInfo(f,:).Plane = str2double(name(11:12));
    FileInfo(f,:).Channel = str2double(name(16));
    FileInfo(f,:).Remainder = name(17:end);
    Channel_Idx(f,1) = FileInfo(f).Channel;
 end

for f = 1:FINISH
    
    time(f,1).ElapsedSeconds = toc;
    clc
    filename = strcat(Folder,srcFiles(f).name);
    progress = ((FINISH-(FINISH-f))/FINISH)*100;
    progress2 = sprintf('On image %d of %d; %0.2f%c complete.',f,FINISH,progress,'%');
    disp(progress2);
    if progress < 10,
        disp('Estimated time remaining will display after 10% of images are analyzed...');
    else
        time(f).AverageSecondsPerLoop = time(f).ElapsedSeconds/(FINISH-(FINISH-f));
        time(f).EstimatedTotalSeconds = time(f).AverageSecondsPerLoop*(FINISH);
        time(f).EstimatedSecondsLeft = time(f).EstimatedTotalSeconds-time(f).ElapsedSeconds;
        time(f).EstimatedMinutesLeft = time(f).EstimatedSecondsLeft/60;
        time(f).EstimatedMinutesElapsed = time(f).ElapsedSeconds/60;
        estimate = sprintf('Run time: %0.2f minutes.',time(f).EstimatedMinutesElapsed);
        estimate2 = sprintf('Estimated time remaining: %0.2f minutes.',time(f).EstimatedMinutesLeft);
        disp(estimate);
        disp(estimate2);
    end
    cd(Folder);
    if f == 1
        mkdir('ImageStacks');          
    else
    end
        
    if f == 1 || size(strfind(FileInfo(f).RCF,FileInfo(f-1).RCF),1) > 0
        
        I = imread(srcFiles(f).name);
   
        if FileInfo(f).Channel == 1
           IMAGE(:,:,FileInfo(f).Plane,1) = I; else end;
        if max(Channel_Idx) > 1 && FileInfo(f).Channel == 2
           IMAGE(:,:,FileInfo(f).Plane,2) = I; else end;
        if max(Channel_Idx) > 2 && FileInfo(f).Channel == 3
           IMAGE(:,:,FileInfo(f).Plane,3) = I; else end;
        if max(Channel_Idx) > 3 && FileInfo(f).Channel == 4
           IMAGE(:,:,FileInfo(f).Plane,4) = I; else end;
    
    else
        cd ./ImageStacks;
        disp('Saving image stack...');
        SAVE = strcat(FileInfo(f-1).RCF,'.ome.tiff');
        bfsave(IMAGE(:,:,:,:),SAVE);
        cd ../;
        I = imread(srcFiles(f).name);
        clearvars IMAGE;
        IMAGE(:,:,FileInfo(f).Plane,1) = I;
    end
    
    if f == FINISH
        cd ./ImageStacks;
        SAVE = strcat(FileInfo(f).RCF,'.ome.tiff');
        bfsave(IMAGE,SAVE);
        cd ../;
    else end    
end

toc