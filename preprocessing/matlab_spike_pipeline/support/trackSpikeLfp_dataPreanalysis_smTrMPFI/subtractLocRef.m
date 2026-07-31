function [] =subtractLocRef(baseName,mvDat)

%-----------------------------------------------
% This fce loads one dat file and subtracts the mean of all traces from all of the ch =
% as if locally referenced and saves the new file as *.datLocRef
% (see Ludwig et al., 2006)
% input: data folder location
%
% structure of a dat file:
% 16bit
% one file for all channels:
% first: all first samples of all channels
% second: all second samples of all channels
% .....
% the last samples of all channels
%-----------------------------------------------
% e.g.: subtractLocRef('A029-20200626-01',1)

fileList = dir([baseName '.dat']);
if isempty(fileList)
    printf('Any dat file with the specified basename was not found - exiting.....');
    return;
end

% load meta file to find out total N of ch
if exist([baseName '.meta'], 'file') == 2
    fid = fopen([baseName '.meta']);
    metaInfo=textscan(fid, '%s');
    metaInfo=metaInfo{1};
    fclose(fid);
else
    disp([baseName '.meta was not found - exiting.']);
    return;
end

% meta file structure
% 3-Aplipex;
% 4- version
% 9 - max amplitude range
% 14 - min amplitude range
% 21 - file length (sec)
% 26 - file size (bytes)
% 29 - file name
% 32 - gain
% 38 - N of recorded ch
% 43 - record start date - DAY of a week
% 44 - record start date - MONTH
% 45 - record start date - DAY-DATE
% 46 - record start date - YEAR
% 51 - recrod start time - hr:min:sec
% 55 - sampling rate
% 61 - sha1 code of the file
% 73 - T end of recording (ms)
% 85 - T start of recording (ms)

if(mvDat ~= 2)
    totalNch = str2double(metaInfo(38));
    samplRate = str2double(metaInfo(55));
else
    totalNch = str2double(metaInfo(36));
    samplRate = str2double(metaInfo(48));
end
        
% open output dat file
    filenameOut = [baseName '-locRef.dat'];
    fidOut = fopen(filenameOut, 'w');
    disp(['Output file:    ' filenameOut]);
    disp('   ');
    
% load dat file, subtract the ref ch and save as *.datLocRef
buffersize = samplRate * 10;
numel=0;
totCh = 64;
plotFig=0;
for nf = 1:length(fileList)
    datafile = fileList(nf).name;
    fid = fopen(datafile);
    totFileSize = fileList(nf).bytes / 2;
    while ~feof(fid),
        [dat,count] = fread(fid,[totalNch,buffersize],'int16');
        numelm = count/totalNch;
        numel = numel+count;
        [nch ndat] = size(dat);        
%         newRef = mean(dat, 1);
        newRef = median(dat(1:totCh,:), 1);
        locRefDat = dat;        
        for n = 1:totCh 
            locRefDat(n,:) = dat(n,:) - newRef;
            if plotFig == 1                
                figure(2); cla;
                plot(dat(n,:), 'k'); hold on;
                plot(newRef, 'b');                
                plot(locRefDat(n,:), 'r');  
                pause;
            end
        end    

        fwrite(fidOut, locRefDat, 'int16');
        fprintf('Perc of the file written: %d \n', (numel/totFileSize)*100);
        per = 1;
        if mod((numel/totFileSize)*100, per) > 1            
            fprintf('Perc of the file written: %d', (numel/totFileSize)*100);
            per = per + 1;
        end
        
    end
    fclose(fid);
    
    if(mvDat >= 1)
        system(['mv ' datafile ' ' datafile '.orig']);
        system(['ln -s ' datafile(1:end-4) '-locRef.dat ' datafile]);
    end
end

disp('done'); 
    
return; 
