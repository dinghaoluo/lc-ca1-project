function CalCCG(path,fileName,onlyRun)
% Calculate CCG for each subsession
% 
% by Yingxue 8/25/2017

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least three arguments are needed for this function.');
        return;
    elseif nargin == 2
        onlyRun = 1;
    elseif nargin > 3
        disp('Too many input arguments');        
        return;
    end
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameCCG = [fileName '_CCG_Run' num2str(onlyRun) '.mat'];
    
    ind = findstr(fileName,'_');
    fileNameBehEphys = [fileName(1:ind(1)-1) '_Behav2PDataLFP.mat'];
    fullPath = [path fileNameBehEphys];
    if(exist(fullPath) == 0)
        disp('The _Behav2PDataLFP file does not exist');
        return;
    end
    load(fullPath,'Track');
    
    GlobalConst2P;
    
    spikes = Track.spks;
    if(onlyRun == 1)
       indSpeed = Track.speed_MMsec <= minSpeed;      
       spikes(indSpeed,:) = 0; 
    end
   
    %%%%%%%%% initialize constants
    fileNameInfo = [fileName '_Info.mat'];
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
    
    % calculate CCG for each subsession
    disp('Calculate CCG for each subsession')
    HalfBins = 1000;   
    CCGSess = cell(1,length(mazeSess));
    for i = 1:length(mazeSess)
        disp(['Session ' num2str(i)]);
        spikesTmp = spikes(Track.mazeSess == mazeSess(i) & Track.lapID > startTrNo,:);
          % total number of bins on each side
        [CCGSess{i}.ccgVal, CCGSess{i}.ccgT] = ...
            xcorr(spikesTmp,HalfBins); 
    end
    
    fullPath = [path fileNameCCG];
    save(fullPath, 'CCGSess','-v7.3');
                   
    clear mydata;

end
