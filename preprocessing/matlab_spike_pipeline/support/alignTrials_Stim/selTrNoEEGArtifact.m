function selTrNoEEGArtifact(path,fileName,mazeSess)
% select trials without EEG artifacts

    fullPath = [path fileName '_alignRun_msess' num2str(mazeSess) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _TrialType1_alignRun file does not exist');
        return;
    end    
    load(fullPath,'trialsRun');
    nTrials = length(trialsRun.goodTrial);
    
    fullPath = [path fileName '_behPar.mat'];
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end  
    load(fullPath,'behPar');
    artFiltThre = 700;
    artThre = 1500;
    
    GlobalConst;
    thetaFiltParam.order = 2;
    thetaFiltParam.band = [4 10];
    [b,a] = butter(thetaFiltParam.order,thetaFiltParam.band/sampleFq*2);
    indTrNoArt = zeros(1,nTrials);
    for i = 1:nTrials
        if(isempty(trialsRun.eeg{i}))
            continue;
        end
        Eegf = filtfilt(b,a,trialsRun.eeg{i});
        indArt = Eegf > artFiltThre; % & trialsRun.eeg{i} > artThre;
        if(sum(indArt) == 0)
            indTrNoArt(i) = 1;
        end 
    end
    behPar.indTrNoArt = indTrNoArt;
    
    fullPath = [path fileName '_behPar.mat'];
    save(fullPath,'behPar','artFiltThre');
end
