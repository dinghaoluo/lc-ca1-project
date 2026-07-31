function plotSpikePhaseFreqRunOnset(path, fileName, onlyRun, mazeSess)
% plot spike phase and theta frequency over time

    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    fullPath = [path fileName '_Info.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'beh');
        
    GlobalConst;
    ind = beh.mazeSess == mazeSess;
    if(isfield(beh,'indStimLap'))
        trialNo = find(beh.indStimLap(ind) == 0);
    else
        trialNo = 1:length(behPar.indTrBadBeh);
    end
    
    indBadBeh = behPar.indTrBadBeh;
    
    goodTrials = intersect(find(indBadBeh == 0),trialNo);
    badTrials = intersect(find(indBadBeh == 1),trialNo);
    disp(['Number of good trials = ' num2str(length(goodTrials))]);
    
    plotSpikePhase_aligned(path, fileName, onlyRun, mazeSess, goodTrials, badTrials);
end

function plotSpikePhase_aligned(path, fileName, onlyRun, mazeSess, trialNoG, trialNoB, neuronNo)
% plot spikes rasters across run segments
% E.G.: plotSpikeRaster_aligned('./','A002-20181005-01_DataStructure_mazeSection1_TrialType',1)

    if(nargin == 4)
        trialNoG = [];
        trialNoB = [];
        neuronNo = [];
    elseif(nargin == 5)
        trialNoB = [];
        neuronNo = [];
    elseif(nargin == 6)
        neuronNo = [];
    end

    %%%%%%%%% load recording file
    fullPath = [path fileName '.mat']; 
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
        return;
    end
    load(fullPath,'cluList');
    if(isempty(neuronNo))
        neuronNo = 1:length(cluList.all);
    end
       
    fullPath = [path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT file does not exist');
        return;
    end
    load(fullPath);
    fullPath = [path fileName '_alignRun_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'trialsRun');
            
    if(isempty(trialNoG))
        trialNoG = 1:size(trialsRunSpikes.Time,2);
        trialNoB = [];
    end
    
    GlobalConst;
    
    trialLenT = 7; %sec
    vLen = length(-nSampBef+1:trialLenT*sampleFq);
    totThetaFreqGoodTr = zeros(length(trialNoG),vLen);
    totThetaAmpGoodTr = zeros(length(trialNoG),vLen);
    totSpeedGoodTr = zeros(length(trialNoG),vLen);
    
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels') 
    figTitle = 'Spikes vs Time';
%             figTitle = 'Spikes vs Dist';
    set(figure(figNew),'OuterPosition',...
        [pos(1) pos(2)-300 pos(3)*1.5 pos(4)*1.7],'Name',figTitle)
        
    subplot(3,2,1)
    hold on;
    n = 1;
    for j = trialNoG   
%         disp(['Trial ' num2str(j)])
        if(isempty(trialsRun.ThetaFreq{j}))
            continue;
        end
        vSamp = [-nSampBef+1:trialsRun.numSamples(j)]/sampleFq;
        vFreq = [trialsRun.ThetaFreqBef{j}' trialsRun.ThetaFreq{j}'];
        plot(vSamp,vFreq,'k.','MarkerSize',3);  
        if(length(vSamp) > vLen)
            totThetaFreqGoodTr(n,:) = vFreq(1:vLen);
        else
            totThetaFreqGoodTr(n,1:length(vSamp)) = vFreq;
        end
        n = n+1;
    end
    set(gca, 'XLim', [-nSampBef/sampleFq trialLenT],'YLim',[-10 40]);
    ylabel('Theta freq (Hz)');
       
    subplot(3,2,2)
    hold on;
    for j = trialNoB   
%         disp(['Trial ' num2str(j)])
        if(isempty(trialsRun.ThetaFreq{j}))
            continue;
        end
        plot([-nSampBef+1:trialsRun.numSamples(j)]/sampleFq,...
             [trialsRun.ThetaFreqBef{j}' trialsRun.ThetaFreq{j}'],'k.','MarkerSize',3);    
    end
    set(gca, 'XLim', [-nSampBef/sampleFq trialLenT],'YLim',[-10 40]);
    
    subplot(3,2,3)
    hold on;
    n = 1;
    for j = trialNoG   
%         disp(['Trial ' num2str(j)])
        if(isempty(trialsRun.ThetaFreq{j}))
            continue;
        end     
        vSamp = [-nSampBef+1:trialsRun.numSamples(j)]/sampleFq;
        vAmp = [trialsRun.ThetaAmpBef{j}' trialsRun.ThetaAmp{j}'];
        plot(vSamp,vAmp,'k.','MarkerSize',3); 
        if(length(vSamp) > vLen)
            totThetaAmpGoodTr(n,:) = vAmp(1:vLen);
        else
            totThetaAmpGoodTr(n,1:length(vSamp)) = vAmp;
        end      
        n = n+1;
    end
    set(gca, 'XLim', [-nSampBef/sampleFq trialLenT]);
    ylabel('Theta amp');
    xlabel('Time (s)');
    
    subplot(3,2,4)
    hold on;
    for j = trialNoB   
%         disp(['Trial ' num2str(j)])
        if(isempty(trialsRun.ThetaFreq{j}))
            continue;
        end
        plot([-nSampBef+1:trialsRun.numSamples(j)]/sampleFq,...
             [trialsRun.ThetaAmpBef{j}' trialsRun.ThetaAmp{j}'],'k.','MarkerSize',3);     
    end
    set(gca, 'XLim', [-nSampBef/sampleFq trialLenT]);
    ylabel('Theta amp');
    xlabel('Time (s)');
    
    subplot(3,2,5)
    hold on;
    n = 1;
    for j = trialNoG   
%         disp(['Trial ' num2str(j)])
        if(isempty(trialsRun.speed_MMsec{j}))
            continue;
        end
        vSamp = [-nSampBef+1:trialsRun.numSamples(j)]/sampleFq;
        speed = [trialsRun.speed_MMsecBef{j}' trialsRun.speed_MMsec{j}'];
        indSpeed = speed < -1;
        speed(indSpeed) = 0;
        plot(vSamp,speed,'k.','MarkerSize',3);     
        if(length(vSamp) > vLen)
            totSpeedGoodTr(n,:) = speed(1:vLen);
        else
            totSpeedGoodTr(n,1:length(vSamp)) = speed;
        end
        n = n+1;
    end
    set(gca, 'XLim', [-nSampBef/sampleFq trialLenT]);
    ylabel('Speed');
    xlabel('Time (s)');
    
    subplot(3,2,6)
    hold on;
    for j = trialNoB   
        disp(['Trial ' num2str(j)])
        if(isempty(trialsRun.speed_MMsec{j}))
            continue;
        end
        vSamp = [-nSampBef+1:trialsRun.numSamples(j)]/sampleFq;
        speed = [trialsRun.speed_MMsecBef{j}' trialsRun.speed_MMsec{j}'];
        indSpeed = speed < -1;
        speed(indSpeed) = 0;
        plot(vSamp,speed,'k.','MarkerSize',3);     
    end
    set(gca, 'XLim', [-nSampBef/sampleFq trialLenT]);
    ylabel('Speed');
    xlabel('Time (s)');
    
end
