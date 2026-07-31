function plotSpikeThetaFreqAmpRunOnsetMeanStd(path, fileName, onlyRun, mazeSess)
% plot spike phase and theta frequency over time

    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    plotSpikePhase_aligned(path, fileName, onlyRun, mazeSess);
end

function plotSpikePhase_aligned(path, fileName, onlyRun, mazeSess, neuronNo)
% plot spikes rasters across run segments
% E.G.: plotSpikeRaster_aligned('./','A002-20181005-01_DataStructure_mazeSection1_TrialType',1)

    if(nargin == 4)
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
    
    fullPath = [path fileName '_Info.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'beh');
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    fullPath = [path fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','trialNoStimCtrl','pulseMeth');
                               
    GlobalConst;
    
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
           
    GlobalConst;
    
    trialLenT = 7; %sec
    meanThetaNonStimGood = calTheta(trialsRun,trialNoNonStimGood',trialLenT,nSampBef,sampleFq);
    meanThetaNonStimBad = calTheta(trialsRun,trialNoNonStimBad',trialLenT,nSampBef,sampleFq);
    meanThetaStim = [];
    for i = 1:length(trialNoStim)
        meanThetaStim = calTheta(trialsRun,trialNoStim{i}',trialLenT,nSampBef,sampleFq);
    end
    t = (-nSampBef+1:trialLenT*sampleFq)/sampleFq;
    
    fullPath = [path fileName '_ThetaAlignRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']; 
    save(fullPath,'meanThetaNonStimGood','meanThetaNonStimBad','meanThetaStim','t');
    
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels') 
    figTitle = 'Theta Freq vs Time mean';
    set(figure(figNew),'OuterPosition',...
        [pos(1) pos(2) 280 210],'Name',figTitle)

    options.handle     = figure(figNew);
    options.color_area = [27 117 187]./255;    % Blue theme
    options.color_line= [ 39 169 225]./255;
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'std';
    options.x_axis = [-nSampBef+1:trialLenT*sampleFq]/sampleFq;
    plot_areaerrorbarSub(meanThetaNonStimGood.totThetaFreqGoodTr,options);
    set(gca,'YLim',[0 20],'XLim',[-nSampBef/sampleFq trialLenT])
    xlabel('Time (s)')
    ylabel('Theta freq (Hz)')
    
    ind = strfind(fileName,'_');
    fileName1 = ['Z:\Yingxue\DataAnalysisRaphi\' ...
        fileName(1:ind(1)) 'ThetaFreqVsTAlignedToRun' num2str(onlyRun)];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    %% 
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels') 
    figTitle = 'Theta Amp vs Time mean';
    set(figure(figNew),'OuterPosition',...
        [pos(1) pos(2) 280 210],'Name',figTitle)

    options.handle     = figure(figNew);
    plot_areaerrorbarSub(meanThetaNonStimGood.totThetaAmpGoodTr,options);
    set(gca,'YLim',[-100 800],'XLim',[-nSampBef/sampleFq trialLenT])
    xlabel('Time (s)')
    ylabel('Theta amp')
    
    fileName1 = ['Z:\Yingxue\DataAnalysisRaphi\' ...
        fileName(1:ind(1)) 'ThetaAmpVsTAlignedToRun' num2str(onlyRun)];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    %%
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels') 
    figTitle = 'Speed vs Time mean';
    set(figure(figNew),'OuterPosition',...
        [pos(1) pos(2) 280 210],'Name',figTitle)

    options.handle     = figure(figNew);
    plot_areaerrorbarSub(meanThetaNonStimGood.totSpeedGoodTr,options);
    set(gca,'YLim',[-100 600],'XLim',[-nSampBef/sampleFq trialLenT])
    xlabel('Time (s)')
    ylabel('Speed')
    
    fileName1 = ['Z:\Yingxue\DataAnalysisRaphi\' ...
        fileName(1:ind(1)) 'SpeedVsTAlignedToRun' num2str(onlyRun)];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')

end

function meanTheta =calTheta(trialsRun,trialNoG,trialLenT,nSampBef,sampleFq)
    vLen = length(-nSampBef+1:trialLenT*sampleFq);
    meanTheta.totThetaFreqGoodTr = zeros(length(trialNoG),vLen);
    meanTheta.totThetaAmpGoodTr = zeros(length(trialNoG),vLen);
    meanTheta.totSpeedGoodTr = zeros(length(trialNoG),vLen);
    
    n = 1;
    for j = trialNoG   
        if(isempty(trialsRun.ThetaFreq{j}))
            continue;
        end
        vSamp = [-nSampBef+1:trialsRun.numSamples(j)];
        vFreq = [trialsRun.ThetaFreqBef{j}' trialsRun.ThetaFreq{j}'];
        if(length(vSamp) > vLen)
            meanTheta.totThetaFreqGoodTr(n,:) = vFreq(1:vLen);
        else
            meanTheta.totThetaFreqGoodTr(n,1:length(vSamp)) = vFreq;
        end
        n = n+1;
    end
    
    n = 1;
    for j = trialNoG  
        if(isempty(trialsRun.ThetaFreq{j}))
            continue;
        end
        vSamp = [-nSampBef+1:trialsRun.numSamples(j)]/sampleFq;
        vAmp = [trialsRun.ThetaAmpBef{j}' trialsRun.ThetaAmp{j}'];
        if(length(vSamp) > vLen)
            meanTheta.totThetaAmpGoodTr(n,:) = vAmp(1:vLen);
        else
            meanTheta.totThetaAmpGoodTr(n,1:length(vSamp)) = vAmp;
        end      
        n = n+1;
    end
    
    n = 1;
    for j = trialNoG   
        if(isempty(trialsRun.ThetaFreq{j}))
            continue;
        end
        vSamp = [-nSampBef+1:trialsRun.numSamples(j)]/sampleFq;
        speed = [trialsRun.speed_MMsecBef{j}' trialsRun.speed_MMsec{j}'];
        indSpeed = speed < -1;
        speed(indSpeed) = 0;
        if(length(vSamp) > vLen)
            meanTheta.totSpeedGoodTr(n,:) = speed(1:vLen);
        else
            meanTheta.totSpeedGoodTr(n,1:length(vSamp)) = speed;
        end
        n = n+1;
    end
    
    [b,a] = butter(2,0.01); % low pass filtered phase histogram
    meanTheta.meanThetaFreqGoodTr = mean(meanTheta.totThetaFreqGoodTr);
    if(length(meanTheta.meanThetaFreqGoodTr) > 1)
        meanTheta.meanThetaFreqFiltGoodTr = filtfilt(b,a,meanTheta.meanThetaFreqGoodTr); 
    else
        meanTheta.meanThetaFreqFiltGoodTr = meanTheta.meanThetaFreqGoodTr;
    end
    meanTheta.stdThetaFreqGoodTr = std(meanTheta.totThetaFreqGoodTr)/sqrt(length(trialNoG));
    meanTheta.meanThetaAmpGoodTr = mean(meanTheta.totThetaAmpGoodTr);
    if(length(meanTheta.meanThetaFreqGoodTr) > 1)
        meanTheta.meanThetaAmpFiltGoodTr = filtfilt(b,a,meanTheta.meanThetaAmpGoodTr); 
    else
        meanTheta.meanThetaAmpFiltGoodTr = meanTheta.meanThetaAmpGoodTr;
    end
    meanTheta.stdThetaAmpGoodTr = std(meanTheta.totThetaAmpGoodTr)/sqrt(length(trialNoG));
    meanTheta.meanSpeedGoodTr = mean(meanTheta.totSpeedGoodTr);
    meanTheta.stdSpeedGoodTr = std(meanTheta.totSpeedGoodTr)/sqrt(length(trialNoG));
end
