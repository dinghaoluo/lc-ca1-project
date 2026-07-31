function plotChangeWithDepth(path, filename)

    load([path '/' filename '_Depth.mat']);
    load([path '/' filename '_ThetaMod_Run1.mat']);
    load([path '/' filename '_ThetaPhaseL_Run1.mat']);
    load([path '/' filename '_Info.mat']);
    load([path '/' filename '.mat'],'cluList');
    load([path '/' filename '_burstAll_THL_Run1.mat']);
    load([path '/' filename '_ThetaMod_Run1.mat']);

    % select pyramidal neurons
    selPyramidal = autoCorr.isPyrneuron == 1 & cluList.firingRate > 0.3;
    indSelPyr = find(selPyramidal == 1);
    
    % theta modulation
    figure
    md.thetaMod = fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
        thetaModSess{1}.thetaMod(selPyramidal));
    plotAdded(md.thetaMod);
    xlabel('depth')
    ylabel('thetaMod')
   
    figure
    histPhaseAll = zeros(1,length(spikeThetaPhaseStruct.posPhase));
    for i = 1:length(indSelPyr)
        histNorm = spikeThetaPhaseStruct.histPhaseFilPerNeuron{indSelPyr(i)}...
            /max(spikeThetaPhaseStruct.histPhaseFilPerNeuron{indSelPyr(i)});
        histPhaseAll = histPhaseAll + histNorm;
    end
    stepPhase = 5;
    histPhase = [stepPhase/2:stepPhase:1080-stepPhase/2]/360*2*pi;
    indPhase = find(histPhase >= 2*pi & histPhase <4*pi);
    histPhase = histPhaseAll(indPhase);
    [~,indMinPhase] = min(histPhase);
    troughPhase = histPhase(indPhase(indMinPhase)) - 2*pi;
    if(troughPhase > pi)
        indShift = meanDire > troughPhase;
        meanDire(indShift) = meanDire(indShift) - 2*pi;
    else
        indShift = meanDire < troughPhase;
        meanDire(indShift) = meanDire(indShift) + 2*pi;
    end
    md.thetaMeanDire1 = fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
        meanDire(selPyramidal));
    plotAdded(md.thetaMeanDire1);
    xlabel('depth')
    ylabel('thetaMeanDire (based on population theta phase)')
    
    figure
    md.thetaMeanResultantLen = fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
        spikeThetaPhaseStruct.meanResultantLen(selPyramidal));
    plotAdded(md.thetaMeanResultantLen);
    xlabel('depth')
    ylabel('thetaMeanResultantLen')
        %%% significant parameter
        
    figure
    md.thetaMeanResultantLen = fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
        spikeThetaPhaseStruct.meanResultantLen(selPyramidal));
    plotAdded(md.thetaMeanResultantLen);
    xlabel('depth')
    ylabel('thetaMeanResultantLen')
        %%% significant parameter

    % ACG spectrum    
    figure
    md.filACGSpectrumThetaPeak = fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
        thetaModSess{1}.filACGSpectrumThetaPeak(selPyramidal));
    plotAdded(md.filACGSpectrumThetaPeak);
    xlabel('depth')
    ylabel('filACGSpectrumThetaPeak')

    figure
    plot(depthNeu.relDepthNeuHDef(selPyramidal),...
        thetaModSess{1}.ACGSpectrumThetaPeak(selPyramidal),'o')
    xlabel('depth')
    ylabel('ACGSpectrumThetaPeak')

    figure
    plot(depthNeu.relDepthNeuHDef(selPyramidal),...
        thetaModSess{1}.filACGSpectrumDeltaPeak(selPyramidal),'o')
    xlabel('depth')
    ylabel('filACGSpectrumDeltaPeak')

    figure
    plot(depthNeu.relDepthNeuHDef(selPyramidal),...
        thetaModSess{1}.filACGSpectrumBetaPeak(selPyramidal),'o')
    xlabel('depth')
    ylabel('filACGSpectrumBetaPeak')

    figure
    plot(depthNeu.relDepthNeuHDef(selPyramidal),...
        thetaModSess{1}.filACGSpectrumGammaPeak(selPyramidal),'o')
    xlabel('depth')
    ylabel('filACGSpectrumGammaPeak')

    figure
    plot(depthNeu.relDepthNeuHDef(selPyramidal),...
        thetaModSess{1}.filACGSpectrumHighGammaPeak(selPyramidal),'o')
    xlabel('depth')
    ylabel('filACGSpectrumHighGammaPeak')

    % burstiness
    figure
    meanDire = burstIsiPerNeuron.meanDire;
    ind = burstIsiPerNeuron.meanDire < 0;
    meanDire(ind) = meanDire(ind) + 2*pi;
    md.burstMeanDire = fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
        meanDire(selPyramidal));
    plotAdded(md.burstMeanDire);
    xlabel('depth')
    ylabel('burstMeanDire')
        %%% significant parameter
        
    figure
    meanResultantLenNonBurst = burstIsiPerNeuron.meanResultantLenNonBurst;
    md.meanResultantLenNonBurst = fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
        meanResultantLenNonBurst(selPyramidal));
    plotAdded(md.meanResultantLenNonBurst);
    xlabel('depth')
    ylabel('meanResultantLenNonBurst')
            
    figure
    meanDireNonBurst = burstIsiPerNeuron.meanDireNonBurst;
    ind = burstIsiPerNeuron.meanDireNonBurst < 0;
    meanDireNonBurst(ind) = meanDireNonBurst(ind) + 2*pi;
    md.nonburstMeanDire = fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
        meanDireNonBurst(selPyramidal));
    plotAdded(md.nonburstMeanDire);
    xlabel('depth')
    ylabel('nonburstMeanDire')
        %%% significant parameter

    figure
    meanResultantLenStart = burstIsiPerNeuron.meanResultantLenStart;
    ind = burstIsiPerNeuron.meanResultantLenStart < 0;
    meanResultantLenStart(ind) = meanResultantLenStart(ind) + 2*pi;
    md.meanResultantLenStart = fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
        meanResultantLenStart(selPyramidal));
    plotAdded(md.meanResultantLenStart);
    xlabel('depth')
    ylabel('meanResultantLenStart')
        %%% significant parameter

    for i = 1:length(burstIsiPerNeuron.meanDire)
        meanNumSpikesPerBurst(i) = mean(burstIsiPerNeuron.numSpPerBurst{i});
        meanISIBurst(i) = mean(burstIsiPerNeuron.isiBurst{i});
        meanPercBurst(i) = ...
            length(burstIsiPerNeuron.numSpPerBurst{i})...
            /burstIsiPerNeuron.numSp(i);
    end
    
    figure
    md.meanNumSpikesPerBurst = fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
        meanNumSpikesPerBurst(selPyramidal));
    plotAdded(md.meanNumSpikesPerBurst);
    xlabel('depth')
    ylabel('meanNumSpikesPerBurst')

    figure
    md.meanISIBurst = fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
        meanISIBurst(selPyramidal));
    plotAdded(md.meanISIBurst);
    xlabel('depth')
    ylabel('meanISIBurst')

    figure
    md.meanPercBurst = fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
        meanPercBurst(selPyramidal));
    plotAdded(md.meanPercBurst);
    xlabel('depth')
    ylabel('meanPercBurst')

    % autocorrelation
    figure
    md.peakTo20ms = fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
        autoCorr.peakTo20ms(selPyramidal));
    plotAdded(md.peakTo20ms);
    xlabel('depth')
    ylabel('peakTo20ms autoCorr')
    
    fileNameDepthPlot = [path '/' filename '_DepthPlots.mat'];
    save(fileNameDepthPlot,'md');