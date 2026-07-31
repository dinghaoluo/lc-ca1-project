function plotThetaModVsDepth(path, filename)
    
    load([path '/' filename '_Depth.mat']);
    load([path '/' filename '_ThetaMod_Run1.mat']);
    load([path '/' filename '_ThetaPhaseL_Run1.mat']);
    load([path '/' filename '_Info.mat']);
    load([path '/' filename '.mat'],'cluList');
    load([path '/' filename '_burstAll_THL_Run1.mat']);
    
    std = 3;
    paramC.gaussFilt = gaussFilter(12*std, std);
    lenGaussKernel = length(paramC.gaussFilt);
    normFactor = sum(paramC.gaussFilt);
    paramC.gaussFilt = paramC.gaussFilt./normFactor;
    
    % select pyramidal neurons
    selPyramidal = autoCorr.isPyrneuron == 1 & cluList.firingRate > 0.3;
%                     & spikeThetaPhaseStruct.pOmnibus < 0.01 ...
%                     & spikeThetaPhaseStruct.pRayleigh < 0.01;
    indSelPyr = find(selPyramidal == 1);
    indSuper = depthNeu.relDepthNeuHDef >= -5;
    indDeep =  depthNeu.relDepthNeuHDef < -5;
    
    thetaPeakVsDepth = ...
        plotSpikePeakThetaPhaseVsDepth(spikeThetaPhaseStruct,depthNeu,selPyramidal);
    
    thetaMeanDireVsDepth = ...
        plotThetaMeanDireVsDepth(spikeThetaPhaseStruct,depthNeu,selPyramidal,indSuper);
    
    % burstiness
    figure;
    meanDire = burstIsiPerNeuron.meanDire;
    ind = burstIsiPerNeuron.meanDire < 0;
    meanDire(ind) = meanDire(ind) + 2*pi;
%     if(thetaModVsDepth.troughPhase > pi)
%         indShift = meanDire > thetaModVsDepth.troughDire;
%         meanDire(indShift) = meanDire(indShift) - 2*pi;
%     else
%         indShift = meanDire < thetaModVsDepth.troughDire;
%         meanDire(indShift) = meanDire(indShift) + 2*pi;
%     end
    thetaModVsDepth.burstMeanDireAdj = meanDire;
    thetaModVsDepth.burstMeanDireFit = ...
        fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
              meanDire(selPyramidal));
    plotAdded(thetaModVsDepth.burstMeanDireFit);
    T = thetaModVsDepth.burstMeanDireFit.anova;
    pValue = T.Variables;
    pValue = pValue(1,end);
    title(['p = ' num2str(pValue)]);
    xlabel('depth')
    ylabel('Mean dire. burst')
        %%% significant parameter
        
    % nonburst spikes
    figure
    meanDireNonBurst = burstIsiPerNeuron.meanDireNonBurst;
    ind = burstIsiPerNeuron.meanDireNonBurst < 0;
    meanDireNonBurst(ind) = meanDireNonBurst(ind) + 2*pi;
%     if(thetaModVsDepth.troughPhase > pi)
%         indShift = meanDireNonBurst > thetaModVsDepth.troughDire;
%         meanDireNonBurst(indShift) = meanDireNonBurst(indShift) - 2*pi;
%     else
%         indShift = meanDireNonBurst < thetaModVsDepth.troughDire;
%         meanDireNonBurst(indShift) = meanDireNonBurst(indShift) + 2*pi;
%     end
    thetaModVsDepth.nonburstMeanDireAdj = meanDireNonBurst;
    thetaModVsDepth.nonburstMeanDireFit = ...
        fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
              meanDireNonBurst(selPyramidal));
    T = thetaModVsDepth.nonburstMeanDireFit.anova;
    pValue = T.Variables;
    pValue = pValue(1,end);
    plotAdded(thetaModVsDepth.nonburstMeanDireFit);
    title(['p = ' num2str(pValue)]);
    xlabel('depth')
    ylabel('Mean dire. nonburst')
        %%% significant parameter
    
    if(~exist('dataAnalysis','dir'))
        mkdir([path '\dataAnalysis']);
    end
    fullName = [path '/dataAnalysis/' filename '_thetaModVsDepth.mat'];
    save(fullName,'thetaModVsDepth');