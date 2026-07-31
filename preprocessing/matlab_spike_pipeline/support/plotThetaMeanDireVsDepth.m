function thetaMeanDireVsDepth = plotThetaMeanDireVsDepth(spikeThetaPhaseStruct,depthNeu,selPyramidal,indSuper)
% plot theta modulation against the depth of the cell in the pyramidal
% layer

    figure
    meanDire = spikeThetaPhaseStruct.meanDire;
    ind = meanDire < 0;
    meanDire(ind) = meanDire(ind) + 2*pi;
    thetaMeanDireVsDepth.meanDire = meanDire;
    
%     % histogram
%     meanDireTmp = [meanDire,meanDire+2*pi,meanDire+4*pi];
%     thetaMeanDireVsDepth.stepPhase = 0:0.01*pi:6*pi;
%     thetaMeanDireVsDepth.histMeanDire = ...
%         hist(meanDireTmp,thetaMeanDireVsDepth.stepPhase);
%     % smooth histogram
%     std = 3;
%     paramC.gaussFilt = gaussFilter(12*std, std);
%     lenGaussKernel = length(paramC.gaussFilt);
%     normFactor = sum(paramC.gaussFilt);
%     paramC.gaussFilt = paramC.gaussFilt./normFactor;
%     thetaMeanDireVsDepth.filMeanDire = ...
%             conv(thetaMeanDireVsDepth.histMeanDire,paramC.gaussFilt);
%     thetaMeanDireVsDepth.filMeanDire = ...
%             thetaMeanDireVsDepth.filMeanDire(floor(lenGaussKernel/2)+1:...
%                         (end-lenGaussKernel+floor(lenGaussKernel/2)+1));
%     indPhase = find(thetaMeanDireVsDepth.stepPhase >= 2*pi ...
%                     & thetaMeanDireVsDepth.stepPhase <4*pi);
%     filMeanDire = thetaMeanDireVsDepth.filMeanDire(indPhase); 
%     % correct theta mean direction of individual cells based on the trough 
%     % of the meandire histogram
%     [~,indMinPhase] = min(filMeanDire);
%     thetaMeanDireVsDepth.troughDire = ...
%         thetaMeanDireVsDepth.stepPhase(indPhase(indMinPhase)) - 2*pi;
% %     if(thetaMeanDireVsDepth.troughDire > pi)
% %         indShift = meanDire > thetaMeanDireVsDepth.troughDire;
% %         meanDire(indShift) = meanDire(indShift) - 2*pi;
% %     else
% %         indShift = meanDire < thetaMeanDireVsDepth.troughDire;
% %         meanDire(indShift) = meanDire(indShift) + 2*pi;
% %     end
    thetaMeanDireVsDepth.meanDire = meanDire;
    thetaMeanDireVsDepth.meanDireSupFit = ...
        fitlm(depthNeu.relDepthNeuHDef(selPyramidal & indSuper),...
        meanDire(selPyramidal & indSuper));

    hold on;
    plot(depthNeu.relDepthNeuHDef(selPyramidal),...
        meanDire(selPyramidal),'rx')
    plotAdded(thetaMeanDireVsDepth.meanDireSupFit);
    T = thetaMeanDireVsDepth.meanDireSupFit.anova;
    pValue = T.Variables;
    pValue = pValue(1,end);
    xlabel('depth')
    ylabel('thetaMeanDire (based on MeanDire)')
    title(['p = ' num2str(pValue)]);
        %%% significant parameter
        
%     % calculate population phase histogram
%     histPhaseAll = zeros(1,length(spikeThetaPhaseStruct.posPhase));
%     for i = 1:length(indSelPyr)
%         histNorm = spikeThetaPhaseStruct.histPhaseFilPerNeuron{indSelPyr(i)}...
%             /max(spikeThetaPhaseStruct.histPhaseFilPerNeuron{indSelPyr(i)});
%         histPhaseAll = histPhaseAll + histNorm;
%     end
%     stepPhase = 5;
%     thetaModVsDepth.stepPhase1 = ...
%         [stepPhase/2:stepPhase:1080-stepPhase/2]/360*2*pi;
%     indPhase = find(thetaModVsDepth.stepPhase1 >= 2*pi ...
%                     & thetaModVsDepth.stepPhase1 <4*pi);
%     thetaModVsDepth.histPhase = histPhaseAll;
%     histPhase = histPhaseAll(indPhase);
%     % estimate the trough of the population phase histogram, and correct
%     % theta mean direction of individual cells correspondingly
%     [~,indMinPhase] = min(histPhase);
%     thetaModVsDepth.troughPhase = ...
%         thetaModVsDepth.stepPhase1(indPhase(indMinPhase)) - 2*pi;
%     if(thetaModVsDepth.troughPhase > pi)
%         indShift = meanDire > thetaModVsDepth.troughPhase;
%         meanDire(indShift) = meanDire(indShift) - 2*pi;
%     else
%         indShift = meanDire < thetaModVsDepth.troughPhase;
%         meanDire(indShift) = meanDire(indShift) + 2*pi;
%     end
%     thetaModVsDepth.meanDireAdj1 = meanDire;
%     thetaModVsDepth.meanDireFit1 = ...
%         fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
%               meanDire(selPyramidal));
%     subplot(2,2,2)
%     plotAdded(thetaModVsDepth.meanDireFit1);
%     T = thetaModVsDepth.meanDireFit1.anova;
%     pValue = T.Variables;
%     pValue = pValue(1,end);
%     title(['p = ' num2str(pValue)]);
%     xlabel('depth')
%     ylabel('thetaMeanDire (based on population theta phase)')
end
