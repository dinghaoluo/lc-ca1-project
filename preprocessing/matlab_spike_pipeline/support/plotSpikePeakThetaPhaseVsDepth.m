function thetaPeakVsDepth = plotSpikePeakThetaPhaseVsDepth(spikeThetaPhaseStruct,depthNeu,selPyramidal)
    
    % spike theta phase peak for each cell
    figure
    maxPhase = spikeThetaPhaseStruct.maxPhaseFilArr/360*2*pi;
    ind = maxPhase < 0;
    maxPhase(ind) = maxPhase(ind) + 2*pi;
    thetaPeakVsDepth.maxPhase = maxPhase;
    thetaPeakVsDepth.maxPhaseFit = fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
        maxPhase(selPyramidal));
    plotAdded(thetaPeakVsDepth.maxPhaseFit);
    T = thetaPeakVsDepth.maxPhaseFit.anova;
    pValue = T.Variables;
    pValue = pValue(1,end);
    xlabel('depth')
    ylabel('thetaMaxPhase')
    title(['p = ' num2str(pValue)]);

    % spike theta phase trough for each cell
%     figure
%     minDire = spikeThetaPhaseStruct.minPhaseFilArr/360*2*pi;
%     ind = minDire < 0;
%     minDire(ind) = minDire(ind) + 2*pi;
%     % histogram
%     minDireTmp = [minDire,minDire+2*pi,minDire+4*pi];
%     thetaPeakVsDepth.stepPhase = 0:0.01*pi:6*pi;
%     thetaPeakVsDepth.histMinDire = hist(minDireTmp,thetaPeakVsDepth.stepPhase);
%     % smooth histogram 
%     thetaPeakVsDepth.filMinDire = ...
%             conv(thetaPeakVsDepth.histMinDire,paramC.gaussFilt);
%     thetaPeakVsDepth.filMinDire = ...
%             thetaPeakVsDepth.filMinDire(floor(lenGaussKernel/2)+1:...
%                         (end-lenGaussKernel+floor(lenGaussKernel/2)+1));
%     indPhase = find(thetaPeakVsDepth.stepPhase >= 2*pi ...
%                     & thetaPeakVsDepth.stepPhase <4*pi);
%     filMinDire = thetaPeakVsDepth.filMinDire(indPhase); 
%     % correct theta mean direction of individual cells based on the trough 
%     % of the meandire histogram
%     [~,indMaxPhase] = max(filMinDire);
%     thetaPeakVsDepth.peakMinDire = ...
%         thetaPeakVsDepth.stepPhase(indPhase(indMaxPhase)) - 2*pi;
%     [~,indMinPhase] = min(filMinDire);
%     thetaPeakVsDepth.troughMinDire = ...
%         thetaPeakVsDepth.stepPhase(indPhase(indMinPhase)) - 2*pi;
%     if(thetaPeakVsDepth.troughMinDire > pi)
%         indShift = minDire > thetaPeakVsDepth.troughMinDire;
%         minDire(indShift) = minDire(indShift) - 2*pi;
%     else
%         indShift = minDire < thetaPeakVsDepth.troughMinDire;
%         minDire(indShift) = minDire(indShift) + 2*pi;
%     end
%     thetaPeakVsDepth.minDireAdj = minDire;
%     thetaPeakVsDepth.minDireFit = fitlm(depthNeu.relDepthNeuHDef(selPyramidal),...
%         minDire(selPyramidal));
%     subplot(2,2,1)
%     plotAdded(thetaPeakVsDepth.minDireFit);
%     T = thetaPeakVsDepth.minDireFit.anova;
%     pValue = T.Variables;
%     pValue = pValue(1,end);
%     xlabel('depth')
%     ylabel('thetaMinDire (based on minPhaseFilArr)')
%     title(['p = ' num2str(pValue)]);
end
