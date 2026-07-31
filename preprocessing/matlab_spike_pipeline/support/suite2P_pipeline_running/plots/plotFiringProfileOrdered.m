function plotFiringProfileOrdered(path,filename)
% plot the firing profile ordered based on their peak
% e.g.: plotFiringProfileOrdered('Z:\Kori\mice-expdata\ANMC114\A114-20220307\A114-20220307-02\','A114-20220307-02_DataStructure_mazeSection1_TrialType1')

    %% plot the firing profile over distance
    fullpath = [path filename '_PeakFR20mm_Run1.mat'];
    load(fullpath,'pFRStructNormT');
    fullpath = [path filename '_convSpikesDist20mm_Run1.mat'];
    load(fullpath,'paramC');
    [~,I] = sort(pFRStructNormT.peakFRInd');
    SORTEDin = pFRStructNormT.avgFRProfileNorm(I,:);
    figure;
    imagesc(paramC.spaceSteps{1},1:length(I),SORTEDin);
    xlabel('Dist (mm)')
    ylabel('Neuron no.')
    savefig([path filename 'FiringProfileDistOrdered.fig'])
    print('-painters','-dpdf',[path filename 'FiringProfileDistOrdered'],'-r600');
    
    %% plot the firing profile aligned to run onset
    fullpath = [path filename '_PeakFRAlignedRun_msess1_Run1.mat'];
    load(fullpath,'pFRNonStimGoodStruct');
    fullpath = [path filename '_convSpikesAlignedRun_msess1_Run1.mat'];
    load(fullpath,'paramC');
    [~,I] = sort(pFRNonStimGoodStruct.peakFRInd');
    SORTEDin = pFRNonStimGoodStruct.avgFRProfileNorm(I,:);
    figure;
    imagesc(0:paramC.timeBin:size(pFRNonStimGoodStruct.avgFRProfileNorm,2)*paramC.timeBin,1:length(I),SORTEDin);
    set(gca,'XLim',[0 6])
    xlabel('Time (s)')
    ylabel('Neuron no.')
    savefig([path filename 'FiringProfileTimeOrdered.fig'])
    print('-painters','-dpdf',[path filename 'FiringProfileTimeOrdered'],'-r600');
end
