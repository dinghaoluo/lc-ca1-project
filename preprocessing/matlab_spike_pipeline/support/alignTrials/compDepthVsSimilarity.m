function compDepthVsSimilarity(path,fileName,onlyRun,mazeSess,intervalT)

    fileNameCorr = [path fileName '_meanSpikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
        
    fileNameCorrSim = [fileName '_meanSpikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
        
    fileNameQualNeurons = [path fileName '_qualityNeuronDepths' '.mat'];
    
    load(fileNameCorr, 'meanCorrTRun');
    load(fileNameCorrSim, 'meanSimTRun');
    load(fileNameQualNeurons, 'QualityNeurons');
    
    clusterID = QualityNeurons.clusterID;
    DepthVsSim.clusterID = clusterID;
    DepthVsSim.meanGoodCorrTRun = meanCorrTRun.meanGood(clusterID);
    DepthVsSim.meanGoodNZCorrTRun = meanCorrTRun.meanGoodNZ(clusterID);
    DepthVsSim.meanGoodSimTRun = meanSimTRun.meanGood(clusterID);
    DepthVsSim.meanGoodNZSimTRun = meanSimTRun.meanGoodNZ(clusterID);
    DepthVsSim.relDepthNeuHDef = QualityNeurons.relDepthNeuHDef;
    DepthVsSim.peakTo20ms = QualityNeurons.peakTo20ms;
    
    fileNameQualNeurons = [path fileName '_neuronDepthsVsSimilarity' '.mat'];
    save(fileNameQualNeurons,'DepthVsSim');
    
    figure;
    subplot(2,2,1)
    h = plot(DepthVsSim.meanGoodCorrTRun,DepthVsSim.relDepthNeuHDef,'o');
    ylabel('Rel. depth');
    xlabel('Mean good corrT Run');
    
    subplot(2,2,2)
    h = plot(DepthVsSim.meanGoodSimTRun,DepthVsSim.relDepthNeuHDef,'o');
    ylabel('Rel. depth');
    xlabel('Mean good simT Run');
    
    subplot(2,2,3)
    h = plot(DepthVsSim.peakTo20ms,DepthVsSim.relDepthNeuHDef,'o');
    ylabel('Rel. depth');
    xlabel('Autocorr. peak to 20 ms');
    
    subplot(2,2,4)
    h = plot(DepthVsSim.peakTo20ms,DepthVsSim.meanGoodCorrTRun,'o');
    ylabel('Mean good corrT Run');
    xlabel('Autocorr. peak to 20 ms');
