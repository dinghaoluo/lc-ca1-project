function ProcessingAligned_Corr(path,fileName,onlyRun,mazeSess)

    GlobalConst;
    
    spaceBin = 20; % mm
    corrIntervalT = 20; %20; % sec
    corrIntervalTMin = -10;
    corrIntervalD = 1800; % mm
    tc = 500; % ms
    cost = 2/(tc/1000*sampleFq);
    methodTheta = 1;
    
    %% single neuron corr between trials
    disp('Single neuron correlation T')
    neuronSpikeCorrT(path,fileName,onlyRun,mazeSess,corrIntervalT,corrIntervalTMin);
    meanNeuronSpikeCorrT(path,fileName,onlyRun,mazeSess,corrIntervalT);
    meanNeuronSpikeCorrTLastToCurTr(path,fileName,onlyRun,mazeSess,corrIntervalT,corrIntervalTMin);
    
%     disp('Single neuron correlation distance')
%     neuronSpikeCorrDist(path,fileName,onlyRun,mazeSess,corrIntervalD);
%     meanNeuronSpikeCorrDist(path,fileName,onlyRun,mazeSess,corrIntervalD);
    
    %% single neuron trial similarity
    disp('Single neuron cosine similarity')
    spikeTrainSimilarityT(path,fileName,onlyRun,mazeSess,corrIntervalT,corrIntervalTMin);    
    meanSpikeTrainSimilarityT(path,fileName,onlyRun,mazeSess,corrIntervalT);
    meanSpikeTrainSimilarityTLastToCurTr(path,fileName,onlyRun,mazeSess,corrIntervalT,corrIntervalTMin);
    
    %% single neuron Victor & Purpura spike time distance
    disp('Single neuron similarity Victor & Purpura (spike time and interval)')
    spikeTrainSimilarityVP(path,fileName,onlyRun,mazeSess,cost,corrIntervalT,corrIntervalTMin);
    meanSpikeTrainSimilarityVP(path,fileName,onlyRun,mazeSess,cost,corrIntervalT);
    meanSpikeTrainSimilarityVPLastToCurTr(path,fileName,onlyRun,mazeSess,cost,corrIntervalT,corrIntervalTMin);
    %% single neuron Victor & Purpura spike time interval
    meanSpikeTrainSimilarityVPI(path,fileName,onlyRun,mazeSess,cost,corrIntervalT);
%     meanSpikeTrainSimilarityVPILastToCurTr(path,fileName,onlyRun,mazeSess,cost,corrIntervalT,corrIntervalTMin);
    
    %% single neuron Victor & Purpura spike time interval
    disp('Single neuron similarity Van Rossum')
    spikeTrainSimilarityVanRossum(path,fileName,onlyRun,mazeSess,tc,corrIntervalT,corrIntervalTMin);
    meanSpikeTrainSimilarityVanRossum(path,fileName,onlyRun,mazeSess,tc,corrIntervalT);
    meanSpikeTrainSimilarityVanRossumLastToCurTr(path,fileName,onlyRun,mazeSess,tc,corrIntervalT,corrIntervalTMin);
    
    % population corrT between trials
    disp('Population vector correlation')
    popSpikeCorrT(path,fileName,onlyRun,mazeSess,corrIntervalT,corrIntervalTMin);
    meanPopSpikeCorrT(path,fileName,onlyRun,mazeSess,corrIntervalT);
    meanPopSpikeCorrTLastToCurTr(path,fileName,onlyRun,mazeSess,corrIntervalT,corrIntervalTMin);
    
    %% population corrDist between trials
    popSimilarityT(path,fileName,onlyRun,mazeSess,corrIntervalT,corrIntervalTMin);
    meanPopSpikeSimT(path,fileName,onlyRun,mazeSess,corrIntervalT);
    meanPopSpikeSimTLastToCurTr(path,fileName,onlyRun,mazeSess,corrIntervalT,corrIntervalTMin);
    
    collectQualityNeuronDepths(path, fileName);
    compDepthVsSimilarity(path,fileName,onlyRun,mazeSess,corrIntervalT);
end
