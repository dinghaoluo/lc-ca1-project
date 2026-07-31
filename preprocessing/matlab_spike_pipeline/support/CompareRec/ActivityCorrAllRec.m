function ActivityCorrAllRec(onlyRun, recSel)

    meanLickAlignedOverDistAllRec;
    pMeanLickAlignedOverDistAllRec;
    
    meanRunSpeedOverDistAllRec(onlyRun);
    
    %% spike correlation of each neuron
    spikeTrainSimilarityAllRec(onlyRun); % AL only
    spikeTrainSimilarityAllRecNoCue(onlyRun); % no cue only
    spikeTrainSimilarityAllRecPassive(onlyRun); % PL only
    spikeTrainSimilarityAllRecALPL(onlyRun); % combine AL and PL
    spikeTrainSimilarityAllRecCompExp(onlyRun); % compare different conditions
    
    %% spike similarity of each neuron
%     spikeTrainSimilarityTAllRec(onlyRun); % AL only
    spikeTrainSimilarityTAllRecCompExp(onlyRun); % compare different conditions
    
    %% spike VP similarity of each neuron
    spikeTrainSimilarityVPAllRec(onlyRun); % AL only
    
    %% spike distance correlation of each neuron
    spikeTrainCorrDistAllRec(onlyRun, recSel); % AL only or combine AL and PL
    spikeTrainCorrDistAllRecPassive(onlyRun); % PL only
    spikeTrainCorrDistAllRecNoCue(onlyRun); % no cue only
    spikeTrainCorrDistAllRecCompExp(onlyRun); % compare different conditions
    
    %% population similarity 
    popSimilarityAllRec(onlyRun,recSel); % AL only or combine AL and PL
    
end
