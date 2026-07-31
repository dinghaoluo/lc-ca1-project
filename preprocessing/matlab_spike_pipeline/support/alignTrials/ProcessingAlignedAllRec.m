function ProcessingAlignedAllRec(onlyRun)
% no cue passive --- condition 1
% 60 cm cue passive --- condition 2
% 0.5s start cue passive --- condition 3
% 0.5s start cue active --- condition 4

    % plot spike train corrD for active licking task (compare Run Cue and Reward alignment)
    spikeTrainCorrDistAllRec(onlyRun);
    
    % plot spike train corrD for passive no cue task (compare Run Cue and Reward alignment)
    spikeTrainCorrDistAllRecNoCue(onlyRun);
    
    % plot spike train corrD for passive licking task (compare Run Cue and Reward alignment)
    spikeTrainCorrDistAllRecPassive(onlyRun);
    
    % plot spike train corrD for passive no cue task (compare Run Cue and Reward alignment)
    spikeTrainCorrDistAllRecPassive(onlyRun);
    
    % plot spike train corrD (compare PassiveNoCue Passive and Active tasks)
    spikeTrainCorrDistAllRecCompExp(onlyRun);
     
    % plot spike train corrT for active licking task (compare Run Cue and Reward alignment)
    spikeTrainSimilarityAllRec(onlyRun);
    
    % plot spike train similarityT for active licking task(compare Run Cue and Reward alignment)
    spikeTrainSimilarityTAllRec(onlyRun);
    
    % plot spike train VP distance for active licking task (compare Run Cue and Reward alignment)
    spikeTrainSimilarityVPAllRec(onlyRun);
    
    % plot spike train corrT for passive no cue task (compare Run Cue and Reward alignment)
    spikeTrainSimilarityAllRecNoCue(onlyRun);
    
    % plot spike train corrT for passive licking task (compare Run Cue and Reward alignment)
    spikeTrainSimilarityAllRecPassive(onlyRun);
    
    % plot spike train corrT for active licking and passive licking task (compare Run Cue and Reward alignment)
    spikeTrainSimilarityAllRecALPL(onlyRun);
    
    % plot spike train corrT (compare PassiveNoCue Passive and Active tasks)
    spikeTrainSimilarityAllRecCompExp(onlyRun);
    
    % plot spike train similarityT (compare PassiveNoCue Passive and Active tasks)
    spikeTrainSimilarityTAllRecCompExp(onlyRun);
        
    % plot population corrT for active licking task (compare Run Cue and Reward alignment)
    popSimilarityAllRec(onlyRun)
    
    % plot population corrT for passive no cue task (compare Run Cue and Reward alignment)
    popSimilarityAllRecNoCue(onlyRun)
    
    % compare population corrT between passive no cue and active licking (aligned to run onset) 
    popSimilarityAllRecCompExp(onlyRun);
end