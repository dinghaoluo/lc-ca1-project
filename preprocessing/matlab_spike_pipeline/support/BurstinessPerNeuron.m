function burstNeuronStruct = BurstinessPerNeuronVR(burstIsiNeuron,spikeThetaPhaseStruct,param)
% test the burstiness of individual neurons
%
% Inputs:
% burstIsiNeuron:            struct with burst isi of individual neurons
%                               (referring to function BurstPerNeuron)
% spikeThetaPhaseStruct:        struct with spike thetaphase (referring to function SpikeThetaPhase)
% param:                   struct containing parameters for the test
%                               'minFractBurst'     % min percentage of burst spikes
%                               'pOmnibusPhPref'    % when pOmnibus < pOmnibusPhPref, the neuron is considered as with phase preference
%                               'pOmnibusNoPhPref'   % when pOmnibus > pOmnibusNoPhPref, the neuron is considered as without phase preference
%                               'minNumSp'           % min number of spikes
% 
% Outputs:
% burstNeuronStruct:         struct with the test results

    if(isempty(burstIsiNeuron) || isempty(spikeThetaPhaseStruct))
        burstNeuronStruct = [];
        return;
    end
    
    numNeurons = length(burstIsiNeuron.fractBurstMean);
    burstNeuronStruct = struct('neuronBurst', zeros(1,numNeurons),... % 1: bursty neurons, 0: others
                              'neuronNonBurst', zeros(1,numNeurons),... % 1: non-bursty neurons, 0: others
                              ...
                              'neuronBurstPhPref',zeros(1,numNeurons),... % 1: bursty neurons and with phase preference, 0: otherwise
                              'neuronBurstNoPhPref', zeros(1,numNeurons), ... % 1: bursty neurons and without phase preference, 0: otherwise
                              'neuronNonBurstPhPref', zeros(1,numNeurons), ... % 1: non-bursty neurons and with phase preference, 0: otherwise
                              'neuronNonBurstNoPhPref', zeros(1,numNeurons),... % 1: non-bursty neurons and without phase preference, 0: otherwise
                              ...
                              'neuronBurstPhPBurstSpPhP',zeros(1,numNeurons),... % 1: bursty neurons with phase preference and the burst spikes are also with phase preference, 0: otherwise
                              'neuronBurstPhPBurstSpNoPhP',zeros(1,numNeurons),... % 1: bursty neurons with phase preference but the burst spikes are without phase preference, 0: otherwise
                              'neuronBurstNoPhPBurstSpPhP',zeros(1,numNeurons),... % 1: bursty neurons without phase preference but the burst spikes are with phase preference, 0: otherwise
                              'neuronBurstNoPhPBurstSpNoPhP',zeros(1,numNeurons),... % 1: bursty neurons without phase preference and the burst spikes are also without phase preference, 0: otherwise
                              ...
                              'neuronBurstPhPNonBurstSpPhP',zeros(1,numNeurons),... % 1: bursty neurons with phase preference and the non-burst spikes are also with phase preference, 0: otherwise
                              'neuronBurstPhPNonBurstSpNoPhP',zeros(1,numNeurons),... % 1: bursty neurons with phase preference but the non-burst spikes are without phase preference, 0: otherwise
                              'neuronBurstNoPhPNonBurstSpPhP',zeros(1,numNeurons),... % 1: bursty neurons without phase preference but the non-burst spikes are with phase preference, 0: otherwise
                              'neuronBurstNoPhPNonBurstSpNoPhP',zeros(1,numNeurons)); % 1: bursty neurons without phase preference and the non-burst spikes are also without phase preference, 0: otherwise
                                                            
    % bursty neurons                         
    indNeuronsBurst = find(burstIsiNeuron.fractBurstMean > param.minFractBurst); 
        % define burstiness as fraction of burst spikes > param.minFractBurst 
    indNeuronsBurstTmp = []; 
        % find the neurons which have at least param.minNumSp spikes
    for i = 1:length(indNeuronsBurst)
        if(length(spikeThetaPhaseStruct.spTimePerNeuron{indNeuronsBurst(i)})...
                > param.minNumSp)
            indNeuronsBurstTmp = [indNeuronsBurstTmp indNeuronsBurst(i)];
        end
    end  
    burstNeuronStruct.neuronBurst(indNeuronsBurstTmp) = 1; 
    
    % bursty neurons with phase preference
    indNeuronsBurstPhPref = indNeuronsBurstTmp(...
        spikeThetaPhaseStruct.pOmnibus(indNeuronsBurstTmp) < param.pOmnibusPhPref);
    burstNeuronStruct.neuronBurstPhPref(indNeuronsBurstPhPref) = 1; 
    % bursty neurons without phase preference
    indNeuronsBurstNoPhPref = indNeuronsBurstTmp(...
        spikeThetaPhaseStruct.pOmnibus(indNeuronsBurstTmp) > param.pOmnibusNoPhPref);
    burstNeuronStruct.neuronBurstNoPhPref(indNeuronsBurstNoPhPref) = 1; 
    
    % bursty neurons with phase preference, and the burst spikes are also with
    % phase preference
    indNeuronsBurstPhPBurstSpPhP = indNeuronsBurstPhPref(...
        burstIsiNeuron.pOmnibus(indNeuronsBurstPhPref) < param.pOmnibusPhPref);
    burstNeuronStruct.neuronBurstPhPBurstSpPhP(indNeuronsBurstPhPBurstSpPhP) = 1;
    % bursty neurons with phase preference, while the burst spikes are
    % without phase preference
    indNeuronsBurstPhPBurstSpNoPhP = indNeuronsBurstPhPref(...
        burstIsiNeuron.pOmnibus(indNeuronsBurstPhPref) > param.pOmnibusNoPhPref);
    burstNeuronStruct.neuronBurstPhPBurstSpNoPhP(indNeuronsBurstPhPBurstSpNoPhP) = 1;
    % bursty neurons without phase preference, while the burst spikes are with
    % phase preference
    indNeuronsBurstNoPhPBurstSpPhP = indNeuronsBurstNoPhPref(...
        burstIsiNeuron.pOmnibus(indNeuronsBurstNoPhPref) < param.pOmnibusPhPref);
    burstNeuronStruct.neuronBurstNoPhPBurstSpPhP(indNeuronsBurstNoPhPBurstSpPhP) = 1;
    % bursty neurons without phase preference, and the burst spikes are
    % also without phase preference
    indNeuronsBurstNoPhPBurstSpNoPhP = indNeuronsBurstNoPhPref(...
        burstIsiNeuron.pOmnibus(indNeuronsBurstNoPhPref) > param.pOmnibusNoPhPref);
    burstNeuronStruct.neuronBurstNoPhPBurstSpNoPhP(indNeuronsBurstNoPhPBurstSpNoPhP) = 1;
    
    % bursty neurons with phase preference, and the non-burst spikes are also with
    % phase preference
    indNeuronsBurstPhPNonBurstSpPhP = indNeuronsBurstPhPref(...
        burstIsiNeuron.pOmnibusNonBurst(indNeuronsBurstPhPref) < param.pOmnibusPhPref);
    burstNeuronStruct.neuronBurstPhPNonBurstSpPhP(indNeuronsBurstPhPNonBurstSpPhP) = 1;
    % bursty neurons with phase preference, while the non-burst spikes are
    % without phase preference
    indNeuronsBurstPhPNonBurstSpNoPhP = indNeuronsBurstPhPref(...
        burstIsiNeuron.pOmnibusNonBurst(indNeuronsBurstPhPref) > param.pOmnibusNoPhPref);
    burstNeuronStruct.neuronBurstPhPNonBurstSpNoPhP(indNeuronsBurstPhPNonBurstSpNoPhP) = 1;
    % bursty neurons without phase preference, while the non-burst spikes are with
    % phase preference
    indNeuronsBurstNoPhPNonBurstSpPhP = indNeuronsBurstNoPhPref(...
        burstIsiNeuron.pOmnibusNonBurst(indNeuronsBurstNoPhPref) < param.pOmnibusPhPref);
    burstNeuronStruct.neuronBurstNoPhPNonBurstSpPhP(indNeuronsBurstNoPhPNonBurstSpPhP) = 1;
    % bursty neurons without phase preference, and the non-burst spikes are
    % also without phase preference
    indNeuronsBurstNoPhPNonBurstSpNoPhP = indNeuronsBurstNoPhPref(...
        burstIsiNeuron.pOmnibusNonBurst(indNeuronsBurstNoPhPref) > param.pOmnibusNoPhPref);
    burstNeuronStruct.neuronBurstNoPhPNonBurstSpNoPhP(indNeuronsBurstNoPhPNonBurstSpNoPhP) = 1;
    
    % non-bursty neurons
    indNeuronsNonBurst = find(burstIsiNeuron.fractBurstMean <= param.maxFractNonBurst); 
                        % define burstiness as fraction of burst spikes <= param.minFractBurst 
    indNeuronsNonBurstTmp = []; % find the neurons which have at least param.minNumSp spikes
    for i = 1:length(indNeuronsNonBurst)
        if(length(spikeThetaPhaseStruct.spTimePerNeuron{indNeuronsNonBurst(i)})...
                > param.minNumSp)
            indNeuronsNonBurstTmp = [indNeuronsNonBurstTmp indNeuronsNonBurst(i)];
        end
    end  
    burstNeuronStruct.neuronNonBurst(indNeuronsNonBurstTmp) = 1; 
    % non-bursty neurons with phase preference
    indNeuronsNonBurstPhPref = indNeuronsNonBurstTmp(...
        spikeThetaPhaseStruct.pOmnibus(indNeuronsNonBurstTmp) < param.pOmnibusPhPref);         
    burstNeuronStruct.neuronNonBurstPhPref(indNeuronsNonBurstPhPref) = 1; 
    % non-bursty neurons without phase preference
    indNeuronsNonBurstNoPhPref = indNeuronsNonBurstTmp(...
        spikeThetaPhaseStruct.pOmnibus(indNeuronsNonBurstTmp) > param.pOmnibusNoPhPref);         
    burstNeuronStruct.neuronNonBurstNoPhPref(indNeuronsNonBurstNoPhPref) = 1; 
