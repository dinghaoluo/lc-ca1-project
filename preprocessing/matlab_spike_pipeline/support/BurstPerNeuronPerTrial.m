function [burstPerNeuPerTr] = BurstPerNeuronPerTrialVR(spikes,theta,dist,lapList,param)
% This function is to extract the burst spikes within each trial for each
% neuron
%
% Inputs:       
% spikes:           spikes per trial
% theta:            spike phase per trial
% dist:             distance on the track for each spike 
% lapList:          list of all valid trials
% param:            a struct including numNeurons, numSamples, and burstIsi
%
% Output:
% burstPerNeuPerTr:    a struct including all the informations about
%                           bursts within each trial (i) for each neuron (j)

    numTrials = length(lapList);
    if(numTrials == 0)
        burstPerNeuPerTr = [];
        return;
    end
    
    burstPerNeuPerTr = struct('lapList',lapList,...
                               'numSpPerTrial', zeros(numTrials,param.numNeurons),... % number of spikes per trial
                               'startTimeBurst', [],... % the start time of the each burst
                               'timeBurst', [],... % time of all the spikes within each burst
                               'numSpPerBurst', [],... % number of spikes per burst                                   
                               'startPhaseBurst', [],... % start phase of each burst
                               'phaseBurst', [],... % phase of all the spikes within each burst
                               'isiBurst',[],... % ISI of burst spikes
                               'fractBurst',zeros(numTrials,param.numNeurons),... % fraction of spikes which belongs to a burst
                               'burstyTrials', zeros(numTrials,param.numNeurons),... % bursty trials for each neuron
                               ...
                               'StartDistBurst', [],... % the x location of where the first spike of a burst occurs for the maze run
                               'DistBurst', [],... % the x location of where the burst occurs for the maze run
                               ...
                               'timeNonBurst',[],... % time of all the non-burst spikes
                               'phaseNonBurst',[],... % phase of all the non-burst spikes
                               'distNonBurst',[]); % distance of all the non-burst spikes

    for tr = 1:numTrials % calculate isi and extract all the bursting spike pairs
        % extract the bursts for each neuron and each trial
        i = lapList(tr);
        for j = 1:param.numNeurons
%             disp(['Trial ' num2str(i) ' Neuron ' num2str(j)]);
            burstPerNeuPerTr.numSpPerTrial(tr,j) = length(spikes{i}{j});
            isiSpikes = diff(spikes{i}{j}); 
                % calculate the isi between spikes
            lenArr = length(isiSpikes);
            spikeLabelArr = zeros(1,lenArr);
            spikeLabelArr(isiSpikes <= param.burstIsi) = 1; 
                % the array used to label each pair of burst spikes
            
            % go through the burst label array and record the starting
            % point of each burst and the number of spikes within each
            % burst
            kStart = 1;
            kStop = kStart;
            tmpStartTimeBurst = [];
            tmpStartPhaseBurst = [];
            tmpNumSpPerBurst = [];
            tmpTimeBurst = [];
            tmpPhaseBurst = [];
            tmpISIBurst = [];
            tmpStartDistBurst = [];
            tmpDistBurst = [];
            while kStart <= lenArr
                if(kStart <= kStop)
                    if(kStop <= lenArr)
                        if(spikeLabelArr(kStop) == 1)
                            kStop = kStop + 1;
                        else
                            if(kStop > kStart)
                                if(isiSpikes(kStart) > param.burstIsi1st)
                                    spikeLabelArr(kStart) = 0;
                                    kStart = kStart + 1;
                                    continue;
                                end
                                tmpTime = spikes{i}{j}(kStart);
                                tmpStartTimeBurst = [tmpStartTimeBurst tmpTime];
                                tmpTimeBurst = [tmpTimeBurst spikes{i}{j}(kStart:kStop)'];
                                tmpNumSpPerBurst = [tmpNumSpPerBurst kStop-kStart+1];
                                tmpISIBurst = [tmpISIBurst isiSpikes(kStart:kStop-1)'];
                                tmpStartPhaseBurst = [tmpStartPhaseBurst theta{i}{j}(kStart)];
                                tmpPhaseBurst = [tmpPhaseBurst theta{i}{j}(kStart:kStop)'];
                                
                                tmpStartDistBurst = [tmpStartDistBurst dist{i}{j}(kStart)];
                                tmpDistBurst = [tmpDistBurst dist{i}{j}(kStart:kStop)'];
                            end
                            kStart = kStop + 1;
                        end
                    else
                        if(isiSpikes(kStart) > param.burstIsi1st)
                            spikeLabelArr(kStart) = 0;
                            kStart = kStart + 1;
                            continue;
                        end
                        tmpTime = spikes{i}{j}(kStart);
                        tmpStartTimeBurst = [tmpStartTimeBurst tmpTime];
                        tmpTimeBurst = [tmpTimeBurst spikes{i}{j}(kStart:lenArr+1)'];
                        tmpNumSpPerBurst = [tmpNumSpPerBurst lenArr-kStart+2];
                        tmpISIBurst = [tmpISIBurst isiSpikes(kStart:lenArr)'];
                        tmpStartPhaseBurst = [tmpStartPhaseBurst theta{i}{j}(kStart)];
                        tmpPhaseBurst = [tmpPhaseBurst theta{i}{j}(kStart:lenArr+1)'];
                        
                        tmpStartDistBurst = [tmpStartDistBurst dist{i}{j}(kStart)];
                        tmpDistBurst = [tmpDistBurst dist{i}{j}(kStart:lenArr+1)'];
                        
                        kStart = lenArr+1;
                    end
                else
                    if(spikeLabelArr(kStart) == 1)
                        kStop = kStart + 1;
                    else
                        kStart = kStart + 1;
                    end
                end
            end
            burstPerNeuPerTr.startTimeBurst{tr,j} = tmpStartTimeBurst; 
            burstPerNeuPerTr.timeBurst{tr,j} = tmpTimeBurst;
            burstPerNeuPerTr.numSpPerBurst{tr,j} = tmpNumSpPerBurst; 
            burstPerNeuPerTr.startPhaseBurst{tr,j} = tmpStartPhaseBurst;
            burstPerNeuPerTr.phaseBurst{tr,j} = tmpPhaseBurst;
            burstPerNeuPerTr.isiBurst{tr,j} = tmpISIBurst;
            if(~isempty(tmpTimeBurst))
                burstPerNeuPerTr.fractBurst(tr,j) = length(tmpTimeBurst)/(lenArr+1);
                burstPerNeuPerTr.burstyTrials(tr,j) = 1;
            end
            burstPerNeuPerTr.StartDistBurst{tr,j} = tmpStartDistBurst;
            burstPerNeuPerTr.DistBurst{tr,j} = tmpDistBurst;
            
            % find all non-burst spikes, changed on 3/28/2021 by Yingxue,
            % removed the first spikes from the bursts beginning from the second
            % burst
            [~,indNonBurst] = setdiff(spikes{i}{j},tmpTimeBurst);
            burstPerNeuPerTr.timeNonBurst{tr,j} = [];
            burstPerNeuPerTr.phaseNonBurst{tr,j} = [];
            if(~isempty(indNonBurst))
                burstPerNeuPerTr.timeNonBurst{tr,j} = ...
                    spikes{i}{j}(indNonBurst)'; 
                
                burstPerNeuPerTr.phaseNonBurst{tr,j} = ...
                    theta{i}{j}(indNonBurst)';
                
                burstPerNeuPerTr.distNonBurst{tr,j} = ...
                    dist{i}{j}(indNonBurst)';
                
            end
            
        end
    end