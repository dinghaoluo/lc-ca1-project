function [pFRStruct] = PeakFR(filteredSpArr,indLapList,numNeurons,numSamples)
% calculate the peak firing rate 
% Inputs:
% filteredSpArr:        the filetered spike profiles containing a one dimensional cell
%                       array, with filteredSpArr{i} the filtered spike train from trial i for all the neurons 
% indLapList:           the array containing the index numbers of the valid trials (Caution: since theta strcuture records the theta traces from all the valid trials in 
%                       the task, indLapList here includes the indices of all to be selected cells in the theta structure, rather than referring back to the 
%                       original trial no.)
% numNeurons:           the total number of neurons
% numSamples:           only consider the samples between 0 and numSamples in the
%                       peak firing rate calculation
%
% Output: pFRStruct: a structure with following fields  
% indLapList: 
% avgFRProfile:         the firing rate profile for each neuron averaged
%                       across all the valid trials
% peakFR:               peak firing rate for each neuron over
%                       all valid trials
% peakFRInd:            index of the peak firing rate for each neuron over
%                       all valid trials
% meanInstFR:           mean instantaneous firing rate for each neuron over all valid trials
% p2MInstRatio:         peak to mean instantaneous firing rate for each
%                       neuron over all valid trials
% nonZeroNumTrials:     real number of trials with non-zero number of
%                       spikes

numTrials = length(indLapList);
if(numTrials == 0)
    pFRStruct = [];
    return;
end

pFRStruct = struct('indLapList',indLapList,...
                   'peakFRPerTrial',zeros(numNeurons,numTrials),... % peak firing rate of each individual trials
                   'peakFRIndPerTrial',-1*ones(numNeurons,numTrials),... % index of the peak of each individual trials
                   'peakFRIndPerTrialAvg',zeros(1,numNeurons),... % averaged peak location
                   'peakFRIndPerTrialStd',zeros(1,numNeurons),... % std of the peak location
                   'nonZeroTrialInd',[],...   % the trial indices with non-zero number of spikes
                   'avgFRProfile',zeros(numNeurons,numSamples),... % the averaged firing profile
                   'peakFR',zeros(1,numNeurons),... % peak firing rate calculated from the averaged firing profile
                   'peakFRInd',zeros(1,numNeurons),... % index of the peak calculated from the averaged firing profile
                   'peakFRStd',zeros(1,numNeurons),... % std of the firing rate at the peak location calculated from the firing profile of non-zero spiking trials
                   'meanInstFR',zeros(1,numNeurons),... % mean instantaneous firing rate of the averaged firing profile
                   'p2MInstRatio',zeros(1,numNeurons),... % peak to mean inst firing rate ratio of the averaged firing profile
                   'nonZeroNumTrials',zeros(1,numNeurons),... % real number of trials
                   'filteredSpArrNonZero',[],... % filtered spike arrays for each neuron with only trials with firings
                   'filteredSpArrAll',[]); % spike arrays for each neuron and every trial (no matter how many spikes there were

for i = 1:numNeurons
    nonZeroTrialInd = [];
    nonZeroTrialIndInd = [];
    pFRStruct.filteredSpArrNonZero{i} = zeros(numTrials,numSamples);
    pFRStruct.filteredSpArrAll{i} = zeros(numTrials,numSamples);
    for j = 1:numTrials
        filteredSpTmp = filteredSpArr{indLapList(j)}(i,:);

        if(~isempty(find(filteredSpTmp ~= 0, 1))) 
            % guarantee that there is at least one spike in the trial
            pFRStruct.nonZeroNumTrials(i) = pFRStruct.nonZeroNumTrials(i) + 1;
            nonZeroTrialInd = [nonZeroTrialInd indLapList(j)];
            nonZeroTrialIndInd = [nonZeroTrialIndInd j];
            
            %%%%% changed on 4/8/2019 by Yingxue to accommodate different
            %%%%% track length
            pFRStruct.avgFRProfile(i,1:length(filteredSpTmp)) = ...
                pFRStruct.avgFRProfile(i,1:length(filteredSpTmp)) ...
                        + filteredSpTmp;
            %%%%%
            pFRStruct.peakFRIndPerTrial(i,j) = ...
                        find(filteredSpTmp == max(filteredSpTmp),1);
            pFRStruct.peakFRPerTrial(i,j) = ...
                        filteredSpTmp(pFRStruct.peakFRIndPerTrial(i,j));
            pFRStruct.filteredSpArrNonZero{i}(pFRStruct.nonZeroNumTrials(i),1:length(filteredSpTmp)) = ...
                        filteredSpTmp;            
        end
        %%%%% changed on 4/8/2019 by Yingxue to accommodate different
        %%%%% track length
        pFRStruct.filteredSpArrAll{i}(j,1:length(filteredSpTmp)) = filteredSpTmp;   
        %%%%%
    end
    
    if(pFRStruct.nonZeroNumTrials(i) ~= 0)
        pFRStruct.avgFRProfile(i,:) = pFRStruct.avgFRProfile(i,:)/numTrials; 
        % average profile over trials with non-zero number of spikes
    else 
        disp(['Neuron ' num2str(i) ' has generated 0 spikes.']);
    end

    pFRStruct.nonZeroTrialInd{i} = nonZeroTrialInd;
    pFRStruct.peakFRIndPerTrialAvg(i) = ...
        mean(pFRStruct.peakFRIndPerTrial(i,nonZeroTrialIndInd));
    pFRStruct.peakFRIndPerTrialStd(i) = ...
        std(pFRStruct.peakFRIndPerTrial(i,nonZeroTrialIndInd));
    [pFRStruct.peakFR(i),indTmp] = max(pFRStruct.avgFRProfile(i,:)); 
        % peak firing rate of each neuron
    pFRStruct.peakFRInd(i) = indTmp(1); % index of the peak firing rate
    peakFRPerTrial = zeros(1,numTrials);
    if(~isempty(pFRStruct.filteredSpArrAll{i}))
        peakFRPerTrial = ...
            pFRStruct.filteredSpArrAll{i}(:,indTmp(1))';
    end
    pFRStruct.peakFRStd(i) = std(peakFRPerTrial(nonZeroTrialIndInd));
    pFRStruct.meanInstFR(i) = mean(pFRStruct.avgFRProfile(i,:)); 
        % mean instantaneous firing rate
    if(pFRStruct.meanInstFR(i) ~= 0)
        pFRStruct.p2MInstRatio(i) = pFRStruct.peakFR(i)/pFRStruct.meanInstFR(i); 
        % peak to mean instantaneous firing rate ratio
    end
end
