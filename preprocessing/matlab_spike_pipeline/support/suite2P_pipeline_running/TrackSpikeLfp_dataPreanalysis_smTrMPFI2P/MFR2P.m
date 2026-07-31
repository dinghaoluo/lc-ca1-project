function [mFRStruct] = MFR2P(spikes,indLapList,numNeurons,lenTrialsT,indRunInLap,startTrNo)
% calculate the mean firing rate 
% Inputs:
% spikes:       spike structure which contains a two dimensional cell
%               structure, with spikes{i,j} representing the spike train 
%               from neuron i in trial j 
% indLapList:   the array containing the index numbers of the valid trials 
%               (Caution: since theta strcuture records the theta traces 
%               from all the valid trials in the task, indLapList here 
%               includes the indices of all to be selected cells in the 
%               theta structure, rather than referring back to the 
%               original trial no.)
% lenTrialsT:   Trial length in sec; changed by Yingxue 09/21/2018
% numNeurons:   the total number of neurons
% timeStep:     the sampling time step
%
% Output: mFRStruct: a structure with following fields   
% indLapList: 
% mFRArr:       mean firing rate array, mFRArr(i,j) is the mean firing
%               rate of neuron i in trial j
% mFR:          mean firing rate for each neuron over all valid trials
% stdFR:        std of the mean firing rate for each neuron over all valid trials

if(nargin == 4)
    indRunInLap = [];
end

numTrials = length(indLapList);
if(numTrials == 0)
    mFRStruct = [];
    return;
end

indLap = find(indLapList == startTrNo);
if(isempty(indLap))
    indLap = 0;
end

mFRStruct = struct('indLapList',indLapList,...
                   'mFRArr',zeros(numNeurons,numTrials),...
                   'mFR',zeros(1,numNeurons),...
                   'stdMFR',zeros(1,numNeurons));
               
%   sum each neuron, and get average
for j = 1:numNeurons 
    for i = 1:numTrials
%         disp(['Neu ' num2str(j) ' Tr ' num2str(i)])
        spPerLap = spikes{indLapList(i)}(:, j);
        if(isempty(indRunInLap))
            totRunSample = 0;
        else
            totRunSample = sum(indRunInLap{indLapList(i)});
        end
        if(totRunSample > 0)
            lenTrial = totRunSample/length(indRunInLap{indLapList(i)})...
                        *lenTrialsT(indLapList(i));
        else
            lenTrial = lenTrialsT(indLapList(i));
        end
            
        numSpikes = sum(spPerLap);
        if(numSpikes ~= 0)
            mFRStruct.mFRArr(j,i) = numSpikes/lenTrial;
        end
    end
end

mFRStruct.mFR = mean(mFRStruct.mFRArr(:,indLap+1:end),2)';
mFRStruct.stdMFR = std(mFRStruct.mFRArr(:,indLap+1:end),0,2)';