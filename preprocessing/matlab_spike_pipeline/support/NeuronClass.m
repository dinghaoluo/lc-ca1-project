function neuronClassP2MFRStruct = NeuronClass(pFRStruct,indNeurons,paramClass)
% classify the neurons into different classes depending on their peak to
% mean firing rate ratio
%
% Inputs:
% pFRStruct:            peak firing rate structure (referring to function PeakFR)
% indNeurons:           indices of neurons
% paramClass:              parameter structure of the field, which should
%                       include:
%           highInstFR:                 lower the peak to mean instantaneous firing
%                                       rate ratio for neurons whose firing rate is larger
%                                       than highInstFR
%           pToMFRRatioMinHighFR:       the peak to mean instantaneous
%                                       firing rate ratio separating neurons with a
%                                       neuron with constant firing and with a peak
%                                       for neurons with high inst FR
%           pToMFRRatioMin1HighFR:      the peak to mean instantaneous
%                                       firing rate ratio separating neurons with a
%                                       peak and neurons that could have a field 
%                                       for neurons with low inst FR
%           pToMFRRatioMin2HighFR:      the peak to mean instantaneous
%                                       firing rate ratio above which neuron is 
%                                       highly likely to have a field for neurons 
%                                       with high inst FR
%           lowerBoundFRFieldNeuronHighFR: the criteria to distinguish a neuron with 
%                                          a field from a neuron with constant firing 
%                                          rate and with a initial bump.
%                                          that is if the neuron inst firing rate
%                                          returns back to mean inst FR *
%                                          lowerBoundFRFieldNeuron for
%                                          neurons with high inst FR
%           pToMFRRatioMinLowFR:        the peak to mean instantaneous
%                                       firing rate ratio separating neurons with a 
%                                       neuron with constant firing and with a peak
%                                       for neurons with low inst FR
%           pToMFRRatioMin1LowFR:       the peak to mean instantaneous firing rate 
%                                       ratio separating neurons with a peak and 
%                                       neurons that could have a field for neurons 
%                                       with low inst FR
%           pToMFRRatioMin2LowFR:       the peak to mean instantaneous
%                                       firing rate ratio above which neuron is 
%                                       highly likely to have a
%                                       field for neurons with low inst FR
%           lowerBoundFRFieldNeuronLowFR:   the criteria to distinguish a neuron 
%                                           with a field from a neuron with constant 
%                                           firing rate and with a initial bump.
%                                           that is if the neuron inst firing rate
%                                           returns back to mean inst FR *
%                                           lowerBoundFRFieldNeuron for neurons with
%                                           low inst FR
%
% Outputs: neuronClassP2MFRStruct which includes the following fields
% neuronConstFiring:        indices of neurons with constant firing rate
% neuronInitPeakConstF:     indices of neurons with constant firing rate
%                           plus a initial peak
% neuronPotentialField:     indices of neurons with potential fields

    if(isempty(pFRStruct))
        neuronClassP2MFRStruct = [];
        return;
    end
    numNeurons = length(indNeurons);
    
    %%%%%%%%% define variables to record different types of neurons
    %%%%%%%%% depending on their firing pattern
    neuronClassP2MFRStruct = struct(...
                'neuronConstFiring', [],... % neurons with constant firing rate
                'neuronInitPeakConstF', [],... % neurons with constant firing rate plus a initial peak
                'neuronPotentialField', []); % neurons with potential fields 
            
    for i = 1:numNeurons
       
        if(pFRStruct.meanInstFR(indNeurons(i)) >= paramClass.highInstFR) 
            % reassign the threshold to a lower value if the mean inst firing
            % rate is high
            pToMFRRatioMin = paramClass.pToMFRRatioMinHighFR;
            pToMFRRatioMin1 = paramClass.pToMFRRatioMin1HighFR;
            pToMFRRatioMin2 = paramClass.pToMFRRatioMin2HighFR;
            lowerBoundFRFieldNeuron = paramClass.lowerBoundFRFieldNeuronHighFR;
        else
            pToMFRRatioMin = paramClass.pToMFRRatioMinLowFR;
            pToMFRRatioMin1 = paramClass.pToMFRRatioMin1LowFR;
            pToMFRRatioMin2 = paramClass.pToMFRRatioMin2LowFR;
            lowerBoundFRFieldNeuron = paramClass.lowerBoundFRFieldNeuronLowFR;
        end
        
        if(pFRStruct.meanInstFR(indNeurons(i)) > paramClass.minInstFR) 
            % all the neurons which could have a field should have a firing 
            % rate > minInstFR
            pToMFRRatio = pFRStruct.peakFR(indNeurons(i))...
                /pFRStruct.meanInstFR(indNeurons(i));
            if(pToMFRRatio < pToMFRRatioMin) 
                % neurons whose peak/mean instantaneous firing rate is smaller
                % than pToMFRRatioMin is classified as neurons with constant firing
                neuronClassP2MFRStruct.neuronConstFiring = ...
                    [neuronClassP2MFRStruct.neuronConstFiring indNeurons(i)];
                disp(['Neuron ' num2str(indNeurons(i)) ' is constant firing']);
                
            elseif(pToMFRRatio >= pToMFRRatioMin && pToMFRRatio < pToMFRRatioMin1)
                % % neurons whose peak/mean instantaneous firing rate falls
                % between pToMFRRatioMin and pToMFRRatioMin1 are classified
                % as neurons with initial peak and constant firing rate
                neuronClassP2MFRStruct.neuronInitPeakConstF = ...
                    [neuronClassP2MFRStruct.neuronInitPeakConstF indNeurons(i)]; 
                disp(['Neuron ' num2str(indNeurons(i)) ' is initial peak firing']);
                
            elseif(pToMFRRatio >= pToMFRRatioMin1 && pToMFRRatio < pToMFRRatioMin2)
                % neurons whose peak/mean instantaneous firing rate falls
                % between pToMFRRatioMin1 and pToMFRRatioMin2 could be
                % either neurons with initial peak and constant firing rate
                % or neurons with a field
                
                if(~isempty(find(pFRStruct.avgFRProfile(indNeurons(i),:) < ...
                        pFRStruct.meanInstFR(indNeurons(i))*lowerBoundFRFieldNeuron, 1)))
                        % for the neurons whose firing rate return back to
                        % very low level, check for a field
                    neuronClassP2MFRStruct.neuronPotentialField = ...
                        [neuronClassP2MFRStruct.neuronPotentialField indNeurons(i)];
                    disp(['Neuron ' num2str(indNeurons(i)) ' has potential fields']);
                else % otherwise, it is classified as a neuron with constant 
                     % firing and a initial peak
                    neuronClassP2MFRStruct.neuronInitPeakConstF = ...
                        [neuronClassP2MFRStruct.neuronInitPeakConstF indNeurons(i)]; 
                    disp(['Neuron ' num2str(indNeurons(i)) ' is initial peak firing']);
                end
                
            else % neurons whose peak/mean instantaneous firing rate is high 
                 % than pToMFRRatioMin2 should have a field
                neuronClassP2MFRStruct.neuronPotentialField = ...
                    [neuronClassP2MFRStruct.neuronPotentialField indNeurons(i)];
                disp(['Neuron ' num2str(indNeurons(i)) ' has potential fields']);
            end
        end
    end
