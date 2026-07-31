function fieldInfo = getFieldInfoIndNeuron(neuronNo,fieldStruct)
% this function is to check whether a neuron has at least a field, and if
% so, extract the field information: start time, stop time, and peak time
% neuronNo:         number of the neuron
% fieldStruct:      field structure (referring to function Fields)
%
% return:
% fieldInfo:        the field information: start time, stop time, and peak
%                   time
    
    fieldInfo = [];
    if(isempty(neuronNo))
        disp('neuronNo should contain only one neuron.');
        return;
    end
    if(isempty(fieldStruct))
       return; 
    end
    indTmp = find(fieldStruct.indNeuron == neuronNo(1));
    numField = length(indTmp);
    if(numField ~= 0)
        for k = 1:numField % label the start, peak and end point of each field
            fieldInfo(k,:) = [fieldStruct.indStartField(indTmp(k)),...
                                 fieldStruct.indStartField(indTmp(k))+fieldStruct.FW(indTmp(k)),...
                                 fieldStruct.indPeakField(indTmp(k))];                    
        end
    end
