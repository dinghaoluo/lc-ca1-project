function neuronClass = getNeuronClassIndNeuron(neuronNo,neuronClassStruct)
% get the class which the neuron belongs to
% neuronNo:             number of the neuron
% neuronClassStruct:    neuron class structure (referring to function NeuronClass)
%
% neuronClass:          return a string labelling the neuron class

neuronClass = [];
if(isempty(neuronNo))
    disp('neuronNo should contain only one neuron.');
    return;
end

if(~isempty(find(neuronClassStruct.neuronConstFiring == neuronNo(1), 1))) 
    neuronClass = 'Const firing';
elseif(~isempty(find(neuronClassStruct.neuronInitPeakConstF == neuronNo(1), 1)))
    neuronClass = 'Initial peak';
elseif(~isempty(find(neuronClassStruct.neuronPotentialField == neuronNo(1), 1)))
    neuronClass = 'Potential field';
end