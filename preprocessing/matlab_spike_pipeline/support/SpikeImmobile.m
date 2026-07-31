function SpikeImmobile(path,fileName)
% extract information for all the spikes during immobile period (speed < certain threshold)
%
% e.g.: SpikeImmobileVR('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1')

    %%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin > 3
        disp('Too many arguments');        
        return;
    end
    
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(isempty(indexFileName))
        fileNameInfo = [fileName '_Info.mat'];
        fileNameIm = [fileName '_im.mat'];
        fileName = [fileName '.mat'];
    else
        fileNameInfo = [fileName(1:indexFileName(end)-1) '_Info.mat'];
        fileNameIm = [fileName(1:indexFileName(end)-1) '_im.mat'];
    end 
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath,'trials');
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath);
    
    GlobalConst;
    
    trialsIm = cell(1,beh.numTrials);
    for i = 1:beh.numTrials
        if(sum(beh.indGoodLap == i)>0)
            trialsIm{i}.spikes = cell(1,rec.numNeurons);
            trialsIm{i}.spikes20kHz = cell(1,rec.numNeurons);
            trialsIm{i}.spikesThetaHil = cell(1,rec.numNeurons);
            trialsIm{i}.spikesThetaLin = cell(1,rec.numNeurons);
            trialsIm{i}.spikesMM = cell(1,rec.numNeurons);
            trialsIm{i}.spikesSpeed = cell(1,rec.numNeurons);
            
            for j = 1:rec.numNeurons
                ind = trials{i}.spikesSpeed{j} <= minSpeed;
                trialsIm{i}.spikes{j} = trials{i}.spikes{j}(ind);
                trialsIm{i}.spikes20kHz{j} = trials{i}.spikes20kHz{j}(ind);
                trialsIm{i}.spikesThetaHil{j} = trials{i}.spikesThetaHil{j}(ind);
                trialsIm{i}.spikesThetaLin{j} = trials{i}.spikesThetaLin{j}(ind);
                trialsIm{i}.spikesMM{j} = trials{i}.spikesMM{j}(ind);
                trialsIm{i}.spikesSpeed{j} = trials{i}.spikesSpeed{j}(ind);
            end
        end
    end
    
    fullPath = [path fileNameIm];
    save(fullPath, 'trialsIm');
end
