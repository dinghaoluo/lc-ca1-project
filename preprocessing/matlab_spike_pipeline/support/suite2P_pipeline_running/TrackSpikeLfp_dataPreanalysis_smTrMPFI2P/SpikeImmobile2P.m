function SpikeImmobile2P(path,fileName)
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
    
    GlobalConst2P;

    trialsIm = cell(1,beh.numTrials);
    for i = 1:beh.numTrials
        if(sum(beh.indGoodLap == i)>0)
            ind = find(trials{i}.speed <= minSpeed);
            trialsIm{i}.ind = ind;
            trialsIm{i}.spikes = trials{i}.spikes(ind,:);
            trialsIm{i}.spikesSM = trials{i}.spikesSM(ind,:);
            trialsIm{i}.F = trials{i}.F(ind,:);
            trialsIm{i}.Fneu = trials{i}.Fneu(ind,:);
            trialsIm{i}.dFF = trials{i}.dFF(ind,:);
            trialsIm{i}.dFFSM = trials{i}.dFFSM(ind,:);
            trialsIm{i}.dFFGF = trials{i}.dFFGF(ind,:); % added on 9/6/2023
            trialsIm{i}.xMM = trials{i}.xMM(ind);
            trialsIm{i}.speed = trials{i}.speed(ind);
        end
    end
        
    fullPath = [path fileNameIm];
    save(fullPath, 'trialsIm','-v7.3');
end
