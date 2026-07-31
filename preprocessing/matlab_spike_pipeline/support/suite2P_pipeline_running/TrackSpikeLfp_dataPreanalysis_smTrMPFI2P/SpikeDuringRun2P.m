function SpikeDuringRun2P(path,fileName)
% extract information for all the spikes during running (speed > certain threshold)
%
% e.g.: SpikeDuringRunVR('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1')

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
        fileNameExt = [fileName '_ext.mat'];
        fileName = [fileName '.mat'];
    else
        fileNameInfo = [fileName(1:indexFileName(end)-1) '_Info.mat'];
        fileNameExt = [fileName(1:indexFileName(end)-1) '_ext.mat'];
    end 
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath,'trials');
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr2P(path,fileName);
    end
    load(fullPath);
    
    GlobalConst2P;

    trialsExt = cell(1,beh.numTrials);
    for i = 1:beh.numTrials
        if(sum(beh.indGoodLap == i)>0)
            ind = find(trials{i}.speed > minSpeed);
            trialsExt{i}.ind = ind;
            trialsExt{i}.spikes = trials{i}.spikes(ind,:);
            trialsExt{i}.spikesSM = trials{i}.spikesSM(ind,:);
            trialsExt{i}.F = trials{i}.F(ind,:);
            trialsExt{i}.Fneu = trials{i}.Fneu(ind,:);
            trialsExt{i}.dFF = trials{i}.dFF(ind,:);
            trialsExt{i}.dFFSM = trials{i}.dFFSM(ind,:);
            trialsExt{i}.dFFGF = trials{i}.dFFGF(ind,:); % added on 9/6/2023
            trialsExt{i}.xMM = trials{i}.xMM(ind);
            trialsExt{i}.speed = trials{i}.speed(ind);
        end
    end
    
    fullPath = [path fileNameExt];
    save(fullPath, 'trialsExt','-v7.3');
end
