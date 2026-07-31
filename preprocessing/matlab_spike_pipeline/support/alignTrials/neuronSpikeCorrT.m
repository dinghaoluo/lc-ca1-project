function neuronSpikeCorrT(path,fileName,onlyRun,mazeSess,intervalT,intervalTMin)
% single neuron level spike correlation across trials
% e.g.: neuronSpikeCorrT('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,20,-10)
    
    if(nargin == 4)
        intervalT = 0;
    end
    
    fullPath = [path fileName '_convSpikesAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _convSpikesAligned_Run file does not exist');
        return;
    end
%     load(fullPath,'filteredSpikeArrayRun','filteredSpikeArrayRew','filteredSpikeArrayCue',...
%             'paramC','timeStep1');
    load(fullPath,'filteredSpikeArrayRun','filteredSpikeArrayRew','filteredSpikeArrayCue',...
            'filteredSpikeArrayRun_LasttoCurTr','filteredSpikeArrayRew_LasttoCurTr','filteredSpikeArrayCue_LasttoCurTr',...
            'paramC','timeStep1');
    
    fileNameCorr = [fileName '_spikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    
    GlobalConst;
    
    if(intervalT > paramC.trialLenT | intervalT == 0)
        intervalT = paramC.trialLenT;
    end
    intervalT1 = floor(intervalT/paramC.timeBin);
    neuronNo = length(filteredSpikeArrayRun);  
    
    indMin = find(timeStep1 >= intervalTMin*sampleFq,1);
    indMax = find(timeStep1 <= intervalT*sampleFq,1,'last');
    indInt = [indMin indMax];
    
    disp('calculating single neuron correlation aligned to run onset');
    [spikeCorrTRun,nonZeroTrRun] = calSpikeCorrT(neuronNo,filteredSpikeArrayRun,intervalT1); 
    disp('calculating single neuron correlation aligned to reward onset');
    [spikeCorrTRew,nonZeroTrRew] = calSpikeCorrT(neuronNo,filteredSpikeArrayRew,intervalT1);
    disp('calculating single neuron correlation aligned to cue onset');
    [spikeCorrTCue,nonZeroTrCue] = calSpikeCorrT(neuronNo,filteredSpikeArrayCue,intervalT1);
    
    disp('calculating single neuron correlation aligned to run onset (all spikes from last reward)');
    [spikeCorrTRun_LasttoCurTr,nonZeroTrRun_LasttoCurTr] = calSpikeCorrT(neuronNo,filteredSpikeArrayRun_LasttoCurTr,indInt); 
    disp('calculating single neuron correlation aligned to reward onset (all spikes from last reward)');
    [spikeCorrTRew_LasttoCurTr,nonZeroTrRew_LasttoCurTr] = calSpikeCorrT(neuronNo,filteredSpikeArrayRew_LasttoCurTr,indInt);
    disp('calculating single neuron correlation aligned to cue onset (all spikes from last reward)');
    [spikeCorrTCue_LasttoCurTr,nonZeroTrCue_LasttoCurTr] = calSpikeCorrT(neuronNo,filteredSpikeArrayCue_LasttoCurTr,indInt);
    
    save([path fileNameCorr], 'spikeCorrTRun','spikeCorrTRew','spikeCorrTCue',...
        'nonZeroTrRun','nonZeroTrRew','nonZeroTrCue',...
        'spikeCorrTRun_LasttoCurTr','spikeCorrTRew_LasttoCurTr','spikeCorrTCue_LasttoCurTr',...
        'nonZeroTrRun_LasttoCurTr','nonZeroTrRew_LasttoCurTr','nonZeroTrCue_LasttoCurTr',...
        'intervalT','indInt');    

%     save([path fileNameCorr], 'spikeCorrTRun','spikeCorrTRew','spikeCorrTCue',...
%         'nonZeroTrRun','nonZeroTrRew','nonZeroTrCue',...
%         'intervalT'); 
end

function [spikeCorrT,nonZeroTr] = calSpikeCorrT(neuronNo,filteredSpikeArray,intervalT)
    parfor n = 1:neuronNo
        spikeArr = filteredSpikeArray{n}';
        if(length(intervalT) == 1)
            spikeArr = spikeArr(1:intervalT,:);
        else
            spikeArr = spikeArr(intervalT(1):intervalT(2),:);
        end
        nonZeroTr{n} = sum(spikeArr) > 0;
        spikesCorrTmp = corr(spikeArr,'Type','Spearman');
%         for i = 1:trialNo
%             if(sum(spikeArr(i,:)) ~= 0)
%                 for j = i+1:trialNo
%                     if(sum(spikeArr(j,:)) ~= 0)
%                         spikesCorrTmp(i,j) = ...
%                             corr(spikeArr(i,:)',spikeArr(j,:)','Type','Spearman');
%                     end
%                 end
%             end            
%         end
        spikeCorrT{n} = spikesCorrTmp;
    end
end