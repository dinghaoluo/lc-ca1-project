function meanNeuronSpikeCorrDist_GoodTr(path,fileName,onlyRun)
% single neuron level mean spike correlation across trials
% meanNeuronSpikeCorrDist_GoodTr('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1)

    fullPath = [path fileName '_spikesCorrDist_GoodTr_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _spikesCorrDist_GoodTr_Run file does not exist');
        return;
    end
    load(fullPath);
                
    fileNameCorr = [fileName '_meanSpikesCorrDist_GoodTr_Run' num2str(onlyRun) '.mat'];
    
    fileNameInfo = [fileName '_Info.mat'];     
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
    
    disp('Calculate mean spike correlation over distance for each session')
    meanCorrDistGoodTr = cell(length(mazeSess),1);
    meanCorrDistBadTr = cell(length(mazeSess),1);
    meanCorrDistOKTr = cell(length(mazeSess),1);
    for i = 1:length(mazeSess) 
        trialNoGood = length(indLapsGoodTr{i});
        nElemGood = (trialNoGood*trialNoGood-trialNoGood)/2;
        trialNoBad = length(indLapsBadTr{i});
        nElemBad = (trialNoBad*trialNoBad-trialNoBad)/2;
        trialNoOK = length(indLapsOKTr{i});
        nElemOK = (trialNoOK*trialNoOK-trialNoOK)/2;
        
        if(trialNoGood > 0)
            neuronNo = size(spikeCorrDistGoodTr{i},2);
            for n = 1:neuronNo            
                corrArr = triu(spikeCorrDistGoodTr{i}{n},1);
                corrArr = corrArr(:); 
                meanCorrDistGoodTr{i}.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElemGood;
                nNonZeroTr = sum(nonZeroGoodTr{i}{n} == 1);
                nElemNonZero = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
                meanCorrDistGoodTr{i}.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZero;
                meanCorrDistGoodTr{i}.nNonZeroTr(n) = nNonZeroTr;
            end
        end

        if(trialNoBad > 0)
            neuronNo = size(spikeCorrDistBadTr{i},2);
            for n = 1:neuronNo 
                corrArr = triu(spikeCorrDistBadTr{i}{n},1);
                corrArr = corrArr(:); 
                meanCorrDistBadTr{i}.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElemBad;
                nNonZeroTr = sum(nonZeroBadTr{i}{n} == 1);
                nElemNonZero = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
                meanCorrDistBadTr{i}.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZero;
                meanCorrDistBadTr{i}.nNonZeroTr(n) = nNonZeroTr;
            end
        end
        
        if(trialNoOK > 0)
            neuronNo = size(spikeCorrDistOKTr{i},2);
            for n = 1:neuronNo 
                corrArr = triu(spikeCorrDistOKTr{i}{n},1);
                corrArr = corrArr(:); 
                meanCorrDistOKTr{i}.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElemOK;
                nNonZeroTr = sum(nonZeroOKTr{i}{n} == 1);
                nElemNonZero = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
                meanCorrDistOKTr{i}.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZero;
                meanCorrDistOKTr{i}.nNonZeroTr(n) = nNonZeroTr;
            end
        end
        
        if(~isempty(meanCorrDistGoodTr{i}) && ~isempty(meanCorrDistBadTr{i}))
            % compare single neuron correlation between good and bad trials --
            % aligned to run
            plotCompCorr(meanCorrDistGoodTr{i}.mean,meanCorrDistBadTr{i}.mean,...
                'Neu Tr CorrDist - good trials','Neu Tr CorrDist - bad trials',['Session ' num2str(i)]);

            % compare single neuron correlation between aligned to run and aligned
            % to cue for good trials
            indNeu = meanCorrDistGoodTr{i}.nNonZeroTr > 10 & meanCorrDistBadTr{i}.nNonZeroTr > 10;
            if(sum(indNeu) ~= 0)
                plotCompCorr(meanCorrDistGoodTr{i}.meanNZ(indNeu),meanCorrDistBadTr{i}.meanNZ(indNeu),...
                    'Neu None-zero Tr CorrDist - good trials','Neu None-zero Tr CorrDist - bad trials',['Session ' num2str(i)]); 
            end
        end
    
    end
    
    save([path fileNameCorr],'meanCorrDistGoodTr','meanCorrDistBadTr','meanCorrDistOKTr');
    
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     plotCompCorr(meanCorrDistRew.meanGood,meanCorrDistRew.meanBad,...
%         'Neu Tr CorrDist Rew - good trials','Neu Tr CorrDist Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     plotCompCorr(meanCorrDistCue.meanGood,meanCorrDistCue.meanBad,...
%         'Neu Tr CorrDist Cue - good trials','Neu Tr CorrDist Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     plotCompCorr(meanCorrDistRun.mean,meanCorrDistRew.mean,...
%         'Neu Tr CorrDist Run','Neu Tr CorrDist Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     plotCompCorr(meanCorrDistRun.mean,meanCorrDistCue.mean,...
%         'Neu Tr CorrDist Run','Neu Tr CorrDist Cue');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward for good trials
%     plotCompCorr(meanCorrDistRun.meanGood,meanCorrDistRew.meanGood,...
%         'Neu Tr CorrDist Run - good trials','Neu Tr CorrDist Rew - good trials');
%    
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue for good trials
%     plotCompCorr(meanCorrDistRun.meanGood,meanCorrDistCue.meanGood,...
%         'Neu Tr CorrDist Run - good trials','Neu Tr CorrDist Cue - good trials');
    
%     %% only consider nonzero trials
%     % compare single neuron correlation between good and bad trials --
%     % aligned to run
%     indNeu = meanCorrDistRun.nGoodNonZeroTr > 10 & meanCorrDistRun.nBadNonZeroTr > 5;
%     plotCompCorr(meanCorrDistRun.meanGoodNZ(indNeu),meanCorrDistRun.meanBadNZ(indNeu),...
%         'Neu None-zero Tr CorrDist Run - good trials','Neu None-zero Tr CorrDist Run - bad trials');
%     
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     indNeu = meanCorrDistRew.nGoodNonZeroTr > 10 & meanCorrDistRew.nBadNonZeroTr > 5;
%     plotCompCorr(meanCorrDistRew.meanGoodNZ(indNeu),meanCorrDistRew.meanBadNZ(indNeu),...
%         'Neu None-zero Tr CorrDist Rew - good trials','Neu None-zero Tr CorrDist Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     indNeu = meanCorrDistCue.nGoodNonZeroTr > 10 & meanCorrDistCue.nBadNonZeroTr > 5;
%     plotCompCorr(meanCorrDistCue.meanGoodNZ(indNeu),meanCorrDistCue.meanBadNZ(indNeu),...
%         'Neu None-zero Tr CorrDist Cue - good trials','Neu None-zero Tr CorrDist Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     indNeu = meanCorrDistRun.nNonZeroTr > 10 & meanCorrDistRew.nNonZeroTr > 10;
%     plotCompCorr(meanCorrDistRun.meanNZ(indNeu),meanCorrDistRew.meanNZ(indNeu),...
%         'Neu None-zero Tr CorrDist Run ','Neu None-zero Tr CorrDist Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     indNeu = meanCorrDistRun.nNonZeroTr > 10 & meanCorrDistCue.nNonZeroTr > 10;
%     plotCompCorr(meanCorrDistRun.meanNZ(indNeu),meanCorrDistCue.meanNZ(indNeu),...
%         'Neu None-zero Tr CorrDist Run','Neu None-zero Tr CorrDist Cue');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward for good trials
%     indNeu = meanCorrDistRun.nGoodNonZeroTr > 10 & meanCorrDistRew.nGoodNonZeroTr > 10;
%     plotCompCorr(meanCorrDistRun.meanGoodNZ(indNeu),meanCorrDistRew.meanGoodNZ(indNeu),...
%         'Neu None-zero Tr CorrDist Run - good trials','Neu None-zero Tr CorrDist Rew - good trials');
%    
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue for good trials
%     indNeu = meanCorrDistRun.nGoodNonZeroTr > 10 & meanCorrDistCue.nGoodNonZeroTr > 10;
%     plotCompCorr(meanCorrDistRun.meanGoodNZ(indNeu),meanCorrDistCue.meanGoodNZ(indNeu),...
%         'Neu None-zero Tr CorrDist Run - good trials','Neu None-zero Tr CorrDist Cue - good trials');    
    
end

function plotCompCorr(x,y,xlab,ylab,ti)
    figure
    plot(x,y,'ro');
    hold on;
    corrAll = [x y];
    maxCorr = max(corrAll);
    minCorr = min(corrAll);
    plot([minCorr maxCorr],[minCorr maxCorr],'k:');
    xlim(gca,[minCorr maxCorr]);
    ylim(gca,[minCorr maxCorr]);
    xlabel(xlab);
    ylabel(ylab);
    title(ti);
end
