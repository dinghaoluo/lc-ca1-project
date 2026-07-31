function stimEffect_NewData( path, fileName )
% find the neurons respond to light stimulation
% stimEffect('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1')

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    end
    
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(isempty(indexFileName))
        fileNameFR = [fileName '_FR_Run1.mat'];
        fileNameStim = [fileName '_stimEff.mat'];
        fileName = [fileName '.mat'];
    else
        fileNameFR = [fileName(1:indexFileName(end)-1) '_FR_Run1.mat'];
        fileNameStim = [fileName(1:indexFileName(end)-1) '_stimEff.mat'];
    end
    
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath);
    
    fullPath = [path fileNameFR];
    if(exist(fullPath) == 0)
        disp('The mean firing rate file does not exist. Please call function "MeanFiringRate" first.');
        return;
    end
    load(fullPath);
    
    fileNameInfo = [fileName(1:end-4) '_Info.mat'];      
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath);
    
    GlobalConst;
    
    numStim = length(stims);
    param.timeStep = 0.01;             
    param.timeStepsBef = (-timeAftStim+param.timeStep/2 : param.timeStep : 0)'*sampleFq;
    param.timeStepsAft = (0:param.timeStep:timeAftStim)'*sampleFq;
    param.factorStd = 6; 
    param.indPeakThre = 15; %% the delay (in ms) of peak response within which the response will be considered stimulation driven
    param.fid = 0.5; % fidelity of stimulation
    param.jit = 5; % jitter of response
    param.percAct = 0.3; % the percentage of trials (one trial is one train of stimulation pulses on one shank) which have stimulation effect
    param.percActNonStim = 0.3; % upper bound of the percentage of trials where the shank is not stimulated, but the cell on the shank has been identified as being activated 
    lenTimStepBef = length(param.timeStepsBef);
    lenTimStepAft = length(param.timeStepsAft);
    
    stimResp = struct('indNeurons',1:numNeurons,...
                     'histStimBef',{cell(numNeurons,1)},...
                     'histStimBefMean',zeros(numNeurons,lenTimStepBef),...
                     'histStimBefStd',zeros(numNeurons,lenTimStepBef),...
                     'histStimAft',{cell(numNeurons,1)},...
                     'meanHistStimBef',zeros(numNeurons,1),...
                     'stdHistStimBef',zeros(numNeurons,1),...
                     'maxHistStimBef',zeros(numNeurons,1),...
                     'spPerPulseBef',{cell(numStim,numNeurons)},...
                     'spPerPulseAft',{cell(numStim,numNeurons)},...
                     'firstSpPerPulse',{cell(numStim,numNeurons)},...
                     'numPulseEff',zeros(numStim,numNeurons),...
                     'fidelity',zeros(numStim,numNeurons),...
                     'jitter',zeros(numStim,numNeurons),...
                     'latency',zeros(numStim,numNeurons),...
                     'indPeakStim',-1*ones(numStim,numNeurons),...
                     'stimShank',zeros(numStim,1),...
                     'stimDuration',zeros(numStim,1),...
                     'indPeak',-1*ones(numStim,numNeurons),...
                     'isActivated',zeros(numNeurons,1),...
                     'isStimulated',zeros(numNeurons,1),...
                     'isNotStimulated',zeros(numNeurons,1),...
                     'isActivatedNotStim',zeros(numNeurons,1),...
                     'indNeuronStim',[]);
    
    %% calculate the baseline firing histogram and firing histogram during stimulation
    for j = 1:numNeurons
        histStimsBef = zeros(numStim,lenTimStepBef);
        histStimsAft = zeros(numStim,lenTimStepAft);
        spikesBef = [];
        spikesAft = [];
        spikesStim = [];
        for i = 1:numStim
            numTrainsInStim = max(stims{i}.indTrainInStim);
            for m = 1:numTrainsInStim
                indTrain = stims{i}.indTrainInStim == m; % find the pulses that belong to the same train 
                pulseNo = stims{i}.indPulseInStim(indTrain); % find the corresponding pulse numbers that belong to the pulse train
                indBef = stims{i}.SpikesPulseNoBef{m,j} == pulseNo(1);
                if(isempty(spikesBef))
                    spikesBef = stims{i}.SpikesBef{m,j}(indBef);
                else
                    spikesBef = [spikesBef; stims{i}.SpikesBef{m,j}(indBef)];
                end
                
                indAft = stims{i}.SpikesPulseNoAft{m,j} == pulseNo(end);
                if(isempty(spikesAft))
                    spikesAft = stims{i}.SpikesAft{m,j}(indAft);
                else
                    spikesAft = [spikesAft; stims{i}.SpikesAft{m,j}(indAft)];
                end
                spikesStim = [spikesStim; stims{i}.Spikes{m,j}];                
            end
            
            histStimsBef(i,:) = hist(spikesBef,param.timeStepsBef);
            spikesAft = [stims{i}.Spikes{j}; repmat(stims{i}.SpikesAft{j},numStim,1)];
            histStimsAft(i,:) = hist(spikesAft,param.timeStepsAft)/numStim;
            if(i == 4 & j == 165)
                figure(3)
                plot([histStimsBef(i,:) histStimsAft(i,:)]);
            end
        end
        stimResp.histStimBef{j} = histStimsBef;
        stimResp.histStimAft{j} = histStimsAft;
        stimResp.histStimBefMean(j,:) = mean(histStimsBef);
        stimResp.histStimBefStd(j,:) = std(histStimsBef);
        stimResp.meanHistStimBef(j) = mean(histStimsBef(:));
        stimResp.maxHistStimBef(j) = max(histStimsBef(:));
        stimResp.stdHistStimBef(j) = std(histStimsBef(:));
    end
    
    %% for each train of light stimulation, find the neurons response to the stimulation, get the index of the first peak response
    for i = 1:numStim
        stimResp.stimShank(i) = 6 - stims{i}.indDiode + 1;
        stimResp.stimDuration(i) = stims{i}.pulsePeriod;
        for j = 1:numNeurons
            numPulsesInStim = length(stims{i}.indPulseInStim);
            numTrainsInStim = max(stims{i}.indTrainInStim);
            stimResp.spPerPulseBef{i,j} = cell(numPulsesInStim,1);
            stimResp.spPerPulseAft{i,j} = cell(numPulsesInStim,1);
            stimResp.firstSpPerPulse{i,j} = -1*ones(numPulsesInStim,1);
            stimResp.numPulseEffPerTrain{i,j} = zeros(1,numTrainsInStim);
            stimResp.fidelPerTrain{i,j} = zeros(1,numTrainsInStim);
            for m = 1:numPulsesInStim
                indTrain = stims{i}.indTrainInStim(m);
                tmp = find(stims{i}.indTrainInStim == indTrain);
                indTrain1st = tmp(1);
                indTrainlast = tmp(end);
                indBef = stims{i}.SpikesPulseNoBef{indTrain,j} == stims{i}.indPulseInStim(indTrain1st);
                stimResp.spPerPulseBef{i,j}{m} = stims{i}.SpikesBef{indTrain,j}(indBef);
                indAft = stims{i}.SpikesPulseNoAft{indTrain,j} == stims{i}.indPulseInStim(indTrainlast);
                indStim = stims{i}.SpikesPulseNo{indTrain,j} == stims{i}.indPulseInStim(m);
                stimResp.spPerPulseAft{i,j}{m} = [stims{i}.Spikes{indTrain,j}(indStim); stims{i}.SpikesAft{indTrain,j}(indAft)];
                indFirst = find(stimResp.spPerPulseAft{i,j}{m} <= (stims{i}.pulsePeriod + param.indPeakThre)/1000*sampleFq); % find the first spike after stimulation, and before the end of stimulation + param.indPeakThre
                if(i == 1 & j == 17)
                    a= 1;
                end
                if(~isempty(indFirst))
                    stimResp.firstSpPerPulse{i,j}(m) = stimResp.spPerPulseAft{i,j}{m}(indFirst(1));
                    stimResp.numPulseEffPerTrain{i,j}(indTrain) = stimResp.numPulseEffPerTrain{i,j}(indTrain) + 1;
                    stimResp.fidelPerTrain{i,j}(indTrain) = stimResp.fidelPerTrain{i,j}(indTrain) + 1/length(tmp);
                    stimResp.numPulseEff(i,j) = stimResp.numPulseEff(i,j) + 1;
                    stimResp.fidelity(i,j) = stimResp.fidelity(i,j) + 1/numPulsesInStim;
                end
            end
            
            indResp = stimResp.firstSpPerPulse{i,j} ~= -1;
            stimResp.jitter(i,j) = std(stimResp.firstSpPerPulse{i,j}(indResp));
            stimResp.latency(i,j) = mean(stimResp.firstSpPerPulse{i,j}(indResp));
            stimResp.jitterPerTrain{i,j} = zeros(1,numTrainsInStim);
            stimResp.latencyPerTrain{i,j} = zeros(1,numTrainsInStim);
            for m = 1:numTrainsInStim
                indResp = stimResp.firstSpPerPulse{i,j} ~= -1 & stims{i}.indTrainInStim == m;
                stimResp.jitterPerTrain{i,j}(m) = std(stimResp.firstSpPerPulse{i,j}(indResp));
                stimResp.latencyPerTrain{i,j}(m) = mean(stimResp.firstSpPerPulse{i,j}(indResp));
            end
                        
            if(cluList.shank(j) == stimResp.stimShank(i))
                stimResp.isStimulated(j) = stimResp.isStimulated(j) + 1;
            else
                stimResp.isNotStimulated(j) = stimResp.isNotStimulated(j) + 1;
            end
            indPeakTmp = find(stimResp.histStimAft{j}(i,:) > stimResp.meanHistStimBef(j) + stimResp.stdHistStimBef(j)*param.factorStd);
            if(~isempty(indPeakTmp))
                [~,indPeak] = max(stimResp.histStimAft{j}(i,indPeakTmp));
                stimResp.indPeakStim(i,j) = indPeakTmp(indPeak(1));
            end
            if(stimResp.fidelity(i,j) >= param.fid & stimResp.indPeakStim(i,j) <= (stims{i}.pulsePeriod + param.indPeakThre)/1000*sampleFq & stimResp.indPeakStim(i,j) >= 0)
                if(cluList.shank(j) == stimResp.stimShank(i))
                    stimResp.isActivated(j) = stimResp.isActivated(j) + 1;
                else
                    stimResp.isActivatedNotStim(j) = stimResp.isActivatedNotStim(j) + 1;
                end
                figure(1)
                plotHistStim([stimResp.histStimBef{j}(i,:) stimResp.histStimAft{j}(i,:)],[param.timeStepsBef; param.timeStepsAft]',...
                    j,stimResp.stimShank(i),stimResp.indPeakStim(i,j),stims{i}.pulsePeriod,sampleFq);
                figure(2)
                plotStimRaster(stimResp.spPerPulseBef{i,j},stimResp.spPerPulseAft{i,j},j,cluList.shank(j),stimResp.stimShank(i),stimResp.stimDuration(i),i,sampleFq,[-timeBefStim timeAftStim]);
%                 a = 1;
%                 pause;
            end
        end
    end
    
    stimResp.indNeuronStim = find(stimResp.isActivated./stimResp.isStimulated > param.percAct & stimResp.isActivatedNotStim./stimResp.isNotStimulated < param.percActNonStim);
    
    save([path fileNameStim], 'stimResp','param');
    clear('stimResp', 'param');
end

function plotStimRaster(spPerPulseBef,spPerPulseAft,neuronNo,shank,shankStim,pulsePeriod,stimNo,sampleFq,xLim)    
    hold off;
    numStimPul = length(spPerPulseBef);
    for i = 1:numStimPul
        h = plot(spPerPulseBef{i}/sampleFq*1000,i*ones(length(spPerPulseBef{i})),'k.');
        hold on;
        set(h,'MarkerSize',8);
        h = plot(spPerPulseAft{i}/sampleFq*1000,i*ones(length(spPerPulseAft{i})),'k.');
        set(h,'MarkerSize',8);
    end
    plot([0 0],[0 numStimPul+1],'r-')
    plot([pulsePeriod pulsePeriod],[0 numStimPul+1],'g.')
    set(gca,'FontSize',8.0,'Box','on','XLim',xLim*1000,'YLim',[0 numStimPul+1]);
    xlabel('Time (ms)');
    ylabel('Pulse no.');
    title(['neu ' num2str(neuronNo) ' ShSt ' num2str(shankStim) ' Sh ' num2str(shank) ' stimNo ' num2str(stimNo)]);
end

function plotHistStim(histStims,timeSteps,neuronNo,shank,indPeak,pulsePeriod,sampleFq)
    hold off;   
    color = [0 0 0];
    maxVal = max(1,max(histStims));
    h = bar(timeSteps/sampleFq*1000,histStims,1);
    set(h,'FaceColor',color,'EdgeColor',color,'LineWidth',1);
    hold on
    plot([0 0],[0 maxVal],'r-')
    plot([pulsePeriod pulsePeriod],[0 maxVal],'g.')
    set(gca,'FontSize',8.0,'Box','on','XLim',[timeSteps(1) timeSteps(end)]/sampleFq*1000,'YLim',[0 maxVal]);
    xlabel('Time (ms)');
    ylabel('Spike count');
    title(['neu ' num2str(neuronNo) ' Sh ' num2str(shank) ' indP ' num2str(indPeak)]);
end
