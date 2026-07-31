function plotdFFAligned2P(path,fileName,mazeSess,cond,intervalTSpInfo,begT,onlyRun)

    %%%%%%%%% initialize constants
    GlobalConst2P;
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    if(cond == 1)
        condStr = 'Run';
    elseif(cond == 2)
        condStr = 'Rew';
    elseif(cond == 3)
        condStr = 'Cue';
    end
    
    if(cond == 1)
        fileNameConv = [fileName '_convSpikesAligned_msess' num2str(mazeSess) '_BefRun' num2str(onlyRun) '.mat'];
        fullPath = [path fileNameConv];
        if(exist(fullPath) == 0)
            disp(['The firing profile file does not exist. Please call ',...
                    'function "ConvSpikeTrain_Aligned" first.']);
            return;
        end
        load(fullPath,'dFFArrayRunOnset','paramC');
        dFFArray = dFFArrayRunOnset;
    elseif(cond == 2)
        fileNameConv = [fileName '_convSpikesAligned_msess' num2str(mazeSess) '_BefRew' num2str(onlyRun) '.mat'];
        fullPath = [path fileNameConv];
        if(exist(fullPath) == 0)
            disp(['The firing profile file does not exist. Please call ',...
                    'function "ConvSpikeTrain_Aligned" first.']);
            return;
        end
        load(fullPath,'dFFArrayReward','paramC');
        dFFArray = dFFArrayReward;
    else
        fileNameConv = [fileName '_convSpikesAligned' condStr '_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
        fullPath = [path fileNameConv];
        if(exist(fullPath) == 0)
            disp(['The firing profile file does not exist. Please call ',...
                    'function "ConvSpikeTrain_Aligned" first.']);
            return;
        end
        if(cond == 3)
            load(fullPath,'dFFArrayCue','paramC');
            dFFArray = dFFArrayCue;
        end
        
        begT = 0;
    end
    
    fileNameFR = [fileName '_FRAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    fullPath = [path fileNameFR];
    if(exist(fullPath) == 0)
        disp(['The mean firing rate file does not exist. Please call ',...
                'function "MeanFiringRateAligned" first.']);
        return;
    end
    load(fullPath,'mFRStructNonStimGood');
    
    fileNamePeakFR = [fileName '_PeakFRAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    fullPath = [path fileNamePeakFR];
    if(exist(fullPath) == 0)
        disp('The peak firing rate aligned to run file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad');
    
    fullPath = [path fileName '.mat']; 
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
        return;
    end
    load(fullPath,'cluList');
    
    %%% all the trials
    neuronNo = 1:length(mFRStructNonStimGood.mFR);
    figNo = 1;
    if(~isempty(dFFArray))
    %% good trials
        count = 0;
        indGoodLap = [trialNoNonStimGood;trialNoNonStimBad];
        indTrialCut = length(trialNoNonStimGood);
        for i = neuronNo 
            if(mFRStructNonStimGood.mFR(i) > 0) %minFR
                count = count + 1;

                if(mod(count-1,16) == 0)
                    if(count > 1)
                        fullpath = [path fileName(1:end-4) '_' condStr '_dFF' num2str(figNo)];
                        print('-painters', '-dpdf', fullpath, '-r600')
                        savefig([fullpath '.fig']);
                        figNo = figNo+1;
                    end

                    [figNew,pos] = CreateFig2P();
                    set(0,'Units','pixels') 
                    figTitle = 'All the trials';
                    set(figure(figNew),'OuterPosition',...
                        [pos(1) pos(2)-300 pos(3)*1.6 pos(4)*1.6],'Name',figTitle)
                end

                subplot(4,4,mod(count-1,16)+1)

                figTitle = ['Neu ' num2str(i) '(' num2str(cluList.localClu(i)-1) ')'];                
                plotFRProfIndNeuronIndTrial(gca,...
                    (dFFArray{i}-min(dFFArray{i},[],2))./(max(dFFArray{i},[],2)-min(dFFArray{i},[],2)),... % changed on 9/7/2023
                    intervalTSpInfo,indGoodLap',figTitle,...
                    -begT:timeStep:timeStep*(size(dFFArray{i},2)-1)-begT);
                if(~isempty(indTrialCut))
                    hold on;
                    plot([-begT timeStep*(size(dFFArray{i},2)-1)-begT],[indTrialCut indTrialCut],'r');
                end
            else
                disp(['Firng rate of neuron ' num2str(i) ' is too low: ' ...
                    num2str(mFRStructNonStimGood.mFR(i)) ' Hz']);
                continue;
            end
        end

        fullpath = [path fileName(1:end-4) '_' condStr '_dFF' num2str(figNo)];
        print('-painters', '-dpdf', fullpath, '-r600')
        savefig([fullpath '.fig']);
    end
end
    
function plotFRProfIndNeuronIndTrial...
            (handle,filteredSpikeArrayT,intervalT,indLaps,figTitle,timeStep)
    numTr = length(indLaps);
    h = imagesc(timeStep,1:numTr,filteredSpikeArrayT(indLaps,:));
%      colormap jet 
    set(gca,'FontSize',8.0,'Box','on','YLim',[0 numTr],'XLim',[min(timeStep) 10]);
    xlabel('Time (s)');
    ylabel('Trial No.');
    title(figTitle);
end