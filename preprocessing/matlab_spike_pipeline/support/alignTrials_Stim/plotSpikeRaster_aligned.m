function plotSpikeRaster_aligned(path, fileName, onlyRun, mazeSess, trialNo1, cond, frProfile, neuronNo)
% plot spikes rasters across run segments
% cond = 1: only the align-to-run plots are extended to before run
% cond = 0: both the align-to-run and align-to-rew plots are extended to
% before the align onset
% added "cond" on 3/4/2022
% E.G.: plotSpikeRaster_aligned('./','A002-20181005-01_DataStructure_mazeSection1_TrialType',1,1,1)

    % edited nargin to +1 due to addition of frProfile, Dinghao 11 Mar 2022
    if(nargin == 5)
        trialNo1 = [];
        cond = 1; %% added by Yingxue on 3/4/2022
        neuronNo = [];
    elseif(nargin == 6)
        cond = 1; %% added by Yingxue on 3/4/2022
        neuronNo = [];
    elseif(nargin == 7)
        neuronNo = [];
    end

    %%%%%%%%% load recording file
    fullPath = [path fileName '.mat']; 
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
        return;
    end
    load(fullPath,'cluList','lap');
    if(isempty(neuronNo))
        neuronNo = 1:length(cluList.all);
    end
       
    fullPath = [path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT file does not exist');
        return;
    end
    load(fullPath);
    if(isempty(trialNo1))
        trialNo1 = 1:size(trialsRunSpikes.Time,2);
    end
    
    % added for loading BefRew fr curve aligned to Rew, Dinghao 7 Mar, 2022
    fullPath = [path fileName '_convSpikesAligned_msess' num2str(mazeSess) '_BefRew' num2str(onlyRun) '.mat'];
    load(fullPath,'timeStepRun','paramC');
    timeStepRunRew = timeStepRun;
    
    fullPath = [path fileName '_convSpikesAligned_msess' num2str(mazeSess) '_BefRun' num2str(onlyRun) '.mat'];

    if(exist(fullPath,'file') == 0)
        disp(['The firing profile file does not exist. Try to run the',...
                    'function again with fileState = 0.']);
    end
    load(fullPath,'timeStepRun','paramC');
    
    GlobalConst;
    
    segLen = 300;
    trialLenT = 10; %sec
    count = 0;
    numFig = 0;
    rasterColour = [0 0.4470 0.7410];
    frColour = [0.8500 0.3250 0.0980];
    greyColour = [0.7 0.7 0.7];
    indTrialCut = find(diff(trialNo1)<0);
    if(~isempty(indTrialCut))
        indTrialCut = indTrialCut + 1;
    end
    for i = neuronNo
        disp(['Neuron ' num2str(i)]);
        count = count + 1;
        
        if(mod(count-1,18) == 0)
            numFig = numFig + 1;
            if(numFig ~= 1)
                fullpath = [path fileName '_Run' num2str(onlyRun) '_Raster' num2str(numFig-1)];
                print('-painters', '-dpdf', fullpath, '-r600')
                savefig([fullpath '.fig']);
            end
            [figNew,pos] = CreateFig();
            set(0,'Units','pixels') 
            figTitle = 'Spikes vs Time';
%             figTitle = 'Spikes vs Dist';
            set(figure(figNew),'OuterPosition',...
                [pos(1) pos(2)-300 pos(3)*1.5 pos(4)*1.7],'Name',figTitle)
        end
        
        subplot(6,3,mod(count-1,18)+1)
        yyaxis left
        hold on;
        for j = 1:length(trialNo1)   
%             disp(['Trial ' num2str(j)])
            if j == 103
                a = 1;
            end
            % spikes over distance
            plot([trialsRunSpikes.TimeBef{i,trialNo1(j)}'/sampleFq ...
                 trialsRunSpikes.Time{i,trialNo1(j)}'/sampleFq],...
                 j*ones(1,length([trialsRunSpikes.TimeBef{i,trialNo1(j)}' ...
                 trialsRunSpikes.Time{i,trialNo1(j)}'])),...
                 '.',...
                 'Color', greyColour,...
                 'MarkerSize',3);
             
%             plot([trialsRunSpikes.TimeBef{i,j}'/sampleFq ...
%                  trialsRunSpikes.Time{i,j}'/sampleFq],...
%                  [trialsRunSpikes.thPhaseInterpSpikeBef{i,j}' ...
%                  trialsRunSpikes.thPhaseInterpSpike{i,j}'],'k.',...
%                  'MarkerSize',3);     
        end
        if(~isempty(indTrialCut))
            plot([0 trialLenT],[indTrialCut indTrialCut],'r');
        end
        if(mod(count-1,3) == 0)
            ylabel('Trial no.');
        end
        
        % added for firing rate profile switch, Dinghao 11 March 2022
        if(frProfile == 1 || frProfile == 3)
            fullPath = [path fileName '_PeakFR_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat'];
            if(exist(fullPath) == 0)
               disp('The peak firing rate aligned to run file does not exist');
               return;
            end
            load(fullPath,'pFRNonStimGoodStruct','pFRNonStimBadStruct');
            yyaxis right
            h = plot(timeStepRun/sampleFq,pFRNonStimGoodStruct.avgFRProfile(i,:),'-', 'Color', frColour);
            set(h,'LineWidth',0.5);
        end
        
        if(mod(count-1,3) == 0)
            ylabel('FR (Hz)');
        end
        set(gca, 'XLim', [-3 trialLenT]);
        
        if(mod(count-1,18) > 14)
            xlabel('Time (s)');
        end
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ')'];
        title(figTitle);
        
        count = count + 1;
        subplot(6,3,mod(count-1,18)+1)
        hold on;
        for j = 1:(length(trialNo1)-1)   
            %% changed by Yingxue on 3/4/2022
            if(cond == 1)
                % spikes over distance
                plot(trialsRewSpikes.Time{i,trialNo1(j)}/sampleFq,...
                     j*ones(1,length(trialsRewSpikes.Time{i,trialNo1(j)})),'k.',...
                     'MarkerSize',3);    
            else
                % spikes over distance + bef, added 3 Mar 2022, Dinghao, for
                % monitoring pre-reward ramps
                plot([trialsRewSpikes.TimeBef{i,trialNo1(j)}'/sampleFq ...
                     trialsRewSpikes.Time{i,trialNo1(j)}'/sampleFq],...
                     j*ones(1,length([trialsRewSpikes.TimeBef{i,trialNo1(j)}' ...
                     trialsRewSpikes.Time{i,trialNo1(j)}'])),...
                     '.',...
                     'Color', greyColour,...
                     'MarkerSize',3);    
            end
        end
        if(~isempty(indTrialCut))
            plot([0 trialLenT],[indTrialCut indTrialCut],'r');
        end
        
        % load reward onset fr struct, 10 March 2022 Dinghao
        % added for firing rate profile switch, Dinghao 11 March 2022
        if(frProfile == 2 || frProfile == 3)
            fullPath = [path fileName '_PeakFR_msess' num2str(mazeSess) '_RewOnset' num2str(onlyRun) '.mat'];
            if(exist(fullPath) == 0)
               disp('The peak firing rate aligned to rew file does not exist');
               return;
            end
            load(fullPath,'pFRNonStimGoodStruct','pFRNonStimBadStruct');
            yyaxis left
            set(gca,'YColor',rasterColour)
            ylabel('Trial no.')
            yyaxis right
            set(gca,'YColor',frColour)
            ylabel('FR (Hz)')
            h = plot(timeStepRunRew/sampleFq,pFRNonStimGoodStruct.avgFRProfile(i,:),'-', 'Color', frColour);
            set(h,'LineWidth',0.5);
        end
        %% changed by Yingxue on 3/4/2022
        if(cond == 1)
            set(gca, 'XLim', [0 trialLenT]);
        else
            set(gca, 'XLim', [-3 trialLenT]);
        end
        %%
        if(mod(count-1,3) == 0)
            ylabel('Trial no.');
        end
        if(mod(count-1,18) > 14)
            xlabel('Time (s)');
        end
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ')'];
        title(figTitle);
        
        count = count + 1;
        subplot(6,3,mod(count-1,18)+1)
        hold on;
        for j = 1:length(trialNo1)
            % spikes over distance
            plot(trialsCueSpikes.Time{i,trialNo1(j)}/sampleFq,...
                 j*ones(1,length(trialsCueSpikes.Time{i,trialNo1(j)})),'k.',...
                 'MarkerSize',3);    
        end
        if(~isempty(indTrialCut))
            plot([0 trialLenT],[indTrialCut indTrialCut],'r');
        end
        set(gca, 'XLim', [0 trialLenT]);
        if(mod(count-1,3) == 0)
            ylabel('Trial no.');
        end
        if(mod(count-1,18) > 14)
            xlabel('Time (s)');
        end
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ')'];
        title(figTitle);
    end
    fullpath = [path fileName '_Run' num2str(onlyRun) '_Raster' num2str(numFig)];
    print('-painters', '-dpdf', fullpath, '-r600')
    savefig([fullpath '.fig']);