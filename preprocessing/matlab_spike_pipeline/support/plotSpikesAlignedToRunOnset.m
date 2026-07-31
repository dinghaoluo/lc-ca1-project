function plotSpikesAlignedToRunOnset(path, fileName, plotType)
% plot spikes across run segments
% E.G.: plotSpikesAlignedToRunOnset('./','A002-20181005-01_DataStructure_mazeSection1_TrialType_RunOnSet')


    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end

    fullPath = [path fileName '.mat']; 
    if(exist(fullPath) == 0)
        disp('The file does not exist');
        return;
    end
    load(fullPath);
    
    indexFileName = findstr(fileName, '_');
    fullPath = [path fileName(1:indexFileName(end)-1) '1.mat'];
    if(exist(fullPath) == 0)
        disp('The file does not exist');
        return;
    end
    load(fullPath,'cluList');
    
    GlobalConst;
    
    numNeurons = length(trials{1}.spikes);
    numSeg = length(trials);
    segLenArr = [];
    for i = 1:numSeg
        segLenArr = [segLenArr trials{i}.Nsamples];
    end
    segLen = max(segLenArr);
    count = 0;
    for i = 1:numNeurons
        count = count + 1;

        if(mod(count-1,16) == 0)
            [figNew,pos] = CreateFig();
            set(0,'Units','pixels') 
            figTitle = 'Spikes vs Time';
%             figTitle = 'Spikes vs Dist';
            set(figure(figNew),'OuterPosition',...
                [pos(1) pos(2) pos(3)*2 pos(4)*2.2],'Name',figTitle)
        end
        subplot(4,4,mod(count-1,16)+1)  
        hold on;
        for j = 1:numSeg
            ind = trials{j}.spikesThetaLin{i} < 0;
            spikeThetaLin = trials{j}.spikesThetaLin{i};
            spikeThetaLin(ind) = spikeThetaLin(ind) + 2*pi;
            spikeThetaLin = [spikeThetaLin spikeThetaLin + 2*pi];
            if(plotType == 1) % phase plot
                plot([trials{j}.spikes{i}*timeStep ...
                     trials{j}.spikes{i}*timeStep],spikeThetaLin,'k.',...
                'MarkerSize',8);
            elseif(plotType == 2) % spikes over time
                plot(trials{j}.spikes{i}*timeStep,...
                     j*ones(1,length(trials{j}.spikes{i})),'k.',...
                     'MarkerSize',8);
            else % spikes over distance
                plot(trials{j}.spikesMMAdj{i},...
                     j*ones(1,length(trials{j}.spikes{i})),'k.',...
                     'MarkerSize',8);         
            end
        end
        if(plotType == 1) % phase plot
            set(gca,'Xlim',[0 segLen*timeStep],...
                'Ylim',[0,4*pi]);
        elseif(plotType == 2) % spikes over time
            set(gca, 'XLim', [0 segLen*timeStep]);
        else % spike over distance
            set(gca, 'XLim', [0 max(trackLen)]);
        end
        
        if(mod(count-1,4) == 0)
            if(plotType == 1) % phase plot
                ylabel('Phase');
            else
                ylabel('Seg. no.');
            end
        end
        if(mod(count-1,16) > 11)
            if(plotType == 3) % spike over distance
                xlabel('Dist (mm)');
            else
                xlabel('Time (s)');
            end
        end
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ')'];
        title(figTitle);
    end