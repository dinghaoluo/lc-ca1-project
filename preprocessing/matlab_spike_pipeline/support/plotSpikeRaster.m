function plotSpikeRaster(path, fileName, plotType, trialNo, neuronNo)
% plot spikes raster across run segments
% E.G.: plotSpikeRaster('./','A002-20181005-01_DataStructure_mazeSection1_TrialType')


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
       
    GlobalConst;
    trackLen = 180;
    
    if(isempty(neuronNo))
        neuronNo = 1:length(trials{1}.spikes);
    end
    if(isempty(trialNo))
        trialNo = 1:length(trials);
    end
    segLenArr = [];
    for i = trialNo
        segLenArr = [segLenArr trials{i}.Nsamples];
    end
    segLen = max(segLenArr);
    count = 0;
    for i = neuronNo
        count = count + 1;

        if(mod(count-1,1) == 0)
            [figNew,pos] = CreateFig();
            set(0,'Units','pixels') 
            figTitle = 'Spikes vs Time';
%             figTitle = 'Spikes vs Dist';
            set(figure(figNew),'OuterPosition',...
                [pos(1) pos(2) 280 280],'Name',figTitle)
        end
%         subplot(4,4,mod(count-1,16)+1)  
        hold on;
        for j = 1:length(trialNo)
            ind = trials{trialNo(j)}.spikesThetaLin{i} < 0;
            spikeThetaLin = trials{trialNo(j)}.spikesThetaLin{i};
            spikeThetaLin(ind) = spikeThetaLin(ind) + 2*pi;
            spikeThetaLin = [spikeThetaLin spikeThetaLin + 2*pi];
            if(plotType == 1) % phase plot
                plot([trials{trialNo(j)}.spikes{i}*timeStep ...
                     trials{trialNo(j)}.spikes{i}*timeStep],spikeThetaLin,'k.',...
                'MarkerSize',3);
            elseif(plotType == 2) % spikes over time
                plot(trials{trialNo(j)}.spikes{i}*timeStep,...
                     j*ones(1,length(trials{trialNo(j)}.spikes{i})),'k.',...
                     'MarkerSize',3);
            else % spikes over distance
                plot(trials{trialNo(j)}.spikesMM{i}/10,...
                     j*ones(1,length(trials{trialNo(j)}.spikes{i})),'k.',...
                     'MarkerSize',3);         
            end
        end
        if(plotType == 1) % phase plot
            set(gca,'Xlim',[0 segLen*timeStep],...
                'Ylim',[0,4*pi]);
        elseif(plotType == 2) % spikes over time
            set(gca, 'XLim', [0 segLen*timeStep]);
        else % spike over distance
            set(gca, 'XLim', [0 max(trackLen)]);
            set(gca, 'YLim', [0 length(trialNo)]);
        end
        
        if(mod(count-1,4) == 0)
            if(plotType == 1) % phase plot
                ylabel('Phase');
            else
                ylabel('Trial. no.');
            end
        end
        if(mod(count-1,16) == 0)
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