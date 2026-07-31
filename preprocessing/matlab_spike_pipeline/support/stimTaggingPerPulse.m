function stimTaggingPerPulse(path, fileName)
% find spikes around each of the tagging stimulation pulses for each neuron 
% and plot the spike rasters

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    end
    
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(isempty(indexFileName))
        fileNameInfo = [fileName '_Info.mat'];
        
        fileNameStim = [fileName '_spikePerTagPulse.mat'];
        fileName = [fileName '.mat'];
    else
        fileNameInfo = [fileName(1:indexFileName(end)-1) '_Info.mat'];
        fileNameStim = [fileName(1:indexFileName(end)-1) '_spikePerTagPulse.mat'];
    end
    
    indexFileName = findstr(fileName,'_');
    fileNameBehEl = [fileName(1:indexFileName(1)) 'BehavElectrDataLFP.mat'];
    
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath,'stims','cluList');
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath);
    
    fullPath = [path fileNameBehEl];
    if(exist(fullPath) == 0)
        disp('BehavElectrDataLFP file does not exist.');
        return;
    end
    load(fullPath,'Stim');
    
    GlobalConst;
    
    numStim = length(stims);
    numNeurons = rec.numNeurons;
    
    spikesStim = cell(numNeurons,length(Stim.indStim));
    indStim = Stim.indStim;
    indPulseInStim = Stim.indPulseInStim;
    for i = 1:numNeurons
        countP = 1;
        for n = 1:numStim
            numPulses = length(stims{n}.startLfpInd);
            for j = 1:numPulses
                spikesStim{i,countP} = [...
                    stims{n}.spikesBef{j,i}' stims{n}.spikes{j,i}' ...
                    stims{n}.spikesAft{j,i}'];
                countP = countP +1;
            end
        end
    end
    
    save([path fileNameStim],'spikesStim','indStim','indPulseInStim');
    
    count = 0;
    dispWindow = 50; % ms
    for i = 1:rec.numNeurons  
        count = count+1;
        if(mod(count-1,16) == 0)
            [figNew,pos] = CreateFig();
            set(0,'Units','pixels') 
            figTitle = 'All the pulses';
            set(figure(figNew),'OuterPosition',...
                [pos(1) pos(2) pos(3)*2 pos(4)*2.2],'Name',figTitle)
        end

        subplot(4,4,mod(count-1,16)+1)

        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
            ' ' num2str(cluList.localClu(i)) ')']; 
        
        spikeCells = spikesStim(i,:);
        plotTaggingNeurons(gca,spikeCells,sampleFq,dispWindow,figTitle);
    end

end

function plotTaggingNeurons(handle,spikesStim,sampFq,dispWindow,figTitle,color)
% plot the theta phase of spikes from individual neuron
% handle:           axis handle of the figure
% spikesStim:       spikes around stimulation

    if(nargin == 5)
        color = [0.5 0.5 0.5];
    end
    
    hold on;
      
    nTrials = length(spikesStim);
    for i = 1:nTrials
        h = plot(spikesStim{i}/sampFq*1000,i*ones(1,length(spikesStim{i})), '.');
        set(h,'LineWidth',2.0,'Color',color);
    end
    h = plot([0 0],[0 nTrials+1],'r-');
    set(h,'LineWidth',2.0);
     
    set(gca,'FontSize',8.0,'Box','on','XLim',[-dispWindow dispWindow],'YLim',...
            [0 nTrials+1]);
    xlabel('Time (ms)');
    ylabel('Trials');
    title(figTitle);
end
