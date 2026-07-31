function stimDispRaster( path, fileName, stimNum, neuronNo)
% display the stimulation result
% stimDisp('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1',3)

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin == 2
        neuronNo = [];
        stimNo = [];
    elseif nargin == 3
        neuronNo = [];
    elseif nargin > 4
        disp('Too many input arguments.');
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
        fileName = [fileName '.mat'];
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
%     
%     fullPath = [path fileNameStim];
%     if(exist(fullPath) == 0)
%         disp('The stimulation effect file does not exist. Please call function "stimEffect" first.');
%         return;
%     end
%     load(fullPath);
    
    fileNameInfo = [fileName(1:end-4) '_Info.mat'];      
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath);
    
    if(isempty(neuronNo))
        neuronNo = rec.numNeurons;
    end
    if(isempty(stimNum))
        numStim = length(stims);
    end
    
    GlobalConst;
    
%     stimNum = 4;
    count = 0;
    for i = 1:length(neuronNo)   
%         if(~isempty(find(indInhNeurons == i, 1)))% inh neuron
%             neuroclass = 'inh';
%         elseif(~isempty(find(indExcNeurons == i, 1))) % exc neuron
%             neuroclass = 'exc';
%         else
%             disp(['Firng rate of neuron ' num2str(i) ' is too low: ' num2str(mFRStruct.mFR(i)) ' Hz']);
%             continue;
%         end   
        
        count = count + 1;
        if(mod(count,16) == 1)
            [figNew,pos] = CreateFig();
            set(0,'Units','pixels') 
            set(figure(figNew),'OuterPosition',[pos(1) pos(2) pos(3)*2 pos(4)*2.2])
        end
        
        subplot(4,4,mod(count-1,16)+1)
        param.timeAftStim = timeAftStim;
        param.timeBefStim = timeBefStim;
        param.timeStep = 0.005;
        param.sampleFq = sampleFq;
        param.pulsePeriod = stims{stimNum(1)}.pulsePeriod/1000;
        plotRaster(stimResp.spPerPulseBef, stimResp.spPerPulseAft, stimNum, neuronNo(i), param, cluList.shank(neuronNo(i)));
        %plotRaster(stims{stimNum}.SpikesBef{i},stims{stimNum}.Spikes{i},stims{stimNum}.SpikesAft{i},param,i,cluList.shank(i));
    end
end

function plotRaster(spPerPulseBef,spPerPulseAft,stimNum,neuronNo,param,shank)
    
    color = [0 0 0];
    count = 0;
    hold on;
    for n = 1:length(stimNum)
        for i = 1:length(spPerPulseBef{stimNum(n),neuronNo})
            count = count + 1;
            spikesTmp = [spPerPulseBef{stimNum(n),neuronNo}{i}; spPerPulseAft{stimNum(n),neuronNo}{i}];
            h = plot(spikesTmp,count*ones(1,length(spikesTmp)),'k.');
            set(h,'MarkerSize',10);
        end
    end
    plot([0 0],[0 count+1],'r-')
    set(gca,'FontSize',8.0,'Box','on','XLim',[-param.timeAftStim param.timeAftStim+param.pulsePeriod]*1000,'YLim',[0 count+1]);
    xlabel('Time (ms)');
    ylabel('Stimulation trials');
    title(['neu ' num2str(neuronNo) ' Sh ' num2str(shank)]);
end
