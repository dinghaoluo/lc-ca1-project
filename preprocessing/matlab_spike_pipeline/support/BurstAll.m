function BurstAll(path,fileName,methodTheta,figureState,onlyRun)
% Find all the bursts which satisfy the criteria (the isi < 6 ms)
% Dependence:  function "ThetaPhaseLR" should be executed first
% path:         the path of the recording file
% fileName:     name of the recording file
% methodTheta:  method used for theta phase estimation
%               0: Hilbert transform
%               1: Linear interpolation
% figureState:  0: figure off
%               1: plot the histogram of the burst time
%               2: plot the burst time and the burst phase during the analysis
% onlyRun:      1: only consider the time period when the animal is running 
%
% Example:
% BurstAllVR('./','xzvr_PR1-20170727-01_DataStructure_mazeSection1_TrialType1',1,0,1)

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin == 2
        methodTheta = 1;
        figureState = 1;
        onlyRun = 1;
    elseif nargin == 3
        figureState = 1;
        onlyRun = 1;
    elseif nargin == 4
        onlyRun = 1;
    elseif nargin > 5
        disp('Too many input arguments.');
        return;
    end
    
    %%%%%%%%% initialize constants
    GlobalConst;
    
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    if(methodTheta == 0)
        th = 'H';
    else
        th = 'L';
    end
    fileNameBurst = [fileName '_burstAll_TH' th '_Run' num2str(onlyRun) ...
                     '.mat'];

    if(methodTheta == 0)
        fileNameThetaPhase = [fileName '_ThetaPhaseH_Run' num2str(onlyRun) ...
                              '.mat'];
    else
        fileNameThetaPhase = [fileName '_ThetaPhaseL_Run' num2str(onlyRun) ...
                              '.mat'];
    end
    fileNameExt = [fileName '_ext.mat'];
    fileName = [fileName '.mat'];

    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
        return;
    end
    load(fullPath);
    
    fullPath = [path fileNameExt];
    if(exist(fullPath) == 0)
        disp(['Extended file does not exist. Please run function' ...
             ' SpikeDuringRun first']);
        return;
    end
    load(fullPath);
    
    fullPath = [path fileNameThetaPhase];
    if(exist(fullPath) == 0)
        disp(['The theta phase file does not exist. Please call function' ...
              ' "ThetaPhase" first.']);
        return;
    end
    load(fullPath);
    
    fileNameInfo = [fileName(1:end-4) '_Info.mat'];    
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
    
    % burst parameters
    paramBurst = struct('numNeurons', rec.numNeurons,...
                        'burstIsi', burstIsi,...
                        'burstIsi1st', burstIsi1st);
    
    % extract fields from trials structure
    if(onlyRun == 0)
        tmp = trials;
    else
        tmp = trialsExt;
    end
    spikes = getRecField(tmp,'spikes',1:length(lapList));
    dist = getRecField(tmp,'spikesMM',1:length(lapList));
    if(methodTheta == 0)
        theta = getRecField(tmp,'spikesThetaHil',1:length(lapList));
    else
        theta = getRecField(tmp,'spikesThetaLin',1:length(lapList));
    end
    
    disp('Calculate bursts for all the good laps');
    %%%%%%%%%%% calculate burst per neuron per trial
    disp('Calculate bursts per neuron per trial')
    burstPerNeuPerTr = ...
        BurstPerNeuronPerTrial(spikes,theta,dist,beh.indGoodLap,paramBurst);     
                                    % struct including the start
                                    % time of each burst, each cell of a
                                    % field containing the bursts from one
                                    % particular neuron and one particular trial
        
    %%%%%%%%%%% calculate bursts per neuron
    disp('Calculate bursts per neuron')
    burstIsiPerNeuron = BurstPerNeuron(burstPerNeuPerTr,beh.indGoodLap,...
                                       rec.numNeurons);  

    %%%%%%%%%% calculate the burstiness of each neuron
    paramBurst2 = struct('minFractBurst', 0.05,... % min percentage of burst spikes for bursty neurons
                        'maxFractNonBurst',0.05,... % max percentage of burst spikes for non-bursty neurons
                        'pOmnibusPhPref', 0.001,... % when pOmnibus < pOmnibusPhPref, the neuron is considered as with phase preference
                        'pOmnibusNoPhPref',0.001,... % when pOmnibus > pOmnibusNoPhPref, the neuron is considered as without phase preference
                        'minNumSp',100); % min number of spikes
    
    % calculate the burstiness of each neuron
    disp('Calculate burstiness of each neuron')
    burstiness = BurstinessPerNeuron(burstIsiPerNeuron,spikeThetaPhaseStruct,...
                                     paramBurst2);
    
    burstPerNeuPerTrSess = [];
    burstIsiPerNeuronSess = [];
    burstinessSess = [];
    if(length(mazeSess)>1)
        disp('Calculate mean firing rate for each session');
        burstPerNeuPerTrSess = cell(length(mazeSess),1);
        burstIsiPerNeuronSess = cell(length(mazeSess),1);
        burstinessSess = cell(length(mazeSess),1);
        for i = 1:length(mazeSess) 
            fprintf('\nSession %d\n',i);
            indLaps = find(beh.mazeSess == mazeSess(i));
            indLaps = intersect(indLaps,beh.indGoodLap);
            disp('Calculate bursts per neuron per trial')
            burstPerNeuPerTrSess{i} = ...
                BurstPerNeuronPerTrial(spikes,theta,dist,indLaps,paramBurst);   
            disp('Calculate bursts per neuron')
            burstIsiPerNeuronSess{i} = BurstPerNeuron(burstPerNeuPerTrSess{i},...
                                       indLaps,rec.numNeurons);
            disp('Calculate burstiness of each neuron')
            paramBurst2.minNumSp = paramBurst2.minNumSp/length(beh.indGoodLap)...
                                    *length(indLaps);
            burstinessSess{i} = BurstinessPerNeuron(burstIsiPerNeuronSess{i},...
                                     spikeThetaPhaseStruct,paramBurst2); 
        end
    end
    %%%%%%%%% save data
    fullPath = [path fileNameBurst];
    save(fullPath, 'burstPerNeuPerTr','burstIsiPerNeuron','burstiness',...
         'burstPerNeuPerTrSess','burstIsiPerNeuronSess','burstinessSess',...
         'paramBurst','paramBurst2');
   
    %%%%%%%%% draw figure
    if(figureState ~= 0)
        set(0,'Units','pixels')
    end
    
    % draw the burst time for individual neurons
    if(figureState == 2)
        exc = autoCorr.isPyrneuron;
        % relation between number of spikes and burstiness of the neuron
        plotLine(burstIsiPerNeuron.numSp(exc),burstIsiPerNeuron.fractBurstMean(exc),...
            'Num spikes vs Burstiness Pyr','Num spikes','Mean fraction of burst');
        
        inter = autoCorr.isInterneuron == 1;
        % relation between number of spikes and burstiness of the neuron
        plotLine(burstIsiPerNeuron.numSp(inter),burstIsiPerNeuron.fractBurstMean(inter),...
            'Num spikes vs Burstiness Int','Num spikes','Mean fraction of burst');
        
        numSess = length(mazeSess);
        ind = exc & burstIsiPerNeuron.pRayleigh < 0.05 & ...
                    burstIsiPerNeuron.pOmnibus < 0.05;
        burst = zeros(numSess,sum(ind));
        burstInt = zeros(numSess,sum(inter));
        meanBurstDiff = zeros(1,numSess-1);
        meanBurstIntDiff = zeros(1,numSess-1);
        ind = exc & burstIsiPerNeuron.pRayleigh < 0.05 & ...
                    burstIsiPerNeuron.pOmnibus < 0.05;
        for i = 1:numSess
            burst(i,:) = burstIsiPerNeuronSess{i}.fractBurstMean(ind);
            burstInt(i,:) = burstIsiPerNeuronSess{i}.fractBurstMean(inter);
            if(i >= 2)
                meanBurstDiff(i-1) = mean(burst(1,:) - burst(i,:));
                stdBurstDiff(i-1) = std(burst(1,:) - burst(i,:)); 
                meanBurstIntDiff(i-1) = mean(burstInt(1,:) - burstInt(i,:));
                stdBurstIntDiff(i-1) = std(burstInt(1,:) - burstInt(i,:)); 
            end
        end
        barPlot(1:length(mazeSess)-1,meanBurstDiff,stdBurstDiff,...
                'Session','Change in burstiness, ref to sess 1','Pyr');
            
        barPlot(1:length(mazeSess)-1,meanBurstIntDiff,stdBurstIntDiff,...
                'Session','Change in burstiness, ref to sess 1','Int');
        
    end
end
    
function plotLoc(x,y,xStart,yStart,titleF)
% plot the location of the burst spikes
% x:                x locations
% y:                y locations
% xStart:           x locations of the start point of each burst
% yStart:           y locations of the start point of each burst
% titleF:           title of the figure

    h = plot(x,y,'.');
    set(h,'LineWidth',2.0,'Color',[0.502 0.502 0.502]);
    hold on 
    h = plot(xStart,yStart,'ro');
    set(h,'LineWidth',2.0)
    set(gca,'FontSize',14.0,'Box','on','XLim',[0 1200],'YLim',[0 1200]);
    xlabel('X')
    ylabel('Y')
    title(titleF);
    hold off
end

function plotPhPref(theta0,rho0,theta1,rho1,titleFig)
% polar plot for neurron phase preferences

    [figNew,pos] = CreateFig();
    set(figure(figNew),'OuterPosition',pos,'Name',titleFig);
    polar(0,1,'.k');
    hold on 
    h = polar(theta0,rho0,'go');
    set(h,'LineWidth',2.0);
    if(~isempty(theta1))
        h = polar(theta1,rho1,'ro');
        set(h,'LineWidth',2.0);
    end
    set(gca,'FontSize',14.0)
end        
    
function plotHist(handle,x1,x2,x3,bins,xLabel,yLabels,titleFig,xLim,fMean)
% handle:           figure handle
% x1,x2,x3:         plot the histogram of either x1 or x1:x3
% bins:             specified bin centers
% xLabel:           x label of the figure
% yLabels:          y labels of each subfigure
% titleFig:         figure title
% xLim:             limits of x-axis
% fMean:            0: do not calculate and display the mean value
%                   1: calculate and display the mean value

    figure(handle);
    newplot;
    clf;
    hist(x1,bins);
    set(gca,'FontSize',14.0,'Box','on','XLim',xLim);
    xlabel(xLabel);
    ylabel(yLabels(1));

    if(fMean == 1)
        neuCount = hist(x1,bins);
        peakCount = max(neuCount);
        meanX1 = mean(x1);
        hold on
        h = plot([meanX1,meanX1],[peakCount*1.15,peakCount*1.05]);  % plot the mean phase
        set(h,'Color',[0 0 0],'LineWidth',3.0);
        h = plot(meanX1,peakCount*1.05,'v');
        set(h,'Color',[0 0 0],'MarkerFaceColor',[0 0 0]);
    end

    if(~isempty(x2) & ~isempty(x3)) 
        set(gca,'Position',[0.15,0.1,0.8,0.28]);

        subplot(3,1,1)
        hist(x2,bins);
        set(gca,'FontSize',14.0,'Box','on','XLim', xLim,'XTick',[],'Position',[0.15,0.38,0.8,0.28]);
        ylabel(yLabels(2));

        subplot(3,1,1)
        hist(x3,bins);
        set(gca,'FontSize',14.0,'Box','on','XLim', xLim,'XTick',[],'Position',[0.15,0.66,0.8,0.28]);
        ylabel(yLabels(3));
    end  

    title(titleFig);
end