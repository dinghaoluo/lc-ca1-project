function BurstAllAlignedRunCtrlOnly(path,fileName,methodTheta,figureState,onlyRun,mazeSess)
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
        mazeSess = 1;
    elseif nargin == 3
        figureState = 1;
        onlyRun = 1;
        mazeSess = 1;
    elseif nargin == 4
        onlyRun = 1;
        mazeSess = 1;
    elseif nargin == 5
        mazeSess = 1;
    elseif nargin > 6
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
    fileNameBurst = [fileName '_burstAllAlignedRun_TH' th '_msess' num2str(mazeSess) '_Run' num2str(onlyRun) ...
                     '.mat'];

    if(methodTheta == 0)
        fileNameThetaPhase = [fileName '_ThetaPhaseHAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    else
        fileNameThetaPhase = [fileName '_ThetaPhaseLAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    end
    fileNameOrig = [fileName '.mat'];
    
    fullPath = [path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT file does not exist');
        return;
    end
    load(fullPath,'trialsRunSpikes');
    
    fullPath = [path fileNameThetaPhase];
    if(exist(fullPath) == 0)
        disp(['The theta phase file does not exist. Please call function' ...
              ' "ThetaPhaseP2PAlignRun" first.']);
        return;
    end
    load(fullPath,'spikeThetaPhaseRunNoStimGood','spikeThetaPhaseRunNoStimBad');
    
    fullPath = [path fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','trialNoStimCtrl','pulseMeth');
    
    fullPath = [path fileName '_Info.mat']; 
    if(exist(fullPath) == 0)
        disp('The behavioral information file does not exist');
        return;
    end
    load(fullPath,'rec');
   
    % burst parameters
    paramBurst = struct('numNeurons', rec.numNeurons,...
                        'burstIsi', burstIsi,...
                        'burstIsi1st', burstIsi1st);
    
    % extract fields from trials structure
    
    spikes = trialsRunSpikes.Time;
    dist = trialsRunSpikes.Dist;
    if(methodTheta == 0)
        theta = trialsRunSpikes.thPhaseHilbSpike;
    else
        theta = trialsRunSpikes.thPhaseInterpSpike;
    end
   
    %%
    disp('Calculate bursts for all the non-stimulated good laps');
    %%%%%%%%%%% calculate burst per neuron per trial
    disp('Calculate bursts per neuron per trial')
    burstPerNeuPerTrNonStimGood = ...
        BurstPerNeuronPerTrialAligned(spikes,theta,dist,trialNoNonStimGood,paramBurst);     
                                    % struct including the start
                                    % time of each burst, each cell of a
                                    % field containing the bursts from one
                                    % particular neuron and one particular trial
        
    %%%%%%%%%%% calculate bursts per neuron
    disp('Calculate bursts per neuron')
    burstIsiPerNeuronNonStimGood = BurstPerNeuron(burstPerNeuPerTrNonStimGood,trialNoNonStimGood,...
                                       rec.numNeurons);  
                                   
    %%
    disp('Calculate bursts for all the non-stimulated bad laps');
    %%%%%%%%%%% calculate burst per neuron per trial
    disp('Calculate bursts per neuron per trial')
    burstPerNeuPerTrNonStimBad = ...
        BurstPerNeuronPerTrialAligned(spikes,theta,dist,trialNoNonStimBad,paramBurst);     
        
    %%%%%%%%%%% calculate bursts per neuron
    disp('Calculate bursts per neuron')
    burstIsiPerNeuronNonStimBad = BurstPerNeuron(burstPerNeuPerTrNonStimBad,trialNoNonStimBad,...
                                       rec.numNeurons);  
                                   
    %% added by Yingxue on 2/15/2021
    burstPerNeuPerTrStim = [];
    burstIsiPerNeuronStim = [];
    burstPerNeuPerTrStimCtrl = [];
    burstIsiPerNeuronStimCtrl = [];
    for i = 1:length(pulseMeth)
        disp('Calculate bursts for all stimulated laps');
        %%%%%%%%%%% calculate burst per neuron per trial
        disp('Calculate bursts per neuron per trial')
        burstPerNeuPerTrStim{i} = ...
            BurstPerNeuronPerTrialAligned(spikes,theta,dist,trialNoStim{i},paramBurst);     

        %%%%%%%%%%% calculate bursts per neuron
        disp('Calculate bursts per neuron')
        burstIsiPerNeuronStim{i} = BurstPerNeuron(burstPerNeuPerTrStim{i},trialNoStim{i},...
                                           rec.numNeurons);  
                                       
        disp('Calculate bursts for all control trials during stimulation');
        %%%%%%%%%%% calculate burst per neuron per trial
        disp('Calculate bursts per neuron per trial')
        burstPerNeuPerTrStimCtrl{i} = ...
            BurstPerNeuronPerTrialAligned(spikes,theta,dist,trialNoStimCtrl{i},paramBurst);     

        %%%%%%%%%%% calculate bursts per neuron
        disp('Calculate bursts per neuron')
        burstIsiPerNeuronStimCtrl{i} = BurstPerNeuron(burstPerNeuPerTrStimCtrl{i},trialNoStimCtrl{i},...
                                           rec.numNeurons);  
    end

    %%%%%%%%% save data
    fullPath = [path fileNameBurst];
    save(fullPath, 'burstPerNeuPerTrNonStimGood','burstIsiPerNeuronNonStimGood',...
         'burstPerNeuPerTrNonStimBad','burstIsiPerNeuronNonStimBad',...
         'burstPerNeuPerTrStim','burstIsiPerNeuronStim',...
         'burstPerNeuPerTrStimCtrl','burstIsiPerNeuronStimCtrl',...
         'paramBurst');

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