function FieldWidthLRAligned(path,fileName,methodTheta,figureState,onlyRun,mazeSess)
% Calculate the field width of individual neurons
% Dependence: function "PeakFiringRate", "MeanFiringRate" and "ThetaPhaseLR" 
% should be executed first function
% (This function classifies the neuron into 4 different classes depending
% on their firing patterns
%  1. neurons with constant firing rate
%  2. neurons with contsant firing rate and with a initial peak
%  3. neurons with fields (single field and double field))
% path:         the path of the recording file
% fileName:     name of the recording file
% spaceBin:     2SD of the Gaussian filter used to obtain the firing rate 
%               profile
% methodTheta:  0: hilbert transform
%               1: linear interpolation
% figureState:  0: figure off
%               1: plot the histogram of field width
%               2: plot the firing rate profile during the analysis before
%               finally plotting the histogram of field width
% paramState:   0: use the param defined in the function
%               1: use the saved param
% onlyRun:      1: only consider the time period when the animal is running 
%
% Example:
% FieldWidthLR('./','A022-20191107-01_DataStructure_mazeSection1_TrialType1',1,2,1,1)
    
    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin == 2
        methodTheta = 1;
        figureState = 0;
        onlyRun = 1;
    elseif nargin == 3
        figureState = 0;
        onlyRun = 1;
    elseif nargin == 4
        onlyRun = 1;
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
        fileNameFW = [fileName '_FieldWidthAligned_' ...
                        'mSess' num2str(mazeSess) '_H_Run' num2str(onlyRun) '.mat'];            
        fileNameThetaPhase = [fileName '_ThetaPhaseHAligned_msess' num2str(mazeSess)...
                        '_Run' num2str(onlyRun) '.mat'];
    else
        fileNameFW = [fileName '_FieldWidthAligned_'...
                        'mSess' num2str(mazeSess) '_L_Run' num2str(onlyRun) '.mat'];   
        fileNameThetaPhase = [fileName '_ThetaPhaseLAligned_msess' num2str(mazeSess)...
                        '_Run' num2str(onlyRun) '.mat'];
    end
    fileNameConv = [fileName '_convSpikesAligned_msess' num2str(mazeSess) ...
                        '_Run' num2str(onlyRun) '.mat'];
    fileNamePeakFR = [fileName '_PeakFRAligned_msess' num2str(mazeSess) ...
                        '_Run' num2str(onlyRun) '.mat'];
    fileNameFR = [fileName '_FR_Run' num2str(onlyRun) '.mat'];
    fileNameSpike = [fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    fileName = [fileName '.mat'];
     
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('The recording file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath,'cluList');
    
    fullPath = [path fileNameFR];
    if(exist(fullPath) == 0)
        disp(['The mean firing rate file does not exist. Please call ',...
                'function "MeanFiringRate" first.']);
        return;
    end
    load(fullPath,'mFRStructSess','mFRStruct');
    
    fullPath = [path fileNameConv];
    if(exist(fullPath) == 0)
        disp(['The firing profile file does not exist. Please call ',...
                'function "ConvSpikeTrain_Aligned" first.']);
        return;
    end
    load(fullPath,'filteredSpikeArrayRun','paramC','timeStep'); 
    numSamples = size(filteredSpikeArrayRun{2},2);
    
    fullPath = [path fileNamePeakFR];
    if(exist(fullPath) == 0)
        disp(['The peak firing rate file does not exist. Please call ',...
                'function "PeakFiringRate_Aligned" first.']);
        return;
    end
    load(fullPath,'pFRNonStimGoodStruct','trialNoNonStimGood');
         
    fullPath = [path fileNameThetaPhase];
    if(exist(fullPath) == 0)
        disp(['The theta phase file does not exist. Please call function '...
                '"ThetaPhaseP2PAlignedRun" first.']);
        return;
    end
    load(fullPath,'spikeThetaPhaseRunNoStimGood');
    
    fullPath = [path fileNameSpike]; 
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT file does not exist');
        return;
    end
    load(fullPath,'trialsRunSpikes');
    
    fullPath = [path fileName(1:end-4) '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    fileNameInfo = [fileName(1:end-4) '_Info.mat'];
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath,'file') == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath);
    
    %%%%%%%%% parameters

    paramClass = struct(...
        'minInstFR', 0.13,... % the minimum instantaneous firing rate 
        'pToMFRRatioMinLowFR', 8,... % 10 the peak to mean instantaneous firing rate ratio separating neurons with a neuron with constant firing and with a peak for neurons with low inst FR
        'pToMFRRatioMin1LowFR', 9,... % 12 the peak to mean instantaneous firing rate ratio separating neurons with a peak and neurons that could have a field for neurons with low inst FR
        'pToMFRRatioMin2LowFR', 12,... %14 the peak to mean instantaneous firing rate ratio above which neuron is highly likely to have a field for neurons with low inst FR
        'lowerBoundFRFieldNeuronLowFR', 0.1,... % the criteria to distinguish a neuron with a field from a neuron with constant firing rate and with a initial bump.
                                       ...% that is if the neuron inst firing rate
                                       ...% returns back to mean inst FR * lowerBoundFRFieldNeuron for neurons with low inst FR
        'highInstFR', 0.3,... %1.7 lower the peak to mean instantaneous firing rate ratio for neurons whose firing rate is larger than highInstFR
        'pToMFRRatioMinHighFR', 5,... % the peak to mean instantaneous firing rate ratio separating neurons with a neuron with constant firing and with a peak for neurons with high inst FR
        'pToMFRRatioMin1HighFR', 5.5,... % 2.8 the peak to mean instantaneous firing rate ratio separating neurons with a peak and neurons that could have a field for neurons with high inst FR
        'pToMFRRatioMin2HighFR', 10,... % 4 the peak to mean instantaneous firing rate ratio above which neuron is highly likely to have a field for neurons with high inst FR
        'lowerBoundFRFieldNeuronHighFR', 0.5); % the criteria to distinguish a neuron with a field from a neuron with constant firing rate and with a initial bump.
                                       % that is if the neuron inst firing rate
                                       % returns back to mean inst FR * lowerBoundFRFieldNeuron for neurons with high inst FR

    paramFW = struct(...
            'highInstFR', 0.3,... %1.7 lower the peak to mean instantaneous firing rate ratio for neurons whose firing rate is larger than highInstFR
            'pToMFRRatioFieldMinLowFR', 9,... %3 the min peak to mean instantaneous firing rate ratio to detect fields in the neurons which should have fields with low inst FR
            'pToMFRRatioFieldMinHighFR', 5.5,... %2.2 the min peak to mean instantaneous firing rate ratio to detect fields in the neurons which should have fields with high inst FR
            'minNumSpikesPerFieldPerTr', 1.2,... %1.7; % for each trial, the min number of spikes a field should contain
            'minNumSpikesPerFieldPerTrHighFR', 1.2,... %1.7; % for each trial, the min number of spikes a field should contain if its inst FR is larger than highInstFR
            ... %'minSpikeDenPerField', 60,... %40; % the min number of spikes per second in a field
            'minNumSpikesPerFieldPerHz', 2,... %0.4 % the min number of spikes per Hz of the peak firing rate
            'minSpikeDenPerFieldPerHzPerTr', 0.5,... %0.005 % the min number of spikes per second per Hz
            'minSpikeDenPerFieldPerHzPerTrHighFR', 0.5,... % 1.0 % the min number of spikes per second per Hz if its inst FR is larger than highInstFR
            'multiPeakFieldDetThr',0.5,... % detect whether there are more than one peak within one field, this is the threshold from peak 
            'overlapFieldWidthMin', 8,... % the min width of the field which might contains two overlapping fields
            'threshPeakRatio', 0.6,... % when the higher peak is located in the second half of the field, this threshold determines the ratio between the first and
                                   ... % second peak, larger than which the valley between peaks is used as the boundary of two fields, otherwise, the peak of the 
                                   ... % first field is used as the boundary
            'percSpikesInFieldsLowFR', 0.65,... % 0.65 percentage of spikes which fall inside the fields low FR
            'percSpikesInFieldsHighFR', 0.6,... % 0.55 percentage of spikes which fall inside the fields high FR
            'minPeakDist', 6,... % the min distance between the peaks of two fields, if there is a double field
            'highThreFieldPeakFR', 0.3,... % when the inst firing rate increases above highThreFieldPeakFR * peak FR, it is considered to be within a field 
            'lowThreFieldPeakFR', 0.13,... % 0.45 when the  inst firing rate increases above lowThreFieldMeanInstFR * peak FR, it is considered to be within a field (double thresholds) 
            'lowerBoundFRFieldNeuron', 0.15,... % 0.45 the criteria to make sure that at least one side of the field would return to 0
            'FieldMaxLenLowFR',4,... % max field length low firing rate cell
            'FieldMaxLenHighFR',4.1,... % max field length high firing rate cell
            'MaxFringeLenLowFR',2,...
            'MaxFringeLenHighFR',2,... % maximum length of the fringe part of the field (counting from inst firing rate hitting the first local minimum below mean firing rate to
            ...                   % it hitting the lowThreFieldMeanInstFR
            'MaxPercTrLength',0.8); % the width of field should not exceed certain percent of the mean trial length

    paramGen = struct(...
        'numSamples', length(timeStep),...
        'timeBin', (timeStep(2)-timeStep(1))/sampleFq,...
        'timeSteps', timeStep,... % actually time in here
        'percActiveTrials', 0.5);
    
    %%%%%%%%% prepare figure
    if(figureState ~= 0)
        set(0,'Units','pixels') 
    end
    
    %%%%%%%%% classify neurons into different catagories according to their
    %%%%%%%%% peak/mean inst ratio

    disp('Field identification ');
    fieldStructSess = [];
    neuronClassStructSess = [];

    % adjust paramFW.MaxFringeLen based on the trial length
    meanTrLen = mean(behPar.numSamplesRun(pFRNonStimGoodStruct.indLapList));
    paramFW.FieldMaxLenHighFR = meanTrLen*paramFW.MaxPercTrLength/sampleFq;
    if(paramFW.FieldMaxLenHighFR> 5.2)
        paramFW.FieldMaxLenHighFR = 5.2;
    end
    if(meanTrLen*paramFW.MaxPercTrLength/sampleFq < paramFW.FieldMaxLenLowFR)
        paramFW.FieldMaxLenLowFR = meanTrLen*paramFW.MaxPercTrLength/sampleFq;
    end
    
    % low pass filter ethe averaged firing rate profile before field
    % identification
    [nNeu,fs] = size(pFRNonStimGoodStruct.avgFRProfile);
    fc = 20;
    [b,a] = butter(4, fc/(fs/2));
    avgFRProfileTmp = zeros(nNeu,fs);
    peakFRTmp = zeros(1,nNeu);
    meanInstFRTmp = zeros(1,nNeu);
    for i = 1:nNeu
        if(sum(pFRNonStimGoodStruct.avgFRProfile(i,:)) > 0)
            avgFRProfileTmp(i,:) = filtfilt(b,a,pFRNonStimGoodStruct.avgFRProfile(i,:));
            [peakFRTmp(i),indTmp] = max(avgFRProfileTmp(i,:)); 
            meanInstFRTmp(i) = mean(avgFRProfileTmp(i,:)); 
        end
    end
    pFRNonStimGoodStructOrig = pFRNonStimGoodStruct;
    pFRNonStimGoodStruct.avgFRProfile = avgFRProfileTmp;
    pFRNonStimGoodStruct.peakFR = peakFRTmp;
    pFRNonStimGoodStruct.meanInstFR = meanInstFRTmp;
    
    neuronClassStructSess = ...
        NeuronClass(pFRNonStimGoodStruct,1:rec.numNeurons,paramClass);

    fieldStructSess = ...
        FieldsTime(pFRNonStimGoodStruct,spikeThetaPhaseRunNoStimGood,...
        trialsRunSpikes.Time,neuronClassStructSess.neuronPotentialField,...
        paramFW,paramGen,figureState);

    fullPath = [path fileNameFW(1:end-4) '.mat'];
    save(fullPath, 'neuronClassStructSess','fieldStructSess'); 
    fullPath = [path fileNameFW(1:end-4) '_param.mat'];
    save(fullPath,'paramGen','paramClass','paramFW');
    
    if(figureState == 2)
       %%% all the trials
        count = 0;
        for i = 1:rec.numNeurons 
            if(mazeSess == 0)
                mazeSess = 1;
            end
            if(isempty(mFRStructSess{mazeSess}))
                mFRStructTmp = mFRStruct;
            else
                mFRStructTmp = mFRStructSess{mazeSess};
            end
            if(mFRStructTmp.mFR(i) > minFR) %minFR
                % get the class information  
                strNumField = '';
                fieldInfoTmp = getFieldInfoIndNeuron(i,fieldStructSess); 
                if(~isempty(fieldInfoTmp))
                    strNumField = [strNumField, ' ', num2str(size(fieldInfoTmp,1))];
                else
                    strNumField = [strNumField, ' 0'];
                end
                        % get the field information
                count = count + 1;

                if(mod(count-1,16) == 0)
%                     if(count ~= 1)
%                         fullpath = [path fileNameFW '_' num2str(count-1)];
%                         print('-painters', '-dpdf', fullpath, '-r600')
%                         savefig([fullpaths '.fig']);
%                     end
                    [figNew,pos] = CreateFig();
                    set(0,'Units','pixels') 
                    figTitle = 'All the trials';
                    set(figure(figNew),'OuterPosition',...
                        [pos(1) pos(2)-500 pos(3)*2 pos(4)*2.2],'Name',figTitle)
                end

                subplot(4,4,mod(count-1,16)+1)

                figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ') NumField' strNumField];                
                                              
               plotFRProfIndNeuronIndTrial(gca,filteredSpikeArrayRun,...
                    paramGen.numSamples,paramGen.timeSteps/sampleFq,i,trialNoNonStimGood,figTitle,[]);
            else
                disp(['Firng rate of neuron ' num2str(i) ' is too low: ' ...
                    num2str(mFRStructTmp.mFR(i)) ' Hz']);
                continue;
            end
        end
    end
        
    clear mydata         
end

function plotFRProfIndNeuronIndTrial...
            (handle,filteredSpikeArrayNormT,numSamples,timeSteps,neuronNo,indLaps,figTitle,indSessBorder)
    numTr = length(indLaps);
    FRProfilePerTrial = filteredSpikeArrayNormT{neuronNo}(indLaps,1:numSamples); 
            %./max(filteredSpikeArrayNormT{indLaps(i)}(neuronNo,:));
%     FRProfilePerTrial = FRProfilePerTrial/max(FRProfilePerTrial(:));
    h = imagesc(timeSteps,1:numTr,FRProfilePerTrial);
    if(~isempty(indSessBorder))
        hold on
        for i = 1:length(indSessBorder)
            h = plot([0 numSamples],indSessBorder(i)*ones(1,2),'r');
            set(h,'LineWidth',1);
        end
    end    
            
    set(gca,'FontSize',8.0,'Box','on','XLim',[0 timeSteps(end)/2],'YLim',[0 numTr]);
    xlabel('Time (s)');
    ylabel('Trial No.');
    title(figTitle);
end

function [thetaDist,thetaPhase] = ...
    extractThetaPhasePerNeuron(spikeThetaPhaseStructSess,indNeuron,spaceMergeBin)
    
    numSess = length(spikeThetaPhaseStructSess);
    if(numSess == 1)
        thetaDist{1} = spikeThetaPhaseStructSess.spDistPerNeuron{indNeuron}...
                        /spaceMergeBin;
        thetaPhase{1} = spikeThetaPhaseStructSess.spPhaseVsTPerNeuron{indNeuron};
    else
        for i = 1:numSess
            thetaDist{i} = spikeThetaPhaseStructSess{i}.spDistPerNeuron{indNeuron}...
                            /spaceMergeBin;
            thetaPhase{i} = spikeThetaPhaseStructSess{i}.spPhaseVsTPerNeuron{indNeuron};
        end
    end
end

function [thetaDist,thetaPhase] = ...
    extractThetaPhasePerNeuronOverTr(spikeThetaPhaseStructSess,indNeuron,spaceMergeBin,...
        indLaps,lapSess)
    
    numSess = length(spikeThetaPhaseStructSess);
    if(numSess == 1)
        thetaDist{1} = [];
        thetaPhse{1} = [];
        for i = indLaps
            thetaDist{1} = [thetaDist{1} ...
                spikeThetaPhaseStructSess.spDistPerTrialPerNeuron{indNeuron,i}...
                            /spaceMergeBin];
            thetaPhase{1} = [thetaPhase{1} ...
                spikeThetaPhaseStructSess.spPhaseVsTPerTrialPerNeuron{indNeuron,i}];
        end
    else
        for sess = 1:numSess
            thetaDist{sess} = [];
            thetaPhase{sess} = [];
            indLapsSess = indLaps(lapSess == sess);
            for i = indLapsSess
                ind = find(spikeThetaPhaseStructSess{sess}.indLapList == i);
                thetaDist{sess} = [thetaDist{sess}; ...
                    spikeThetaPhaseStructSess{sess}.spDistPerTrialPerNeuron{indNeuron,ind}...
                                /spaceMergeBin];
                thetaPhase{sess} = [thetaPhase{sess}; ...
                    spikeThetaPhaseStructSess{sess}.spPhaseVsTPerTrialPerNeuron{indNeuron,ind}];
            end
        end
    end
end

function plotThetaPhaseSessions(handle,thetaDist,thetaPhase,distInterval,...
    figTitle,color)
% plot the theta phase of spikes from individual neuron
% handle:           axis handle of the figure
% thetaDist:        distance array of spikes
% thetaPhase:       theta phase array of spikes
% distInterval:     time interval of the plot
% timeStep:         the time step 

    if(nargin == 5)
        color = [0.5 0.5 0.5];
    end
    if(length(distInterval) == 1)
        distInterval = [0 distInterval];
    end
    
    hold on;
      
    numSess = length(thetaDist);
    totalPhase = numSess*(360+30);   % add 20 degree gap between sessions   
    for i = 1:numSess
        phaseAcc = (numSess-i)*(360+30);
        
        [thetaDistDouble,thetaSpikesDouble] = ...
        getMultiCycles(thetaDist{i},thetaPhase{i},1);
        
        thetaSpikesDouble = thetaSpikesDouble + phaseAcc;

        h = plot(thetaDistDouble,thetaSpikesDouble, '.');
        set(h,'LineWidth',2.0,'Color',color);
        h = plot(distInterval, [phaseAcc phaseAcc],'r-');
        set(h,'LineWidth',0.2);
    end
     
    set(gca,'FontSize',8.0,'Box','on','XLim',distInterval,'YLim',...
            [0 totalPhase]);
    xlabel('Dist (mm)');
    ylabel('Theta phase (deg)');
    title(figTitle);
end
