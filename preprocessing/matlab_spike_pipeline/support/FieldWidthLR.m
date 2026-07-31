function FieldWidthLR(path,fileName,spaceBin,methodTheta,figureState,paramState,onlyRun)
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
% FieldWidthLR('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1',1,1,2,1,1)
    
    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin == 2
        methodTheta = 1;
        figureState = 0;
        spaceBin = 20; % mm
        paramState = 0;
        onlyRun = 1;
    elseif nargin == 3
        methodTheta = 1;
        figureState = 0;
        paramState = 0;
        onlyRun = 1;
    elseif nargin == 4
        methodTheta = 1;
        paramState = 0;
        onlyRun = 1;
    elseif nargin == 5
        paramState = 0;
        onlyRun = 1;
    elseif nargin == 6
        onlyRun = 1;
    elseif nargin > 7
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
        fileNameFW = [fileName '_FieldWidthLR_' num2str(spaceBin) ...
                        'mm_H_Run' num2str(onlyRun) '.mat'];            
        fileNameThetaPhase = [fileName '_ThetaPhaseH_Run' ...
                        num2str(onlyRun) '.mat'];
    else
        fileNameFW = [fileName '_FieldWidthLR_' num2str(spaceBin)...
                        'mm_L_Run' num2str(onlyRun) '.mat'];
        fileNameThetaPhase = [fileName '_ThetaPhaseL_Run' ...
                        num2str(onlyRun) '.mat'];
    end
    fileNameConv = [fileName '_convSpikesDist' num2str(spaceBin)...
                        'mm_Run' num2str(onlyRun) '.mat'];
    fileNamePeakFR = [fileName '_PeakFR' num2str(spaceBin)...
                        'mm_Run' num2str(onlyRun) '.mat'];
    fileNameFR = [fileName '_FR_Run' num2str(onlyRun) '.mat'];
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
                'function "MeanFiringRateVR" first.']);
        return;
    end
    load(fullPath,'mFRStruct');
    
    fullPath = [path fileNameConv];
    if(exist(fullPath) == 0)
        disp(['The firing profile file does not exist. Please call ',...
                'function "ConvSpikeTrainDistParVR" first.']);
        return;
    end
    load(fullPath,'filteredSpikeArrayNormT','paramC');
    
    fullPath = [path fileNamePeakFR];
    if(exist(fullPath) == 0)
        disp(['The peak firing rate file does not exist. Please call ',...
                'function "PeakFiringRateVR" first.']);
        return;
    end
    load(fullPath);
         
    fullPath = [path fileNameThetaPhase];
    if(exist(fullPath) == 0)
        disp(['The theta phase file does not exist. Please call function '...
                '"ThetaPhaseP2PVR" first.']);
        return;
    end
    load(fullPath);
    
    fileNameInfo = [fileName(1:end-4) '_Info.mat'];
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath,'file') == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
    numSamples = zeros(1,length(unique(beh.trackLen)));
    for i = 1:length(numSamples)
        numSamples(i) = length(paramC.spaceSteps{i});
    end
    [maxNumSamples,ind] = max(numSamples);
    spaceSteps = paramC.spaceSteps{ind};
    
    %%%%%%%%% parameters
    flag = 0;
    if(paramState == 1)
        fullPath = [path fileNameFW(1:end-4) '_param.mat'];
        if(exist(fullPath,'file') == 2)
            load(fullPath);
            flag = 1;
        end
    end
    
    if(flag == 0)
        paramClass = struct(...
            'minInstFR', 0.3,... % the minimum instantaneous firing rate 
            'pToMFRRatioMinLowFR', 2.5,... % 2.8 the peak to mean instantaneous firing rate ratio separating neurons with a neuron with constant firing and with a peak for neurons with low inst FR
            'pToMFRRatioMin1LowFR', 3,... % 3 the peak to mean instantaneous firing rate ratio separating neurons with a peak and neurons that could have a field for neurons with low inst FR
            'pToMFRRatioMin2LowFR', 4,... % the peak to mean instantaneous firing rate ratio above which neuron is highly likely to have a field for neurons with low inst FR
            'lowerBoundFRFieldNeuronLowFR', 0.1,... % the criteria to distinguish a neuron with a field from a neuron with constant firing rate and with a initial bump.
                                           ...% that is if the neuron inst firing rate
                                           ...% returns back to mean inst FR * lowerBoundFRFieldNeuron for neurons with low inst FR
            'highInstFR', 2,... %1.7 lower the peak to mean instantaneous firing rate ratio for neurons whose firing rate is larger than highInstFR
            'pToMFRRatioMinHighFR', 2,... % the peak to mean instantaneous firing rate ratio separating neurons with a neuron with constant firing and with a peak for neurons with high inst FR
            'pToMFRRatioMin1HighFR', 2.5,... % 2.8 the peak to mean instantaneous firing rate ratio separating neurons with a peak and neurons that could have a field for neurons with high inst FR
            'pToMFRRatioMin2HighFR', 2.7,... % 4 the peak to mean instantaneous firing rate ratio above which neuron is highly likely to have a field for neurons with high inst FR
            'lowerBoundFRFieldNeuronHighFR', 0.3); % the criteria to distinguish a neuron with a field from a neuron with constant firing rate and with a initial bump.
                                           % that is if the neuron inst firing rate
                                           % returns back to mean inst FR * lowerBoundFRFieldNeuron for neurons with high inst FR

        paramFW = struct(...
                'highInstFR', 2,... %1.7 lower the peak to mean instantaneous firing rate ratio for neurons whose firing rate is larger than highInstFR
                'pToMFRRatioFieldMinLowFR', 3.5,... %3 the min peak to mean instantaneous firing rate ratio to detect fields in the neurons which should have fields with low inst FR
                'pToMFRRatioFieldMinHighFR', 2.5,... %2.2 the min peak to mean instantaneous firing rate ratio to detect fields in the neurons which should have fields with high inst FR
                'minNumSpikesPerFieldPerTr', 3.0,... %1.7; % for each trial, the min number of spikes a field should contain
                'minNumSpikesPerFieldPerTrHighFR', 4.5,... %1.7; % for each trial, the min number of spikes a field should contain if its inst FR is larger than highInstFR
                ... %'minSpikeDenPerField', 60,... %40; % the min number of spikes per second in a field
                'minNumSpikesPerFieldPerHz', 0.5,... %0.4 % the min number of spikes per Hz of the peak firing rate
                'minSpikeDenPerFieldPerHzPerTr', 0.006,... %0.005 % the min number of spikes per second per Hz
                'minSpikeDenPerFieldPerHzPerTrHighFR', 0.01,... % 1.0 % the min number of spikes per second per Hz if its inst FR is larger than highInstFR
                'overlapFieldWidthMin', 19*meanSpeed,... % the min width of the field which might contains two overlapping fields
                'threshPeakRatio', 0.6,... % when the higher peak is located in the second half of the field, this threshold determines the ratio between the first and
                                       ... % second peak, larger than which the valley between peaks is used as the boundary of two fields, otherwise, the peak of the 
                                       ... % first field is used as the boundary
                'percSpikesInFields', 0.40,... % 0.55 percentage of spikes which fall inside the fields
                'percSpikesInFieldsHigh',0.40,... ... % percentage of spikes which fall inside the fields for high firing neurons
                'percSpikesInFieldsPerUnit',0.00085,... % percentage of spikes falling inside the fields per unit
                'percSpikesInFieldsPerUnitHigh',0.00065,... % percentage of spikes falling inside the fields per unit for high firing neurons
                'maxPercSpikesInFields', 0.6,... % as long as the percentage of spikes exceed this value, ignore the value calculated using percSpikesInFieldsPerUnit 
                'maxPercSpikesInFieldsHigh', 0.7,... % as long as the percentage of spikes exceed this value, ignore the value calculated using percSpikesInFieldsPerUnitHigh
                'minPeakDist', 1*meanSpeed,... % the min distance between the peaks of two fields, if there is a double field
                'highThreFieldMeanInstFR', 0.9,... % when the mean inst firing rate increases above highThreFieldMeanInstFR * mean inst FR, it is considered to be within a field 
                'lowThreFieldMeanInstFR', 0.10,... % 0.45 when the mean inst firing rate increases above lowThreFieldMeanInstFR * mean inst FR, it is considered to be within a field (double thresholds) 
                'lowerBoundFRFieldNeuron', 0.25,... % 0.25 the criteria to make sure that at least one side of the field would return to 0
                'ReboundCheckRegion', 450,... % there is no clear peak with ReboundCheckRegion (mm) region before or after the field
                'ReboundHeight',0.4,... % rebound height should be smaller than ReboundHeight*peakFieldFR
                'FieldMaxLen',14*meanSpeed,... % max field length
                'MaxFringeLen',7*meanSpeed); % maximum length of the fringe part of the field (counting from inst firing rate hitting the first local minimum below mean firing rate to
                                   % it hitting the lowThreFieldMeanInstFR
        
        paramGen = struct(...
            'numSamples', maxNumSamples,...
            'spaceMergeBin', spaceMergeBin,...
            'spaceSteps',spaceSteps,...
            'percActiveTrials', 0.5);
    end
    
    %%%%%%%%% prepare figure
    if(figureState ~= 0)
        set(0,'Units','pixels') 
    end
    
    %%%%%%%%% classify neurons into different catagories according to their
    %%%%%%%%% peak/mean inst ratio
%     neuronClassStruct = ...
%         NeuronClass(pFRStructNormT,1:rec.numNeurons,paramClass);
%             
%     fullPath = [path fileNameFW];
%     save(fullPath, 'neuronClassStruct','-v7.3'); 
%     
%     %%%%%%%%% extract episode fields from neurons with potential fields
%     disp('Identify fields (all the trials)');
%     fieldStruct = Fields(pFRStructNormT,spikeThetaPhaseStruct,...
%         neuronClassStruct.neuronPotentialField,paramFW,paramGen,figureState);
%                 % extract episode fields for exc neurons with potential fields
                
    disp('Field identification for each session');
    fieldStructSess = cell(1,length(mazeSess));
    neuronClassStructSess = cell(1,length(mazeSess));
    if(length(mazeSess) > 1)
        for i = 1:length(mazeSess)
            disp(['Session' num2str(i)]);
            neuronClassStructSess{i} = ...
                NeuronClass(pFRStructNormTSess{i},1:rec.numNeurons,paramClass);
    
            fieldStructSess{i} = ...
                Fields(pFRStructNormTSess{i},spikeThetaPhaseStructSess{i},...
                neuronClassStructSess{i}.neuronPotentialField,...
                paramFW,paramGen,figureState);
                % extract episode fields for exc neurons with potential fields
        end
    else
        disp(['Session' num2str(i)]);
            neuronClassStructSess{1} = ...
                NeuronClass(pFRStructNormT,1:rec.numNeurons,paramClass);
    
            fieldStructSess{1} = ...
                Fields(pFRStructNormT,spikeThetaPhaseStruct,...
                neuronClassStructSess{1}.neuronPotentialField,...
                paramFW,paramGen,figureState);
    end
    fullPath = [path fileNameFW];
    save(fullPath, 'neuronClassStructSess','fieldStructSess'); 
%     save(fullPath, 'fieldStruct','neuronClassStructSess','fieldStructSess','-append'); 
    
    if(figureState == 2)
       %%% all the trials
        count = 0;
%         indNeurons = [3 21 48 53 54 80 95 126 144 165 206 216 226];
                %A011-20190218-01
%         indNeurons = [22 30 49 58 62 66 68 69];
                %A009-20190111-02 
%         for i = indNeurons
        for i = 1:rec.numNeurons  
            if(mFRStruct.mFR(i) > 0.05) %minFR
                % left trials
%                 neuronClass = getNeuronClassIndNeuron(i,neuronClassStruct); 
                        % get the class information  
                strNumField = '';
                for j = 1:length(mazeSess)
                    fieldInfoTmp{j} = getFieldInfoIndNeuron(i,fieldStructSess{j}); 
                    if(~isempty(fieldInfoTmp{j}))
                        strNumField = [strNumField, ' ', num2str(size(fieldInfoTmp,1))];
                    else
                        strNumField = [strNumField, ' 0'];
                    end
                end
                        % get the field information
                count = count + 1;

                if(mod(count-1,16) == 0)
                    [figNew,pos] = CreateFig();
                    set(0,'Units','pixels') 
                    figTitle = 'All the trials';
                    set(figure(figNew),'OuterPosition',...
                        [pos(1) pos(2) pos(3)*2 pos(4)*2.2],'Name',figTitle)
                end

                subplot(4,4,mod(count-1,16)+1)

                figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ') NumField' strNumField];                

%                 if(isempty(spikeThetaPhaseStructSess{1}))
%                     [thetaDist,thetaPhase] = extractThetaPhasePerNeuron(...
%                         spikeThetaPhaseStruct,i,spaceMergeBin);
%                 else
%                     [thetaDist,thetaPhase] = extractThetaPhasePerNeuron(...
%                         spikeThetaPhaseStructSess,i,spaceMergeBin);
%                 end
                
%                 indGoodLap = find(beh.pauseWithinTrial < 1 & ...
%                     beh.pauseBefTrial > 1);
%                 indGoodLap1 = intersect(indGoodLap, beh.indGoodLap);
%                 goodLapSess = beh.mazeSess(indGoodLap1);
%                 indSessBorder = find(diff(goodLapSess) == 1);  
% 
%                 if(isempty(spikeThetaPhaseStructSess{1}))
%                     [thetaDist,thetaPhase] = ...
%                             extractThetaPhasePerNeuronOverTr(...
%                             spikeThetaPhaseStruct,i,spaceMergeBin,...
%                             indGoodLap1,goodLapSess);
%                 else
%                     [thetaDist,thetaPhase] = ...
%                             extractThetaPhasePerNeuronOverTr(...
%                             spikeThetaPhaseStructSess,i,spaceMergeBin,...
%                             indGoodLap1,goodLapSess);
%                 end
                
%                 plotThetaPhaseSessions(gca,thetaDist,thetaPhase,...
%                     max(beh.trackLen)/spaceMergeBin,figTitle);
%                     
%                 indGoodLap1 = ...
%                     [3 4 10 11 12 13 14 16 17 18 20 21 22 23 24 25 27 28 29 ...
%                     45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 ...
%                     64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 ...
%                     85 93 94 97 98 121 122 123 124 125 126 127 128 129 130 ...
%                     131 132 133 134 135 142 143 145 146 150 151 152 153 154 ...
%                     155 156 157 158 159 160 161 162 163 164 165 166 167 168 ...
%                     169 170 171 172 173 174 175 176 177 178 179 180 181 182 ...
%                     183 184 185 186 189 190 195 196 198 199]; %A011-20190215-01
%                 indGoodLap10 = ...
%                     [3 4 5 15 16 18 19 21 22 23 24 26 28 32 33 34 48 49 51 56 ...
%                     57 ];%A011-20190218-01
%                 indGoodLap1 = ...
%                       [82 84 85 86 87 88 89 90 91 93 94 96 97 98 99 100 101 ...
%                     102 104 105 106 107 108 109 110 111 112 117 118 119 120 ...
%                     121 122 123 124 128 129 130 131 132 133 134 135 136 137 ...
%                     138 139 140 141 142 143 144 145 146 147 148 149 150 153 ...
%                     154 155 156 157 158 160 161];%A011-20190218-01
%                 indGoodLap1 = setdiff(beh.indCorrLap,[indGoodLap1 indGoodLap10]);
                            %A011-20190218-01
%                 indGoodLap1 = 1:100; %beh.indGoodLap; %A011-20190219-01
                indGoodLap1 = beh.indGoodLap; %A009-20190219-01
                
%                 indGoodLap1 = [find(beh.indStimLap(beh.indGoodLap) ~= 1 ),...
%                     find(beh.indStimLap(beh.indGoodLap) == 1)];
%                 indGoodLap1 = beh.indGoodLap(indGoodLap1);
%                 indGoodLap1 = 1:60; %A012-20190225-01
                goodLapSess = beh.mazeSess(indGoodLap1);
                
                indSessBorder = find(diff(goodLapSess) > 0);  
                plotFRProfIndNeuronIndTrial(gca,filteredSpikeArrayNormT,...
                    numSamples(1),i,indGoodLap1,figTitle,indSessBorder);
%                 plotFRProfIndNeuronIndTrial(gca,filteredSpikeArrayNormT,...
%                     1801,i,indGoodLap1,figTitle,indSessBorder);
                
%                   goodLapSess = beh.mazeSess(beh.indGoodLap);
%                   indSessBorder = find(diff(goodLapSess) == 1); 
%                   plotFRProfIndNeuronIndTrial(gca,filteredSpikeArrayNormT,i,...
%                     beh.indGoodLap,figTitle,indSessBorder);
                
%                 plotFRProfIndNeuronIndTrial(gca,filteredSpikeArrayNormT,i,...
%                     1:50,figTitle);
                                
    %                     ax = axes('Position',get(gca,'Position'),...
    %                        'XAxisLocation','top',...
    %                        'YAxisLocation','right',...
    %                        'XTickLabel',[],...
    %                        'Xlim',[0 numSamples*timeStep],...
    %                        'Color','none',...
    %                        'XColor','k','YColor','k',...
    %                        'FontSize',14);     
    %                     figTitle = ['Neuron ' num2str(i) ' ' neuronClass];
    %                     plotFRProfileIndNeuron(ax,pFRStructLe.FRArr(i,:),...
    %                       pFRStructLe.meanInstFRArr(i),fieldInfoTmp,figTitle,...
    %                       numSamples,timeStep);
            else
                disp(['Firng rate of neuron ' num2str(i) ' is too low: ' ...
                    num2str(mFRStruct.mFR(i)) ' Hz']);
                continue;
            end
        end
    end
        
    clear mydata
    
    %%%%%%%%% save data to file
     
    fullPath = [path fileNameFW(1:end-4) '_param.mat'];
    save(fullPath,'paramGen','paramClass','paramFW');
        
    clear mydata;
    
    %%%%%%%%% draw figure 
%     if(figureState == 1)
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% all the trials inh
%         plotFWData(fieldInhStructLe, 'L', 'inh', timeStep, '');
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% all the trials exc
%         plotFWData(fieldExcStructLe, 'L', 'exc', timeStep, '');
%             
end

function plotFRProfIndNeuronIndTrial...
            (handle,filteredSpikeArrayNormT,numSamples,neuronNo,indLaps,figTitle,indSessBorder)
    numTr = length(indLaps);
    FRProfilePerTrial = zeros(numTr,numSamples);
    for i = 1:numTr
        if(isempty(filteredSpikeArrayNormT{i}))
            continue;
        end
%         disp(['Tr no. ' num2str(i)]);
        szArr = size(filteredSpikeArrayNormT{i},2);
        if(szArr > numSamples)
            szArr = numSamples;
        end
%         szArr = size(filteredSpikeArrayNormT{indLaps(i)},2);
        FRProfilePerTrial(i,1:szArr) = filteredSpikeArrayNormT{indLaps(i)}(neuronNo,1:numSamples); 
            %./max(filteredSpikeArrayNormT{indLaps(i)}(neuronNo,:));
    end
%     FRProfilePerTrial = FRProfilePerTrial/max(FRProfilePerTrial(:));
    h = imagesc(0:numSamples-1,1:numTr,FRProfilePerTrial);
    if(~isempty(indSessBorder))
        hold on
        for i = 1:length(indSessBorder)
            h = plot([0 numSamples],indSessBorder(i)*ones(1,2),'r');
            set(h,'LineWidth',1);
        end
    end    
            
    set(gca,'FontSize',8.0,'Box','on','XLim',[0 numSamples-1],'YLim',[0 numTr]);
    xlabel('Dist (mm)');
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
