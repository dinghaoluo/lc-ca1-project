function PyrInitPeakAllStimPVRec(methodKMean)
%% only consider stimulation recordings
    
    methodTheta = 1;
    minFR = 0.15;
    maxFR = 7;
    if(nargin == 0)
        methodKMean = 2;
    end
    
    RecordingList;
    pathAnal0 = 'Z:\Yingxue\DataAnalysisRaphi\Pyramidal\';
    pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\PyramidalALPVStim\';
        
    sampleFq = 1250;
    
    load([pathAnal0 'autoCorrPyrAllRec.mat']);
    if(exist([pathAnal 'initPeakPyrAllRecStim.mat']))
        load([pathAnal 'initPeakPyrAllRecStim.mat']);
    end
    
    anmNoInact = [41, 42, 46];
    anmNoAct = [22 23 37 39 44 24 25 40];
    pulseMethod{1} = [2 4]; % inactivation
    pulseMethod{2} = [2 3]; % activation
    
    nNeuWithFieldAligned = [modPyrNoCue.nNeuWithFieldAligned modPyrAL.nNeuWithFieldAligned ...
        modPyrPL.nNeuWithFieldAligned];
    isNeuWithFieldAligned = [modPyrNoCue.isNeuWithFieldAligned modPyrAL.isNeuWithFieldAligned modPyrPL.isNeuWithFieldAligned];
        
%     disp('Active licking - peak firing rate')
%     modPyr1AL = accumPyr3(listRecordingsActiveLickPath,...
%         listRecordingsActiveLickFileName,mazeSessionActiveLick,autoCorrPyrAll,...
%         nNeuWithFieldAligned,isNeuWithFieldAligned,...
%         anmNoInact,anmNoAct,minFR,maxFR,2,sampleFq);
%     
%     save([pathAnal 'initPeakPyrAllRecStim.mat'],'modPyr1AL'); 
    
    actOrInActStr{1} = 'act';
    actOrInActStr{2} = 'inact';
    
    % activation or inactivation
    cond = 0;
    for i = 1: 2
        % different pulse methods
        for j = 1:length(pulseMethod{i})
            cond = cond + 1;
            indRec = find(modPyr1AL.actOrInactPerRec == i & modPyr1AL.pulseMethPerRec == pulseMethod{i}(j));
            numField{cond}.actOrInact = i;
            numField{cond}.pulseMethod = pulseMethod{i}(j);
            numField{cond}.indRec = modPyr1AL.indRecPerRec(indRec);
            numField{cond}.nNeuWithFieldAlignedPerRec = modPyr1AL.nNeuWithFieldAlignedPerRec(indRec);
            numField{cond}.nNeuWithFieldAlignedPerRecGoodNonStim = modPyr1AL.nNeuWithFieldAlignedPerRecGoodNonStim(indRec);
            numField{cond}.nNeuWithFieldAlignedPerRecStim = modPyr1AL.nNeuWithFieldAlignedPerRecStim(indRec);
            numField{cond}.nNeuWithFieldAlignedPerRecStimCtrl = modPyr1AL.nNeuWithFieldAlignedPerRecStimCtrl(indRec);
            
            numField{cond}.percNeuWithFieldAlignedPerRec = modPyr1AL.nNeuWithFieldAlignedPerRec(indRec)...
                ./modPyr1AL.nNeuPerRec(indRec);
            numField{cond}.percNeuWithFieldAlignedPerRecGoodNonStim = modPyr1AL.nNeuWithFieldAlignedPerRecGoodNonStim(indRec)...
                ./modPyr1AL.nNeuPerRec(indRec);
            numField{cond}.percNeuWithFieldAlignedPerRecStim = modPyr1AL.nNeuWithFieldAlignedPerRecStim(indRec)...
                ./modPyr1AL.nNeuPerRec(indRec);
            numField{cond}.percNeuWithFieldAlignedPerRecStimCtrl = modPyr1AL.nNeuWithFieldAlignedPerRecStimCtrl(indRec)...
                ./modPyr1AL.nNeuPerRec(indRec);
            
            numField{cond}.pRSNumFieldGoodNonStimVsStim = ranksum(...
                numField{cond}.nNeuWithFieldAlignedPerRecGoodNonStim, numField{cond}.nNeuWithFieldAlignedPerRecStim);
            numField{cond}.pRSNumFieldStimVsStimCtrl = ranksum(...
                numField{cond}.nNeuWithFieldAlignedPerRecStim, numField{cond}.nNeuWithFieldAlignedPerRecStimCtrl);
            
            numField{cond}.pRSPercFieldGoodNonStimVsStim = ranksum(...
                numField{cond}.percNeuWithFieldAlignedPerRecGoodNonStim,numField{cond}.percNeuWithFieldAlignedPerRecStim);
            numField{cond}.pRSPercFieldStimVsStimCtrl = ranksum(...
                numField{cond}.percNeuWithFieldAlignedPerRecStim,numField{cond}.percNeuWithFieldAlignedPerRecStimCtrl);
    
            % field good non-stim vs. stim
            plotBars(numField{cond}.nNeuWithFieldAlignedPerRecGoodNonStim,...
                numField{cond}.nNeuWithFieldAlignedPerRecStim,...
                [mean(numField{cond}.nNeuWithFieldAlignedPerRecGoodNonStim),...
                mean(numField{cond}.nNeuWithFieldAlignedPerRecStim)],...
                [std(numField{cond}.nNeuWithFieldAlignedPerRecGoodNonStim)/sqrt(length(numField{cond}.nNeuWithFieldAlignedPerRecGoodNonStim)),...
                std(numField{cond}.nNeuWithFieldAlignedPerRecStim)/sqrt(length(numField{cond}.nNeuWithFieldAlignedPerRecStim))],...
                '','No. fields', ['p=' num2str(numField{cond}.pRSNumFieldGoodNonStimVsStim)],pathAnal,...
                ['NumFieldAlignedGoodNoStimVsStim_' actOrInActStr{i} '_P' num2str(pulseMethod{i}(j))]);
            
            % field stim ctrl vs stim
            plotBars(numField{cond}.nNeuWithFieldAlignedPerRecStimCtrl,...
                numField{cond}.nNeuWithFieldAlignedPerRecStim,...
                [mean(numField{cond}.nNeuWithFieldAlignedPerRecStimCtrl),...
                mean(numField{cond}.nNeuWithFieldAlignedPerRecStim)],...
                [std(numField{cond}.nNeuWithFieldAlignedPerRecStimCtrl)/sqrt(length(numField{cond}.nNeuWithFieldAlignedPerRecStimCtrl)),...
                std(numField{cond}.nNeuWithFieldAlignedPerRecStim)/sqrt(length(numField{cond}.nNeuWithFieldAlignedPerRecStim))],...
                '','No. fields', ['p=' num2str(numField{cond}.pRSNumFieldStimVsStimCtrl)],pathAnal,...
                ['NumFieldAlignedStimCtrlVsStim_' actOrInActStr{i} '_P' num2str(pulseMethod{i}(j))]);
            
            % field good non-stim vs. stim
            plotBars(numField{cond}.percNeuWithFieldAlignedPerRecGoodNonStim,...
                numField{cond}.percNeuWithFieldAlignedPerRecStim,...
                [mean(numField{cond}.percNeuWithFieldAlignedPerRecGoodNonStim),...
                mean(numField{cond}.percNeuWithFieldAlignedPerRecStim)],...
                [std(numField{cond}.percNeuWithFieldAlignedPerRecGoodNonStim)/sqrt(length(numField{cond}.percNeuWithFieldAlignedPerRecGoodNonStim)),...
                std(numField{cond}.percNeuWithFieldAlignedPerRecStim)/sqrt(length(numField{cond}.percNeuWithFieldAlignedPerRecStim))],...
                '','No. fields', ['p=' num2str(numField{cond}.pRSPercFieldGoodNonStimVsStim)],pathAnal,...
                ['PercFieldAlignedGoodNoStimVsStim_' actOrInActStr{i} '_P' num2str(pulseMethod{i}(j))]);
            
            % field stim ctrl vs stim
            plotBars(numField{cond}.percNeuWithFieldAlignedPerRecStimCtrl,...
                numField{cond}.percNeuWithFieldAlignedPerRecStim,...
                [mean(numField{cond}.percNeuWithFieldAlignedPerRecStimCtrl),...
                mean(numField{cond}.percNeuWithFieldAlignedPerRecStim)],...
                [std(numField{cond}.percNeuWithFieldAlignedPerRecStimCtrl)/sqrt(length(numField{cond}.percNeuWithFieldAlignedPerRecStimCtrl)),...
                std(numField{cond}.percNeuWithFieldAlignedPerRecStim)/sqrt(length(numField{cond}.percNeuWithFieldAlignedPerRecStim))],...
                '','No. fields', ['p=' num2str(numField{cond}.pRSPercFieldStimVsStimCtrl)],pathAnal,...
                ['PercFieldAlignedStimCtrlVsStim_' actOrInActStr{i} '_P' num2str(pulseMethod{i}(j))]);
        end
    end
    save([pathAnal 'initPeakPyrAllRecStim.mat'],'numField','-append');
    
    avgFRProfile = modPyr1AL.avgFRProfile; 
    avgFRProfileStim = modPyr1AL.avgFRProfileStim;
    avgFRProfileStimCtrl = modPyr1AL.avgFRProfileStimCtrl;
    
    if(methodKMean == 1)
        idxC = modPyr1AL.idxC1; 
    elseif(methodKMean == 2)
        idxC = modPyr1AL.idxC2;
    elseif(methodKMean == 3)
        idxC = modPyr1AL.idxC3;
    end
    
    FRProfileMeanAll = accumMean(avgFRProfile,modPyr1AL.timeStepRun);
    
    cond = 0;
    for i = 1: 2
        for j = 1:length(pulseMethod{i})
            cond = cond+1;
            ind = find(modPyr1AL.actOrInact == i & modPyr1AL.pulseMeth == pulseMethod{i}(j));
            FRProfile{cond}.actOrInact = i;
            FRProfile{cond}.pulseMethod = pulseMethod{i}(j);
            FRProfile{cond}.ind = ind;
            FRProfile{cond}.indRec = modPyr1AL.indRec(ind);
            FRProfile{cond}.indNeu = modPyr1AL.indNeu(ind);
            
            FRProfile{cond}.isNeuWithFieldAligned = modPyr1AL.isNeuWithFieldAligned(ind);
            FRProfile{cond}.isNeuWithFieldAlignedGoodNonStim = modPyr1AL.isNeuWithFieldAlignedGoodNonStim(ind);
            FRProfile{cond}.isNeuWithFieldAlignedStim = modPyr1AL.isNeuWithFieldAlignedStim(ind);
            FRProfile{cond}.isNeuWithFieldAlignedStimCtrl = modPyr1AL.isNeuWithFieldAlignedStimCtrl(ind);
            
            FRProfileMean{cond} = accumMean(avgFRProfile(ind,:),modPyr1AL.timeStepRun);

            FRProfileMeanStim{cond} = accumMean(avgFRProfileStim(ind,:),modPyr1AL.timeStepRun);

            FRProfileMeanStimCtrl{cond} = accumMean(avgFRProfileStimCtrl(ind,:),modPyr1AL.timeStepRun);

            % compare good non-stim and stim trials
            FRProfileMeanStatGoodNonStimVsStim{cond} = accumMeanStatCGoodBad(FRProfileMean{cond},FRProfileMeanStim{cond},idxC(ind),idxC(ind));

            % compare stim ctrl and stim trials
            FRProfileMeanStatStimCtrlVsStim{cond} = accumMeanStatCGoodBad(FRProfileMean{cond},FRProfileMeanStimCtrl{cond},idxC(ind),idxC(ind));
            
        end
    end
    
    save([pathAnal 'initPeakPyrAllRecStim_km' num2str(methodKMean) '.mat'],'FRProfileMeanAll',...
        'FRProfile','FRProfileMean','FRProfileMeanStim',...
        'FRProfileMeanStimCtrl','FRProfileMeanStatGoodNonStimVsStim','FRProfileMeanStatStimCtrlVsStim',...
        'anmNoInact','anmNoAct','pulseMethod'); 
    
end

function modPyr1 = accumPyr3(paths,filenames,mazeSess,autoCorrPyrAll,...
            nNeuWithFieldAligned,isNeuWithFieldAligned,...
            anmNoInact,anmNoAct,minFR,maxFR,task,sampleFq)
    numRec = size(paths,1);
    modPyr1 = struct('actOrInact',[],... % is it activation or inactivation
                              'task',[],... % no cue - 1, AL - 2, PL - 3
                              'indRec',[],... % recording index
                              'indNeu',[],... % neuron indices trials
                              'pulseMeth',[], ... pulse method
                              ...
                              'idxC1',[],...  % cluster no. Good trials
                              'idxC2',[],...  % cluster no. Good trials
                              'idxC3',[],...  % cluster no. Good trials
                              ...
                              'relDepthNeuHDef',[],... % depth
                              'nNeuWithFieldAligned',[],... % number of neurons with fields after aligning to run onset
                              'nNeuWithFieldAlignedGoodNonStim',[],... % number of neurons with fields after aligning to run onset, non-stim good
                              'nNeuWithFieldAlignedStim',[],... % number of neurons with fields after aligning to run onset, stim trials
                              'nNeuWithFieldAlignedStimCtrl',[],... % number of neurons with fields after aligning to run onset, stim ctrl trials
                              'isNeuWithFieldAligned',[],... % whether this neuron has a field
                              'isNeuWithFieldAlignedGoodNonStim',[],... % whether this neuron has a field, good non-stim trials
                              'isNeuWithFieldAlignedStim',[],... % whether this neuron has a field, stim trials
                              'isNeuWithFieldAlignedStimCtrl',[],... % whether this neuron has a field, stim ctrl trials
                              ...
                              'timeStepRun',[],...
                              'avgFRProfile',[],...% average firing rate profile good trials
                              'avgFRProfileStim',[],... % average firing rate profile stim trials
                              'avgFRProfileStimCtrl',[],... % average firing rate profile stim ctrl trials
                              ...
                              'nNeuPerRec',[],... % number of pyramidal neurons per recording
                              'nNeuWithFieldAlignedPerRec',[],... % number of field for each recording after aligning to run onset
                              'nNeuWithFieldAlignedPerRecGoodNonStim',[],... % number of field for each recording after aligning to run onset, good non-stim 
                              'nNeuWithFieldAlignedPerRecStim',[],... % number of field for each recording after aligning to run onset, stim trials
                              'nNeuWithFieldAlignedPerRecStimCtrl',[],... % number of field for each recording after aligning to run onset, stim ctrl trials
                              'taskPerRec',[],...
                              'indRecPerRec',[],...
                              'actOrInactPerRec',[],...
                              'pulseMethPerRec',[]);

    totRec = 0;
    for i = 1:numRec
        actOrInact = 0;
        anmNo = str2num(filenames(i,2:4));
        % does the recording belong to inactivation
        indTmp = find(anmNoInact == anmNo);
        if(~isempty(indTmp))
            actOrInact = 1;
        end
        % does the recording belong to activation
        indTmp1 = find(anmNoAct == anmNo);
        if(~isempty(indTmp1))
            actOrInact = 2;
        end
        % continue if the recording is not an activation or inactivation
        % recording
        if(isempty(indTmp) && isempty(indTmp1))
            continue;
        end
        
        fullPath = [paths(i,:) filenames(i,:) '.mat'];
        if(exist(fullPath) == 0)
            disp('File does not exist.');
            return;
        end
        load(fullPath,'cluList'); 
        
        fileNameInfo = [filenames(i,:) '_Info.mat'];
        fullPath = [paths(i,:) fileNameInfo];
        if(exist(fullPath) == 0)
            disp('_Info.mat file does not exist.');
            return;
        end
        load(fullPath,'autoCorr','beh'); 
        
        fileNamePeakFR = [filenames(i,:) '_PeakFR_msess' num2str(mazeSess(i)) ...
                        '_RunOnset0.mat'];
        fullPath = [paths(i,:) fileNamePeakFR];
        if(exist(fullPath) == 0)
            disp(['The peak firing rate file does not exist. Please call ',...
                    'function "PeakFiringRate_Aligned" first.']);
            return;
        end
        load(fullPath,'pFRNonStimGoodStruct','pFRStimStruct','pFRStimCtrlStruct');
        
        fullPath = [paths(i,:) filenames(i,:) '_PeakFRAligned_msess' num2str(mazeSess(i)) '_Run1.mat'];    
        if(exist(fullPath) == 0)
            disp('The _PeakFRAligned file does not exist');
            return;
        end
        load(fullPath,'trialNoNonStimGood','trialNoStim','trialNoStimCtrl','pulseMeth');
        
        fileNameConv = [filenames(i,:) '_convSpikesAligned_msess' num2str(mazeSess(i)) '_BefRun0.mat'];
        fullPath = [paths(i,:) fileNameConv];
        if(exist(fullPath) == 0)
            disp(['The convSpikesAligned file does not exist. Please call ',...
                    'function "ConvSpikeTrain_AlignedRunOnset" first.']);
            return;
        end
        load(fullPath,'timeStepRun');
        
        fullPathFR = [filenames(i,:) '_FR_Run1.mat'];
        fullPath = [paths(i,:) fullPathFR];
        if(exist(fullPath) == 0)
            disp('_FR_Run.mat file does not exist.');
            return;
        end
        load(fullPath,'mFRStruct','mFRStructSess'); 
        if(length(beh.mazeSessAll) > 1)
            mFR = mFRStructSess{mazeSess(i)};
        else
            mFR = mFRStruct;
        end
        
        fileNameFW = [filenames(i,:) '_FieldSpCorrAligned_Run' num2str(mazeSess(i)) ...
                        '_Run1.mat'];
        fullPath = [paths(i,:) fileNameFW];
        if(exist(fullPath) == 0)
            disp(['The field detection file does not exist. Please call ',...
                    'function "FieldDetectionAligned" first.']);
            return;
        end
        load(fullPath,'fieldSpCorrSessNonStimGood'); 
        fieldSpCorrSessAll = fieldSpCorrSessNonStimGood;
        
        fileNameFWStim = [filenames(i,:) '_FieldSpCorrAlignedStim_Run' num2str(mazeSess(i)) ...
                        '_Run1.mat'];
        fullPath = [paths(i,:) fileNameFWStim];
        if(exist(fullPath) == 0)
            disp(['The field detection file does not exist. Please call ',...
                    'function "FieldDetectionAlignedStim" first.']);
            return;
        end
        load(fullPath,'fieldSpCorrSessNonStimGood','fieldSpCorrSessStim','fieldSpCorrSessStimCtrl','paramF'); 
        
        indNeu = mFR.mFR > minFR & mFR.mFR < maxFR &...
                    autoCorr.isPyrneuron == 1;
        
        indTmp = find(autoCorrPyrAll.task == task & autoCorrPyrAll.indRec == i);
        if(length(indTmp) ~= sum(indNeu))
            disp(['the number of neurons in recording task = ' num2str(task) ' rec. no. = ' num2str(indRec)...
                    'does not match that in the autoCorrPyrAll struct.']);
        end
        
        if(length(pFRNonStimGoodStruct.indLapList) < paramF.minNumTr)
            disp([filenames(i,:) ' only has ' num2str(length(pFRNonStimGoodStruct.indLapList)) ...
                ' good trials.']);
            disp(['No. pyramidal neurons in this recording is ' num2str(sum(indNeu))]);
            
            continue;
        end
        
        for n = 1:length(pulseMeth)
            if(pulseMeth(n) ~= 2 && pulseMeth(n) ~= 3 && pulseMeth(n) ~= 4) % check the pulse method number
                disp(['Wrong pulse method = ' num2str(pulseMeth(n)), 'continue ...']);
                continue;
            end
            totRec = totRec+ 1;
            modPyr1.taskPerRec(totRec) = task;
            modPyr1.actOrInactPerRec(totRec) = actOrInact;
            modPyr1.pulseMethPerRec(totRec) = pulseMeth(n);
            modPyr1.indRecPerRec(totRec) = i;            
            modPyr1.task = [modPyr1.task task*ones(1,sum(indNeu))];
            modPyr1.indRec = [modPyr1.indRec i*ones(1,sum(indNeu))];
            modPyr1.indNeu = [modPyr1.indNeu find(indNeu == 1)]; 
            modPyr1.actOrInact = [modPyr1.actOrInact actOrInact*ones(1,sum(indNeu))];
            modPyr1.pulseMeth = [modPyr1.pulseMeth pulseMeth(n)*ones(1,sum(indNeu))];
            modPyr1.avgFRProfile = [modPyr1.avgFRProfile; pFRNonStimGoodStruct.avgFRProfile(indNeu,:)];
            modPyr1.avgFRProfileStim = [modPyr1.avgFRProfileStim; pFRStimStruct{n}.avgFRProfile(indNeu,:)];
            modPyr1.avgFRProfileStimCtrl = [modPyr1.avgFRProfileStimCtrl; pFRStimCtrlStruct{n}.avgFRProfile(indNeu,:)];
            
            if(length(indTmp) == sum(indNeu))
                modPyr1.idxC1 = [modPyr1.idxC1 autoCorrPyrAll.idxC1(indTmp)'];
                modPyr1.idxC2 = [modPyr1.idxC2 autoCorrPyrAll.idxC2(indTmp)'];
                modPyr1.idxC3 = [modPyr1.idxC3 autoCorrPyrAll.idxC3(indTmp)];
                modPyr1.relDepthNeuHDef = [modPyr1.relDepthNeuHDef autoCorrPyrAll.relDepthNeuHDef(indTmp)];
                modPyr1.nNeuPerRec(totRec) = sum(indNeu == 1);
                
                modPyr1.nNeuWithFieldAlignedPerRec(totRec) = unique(nNeuWithFieldAligned(indTmp));  
                modPyr1.nNeuWithFieldAligned = [modPyr1.nNeuWithFieldAligned nNeuWithFieldAligned(indTmp)];
                modPyr1.isNeuWithFieldAligned = [modPyr1.isNeuWithFieldAligned isNeuWithFieldAligned(indTmp)];
                
                modPyr1.nNeuWithFieldAlignedPerRecGoodNonStim(totRec) = 0;
                if(~isempty(fieldSpCorrSessNonStimGood))
                    modPyr1.nNeuWithFieldAlignedPerRecGoodNonStim(totRec) = length(unique(fieldSpCorrSessNonStimGood.indNeuron));
                end
                modPyr1.nNeuWithFieldAlignedGoodNonStim = [modPyr1.nNeuWithFieldAlignedGoodNonStim...
                    modPyr1.nNeuWithFieldAlignedPerRecGoodNonStim(totRec)*ones(1,sum(indNeu))];
                arrTmp = zeros(1,length(indNeu));
                if(~isempty(fieldSpCorrSessNonStimGood))
                    arrTmp(fieldSpCorrSessNonStimGood.indNeuron) = 1;
                end
                arrTmp = arrTmp(indNeu);
                modPyr1.isNeuWithFieldAlignedGoodNonStim = [modPyr1.isNeuWithFieldAlignedGoodNonStim...
                    arrTmp]; % number of neurons labelled here is <= nNeuWithField, because of the nNeuWithField consider all the neurons in the recording  
                
                modPyr1.nNeuWithFieldAlignedPerRecStim(totRec) = 0;
                if(~isempty(fieldSpCorrSessStim{n}))
                    modPyr1.nNeuWithFieldAlignedPerRecStim(totRec) = length(unique(fieldSpCorrSessStim{n}.indNeuron));
                end
                modPyr1.nNeuWithFieldAlignedStim = [modPyr1.nNeuWithFieldAlignedStim...
                    modPyr1.nNeuWithFieldAlignedPerRecStim(totRec)*ones(1,sum(indNeu))];
                arrTmp = zeros(1,length(indNeu));
                if(~isempty(fieldSpCorrSessStim{n}))
                    arrTmp(fieldSpCorrSessStim{n}.indNeuron) = 1;
                end
                arrTmp = arrTmp(indNeu);
                modPyr1.isNeuWithFieldAlignedStim = [modPyr1.isNeuWithFieldAlignedStim...
                    arrTmp];                
                
                modPyr1.nNeuWithFieldAlignedPerRecStimCtrl(totRec) = 0;
                if(~isempty(fieldSpCorrSessStimCtrl{n}))
                    modPyr1.nNeuWithFieldAlignedPerRecStimCtrl(totRec) = length(unique(fieldSpCorrSessStimCtrl{n}.indNeuron));
                end
                modPyr1.nNeuWithFieldAlignedStimCtrl = [modPyr1.nNeuWithFieldAlignedStimCtrl...
                    modPyr1.nNeuWithFieldAlignedPerRecStimCtrl(totRec)*ones(1,sum(indNeu))];
                arrTmp = zeros(1,length(indNeu));
                if(~isempty(fieldSpCorrSessStimCtrl{n}))
                    arrTmp(fieldSpCorrSessStimCtrl{n}.indNeuron) = 1;
                end
                arrTmp = arrTmp(indNeu);
                modPyr1.isNeuWithFieldAlignedStimCtrl = [modPyr1.isNeuWithFieldAlignedStimCtrl...
                    arrTmp];
            end
       
        end
        modPyr1.timeStepRun = timeStepRun/sampleFq;
        
    end
end

function FRProfileMean = accumMean(avgFRProfile,timeStep)
       
    % baseline
    indFRBaseline = find(timeStep >= -3 & timeStep < -2);
    FRProfileMean.indFRBaseline = indFRBaseline;
    FRProfileMean.meanAvgFRProfileBaseline =  mean(avgFRProfile(:,indFRBaseline),2);
    
    % -1.5- -0.5 sec
    indFRBefRun = find(timeStep >= -1.5 & timeStep < -0.5);  
    FRProfileMean.indFRBefRun = indFRBefRun;
    FRProfileMean.meanAvgFRProfileBefRun = mean(avgFRProfile(:,indFRBefRun),2);
    
    % 0.5-1.5 sec
    indFR0to1 = find(timeStep >= 0.5 & timeStep < 1.5);  
    FRProfileMean.indFR0to1 = indFR0to1;
    FRProfileMean.meanAvgFRProfile0to1 = mean(avgFRProfile(:,indFR0to1),2);
    
    % 3-5 sec
    indFR3to5 = find(timeStep >= 3 & timeStep < 5);  
    FRProfileMean.indFR3to5 = indFR3to5;
    FRProfileMean.meanAvgFRProfile3to5 = mean(avgFRProfile(:,indFR3to5),2);
    
    % perc change from 0.5-1.5 s to baseline
    FRProfileMean.percChange0to1VsBL = FRProfileMean.meanAvgFRProfile0to1...
        ./FRProfileMean.meanAvgFRProfileBaseline;
    
    % perc change -1.5- -0.5 s to baseline
    FRProfileMean.percChangeBefRunVsBL = FRProfileMean.meanAvgFRProfileBefRun...
        ./FRProfileMean.meanAvgFRProfileBaseline;
    
    % perc change 0.5 to 1.5 s to -1.5- -0.5 s 
    FRProfileMean.percChangeBefRunVs0to1 = FRProfileMean.meanAvgFRProfile0to1...
        ./FRProfileMean.meanAvgFRProfileBefRun;
    
    % perc change from 0.5-1.5 s to 3-5s
    FRProfileMean.percChange0to1Vs3to5 = FRProfileMean.meanAvgFRProfile0to1...
        ./FRProfileMean.meanAvgFRProfile3to5;
    
    % perc change -1.5- -0.5 s to 3-5s
    FRProfileMean.percChangeBefRunVs3to5 = FRProfileMean.meanAvgFRProfileBefRun...
        ./FRProfileMean.meanAvgFRProfile3to5;
end

function FRProfileMeanStatC = accumMeanStatC(FRProfileMean,idxC,nNeuWithField,isNeuWithField)
    
    numC = max(idxC);
    for i = 1:numC
        idxCI = idxC == i;
        
        %% recordings with field vs without field
        indCurCField = idxC == i & nNeuWithField >= 2;
        indCurCNoField = idxC == i & nNeuWithField < 1;
        FRProfileMeanStatC.pRS0to1VsBL(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(idxCI),...
                    FRProfileMean.meanAvgFRProfile0to1(idxCI));
        FRProfileMeanStatC.pRS0to1VsBLField(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(indCurCField),...
                    FRProfileMean.meanAvgFRProfile0to1(indCurCField));
        FRProfileMeanStatC.pRS0to1VsBLNoField(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfile0to1(indCurCNoField));        
                            
        FRProfileMeanStatC.pRSBefRunVsBL(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(idxCI),...
                    FRProfileMean.meanAvgFRProfileBefRun(idxCI));
        FRProfileMeanStatC.pRSBefRunVsBLField(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(indCurCField),...
                    FRProfileMean.meanAvgFRProfileBefRun(indCurCField));
        FRProfileMeanStatC.pRSBefRunVsBLNoField(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfileBefRun(indCurCNoField));
                
        FRProfileMeanStatC.pRS3to5VsBL(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(idxCI),...
                    FRProfileMean.meanAvgFRProfile3to5(idxCI));
        FRProfileMeanStatC.pRS3to5VsBLField(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(indCurCField),...
                    FRProfileMean.meanAvgFRProfile3to5(indCurCField));
        FRProfileMeanStatC.pRS3to5VsBLNoField(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfile3to5(indCurCNoField));
                
        FRProfileMeanStatC.pRSBefRunVs0to1(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(idxCI),...
                    FRProfileMean.meanAvgFRProfileBefRun(idxCI));
        FRProfileMeanStatC.pRSBefRunVs0to1Field(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(indCurCField),...
                    FRProfileMean.meanAvgFRProfileBefRun(indCurCField));
        FRProfileMeanStatC.pRSBefRunVs0to1NoField(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfileBefRun(indCurCNoField));
                
        FRProfileMeanStatC.pRS3to5Vs0to1(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(idxCI),...
                    FRProfileMean.meanAvgFRProfile3to5(idxCI));
        FRProfileMeanStatC.pRS3to5Vs0to1Field(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(indCurCField),...
                    FRProfileMean.meanAvgFRProfile3to5(indCurCField));
        FRProfileMeanStatC.pRS3to5Vs0to1NoField(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfile3to5(indCurCNoField));
                
        FRProfileMeanStatC.pRS3to5VsBefRun(i) = ranksum(FRProfileMean.meanAvgFRProfileBefRun(idxCI),...
                    FRProfileMean.meanAvgFRProfile3to5(idxCI));
        FRProfileMeanStatC.pRS3to5VsBefRunField(i) = ranksum(FRProfileMean.meanAvgFRProfileBefRun(indCurCField),...
                    FRProfileMean.meanAvgFRProfile3to5(indCurCField));
        FRProfileMeanStatC.pRS3to5VsBefRunNoField(i) = ranksum(FRProfileMean.meanAvgFRProfileBefRun(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfile3to5(indCurCNoField));
                
        FRProfileMeanStatC.pTTPercChange0to1VsBL(i) = ttest(FRProfileMean.percChange0to1VsBL(idxCI));
        FRProfileMeanStatC.pRSPercChange0to1VsBLFieldVsNoField(i) = ...
            ranksum(FRProfileMean.percChange0to1VsBL(indCurCField),...
            FRProfileMean.percChange0to1VsBL(indCurCNoField));
        
        FRProfileMeanStatC.pTTPercChangeBefRunVsBL(i) = ttest(FRProfileMean.percChangeBefRunVsBL(idxCI));
        FRProfileMeanStatC.pRSPercChangeBefRunVsBLFieldVsNoField(i) = ...
            ranksum(FRProfileMean.percChangeBefRunVsBL(indCurCField),...
            FRProfileMean.percChangeBefRunVsBL(indCurCNoField));
        
        FRProfileMeanStatC.pTTPercChangeBefRunVs0to1(i) = ttest(FRProfileMean.percChangeBefRunVs0to1(idxCI));
        FRProfileMeanStatC.pRSPercChangeBefRunVs0to1FieldVsNoField(i) = ...
            ranksum(FRProfileMean.percChangeBefRunVs0to1(indCurCField),...
            FRProfileMean.percChangeBefRunVs0to1(indCurCNoField));
        
        FRProfileMeanStatC.pTTPercChange0to1Vs3to5(i) = ttest(FRProfileMean.percChange0to1Vs3to5(idxCI));
        FRProfileMeanStatC.pRSPercChange0to1Vs3to5FieldVsNoField(i) = ...
            ranksum(FRProfileMean.percChange0to1Vs3to5(indCurCField),...
            FRProfileMean.percChange0to1Vs3to5(indCurCNoField));
        
        FRProfileMeanStatC.pTTPercChangeBefRunVs3to5(i) = ttest(FRProfileMean.percChangeBefRunVs3to5(idxCI));
        FRProfileMeanStatC.pRSPercChangeBefRunVs3to5FieldVsNoField(i) = ...
            ranksum(FRProfileMean.percChangeBefRunVs3to5(indCurCField),...
            FRProfileMean.percChangeBefRunVs3to5(indCurCNoField));
        
        %% neurons with field vs without field
        indCurCField = idxC == i & isNeuWithField == 1;
        indCurCNoField = idxC == i & isNeuWithField == 0;
        FRProfileMeanStatC.pRS0to1VsBLFieldNeu(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(indCurCField),...
                    FRProfileMean.meanAvgFRProfile0to1(indCurCField));
        FRProfileMeanStatC.pRS0to1VsBLNoFieldNeu(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfile0to1(indCurCNoField));        
                            
        FRProfileMeanStatC.pRSBefRunVsBLFieldNeu(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(indCurCField),...
                    FRProfileMean.meanAvgFRProfileBefRun(indCurCField));
        FRProfileMeanStatC.pRSBefRunVsBLNoFieldNeu(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfileBefRun(indCurCNoField));
                
        FRProfileMeanStatC.pRSBefRunVs0to1FieldNeu(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(indCurCField),...
                    FRProfileMean.meanAvgFRProfileBefRun(indCurCField));
        FRProfileMeanStatC.pRSBefRunVs0to1NoFieldNeu(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfileBefRun(indCurCNoField));
                
        FRProfileMeanStatC.pRS3to5VsBLFieldNeu(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(indCurCField),...
                    FRProfileMean.meanAvgFRProfile3to5(indCurCField));
        FRProfileMeanStatC.pRS3to5VsBLNoFieldNeu(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfile3to5(indCurCNoField));
                
        FRProfileMeanStatC.pRS3to5Vs0to1FieldNeu(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(indCurCField),...
                    FRProfileMean.meanAvgFRProfile3to5(indCurCField));
        FRProfileMeanStatC.pRS3to5Vs0to1NoFieldNeu(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfile3to5(indCurCNoField));
                
        FRProfileMeanStatC.pRS3to5VsBefRunFieldNeu(i) = ranksum(FRProfileMean.meanAvgFRProfileBefRun(indCurCField),...
                    FRProfileMean.meanAvgFRProfile3to5(indCurCField));
        FRProfileMeanStatC.pRS3to5VsBefRunNoFieldNeu(i) = ranksum(FRProfileMean.meanAvgFRProfileBefRun(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfile3to5(indCurCNoField));
                
        FRProfileMeanStatC.pRSPercChange0to1VsBLFieldNeuVsNoFieldNeu(i) = ...
            ranksum(FRProfileMean.percChange0to1VsBL(indCurCField),...
            FRProfileMean.percChange0to1VsBL(indCurCNoField));
        
        FRProfileMeanStatC.pRSPercChangeBefRunVsBLFieldNeuVsNoFieldNeu(i) = ...
            ranksum(FRProfileMean.percChangeBefRunVsBL(indCurCField),...
            FRProfileMean.percChangeBefRunVsBL(indCurCNoField));
        
        FRProfileMeanStatC.pRSPercChangeBefRunVs0to1FieldNeuVsNoFieldNeu(i) = ...
            ranksum(FRProfileMean.percChangeBefRunVs0to1(indCurCField),...
            FRProfileMean.percChangeBefRunVs0to1(indCurCNoField));
        
        FRProfileMeanStatC.pRSPercChange0to1Vs3to5FieldNeuVsNoFieldNeu(i) = ...
            ranksum(FRProfileMean.percChange0to1Vs3to5(indCurCField),...
            FRProfileMean.percChange0to1Vs3to5(indCurCNoField));
        
        FRProfileMeanStatC.pRSPercChangeBefRunVs3to5FieldNeuVsNoFieldNeu(i) = ...
            ranksum(FRProfileMean.percChangeBefRunVs3to5(indCurCField),...
            FRProfileMean.percChangeBefRunVs3to5(indCurCNoField));
    end
end

function FRProfileMeanStatC = accumMeanStatCCmp(FRProfileMean,idxC)
    
    idxC1 = idxC == 1;
    idxC2 = idxC == 2;
        
    FRProfileMeanStatC.pRSPercChange0to1VsBLC = ...
        ranksum(FRProfileMean.percChange0to1VsBL(idxC1),...
        FRProfileMean.percChange0to1VsBL(idxC2));

    FRProfileMeanStatC.pRSPercChangeBefRunVsBLC = ...
        ranksum(FRProfileMean.percChangeBefRunVsBL(idxC1),...
        FRProfileMean.percChangeBefRunVsBL(idxC2));
    
    FRProfileMeanStatC.pRSPercChangeBefRunVs0to1C = ...
        ranksum(FRProfileMean.percChangeBefRunVs0to1(idxC1),...
        FRProfileMean.percChangeBefRunVs0to1(idxC2));

    FRProfileMeanStatC.pRSPercChange0to1Vs3to5C = ...
        ranksum(FRProfileMean.percChange0to1Vs3to5(idxC1),...
        FRProfileMean.percChange0to1Vs3to5(idxC2));

    FRProfileMeanStatC.pRSPercChangeBefRunVs3to5C = ...
        ranksum(FRProfileMean.percChangeBefRunVs3to5(idxC1),...
        FRProfileMean.percChangeBefRunVs3to5(idxC2));

end

function FRProfileMeanStatC = accumMeanStatCGoodBad(FRProfileMean,FRProfileMeanBad,idxC,idxCBad)
    
    numC = max(idxC);
    for i = 1:numC
        idxCI = idxC == i;
        idxCIBad = idxCBad == i;
        FRProfileMeanStatC.pRSBL(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(idxCI),...
                    FRProfileMeanBad.meanAvgFRProfileBaseline(idxCIBad));
                                
        FRProfileMeanStatC.pRSBefRun(i) = ranksum(FRProfileMean.meanAvgFRProfileBefRun(idxCI),...
                    FRProfileMeanBad.meanAvgFRProfileBefRun(idxCIBad));
               
        FRProfileMeanStatC.pRS3to5(i) = ranksum(FRProfileMean.meanAvgFRProfile3to5(idxCI),...
                    FRProfileMeanBad.meanAvgFRProfile3to5(idxCIBad));
               
        FRProfileMeanStatC.pRS0to1(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(idxCI),...
                    FRProfileMeanBad.meanAvgFRProfile0to1(idxCIBad));
                
        % perc change from 0.5-1.5 s to baseline
        FRProfileMeanStatC.pRSPercChange0to1VsBL(i) = ranksum(FRProfileMean.percChange0to1VsBL(idxCI),...
                    FRProfileMeanBad.percChange0to1VsBL(idxCIBad));

        % perc change -1.5- -0.5 s to baseline
        FRProfileMeanStatC.pRSPercChangeBefRunVsBL(i) = ranksum(FRProfileMean.percChangeBefRunVsBL(idxCI),...
                    FRProfileMeanBad.percChangeBefRunVsBL(idxCIBad));

        % perc change 0.5-1.5 s to -1.5- -0.5 s 
        FRProfileMeanStatC.pRSPercChangeBefRunVs0to1(i) = ranksum(FRProfileMean.percChangeBefRunVs0to1(idxCI),...
                    FRProfileMeanBad.percChangeBefRunVs0to1(idxCIBad));

        % perc change from 0.5-1.5 s to 3-5s
        FRProfileMeanStatC.pRSPercChange0to1Vs3to5(i) = ranksum(FRProfileMean.percChange0to1Vs3to5(idxCI),...
                    FRProfileMeanBad.percChange0to1Vs3to5(idxCIBad));

        % perc change -1.5- -0.5 s to 3-5s
        FRProfileMeanStatC.pRSPercChangeBefRunVs3to5(i) = ranksum(FRProfileMean.percChangeBefRunVs3to5(idxCI),...
                    FRProfileMeanBad.percChangeBefRunVs3to5(idxCIBad));
        
    end
    
    FRProfileMeanStatC.pRSBLAll = ranksum(FRProfileMean.meanAvgFRProfileBaseline,...
                    FRProfileMeanBad.meanAvgFRProfileBaseline);
                                
    FRProfileMeanStatC.pRSBefRunAll = ranksum(FRProfileMean.meanAvgFRProfileBefRun,...
                FRProfileMeanBad.meanAvgFRProfileBefRun);

    FRProfileMeanStatC.pRS3to5All = ranksum(FRProfileMean.meanAvgFRProfile3to5,...
                FRProfileMeanBad.meanAvgFRProfile3to5);

    FRProfileMeanStatC.pRS0to1All = ranksum(FRProfileMean.meanAvgFRProfile0to1,...
                FRProfileMeanBad.meanAvgFRProfile0to1);

    % perc change from 0.5-1.5 s to baseline
    FRProfileMeanStatC.pRSPercChange0to1VsBLAll = ranksum(FRProfileMean.percChange0to1VsBL,...
                FRProfileMeanBad.percChange0to1VsBL);

    % perc change -1.5- -0.5 s to baseline
    FRProfileMeanStatC.pRSPercChangeBefRunVsBLAll = ranksum(FRProfileMean.percChangeBefRunVsBL,...
                FRProfileMeanBad.percChangeBefRunVsBL);

    % perc change 0.5-1.5 s to -1.5- -0.5 s 
    FRProfileMeanStatC.pRSPercChangeBefRunVs0to1All = ranksum(FRProfileMean.percChangeBefRunVs0to1,...
                FRProfileMeanBad.percChangeBefRunVs0to1);

    % perc change from 0.5-1.5 s to 3-5s
    FRProfileMeanStatC.pRSPercChange0to1Vs3to5All = ranksum(FRProfileMean.percChange0to1Vs3to5,...
                FRProfileMeanBad.percChange0to1Vs3to5);

    % perc change -1.5- -0.5 s to 3-5s
    FRProfileMeanStatC.pRSPercChangeBefRunVs3to5All = ranksum(FRProfileMean.percChangeBefRunVs3to5,...
                FRProfileMeanBad.percChangeBefRunVs3to5);
end

function plotBoxPlot(x1,x2,yl,fn,pathAnal,ylimit,p,colorSel)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 200 400])
    if(colorSel == 0)
        colorArr = [163 207 98;...
                234 131 114]/255;
    elseif(colorSel == 1)            
        colorArr = [234 131 114;...
                116 53 61]/255;
    else        
        colorArr = [163 207 98;... 
            63 79 37]/255;
    end
    x = [x1;x2];
    g = [repmat({'C1'},length(x1),1);...
        repmat({'C2'},length(x2),1)];
    boxplot(x,g,'Notch','on','Widths',0.3,'Symbol','');
    h = findobj(gca,'Tag','Box');
    for j = 1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colorArr(j,:),'FaceAlpha',0.5);
    end
    if(~isempty(ylimit))
        set(gca,'YLim',ylimit);
    end
    ylabel(yl);
    title(['p = ' num2str(p)]);
    
    fileName1 = [pathAnal fn];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotAvgFRProfile(timeStepRun,avgFRProfile,yl,fileName,pathAnal,ylimit)
    options.handle     = figure;
    set(options.handle,'OuterPosition',...
        [500 500 280 280])
    options.color_area = [27 117 187]./255;    % Blue theme
    options.color_line = [ 39 169 225]./255;
    options.alpha      = 0.5;
    options.line_width = 0.5;
    options.error      = 'sem';
    options.x_axis = timeStepRun;
    plot_areaerrorbar(avgFRProfile,options);
    hold on;
    h = plot([0 0],[min(mean(avgFRProfile)-std(avgFRProfile)/sqrt(size(avgFRProfile,1)))*0.95 ...
        max(mean(avgFRProfile)+std(avgFRProfile)/sqrt(size(avgFRProfile,1)))*1.05],'r-');
    set(h,'LineWidth',1)
    set(gca,'XLim',[timeStepRun(1) 7]);
    if(~isempty(ylimit))
        set(gca,'YLim',ylimit);
    else
        set(gca,'YLim',[min(mean(avgFRProfile)-std(avgFRProfile)/sqrt(size(avgFRProfile,1)))*0.95 ...
        max(mean(avgFRProfile)+std(avgFRProfile)/sqrt(size(avgFRProfile,1)))*1.05]);
    end
    xlabel('Time (s)')
    ylabel(yl)
    
    fileName1 = [pathAnal fileName];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
        
end

function plotAvgFRProfileCmp(timeStepRun,avgFRProfilex,avgFRProfiley, yl,fileName,pathAnal,ylimit)
    options.handle     = figure;
    set(options.handle,'OuterPosition',...
        [500 500 280 280])
    options.color_areaX = [27 117 187]./255;    % Blue theme
    options.color_lineX = [ 39 169 225]./255;
    options.color_areaY = [187 189 192]./255;    % Orange theme
    options.color_lineY = [167 169  171]./255;
    options.alpha      = 0.5;
    options.line_width = 0.5;
    options.error      = 'sem';
    options.x_axisX = timeStepRun;
    options.x_axisY = timeStepRun;
    plot_areaerrorbarXY(avgFRProfilex, avgFRProfiley,...
        options);
    hold on;
    minX = min(mean(avgFRProfilex)-std(avgFRProfilex)/sqrt(size(avgFRProfilex,1)));
    minY = min(mean(avgFRProfiley)-std(avgFRProfiley)/sqrt(size(avgFRProfiley,1)));
    maxX = max(mean(avgFRProfilex)+std(avgFRProfilex)/sqrt(size(avgFRProfilex,1)));
    maxY = max(mean(avgFRProfiley)+std(avgFRProfiley)/sqrt(size(avgFRProfiley,1)));
    if(~isempty(ylimit))
        h = plot([0 0],ylimit,'r-');
    else
        h = plot([0 0],[min([minX minY])*0.95 ...
            max([maxX maxY])*1.05],'r-');
    end
    set(h,'LineWidth',1)
    set(gca,'XLim',[-1 4]);
%     set(gca,'XLim',[timeStepRun(1) 7]);
    if(~isempty(ylimit))
        set(gca,'YLim',ylimit);
    else
        set(gca,'YLim',[min([minX minY])*0.95 ...
        max([maxX maxY])*1.05]);
    end
    xlabel('Time (s)')
    ylabel(yl)
    
    fileName1 = [pathAnal fileName];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
        
end

function plotBars(data1,data2,mean,std,x1,y1,t,pathAnal,fn)
    fig = figure;
    set(0,'Units','pixels') 
    set(figure(fig),'OuterPosition',...
        [500 500 210 280])
    fig.Renderer = 'Painters';
    
    h = bar([1,2],mean,0.5);
    set(h,'EdgeColor',[0.3 0.3 0.3],'FaceColor',[187 189 192]/255);
    hold
    
    h = errorbar([1,2],mean,std);
    set(h,'Marker','.','MarkerSize',0.1,'Color',[0 0 0],'LineStyle','none')
    
    h = plot(1+0.15*rand(1,length(data1)),data1,'o');
    set(h,'MarkerSize',3,'Color',[167 169 171]/255);
    
    h = plot(2+0.15*rand(1,length(data2)),data2,'o');
    set(h,'MarkerSize',3,'Color',[27 117 187]/255);
    set(gca,'XLim',[0.5 2.5],'YLim',[0 max([data1 data2])+0.01]);
    
    xlabel(x1);
    ylabel(y1);
    title(t);
    
    fileName1 = [pathAnal fn];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end
