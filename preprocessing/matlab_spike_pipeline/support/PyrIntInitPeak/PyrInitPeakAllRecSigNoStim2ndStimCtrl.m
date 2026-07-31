function PyrInitPeakAllRecSigNoStim2ndStimCtrl()
            
    RecordingListPyrIntSST;  % include all the early NoCue (no blackout), PL and AL (up to A057)
    pathAnal0 = ['Z:\Yingxue\Draft\SST\Pyramidal\'];
    
    if(exist([pathAnal0 'initPeakPyrAllRecSigNoStim2ndStimCtrl.mat']))
        load([pathAnal0 'initPeakPyrAllRecSigNoStim2ndStimCtrl.mat']);
    end
    pathAnalPeak0 = ['Z:\Yingxue\Draft\SST\PyramidalInitPeakALPL\2\'];
    if(exist([pathAnalPeak0 'initPeakPyrAllRec_km2.mat']))
        load([pathAnalPeak0 'initPeakPyrAllRec_km2.mat'],'FRProfileMean');
    end
    
    GlobalConstFqSST;
    
    if(exist('modPyr1SigAL') == 0)
        disp('No cue - peak firing rate')
        modPyr1SigNoCue = ...
            pyrInitSig(listRecordingsNoCuePath,...
                listRecordingsNoCueFileName,mazeSessionNoCue,...
                FRProfileMean.indFR0to1,FRProfileMean.indFRBefRun,...
                1,minFR,maxFR,p,numShuffle);
            
        fullpath = [pathAnal0 'initPeakPyrAllRecSigNoStim2ndStimCtrl.mat'];
        save(fullpath, 'modPyr1SigNoCue','-v7.3');
            
        disp('Active licking - peak firing rate')
        modPyr1SigAL = ...
            pyrInitSig(listRecordingsActiveLickPath,...
            listRecordingsActiveLickFileName,mazeSessionActiveLick,...
                FRProfileMean.indFR0to1,FRProfileMean.indFRBefRun,...
                2,minFR,maxFR,p,numShuffle);
        
        save(fullpath, 'modPyr1SigAL','-append');
        
        disp('Passive licking - peak firing rate')
        modPyr1SigPL = ...
            pyrInitSig(listRecordingsPassiveLickPath,...
                listRecordingsPassiveLickFileName,mazeSessionPassiveLick,...
                FRProfileMean.indFR0to1,FRProfileMean.indFRBefRun,...
                3,minFR,maxFR,p,numShuffle);
        
        save(fullpath, 'modPyr1SigPL','-append');        
    end
end

function modPyr1Sig = ...
    pyrInitSig(paths,filenames,mazeSess,indAft,indBef,...
    task,minFR,maxFR,p,numShuffle)
   
    modPyr1Sig.task = [];
    modPyr1Sig.indRec = [];
    modPyr1Sig.indNeu = []; 
    modPyr1Sig.isPeakNeuArrAll = [];
    modPyr1Sig.ratioAftBefShufArrAll = [];
    modPyr1Sig.trialNoCtrlGood = [];
    modPyr1Sig.trialNo2ndStimCtrl = [];
    modPyr1Sig.trialNo2ndStimCtrlGood = [];
    
    numRec = size(paths,1);
    for i = 1:numRec
        disp(['rec' num2str(i)])
        fileNamePeakFR = [filenames(i,:) '_PeakFR_msess' num2str(mazeSess(i)) ...
                        '_RunOnset0.mat'];
        fullPath = [paths(i,:) fileNamePeakFR];
        if(exist(fullPath) == 0)
            disp(['The peak firing rate file does not exist. Please call ',...
                    'function "PeakFiringRate_Aligned" first.']);
            return;
        end
        load(fullPath,'pFRNonStimGoodStruct');
        
        fileNameConv = [filenames(i,:) '_convSpikesAligned_msess' num2str(mazeSess(i)) '_BefRun0.mat'];
        fullPath = [paths(i,:) fileNameConv];
        if(exist(fullPath,'file') == 0)
            disp(['The firing profile file does not exist. Try to run the',...
                        ' "ConvSpikeTrain_AlignedRunOnset" function first.']);
            return;
        end
        load(fullPath);
        
        fullPathFR = [filenames(i,:) '_FR_Run1.mat'];
        fullPath = [paths(i,:) fullPathFR];
        if(exist(fullPath) == 0)
            disp('_FR_Run.mat file does not exist.');
            return;
        end
        load(fullPath,'mFRStruct','mFRStructSess'); 
        if(length(mFRStructSess) > 1)
            mFR = mFRStructSess{mazeSess(i)};
        else
            mFR = mFRStruct;
        end
        
        fullPath = [paths(i,:) filenames(i,:) '_behPar_msess' num2str(mazeSess(i)) '.mat']; 
        if(exist(fullPath) == 0)
            disp('The _behPar file does not exist');
            return;
        end
        load(fullPath);

        fullPath = [paths(i,:) filenames(i,:) '_PeakFRAligned_msess' num2str(mazeSess(i))...
            '_Run1.mat'];    
        if(exist(fullPath) == 0)
            disp('The _PeakFRAligned file does not exist');
            return;
        end
        load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','pulseMeth');
        
        fileNameInfo = [filenames(i,:) '_Info.mat'];
        fullPath = [paths(i,:) fileNameInfo];
        if(exist(fullPath) == 0)
            disp('_Info.mat file does not exist.');
            return;
        end
        load(fullPath,'autoCorr','beh'); 
        
        fileNameFW = [filenames(i,:) '_FieldSpCorrAligned_Run' num2str(mazeSess(i)) ...
                        '_Run1.mat'];
        fullPath = [paths(i,:) fileNameFW];
        if(exist(fullPath) == 0)
            disp(['The field detection file does not exist. Please call ',...
                    'function "FieldDetection_GoodTr" first.']);
            return;
        end
        load(fullPath,'paramF'); 
                
        indNeu = find(mFR.mFR > minFR & mFR.mFR < maxFR &...
                    autoCorr.isPyrneuron == 1);

        trialNoStimCtrlGood = [];
        trialNo2ndStimCtrl = [];
        trialNoStimGood = [];
        if(~isempty(trialNoStim))
            for j = 1:length(trialNoStim)
                ind = find(beh.mazeSess == mazeSess(i));
                indGood = find(behPar.indTrBadBeh == 0);
                indGood = ind(intersect(ind,indGood));
                stimInterval = diff(trialNoStim{j});
                interval2 = sum(stimInterval == 2);
                interval3 = sum(stimInterval == 3);
                if(pulseMeth(j) == 4)
                    if(interval2>interval3)
                        trialNoStimCtrl = trialNoStim{j}-1;
                    else
                        trialNoStimCtrl = trialNoStim{j}-2;
                    end
                else
                    trialNoStimCtrl = trialNoStim{j}-1;
                end
                trialNoStimTmp = ismember(trialNoStimCtrl,trialNoStim{j});
                trialNoStimCtrl = trialNoStimCtrl(trialNoStimCtrl > 0 & trialNoStimTmp == 0);
                trialNoStimCtrlGood = [trialNoStimCtrlGood; intersect(trialNoStimCtrl,indGood)];
                trialNo2ndStimCtrl = [trialNo2ndStimCtrl; trialNoStimCtrl];
                trialNoStimGood{j} = intersect(trialNoStim{j},indGood);
            end
        end
        trialNoNonStimGood = [trialNoNonStimGood; trialNoStimCtrlGood];
        
        if(length(trialNoNonStimGood) >= paramF.minNumTr)
        
            % median trial length
            medTrLen = prctile(behPar.numSamplesRun,75)/1250;
            indTime = find(timeStepRun/1250 >= medTrLen,1);
            if(isempty(indTime))
                indTime = length(timeStepRun);
            end

            isPeakNeuArr = zeros(length(indNeu),length(p));
            ratioAftBefShufArr = zeros(numShuffle,length(indNeu));

            % detect neurons with peaks based on all the trials in the session
            for n = 1:length(indNeu)
%                 disp(['Neuron ' num2str(n)])
                [isPeakNeuArr(n,:),ratioAftBefShufArr(:,n)] = neuPeakDetection(...
                    filteredSpikeArrayRunOnSet{indNeu(n)}(trialNoNonStimGood,1:indTime),...
                    indAft,indBef,p,numShuffle);
            end

            modPyr1Sig.task = [modPyr1Sig.task task*ones(1,length(indNeu))];
            modPyr1Sig.indRec = [modPyr1Sig.indRec i*ones(1,length(indNeu))];
            modPyr1Sig.indNeu = [modPyr1Sig.indNeu indNeu]; 
            modPyr1Sig.isPeakNeuArrAll = [modPyr1Sig.isPeakNeuArrAll; isPeakNeuArr];
            modPyr1Sig.ratioAftBefShufArrAll = [modPyr1Sig.ratioAftBefShufArrAll ratioAftBefShufArr];
            modPyr1Sig.trialNoCtrlGood{i} = trialNoNonStimGood;
            modPyr1Sig.trialNo2ndStimCtrlGood{i} = trialNoStimCtrlGood;
            modPyr1Sig.trialNo2ndStimCtrl{i} = trialNo2ndStimCtrl;
            modPyr1Sig.trialNoStimGood{i} = trialNoStimGood;
        else
            disp(['Recording ' filenames(i,:) ' has ' num2str(length(trialNoNonStimGood)) ' good trials']);
        end
    end
    
end