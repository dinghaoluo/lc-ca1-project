function PyrInitPeakAllRecSig()
            
    RecordingList_dinghao;
    pathAnal = 'Z:\Dinghao\code_dinghao\HPC_PyrUp_PyrDown';
    
    if(exist([pathAnal 'initPeakPyrAllRecSig.mat']))
        load([pathAnal 'initPeakPyrAllRecSig.mat']);
    end
    
    indFR0to1 = 601:1000;
    indFRBefRun = 1401:1800;
    
    GlobalConstFqSST;
    
    if(exist('modPyr1SigNoCue') == 0)
%         disp('No cue - peak firing rate')
%         modPyr1SigNoCue = ...
%             pyrInitSig(listRecordingsNoCuePath,...
%                 listRecordingsNoCueFileName,mazeSessionNoCue,...
%                 FRProfileMean.indFR0to1,FRProfileMean.indFRBefRun,...
%                 1,minFR,maxFR,p,numShuffle);
%         
%         fullpath = [pathAnal 'initPeakPyrAllRecSig.mat'];
%         save(fullpath, 'modPyr1SigNoCue','-v7.3');
        
        disp('Active licking - peak firing rate')
        modPyr1SigAL = ...
            pyrInitSig(listRecordingsActiveLickPathHPCLCOpt,...
            listRecordingsActiveLickFileNameHPCLCOpt,mazeSessionActiveLickHPCLCOpt,...
                indFR0to1,indFRBefRun,...
                2,minFR,maxFR,p,numShuffle);
        
        save(fullpath, 'modPyr1SigAL','-append');
        
%         disp('Passive licking - peak firing rate')
%         modPyr1SigPL = ...
%             pyrInitSig(listRecordingsPassiveLickPath,...
%                 listRecordingsPassiveLickFileName,mazeSessionPassiveLick,...
%                 FRProfileMean.indFR0to1,FRProfileMean.indFRBefRun,...
%                 3,minFR,maxFR,p,numShuffle);
        
%         save(fullpath, 'modPyr1SigPL','-append');
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
    
    numRec = size(paths,1);
    for i = 1:numRec
        disp(['rec' num2str(i)])
        fileNamePeakFR = [filenames(i,:) '_PeakFR_msess' num2str(mazeSess(i)) ...
                        '_RunOnset1.mat'];
        fullPath = [paths(i,:) fileNamePeakFR];
        if(exist(fullPath) == 0)
            disp(['The peak firing rate file does not exist. Please call ',...
                    'function "PeakFiringRate_Aligned" first.']);
            return;
        end
        load(fullPath,'pFRNonStimGoodStruct');
        
        fileNameConv = [filenames(i,:) '_convSpikesAligned_msess' num2str(mazeSess(i)) '_BefRun1.mat'];
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
        load(fullPath,'trialNoNonStimGood','trialNoNonStimBad');
        
        fileNameInfo = [filenames(i,:) '_Info.mat'];
        fullPath = [paths(i,:) fileNameInfo];
        if(exist(fullPath) == 0)
            disp('_Info.mat file does not exist.');
            return;
        end
        load(fullPath,'rec');  % changed Dinghao 11 Dec 2023
        
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
                      rec.isIntern ~= 1);
        
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
        else
            disp(['Recording ' filenames(i,:) ' has ' num2str(length(trialNoNonStimGood)) ' good trials']);
        end
    end
    
end