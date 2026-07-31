function InterneuronInitPeakSpeedAllRec()

    methodTheta = 1;
    minFRInt = 3;
        
    RecordingList;
    pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\interneuron\';
    
    sampleFq = 1250;
    nSampBef = 3*sampleFq; % used in align to running onset
    
    load([pathAnal 'autoCorrIntAllRec.mat']);
    if(exist([pathAnal 'initPeakIntSpeedAllRec.mat']))
        load([pathAnal 'initPeakIntSpeedAllRec.mat']);
    end
    if(exist([pathAnal 'initPeakIntAllRec.mat']))
        load([pathAnal 'initPeakIntAllRec.mat']);
    end
    
%     disp('No cue - firing rate vs. speed')
%     modInt2NoCue = accumInterneurons3(listRecordingsNoCuePath,...
%         listRecordingsNoCueFileName,mazeSessionNoCue,minFRInt,1,sampleFq,nSampBef);
%     
%     disp('Active licking - firing rate vs. speed')
%     modInt2AL = accumInterneurons3(listRecordingsActiveLickPath,...
%         listRecordingsActiveLickFileName,mazeSessionActiveLick,minFRInt,2,sampleFq,nSampBef);
%     
%     disp('Passive licking - firing rate vs. speed')
%     modInt2PL = accumInterneurons3(listRecordingsPassiveLickPath,...
%         listRecordingsPassiveLickFileName,mazeSessionPassiveLick,minFRInt,3,sampleFq,nSampBef);
%     
%     save([pathAnal 'initPeakIntSpeedAllRec.mat'],'modInt2NoCue','modInt2AL','modInt2PL','-v7.3'); 
    
    modInt2.corrInstFRRunSpeedGood = [modInt2NoCue.corrInstFRRunSpeedGood,...
        modInt2AL.corrInstFRRunSpeedGood,modInt2PL.corrInstFRRunSpeedGood];
    modInt2.corrInstFRRunSpeedBLGood = [modInt2NoCue.corrInstFRRunSpeedBLGood,...
        modInt2AL.corrInstFRRunSpeedBLGood,modInt2PL.corrInstFRRunSpeedBLGood];
    modInt2.corrInstFRRunSpeedBefRunGood = [modInt2NoCue.corrInstFRRunSpeedBefRunGood,...
        modInt2AL.corrInstFRRunSpeedBefRunGood,modInt2PL.corrInstFRRunSpeedBefRunGood];
    modInt2.corrInstFRRunSpeed0to1Good = [modInt2NoCue.corrInstFRRunSpeed0to1Good,...
        modInt2AL.corrInstFRRunSpeed0to1Good,modInt2PL.corrInstFRRunSpeed0to1Good];
    modInt2.corrInstFRRunSpeed3to5Good = [modInt2NoCue.corrInstFRRunSpeed3to5Good,...
        modInt2AL.corrInstFRRunSpeed3to5Good,modInt2PL.corrInstFRRunSpeed3to5Good];
    modInt2.corrStopTimeRunSpeedGood = [modInt2NoCue.corrStopTimeRunSpeedGood,...
        modInt2AL.corrStopTimeRunSpeedGood,modInt2PL.corrStopTimeRunSpeedGood];
    modInt2.corrStopTimeInstFRBefRunGood = [modInt2NoCue.corrStopTimeInstFRBefRunGood,...
        modInt2AL.corrStopTimeInstFRBefRunGood,modInt2PL.corrStopTimeInstFRBefRunGood];
    modInt2.corrStopTimeInstFR0to1Good = [modInt2NoCue.corrStopTimeInstFR0to1Good,...
        modInt2AL.corrStopTimeInstFR0to1Good,modInt2PL.corrStopTimeInstFR0to1Good];
    
    modInt2.pCorrInstFRRunSpeedGood = [modInt2NoCue.pCorrInstFRRunSpeedGood,...
        modInt2AL.pCorrInstFRRunSpeedGood,modInt2PL.pCorrInstFRRunSpeedGood];
    modInt2.pCorrInstFRRunSpeedBLGood = [modInt2NoCue.pCorrInstFRRunSpeedBLGood,...
        modInt2AL.pCorrInstFRRunSpeedBLGood,modInt2PL.pCorrInstFRRunSpeedBLGood];
    modInt2.pCorrInstFRRunSpeedBefRunGood = [modInt2NoCue.pCorrInstFRRunSpeedBefRunGood,...
        modInt2AL.pCorrInstFRRunSpeedBefRunGood,modInt2PL.pCorrInstFRRunSpeedBefRunGood];
    modInt2.pCorrInstFRRunSpeed0to1Good = [modInt2NoCue.pCorrInstFRRunSpeed0to1Good,...
        modInt2AL.pCorrInstFRRunSpeed0to1Good,modInt2PL.pCorrInstFRRunSpeed0to1Good];
    modInt2.pCorrInstFRRunSpeed3to5Good = [modInt2NoCue.pCorrInstFRRunSpeed3to5Good,...
        modInt2AL.pCorrInstFRRunSpeed3to5Good,modInt2PL.pCorrInstFRRunSpeed3to5Good];
    modInt2.pCorrStopTimeRunSpeedGood = [modInt2NoCue.pCorrStopTimeRunSpeedGood,...
        modInt2AL.pCorrStopTimeRunSpeedGood,modInt2PL.pCorrStopTimeRunSpeedGood];
    modInt2.pCorrStopTimeInstFRBefRunGood = [modInt2NoCue.pCorrStopTimeInstFRBefRunGood,...
        modInt2AL.pCorrStopTimeInstFRBefRunGood,modInt2PL.pCorrStopTimeInstFRBefRunGood];
    modInt2.pCorrStopTimeInstFR0to1Good = [modInt2NoCue.pCorrStopTimeInstFR0to1Good,...
        modInt2AL.pCorrStopTimeInstFR0to1Good,modInt2PL.pCorrStopTimeInstFR0to1Good];
    
    modInt2.lmSlopeInstFRRunSpeedGood = [modInt2NoCue.lmSlopeInstFRRunSpeedGood,...
        modInt2AL.lmSlopeInstFRRunSpeedGood,modInt2PL.lmSlopeInstFRRunSpeedGood];
    modInt2.lmSlopeInstFRRunSpeedBefRunGood = [modInt2NoCue.lmSlopeInstFRRunSpeedBefRunGood,...
        modInt2AL.lmSlopeInstFRRunSpeedBefRunGood,modInt2PL.lmSlopeInstFRRunSpeedBefRunGood];
    modInt2.lmSlopeInstFRRunSpeed0to1Good = [modInt2NoCue.lmSlopeInstFRRunSpeed0to1Good,...
        modInt2AL.lmSlopeInstFRRunSpeed0to1Good,modInt2PL.lmSlopeInstFRRunSpeed0to1Good];
    modInt2.lmSlopeInstFRRunSpeed3to5Good = [modInt2NoCue.lmSlopeInstFRRunSpeed3to5Good,...
        modInt2AL.lmSlopeInstFRRunSpeed3to5Good,modInt2PL.lmSlopeInstFRRunSpeed3to5Good];
    
    modInt2.lmPInstFRRunSpeedGood = [modInt2NoCue.lmPInstFRRunSpeedGood,...
        modInt2AL.lmPInstFRRunSpeedGood,modInt2PL.lmPInstFRRunSpeedGood];
    modInt2.lmPInstFRRunSpeedBefRunGood = [modInt2NoCue.lmPInstFRRunSpeedBefRunGood,...
        modInt2AL.lmPInstFRRunSpeedBefRunGood,modInt2PL.lmPInstFRRunSpeedBefRunGood];
    modInt2.lmPInstFRRunSpeed0to1Good = [modInt2NoCue.lmPInstFRRunSpeed0to1Good,...
        modInt2AL.lmPInstFRRunSpeed0to1Good,modInt2PL.lmPInstFRRunSpeed0to1Good];
    modInt2.lmPInstFRRunSpeed3to5Good = [modInt2NoCue.lmPInstFRRunSpeed3to5Good,...
        modInt2AL.lmPInstFRRunSpeed3to5Good,modInt2PL.lmPInstFRRunSpeed3to5Good];
    
    modInt2.lmRInstFRRunSpeedGood = [modInt2NoCue.lmRInstFRRunSpeedGood,...
        modInt2AL.lmRInstFRRunSpeedGood,modInt2PL.lmRInstFRRunSpeedGood];
    modInt2.lmRInstFRRunSpeedBefRunGood = [modInt2NoCue.lmRInstFRRunSpeedBefRunGood,...
        modInt2AL.lmRInstFRRunSpeedBefRunGood,modInt2PL.lmRInstFRRunSpeedBefRunGood];
    modInt2.lmRInstFRRunSpeed0to1Good = [modInt2NoCue.lmRInstFRRunSpeed0to1Good,...
        modInt2AL.lmRInstFRRunSpeed0to1Good,modInt2PL.lmRInstFRRunSpeed0to1Good];
    modInt2.lmRInstFRRunSpeed3to5Good = [modInt2NoCue.lmRInstFRRunSpeed3to5Good,...
        modInt2AL.lmRInstFRRunSpeed3to5Good,modInt2PL.lmRInstFRRunSpeed3to5Good];
    
    nNeuWithField = [modIntNoCue.nNeuWithField modIntAL.nNeuWithField modIntPL.nNeuWithField];
    idxC2 = [modInt1NoCue.idxC2 modInt1AL.idxC2 modInt1PL.idxC2];
    for i = 1:max(idxC2)
%         indCur = idxC1' == i & nNeuWithField >= 2;
        indCur = idxC2' == i;
        corrStopTimeRunSpeedGoodC = modInt2.corrStopTimeRunSpeedGood(indCur);
        corrStopTimeInstFRBefRunGoodC = modInt2.corrStopTimeInstFRBefRunGood(indCur);
        corrStopTimeInstFR0to1GoodC = modInt2.corrStopTimeInstFR0to1Good(indCur);
        
        pCorrStopTimeRunSpeedGoodC = modInt2.pCorrStopTimeRunSpeedGood(indCur);
        pCorrStopTimeInstFRBefRunGoodC = modInt2.pCorrStopTimeInstFRBefRunGood(indCur);
        pCorrStopTimeInstFR0to1GoodC = modInt2.pCorrStopTimeInstFR0to1Good(indCur);
        
        meanCorrStopTimeRunSpeedGoodC(i) = mean(corrStopTimeRunSpeedGoodC);
        meanCorrStopTimeInstFRBefRunGoodC(i) = mean(corrStopTimeInstFRBefRunGoodC);
        meanCorrStopTimeInstFR0to1GoodC(i) = mean(corrStopTimeInstFR0to1GoodC);
        
        percSigCorrStopTimeRunSpeedGoodC(i) = sum(pCorrStopTimeRunSpeedGoodC < 0.05)/sum(indCur);
        percSigCorrStopTimeInstFRBefRunGoodC(i) = sum(pCorrStopTimeInstFRBefRunGoodC < 0.05)/sum(indCur);
        percSigCorrStopTimeInstFR0to1GoodC(i) = sum(pCorrStopTimeInstFR0to1GoodC < 0.05)/sum(indCur);
        
        corrInstFRRunSpeedGoodC = modInt2.corrInstFRRunSpeedGood(indCur);
        corrInstFRRunSpeedBLGoodC = modInt2.corrInstFRRunSpeedBLGood(indCur);
        corrInstFRRunSpeedBefRunGoodC = modInt2.corrInstFRRunSpeedBefRunGood(indCur);
        corrInstFRRunSpeed0to1GoodC = modInt2.corrInstFRRunSpeed0to1Good(indCur);
        corrInstFRRunSpeed3to5GoodC = modInt2.corrInstFRRunSpeed3to5Good(indCur);
        
        pCorrInstFRRunSpeedGoodC = modInt2.pCorrInstFRRunSpeedGood(indCur);
        pCorrInstFRRunSpeedBLGoodC = modInt2.pCorrInstFRRunSpeedBLGood(indCur);
        pCorrInstFRRunSpeedBefRunGoodC = modInt2.pCorrInstFRRunSpeedBefRunGood(indCur);
        pCorrInstFRRunSpeed0to1GoodC = modInt2.pCorrInstFRRunSpeed0to1Good(indCur);
        pCorrInstFRRunSpeed3to5GoodC = modInt2.pCorrInstFRRunSpeed3to5Good(indCur);
        
        meanCorrInstFRRunSpeedGoodC(i) = mean(corrInstFRRunSpeedGoodC);
        meanCorrInstFRRunSpeedBLGoodC(i) = mean(corrInstFRRunSpeedBLGoodC);
        meanCorrInstFRRunSpeedBefRunGoodC(i) = mean(corrInstFRRunSpeedBefRunGoodC);
        meanCorrInstFRRunSpeed0to1GoodC(i) = mean(corrInstFRRunSpeed0to1GoodC);
        meanCorrInstFRRunSpeed3to5GoodC(i) = mean(corrInstFRRunSpeed3to5GoodC);
        
        percSigCorrInstFRRunSpeedGoodC(i) = sum(pCorrInstFRRunSpeedGoodC < 0.05)/sum(indCur);
        percSigCorrInstFRRunSpeedBLGoodC(i) = sum(pCorrInstFRRunSpeedBLGoodC < 0.05)/sum(indCur);
        percSigCorrInstFRRunSpeedBefRunGoodC(i) = sum(pCorrInstFRRunSpeedBefRunGoodC < 0.05)/sum(indCur);
        percSigCorrInstFRRunSpeed0to1GoodC(i) = sum(pCorrInstFRRunSpeed0to1GoodC < 0.05)/sum(indCur);
        percSigCorrInstFRRunSpeed3to5GoodC(i) = sum(pCorrInstFRRunSpeed3to5GoodC < 0.05)/sum(indCur);
        
        lmSlopeInstFRRunSpeedGoodC = modInt2.lmSlopeInstFRRunSpeedGood(indCur);
        lmSlopeInstFRRunSpeedBefRunGoodC = modInt2.lmSlopeInstFRRunSpeedBefRunGood(indCur);
        lmSlopeInstFRRunSpeed0to1GoodC = modInt2.lmSlopeInstFRRunSpeed0to1Good(indCur);
        lmSlopeInstFRRunSpeed3to5GoodC = modInt2.lmSlopeInstFRRunSpeed3to5Good(indCur);
        
        lmPInstFRRunSpeedGoodC = modInt2.lmPInstFRRunSpeedGood(indCur);
        lmPInstFRRunSpeedBefRunGoodC = modInt2.lmPInstFRRunSpeedBefRunGood(indCur);
        lmPInstFRRunSpeed0to1GoodC = modInt2.lmPInstFRRunSpeed0to1Good(indCur);
        lmPInstFRRunSpeed3to5GoodC = modInt2.lmPInstFRRunSpeed3to5Good(indCur);
        
        lmRInstFRRunSpeedGoodC = modInt2.lmRInstFRRunSpeedGood(indCur);
        lmRInstFRRunSpeedBefRunGoodC = modInt2.lmRInstFRRunSpeedBefRunGood(indCur);
        lmRInstFRRunSpeed0to1GoodC = modInt2.lmRInstFRRunSpeed0to1Good(indCur);
        lmRInstFRRunSpeed3to5GoodC = modInt2.lmRInstFRRunSpeed3to5Good(indCur);
    end
        
end

function modInt2 = accumInterneurons3(paths,filenames,mazeSess,minFRInt,task,sampleFq,nSampBef)
    numRec = size(paths,1);
    modInt2 = struct('task',[],... % no cue - 1, AL - 2, PL - 3
                              'indRec',[],... % recording index
                              'indNeu',[],... % neuron indices
                              'stopTimePerTrGood',[],... % stop time for good trials
                              'stopTimePerTrBad',[],... % stop time for bad trials
                              ...
                              'timeStepRun',[],...                              
                              'filteredSpArrGood',[],... % firing profile per trial for good trials
                              'meanInstFRTrGood',[],... % mean inst FR per trial for good trials
                              'meanInstFRBLTrGood',[],... % mean inst FR per trial during baseline period
                              'meanInstFRBefRunTrGood',[],... % mean inst FR per trial before run 
                              'meanInstFR0to1TrGood',[],... % mean inst FR per trial between 0-1 s after running onset
                              'meanInstFR3to5TrGood',[],... % mean inst FR per trial between 3-5 s after running onset
                              ...
                              'filteredSpArrBad',[],... % firing profile per trial for bad trials
                              'meanInstFRTrBad',[],... % mean inst FR per trial for bad trials
                              'meanInstFRBLTrBad',[],... % mean inst FR per trial during baseline period
                              'meanInstFRBefRunTrBad',[],... % mean inst FR per trial before run 
                              'meanInstFR0to1TrBad',[],... % mean inst FR per trial between 0-1 s after running onset
                              'meanInstFR3to5TrBad',[],... % mean inst FR per trial between 3-5 s after running onset
                              ...
                              'timeStepRunSpeed',[],...
                              'runSpeedGood',[],... % running speed profile for good trials
                              'meanRunSpeedTrGood',[], ... % mean running speed per trial
                              'meanRunSpeedBLTrGood',[], ... % mean running speed per trial during baseline period                              
                              'meanRunSpeedBefRunTrGood',[], ... % mean running speed per trial before run
                              'meanRunSpeed0to1TrGood',[], ... % mean running speed per trial between 0-1 s after running onset
                              'meanRunSpeed3to5TrGood',[], ... % mean running speed per trial between 3-5 s after running onset
                              ...
                              'runSpeedBad',[],... % running speed profile for bad trials
                              'meanRunSpeedTrBad',[], ... % mean running speed per trial
                              'meanRunSpeedBLTrBad',[], ... % mean running speed per trial during baseline period                              
                              'meanRunSpeedBefRunTrBad',[], ... % mean running speed per trial before run
                              'meanRunSpeed0to1TrBad',[], ... % mean running speed per trial between 0-1 s after running onset
                              'meanRunSpeed3to5TrBad',[]); % mean running speed per trial between 3-5 s after running onset
                             
    totExcNeu = 0;
    curNoNeu = 1;
    curNoNeuBad = 1;
    for i = 1:numRec
        disp(i);
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
        
        fileNameBeh = [filenames(i,:) '_behPar_msess' num2str(mazeSess(i)) '.mat'];
        fullPath = [paths(i,:) fileNameBeh];
        if(exist(fullPath) == 0)
            disp('_behPar.mat file does not exist.');
            return;
        end
        load(fullPath,'behPar');
        
        fileNamePeakFR = [filenames(i,:) '_PeakFR_msess' num2str(mazeSess(i)) ...
                        '_RunOnset0.mat'];
        fullPath = [paths(i,:) fileNamePeakFR];
        if(exist(fullPath) == 0)
            disp(['The peak firing rate file does not exist. Please call ',...
                    'function "PeakFiringRate_Aligned" first.']);
            return;
        end
        load(fullPath,'pFRNonStimGoodStruct','pFRNonStimBadStruct');
        
        fullPath = [paths(i,:) filenames(i,:) '_PeakFRAligned_msess' num2str(mazeSess(i)) '_Run0.mat'];    
        if(exist(fullPath) == 0)
            disp('The _PeakFRAligned file does not exist');
            return;
        end
        load(fullPath,'trialNoNonStimGood','trialNoNonStimBad');
        
        fileNameConv = [filenames(i,:) '_convSpikesAligned_msess' num2str(mazeSess(i)) '_BefRun0.mat'];
        fullPath = [paths(i,:) fileNameConv];
        if(exist(fullPath) == 0)
            disp(['The convSpikesAligned file does not exist. Please call ',...
                    'function "ConvSpikeTrain_AlignedRunOnset" first.']);
            return;
        end
        load(fullPath,'timeStepRun');
        
        fileNameThetaPh = [filenames(i,:) '_thetaPhaseOverTimeligned_msess' num2str(mazeSess(i)) '_Run0.mat'];
        fullPath = [paths(i,:) fileNameThetaPh];
        if(exist(fullPath) == 0)
            disp(['The thetaPhaseOverTimeligned file does not exist. Please call ',...
                    'function "thetaPhaseOverTimeAligned" first.']);
            return;
        end
        load(fullPath,'runSpeedNonStimGood','runSpeedNonStimBad','paramC');
        
        fileNameFR = [filenames(i,:) '_FR_Run1.mat'];
        fullPath = [paths(i,:) fileNameFR];
        if(exist(fullPath) == 0)
            disp('_FR file does not exist.');
            return;
        end
        load(fullPath,'mFRStruct','mFRStructSess');
        if(length(beh.mazeSessAll) > 1)
            mFR = mFRStructSess{mazeSess(i)};
        else
            mFR = mFRStruct;
        end
        
        fileNameFW = [filenames(i,:) '_FieldSpCorr_GoodTr_Run1.mat'];
        fullPath = [paths(i,:) fileNameFW];
        if(exist(fullPath) == 0)
            disp(['The field detection file does not exist. Please call ',...
                    'function "FieldDetection_GoodTr" first.']);
            return;
        end
        load(fullPath,'paramF'); 
                       
%         indNeu = autoCorr.isInterneuron == 1 & cluList.firingRate > minFRInt;
        indNeu = autoCorr.isInterneuron == 1 & mFR.mFR > minFRInt;
                
        if(length(pFRNonStimGoodStruct.indLapList) >= paramF.minNumTr)
            modInt2.task = [modInt2.task task*ones(1,sum(indNeu))];
            modInt2.indRec = [modInt2.indRec i*ones(1,sum(indNeu))];
            indNeu1 = find(indNeu == 1);
            modInt2.indNeu = [modInt2.indNeu indNeu1]; 
        
            modInt2.timeStepRun = timeStepRun/sampleFq;  
            modInt2.timeStepRunSpeed = paramC.timeIncBef;
            
            for j = 1:length(indNeu1)           
                timeStep = modInt2.timeStepRun; 
                indLapsGood = pFRNonStimGoodStruct.indLapList~=1;
                filteredSpArr = pFRNonStimGoodStruct.filteredSpArrAll{indNeu1(j)}(indLapsGood,:);
                trialNoGood = trialNoNonStimGood(trialNoNonStimGood~=1);
                modInt2.stopTimePerTrGood{curNoNeu} = behPar.rewardToRun(trialNoGood)';
                modInt2.filteredSpArrGood{curNoNeu} = filteredSpArr;
                modInt2.meanInstFRTrGood{curNoNeu} = mean(filteredSpArr,2);
                indFRBaseline = find(timeStep >= -3 & timeStep < -2);
                modInt2.meanInstFRBLTrGood{curNoNeu} = mean(filteredSpArr(:,indFRBaseline),2);
                indFRBefRun = find(timeStep >= -0.5 & timeStep < 0);
                modInt2.meanInstFRBefRunTrGood{curNoNeu} = mean(filteredSpArr(:,indFRBefRun),2);
                indFR0to1 = find(timeStep >= 0 & timeStep < 1);
                modInt2.meanInstFR0to1TrGood{curNoNeu} = mean(filteredSpArr(:,indFR0to1),2);
                indFR3to5 = find(timeStep >= 3 & timeStep < 5);
                modInt2.meanInstFR3to5TrGood{curNoNeu} = mean(filteredSpArr(:,indFR3to5),2);

                timeStep = modInt2.timeStepRunSpeed;
                runSpeedNonStimGood1 = runSpeedNonStimGood(indLapsGood,:);
                modInt2.runSpeedGood{curNoNeu} = runSpeedNonStimGood1;
                modInt2.meanRunSpeedTrGood{curNoNeu} = mean(runSpeedNonStimGood1(:,nSampBef+1:end),2);
                indFRBaseline = find(timeStep >= -3 & timeStep < -2);
                modInt2.meanRunSpeedBLTrGood{curNoNeu} = mean(runSpeedNonStimGood1(:,indFRBaseline),2);
                indFRBefRun = find(timeStep >= -0.5 & timeStep < 0);
                modInt2.meanRunSpeedBefRunTrGood{curNoNeu} = mean(runSpeedNonStimGood1(:,indFRBefRun),2);
                indFR0to1 = find(timeStep >= 0 & timeStep < 1);
                modInt2.meanRunSpeed0to1TrGood{curNoNeu} = mean(runSpeedNonStimGood1(:,indFR0to1),2);
                indFR3to5 = find(timeStep >= 3 & timeStep < 5);
                modInt2.meanRunSpeed3to5TrGood{curNoNeu} = mean(runSpeedNonStimGood1(:,indFR3to5),2);
                
                [modInt2.corrStopTimeRunSpeedGood(curNoNeu),modInt2.pCorrStopTimeRunSpeedGood(curNoNeu)] = ...
                    corr(modInt2.stopTimePerTrGood{curNoNeu},modInt2.meanRunSpeedTrGood{curNoNeu});
                [modInt2.corrStopTimeRunSpeed0to1Good(curNoNeu),modInt2.pCorrStopTimeRunSpeed0to1Good(curNoNeu)] = ...
                    corr(modInt2.stopTimePerTrGood{curNoNeu},modInt2.meanRunSpeed0to1TrGood{curNoNeu});
                [modInt2.corrStopTimeRunSpeed3to5Good(curNoNeu),modInt2.pCorrStopTimeRunSpeed3to5Good(curNoNeu)] = ...
                    corr(modInt2.stopTimePerTrGood{curNoNeu},modInt2.meanRunSpeed3to5TrGood{curNoNeu});
                [modInt2.corrStopTimeInstFRBefRunGood(curNoNeu),modInt2.pCorrStopTimeInstFRBefRunGood(curNoNeu)] = ...
                    corr(modInt2.stopTimePerTrGood{curNoNeu},modInt2.meanInstFRBefRunTrGood{curNoNeu});
                [modInt2.corrStopTimeInstFR0to1Good(curNoNeu),modInt2.pCorrStopTimeInstFR0to1Good(curNoNeu)] = ...
                    corr(modInt2.stopTimePerTrGood{curNoNeu},modInt2.meanInstFR0to1TrGood{curNoNeu});

                [modInt2.corrInstFRRunSpeedGood(curNoNeu),modInt2.pCorrInstFRRunSpeedGood(curNoNeu)] = ...
                    corr(modInt2.meanInstFRTrGood{curNoNeu},modInt2.meanRunSpeedTrGood{curNoNeu});
                [modInt2.corrInstFRRunSpeedBLGood(curNoNeu),modInt2.pCorrInstFRRunSpeedBLGood(curNoNeu)] = ...
                    corr(modInt2.meanRunSpeedBLTrGood{curNoNeu},modInt2.meanInstFRBLTrGood{curNoNeu});
                [modInt2.corrInstFRRunSpeedBefRunGood(curNoNeu),modInt2.pCorrInstFRRunSpeedBefRunGood(curNoNeu)] = ...
                    corr(modInt2.meanRunSpeedBefRunTrGood{curNoNeu},modInt2.meanInstFRBefRunTrGood{curNoNeu});
                [modInt2.corrInstFRRunSpeed0to1Good(curNoNeu),modInt2.pCorrInstFRRunSpeed0to1Good(curNoNeu)] = ...
                    corr(modInt2.meanRunSpeed0to1TrGood{curNoNeu},modInt2.meanInstFR0to1TrGood{curNoNeu});
                [modInt2.corrInstFRRunSpeed3to5Good(curNoNeu),modInt2.pCorrInstFRRunSpeed3to5Good(curNoNeu)] = ...
                    corr(modInt2.meanRunSpeed3to5TrGood{curNoNeu},modInt2.meanInstFR3to5TrGood{curNoNeu});

                md = fitlm(modInt2.meanRunSpeedTrGood{curNoNeu},modInt2.meanInstFRTrGood{curNoNeu},'linear');
                modInt2.lmSlopeInstFRRunSpeedGood(curNoNeu) = md.Coefficients.Estimate(2);
                modInt2.lmPInstFRRunSpeedGood(curNoNeu) = md.Coefficients.pValue(2);
                modInt2.lmRInstFRRunSpeedGood(curNoNeu) = md.Rsquared.Ordinary;

                md = fitlm(modInt2.meanRunSpeedBefRunTrGood{curNoNeu},modInt2.meanInstFRBefRunTrGood{curNoNeu},'linear');
                modInt2.lmSlopeInstFRRunSpeedBefRunGood(curNoNeu) = md.Coefficients.Estimate(2);
                modInt2.lmPInstFRRunSpeedBefRunGood(curNoNeu) = md.Coefficients.pValue(2);
                modInt2.lmRInstFRRunSpeedBefRunGood(curNoNeu) = md.Rsquared.Ordinary;

                md = fitlm(modInt2.meanRunSpeed0to1TrGood{curNoNeu},modInt2.meanInstFR0to1TrGood{curNoNeu},'linear');
                modInt2.lmSlopeInstFRRunSpeed0to1Good(curNoNeu) = md.Coefficients.Estimate(2);
                modInt2.lmPInstFRRunSpeed0to1Good(curNoNeu) = md.Coefficients.pValue(2);
                modInt2.lmRInstFRRunSpeed0to1Good(curNoNeu) = md.Rsquared.Ordinary;

                md = fitlm(modInt2.meanRunSpeed3to5TrGood{curNoNeu},modInt2.meanInstFR0to1TrGood{curNoNeu},'linear');
                modInt2.lmSlopeInstFRRunSpeed3to5Good(curNoNeu) = md.Coefficients.Estimate(2);
                modInt2.lmPInstFRRunSpeed3to5Good(curNoNeu) = md.Coefficients.pValue(2);
                modInt2.lmRInstFRRunSpeed3to5Good(curNoNeu) = md.Rsquared.Ordinary;
                curNoNeu = curNoNeu+1;
            end
        else
            disp([filenames(i,:) ' only has ' num2str(length(pFRNonStimGoodStruct.indLapList)) ...
                ' good trials.']);
            disp(['No. interneurons in this recording is ' num2str(sum(indNeu))]);
        end  

        if(~isempty(pFRNonStimBadStruct) && length(pFRNonStimBadStruct.indLapList) >= paramF.minNumTr)
            indNeu1 = find(indNeu == 1);
            for j = 1:length(indNeu1)
                indLapsBad = [];
                if(isempty(pFRNonStimBadStruct))
                    modInt2.filteredSpArrBad{curNoNeuBad} = [];
                    modInt2.stopTimePerTrBad{curNoNeuBad} = [];
                    modInt2.meanInstFRTrBad{curNoNeuBad} = [];
                    modInt2.meanInstFRBLTrBad{curNoNeuBad} = [];
                    modInt2.meanInstFRBefRunTrBad{curNoNeuBad} = [];
                    modInt2.meanInstFR0to1TrBad{curNoNeuBad} = [];
                    modInt2.meanInstFR3to5TrBad{curNoNeuBad} = [];
                else
                    if(i == 58)
                        a = 1;
                    end
                    timeStep = modInt2.timeStepRun; 
                    indLapsBad = pFRNonStimBadStruct.indLapList~=1;
                    filteredSpArr = pFRNonStimBadStruct.filteredSpArrAll{indNeu1(j)}(indLapsBad,:);
                    modInt2.filteredSpArrBad{curNoNeuBad} = filteredSpArr;
                    trialNoBad = trialNoNonStimBad(trialNoNonStimBad~=1);
                    modInt2.stopTimePerTrBad{curNoNeuBad} = behPar.rewardToRun(trialNoBad)';
                    modInt2.meanInstFRTrBad{curNoNeuBad} = mean(filteredSpArr,2);                    
                    indFRBaseline = find(timeStep >= -3 & timeStep < -2);
                    modInt2.meanInstFRBLTrBad{curNoNeuBad} = mean(filteredSpArr(:,indFRBaseline),2);
                    indFRBefRun = find(timeStep >= -0.5 & timeStep < 0);
                    modInt2.meanInstFRBefRunTrBad{curNoNeuBad} = mean(filteredSpArr(:,indFRBefRun),2);
                    indFR0to1 = find(timeStep >= 0 & timeStep < 1);
                    modInt2.meanInstFR0to1TrBad{curNoNeuBad} = mean(filteredSpArr(:,indFR0to1),2);
                    indFR3to5 = find(timeStep >= 3 & timeStep < 5);
                    modInt2.meanInstFR3to5TrBad{curNoNeuBad} = mean(filteredSpArr(:,indFR3to5),2);
                end
                
                runSpeedNonStimBad1 = runSpeedNonStimBad(indLapsBad,:);
                modInt2.runSpeedBad{curNoNeuBad} = runSpeedNonStimBad1;
                modInt2.meanRunSpeedTrBad{curNoNeuBad} = mean(runSpeedNonStimBad1,2);
                modInt2.meanRunSpeedBLTrBad{curNoNeuBad} = mean(runSpeedNonStimBad1(:,indFRBaseline),2);                
                modInt2.meanRunSpeedBefRunTrBad{curNoNeuBad} = mean(runSpeedNonStimBad1(:,indFRBefRun),2);
                modInt2.meanRunSpeed0to1TrBad{curNoNeuBad} = mean(runSpeedNonStimBad1(:,indFR0to1),2);
                modInt2.meanRunSpeed3to5TrBad{curNoNeuBad} = mean(runSpeedNonStimBad1(:,indFR3to5),2);
                
                curNoNeuBad = curNoNeuBad+1;
            end
        else
            totExcNeu = totExcNeu + sum(indNeu);
            if(isempty(pFRNonStimBadStruct))
                lenTr = 0;
            else
                lenTr = length(pFRNonStimBadStruct.indLapList);
            end
            disp([filenames(i,:) ' only has ' num2str(lenTr) ' bad trials.']);
            disp(['No. interneurons in this recording is ' num2str(sum(indNeu)) ...
                ', total number of excluded neurons is ' num2str(totExcNeu)]);
        end  
    end
end
