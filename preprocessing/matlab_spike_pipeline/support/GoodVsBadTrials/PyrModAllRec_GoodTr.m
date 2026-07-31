function PyrModAllRec_GoodTr(onlyRun,methodKMean)

    methodTheta = 1;
    minFR = 0.15;
    maxFR = 7;
    spaceBin = 20;
    if(nargin == 1)
        methodKMean = 2; % which kmean method is used
    end
    
    RecordingList;
    pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\Pyramidal\';
    
    if(exist([pathAnal 'autoCorrPyrAllRec.mat']))
        load([pathAnal 'autoCorrPyrAllRec.mat']);
    end
    
    if(exist([pathAnal 'autoCorrPyrAllRec_GoodTr.mat']))
        load([pathAnal 'autoCorrPyrAllRec_GoodTr.mat']);
    end
    
    %% pyramidal neurons in no cue passive task
%     disp('No cue')
%     modPyrNoCue_GoodTr = accumPyrNeurons1(listRecordingsNoCuePath,...
%         listRecordingsNoCueFileName,mazeSessionNoCue,minFR,maxFR,1,methodTheta,onlyRun,1,spaceBin);
%     modPyrNoCue_BadTr = accumPyrNeurons1(listRecordingsNoCuePath,...
%         listRecordingsNoCueFileName,mazeSessionNoCue,minFR,maxFR,1,methodTheta,onlyRun,0,spaceBin);
%     
%     disp('Active licking')
%     modPyrAL_GoodTr = accumPyrNeurons1(listRecordingsActiveLickPath,...
%         listRecordingsActiveLickFileName,mazeSessionActiveLick,minFR,maxFR,2,methodTheta,onlyRun,1,spaceBin);
%     modPyrAL_BadTr = accumPyrNeurons1(listRecordingsActiveLickPath,...
%         listRecordingsActiveLickFileName,mazeSessionActiveLick,minFR,maxFR,2,methodTheta,onlyRun,0,spaceBin);
%     
%     disp('Passive licking')
%     modPyrPL_GoodTr = accumPyrNeurons1(listRecordingsPassiveLickPath,...
%         listRecordingsPassiveLickFileName,mazeSessionPassiveLick,minFR,maxFR,3,methodTheta,onlyRun,1,spaceBin);
%     modPyrPL_BadTr = accumPyrNeurons1(listRecordingsPassiveLickPath,...
%         listRecordingsPassiveLickFileName,mazeSessionPassiveLick,minFR,maxFR,3,methodTheta,onlyRun,0,spaceBin);
%     
%     if(exist([pathAnal 'autoCorrPyrAllRec_GoodTr.mat']))
%         save([pathAnal 'autoCorrPyrAllRec_GoodTr.mat'],'modPyrNoCue_GoodTr','modPyrAL_GoodTr',...
%             'modPyrPL_GoodTr','modPyrNoCue_BadTr','modPyrAL_BadTr','modPyrPL_BadTr','-append'); 
%     else
%         save([pathAnal 'autoCorrPyrAllRec_GoodTr.mat'],'modPyrNoCue_GoodTr','modPyrAL_GoodTr',...
%             'modPyrPL_GoodTr','modPyrNoCue_BadTr','modPyrAL_BadTr','modPyrPL_BadTr');
%     end
    
    mod_GoodTr = modPyrAllCondi(modPyrNoCue_GoodTr, modPyrPL_GoodTr, modPyrAL_GoodTr);
    mod_BadTr = modPyrAllCondi(modPyrNoCue_BadTr, modPyrPL_BadTr, modPyrAL_BadTr);
    
    if(methodKMean == 1)
        idxC = [autoCorrPyrAll.idxC1(autoCorrPyrAll.task == autoCorrPyrNoCue.task(1))' ...
            autoCorrPyrAll.idxC1(autoCorrPyrAll.task == autoCorrPyrAL.task(1))' ...
            autoCorrPyrAll.idxC1(autoCorrPyrAll.task == autoCorrPyrPL.task(1))'];
    elseif(methodKMean == 2)
        idxC = [autoCorrPyrAll.idxC2(autoCorrPyrAll.task == autoCorrPyrNoCue.task(1))' ...
            autoCorrPyrAll.idxC2(autoCorrPyrAll.task == autoCorrPyrAL.task(1))' ...
            autoCorrPyrAll.idxC2(autoCorrPyrAll.task == autoCorrPyrPL.task(1))'];
    elseif(methodKMean == 3)
        idxC = [autoCorrPyrAll.idxC3(autoCorrPyrAll.task == autoCorrPyrNoCue.task(1)) ...
            autoCorrPyrAll.idxC3(autoCorrPyrAll.task == autoCorrPyrAL.task(1)) ...
            autoCorrPyrAll.idxC3(autoCorrPyrAll.task == autoCorrPyrPL.task(1))];
    end
    mod_GoodTr.idxC = idxC;
    mod_BadTr.idxC = idxC;

    indFGood = mod_GoodTr.isNeuWithField == 1;
    plotCompXY(mod_BadTr.meanCorrDistNZ(indFGood),mod_GoodTr.meanCorrDistNZ(indFGood),...
        'Pyr_GoodVsBadTrMeanCorrDist_GoodTrField',pathAnal,'MeanCorrDist BadTr','MeanCorrDist GooddTr','GoodTrField');
    plotCompXY(mod_BadTr.adaptSpatialInfo(indFGood),mod_GoodTr.adaptSpatialInfo(indFGood),...
        'Pyr_GoodVsBadTrSpInfo_GoodTrField',pathAnal,'Spatial info. BadTr','Spatial info. GooddTr','GoodTrField');
    
    indFBad = mod_BadTr.isNeuWithField == 1;
    plotCompXY(mod_BadTr.meanCorrDistNZ(indFBad),mod_GoodTr.meanCorrDistNZ(indFBad),...
        'Pyr_GoodVsBadTrMeanCorrDist_BadTrField',pathAnal,'MeanCorrDist BadTr','MeanCorrDist GooddTr','BadTrField');
    plotCompXY(mod_BadTr.adaptSpatialInfo(indFBad),mod_GoodTr.adaptSpatialInfo(indFBad),...
        'Pyr_GoodVsBadTrSpInfo_BadTrField',pathAnal,'Spatial info. BadTr','Spatial info. GooddTr','BadTrField');
    
    plotCompXY(mod_BadTr.meanCorrDistNZ,mod_GoodTr.meanCorrDistNZ,...
        'Pyr_GoodVsBadTrMeanCorrDist',pathAnal,'MeanCorrDist BadTr','MeanCorrDist GooddTr','All');
    plotCompXY(mod_BadTr.adaptSpatialInfo,mod_GoodTr.adaptSpatialInfo,...
        'Pyr_GoodVsBadTrSpInfo',pathAnal,'Spatial info. BadTr','Spatial info. GooddTr','All');
   
    plotBurstVsTheta(mod_GoodTr.burstMeanDire,mod_GoodTr.phaseMeanDire,mod_GoodTr.fractBurst,pathAnal,'_GoodTr');
    plotBurstVsTheta(mod_BadTr.burstMeanDire,mod_BadTr.phaseMeanDire,mod_BadTr.fractBurst,pathAnal,'_BadTr');
    
    modPyrStatsGoodVsBad = modPyrStats_GoodVsBad(mod_GoodTr,mod_BadTr);

    save([pathAnal 'autoCorrPyrAllRec_GoodTr.mat'],'modPyrStatsGoodVsBad','-append'); 
    
        
    %% compare good vs bad trials
    colorSel = 0;
    
    plotBoxPlot(mod_GoodTr.meanInstFR,mod_BadTr.meanInstFR,...
         'Mean inst. FR (Hz)','Pyr_MeanInstFRGoodVsBadTrBox',...
        pathAnal,[],modPyrStatsGoodVsBad.pRSMeanInstFR,colorSel);
       
    plotBoxPlot(mod_GoodTr.phaseMeanDire/pi*180,...
        mod_BadTr.phaseMeanDire/pi*180,'Mean theta phase (deg.)','Pyr_ThetaMeanGoodVsBadTrBox',...
        pathAnal,[],modPyrStatsGoodVsBad.pWWPhaseMeanDire,colorSel);
    
    plotBoxPlot(mod_GoodTr.phaseMeanDire/pi*180,...
        mod_BadTr.phaseMeanDire/pi*180,'Mean theta phase (deg.)','Pyr_ThetaMeanGoodVsBadTrBox',...
        pathAnal,[],modPyrStatsGoodVsBad.pWWPhaseMeanDire,colorSel);
    
    plotBoxPlot(mod_GoodTr.numSpPerBurstMean,...
        mod_BadTr.numSpPerBurstMean,'Num. spikes per burst','Pyr_NumSpPerBurstMeanGoodVsBadTrBox',...
        pathAnal,[1.5 3.5],modPyrStatsGoodVsBad.pRSNumSpPerBurstMean,colorSel);
    
    plotBoxPlot(mod_GoodTr.minPhaseArr,...
        mod_BadTr.minPhaseArr,'Min theta phase','Pyr_MinThetaPhaseGoodVsBadTrBox',...
        pathAnal,[],modPyrStatsGoodVsBad.pWWMinPhase,colorSel);
    
    plotBoxPlot(mod_GoodTr.meanCorrDistNZ,...
        mod_BadTr.meanCorrDistNZ,'Mean corr. dist. NZ','Pyr_MeanCorrDistNZGoodVsBadTrBox',...
        pathAnal,[-0.1 0.2],modPyrStatsGoodVsBad.pRSMeanCorrDistNZ,colorSel);
    
    plotBoxPlot(mod_GoodTr.adaptSpatialInfo,...
        mod_BadTr.adaptSpatialInfo,'Adapt spatial info. (bit)','Pyr_AdaptSpInfoGoodVsBadTrBox',...
        pathAnal,[0 3],modPyrStatsGoodVsBad.pRSAdaptSpatialInfo,colorSel);
    
    plotBoxPlot(mod_GoodTr.sparsity,...
        mod_BadTr.sparsity,'Sparsity','Pyr_SparsityGoodVsBadTrBox',...
        pathAnal,[0 3],modPyrStatsGoodVsBad.pRSSparsity,colorSel);
    
    plotPolarPlot(mod_GoodTr.phaseMeanDire,...
        mod_BadTr.phaseMeanDire,'Mean theta phase direction',...
        'Pyr_ThetaMeanGoodVsBadTrPolar',pathAnal,modPyrStatsGoodVsBad.pWWPhaseMeanDire);
    
end

function modPyr = accumPyrNeurons1(paths,filenames,mazeSess,minFR,maxFR,task,methodTheta,onlyRun,goodTr,spaceBin)

    %% Pyrs in no cue passive task
    numRec = size(paths,1);
    modPyr = struct('task',[],... % no cue - 1, AL - 2, PL - 3
                              'indRec',[],... % recording index
                              'indNeu',[],... % neuron indices
                              'thetaFreqHMean',[],... % theta frequency hilbert
                              ...
                              'thetaMod',[],... % theta modulation
                              'trough',[],... % % ACG first trough
                              'peak',[],... % % ACG first peak
                              'thetaModInd',[],... % theta modulation index
                              'troughT3',[],... % % time of ACG first trough
                              'peakT3',[],... % % time of ACG first peak
                              'thetaModInd3',[],... % theta modulation index (method 3)
                              'thetaAsym3',[],... % theta asymmetry (method 3)
                              'thetaModFreq3',[],... % theta modulation frequency (method 3)
                              ...
                              'mFR',[],... % mean firing rate
                              'meanInstFR',[], ... % mean instantaneous firing rate
                              ...
                              'nNeuWithField',[],... % number of neurons with fields in the recording
                              'isNeuWithField',[],... % does the neuron have field(s)
                              'fieldWidth',[],... % field width
                              'indStartField',[],... % start index of a field
                              'indPeakField',[],... % peak index of a field
                              ...
                              'fractBurst',[],... % fraction of spikes which belongs to a burst over all the spikes across all the trials
                              'burstMeanDire',[],... % the mean phase direction
                              'burstMeanResultantLen',[],... % the mean resultant length of the mean phase direction
                              'nonBurstMeanDire',[],... % the mean phase direction of non-burst spikes
                              'nonBurstMeanResultantLen',[],... % the mean resultant length of the mean phase direction of non-burst spikes
                              'burstMeanDireStart',[],... % the mean phase direction of the first spike of a burst
                              'burstMeanResultantLenStart',[],... % the mean resultant length of the mean phase direction of the first spike of a burst
                              'numSpPerBurstMean',[],... % mean number of spikes per burst
                              ...  
                              'trialLenMean',[],... % mean trial length
                              'minPhaseFilH',[],... % the phase which fires the least number of spikes (hilbert)
                              'maxPhaseFilH',[],... % the phase which fires the largest number of spikes (hilbert)
                              'thetaModHistH',[],... % theta modulation calculated based on theta phase histogram (hilbert)
                              'phaseMeanDireH',[],... % the mean phase direction (hilbert)
                              'phaseMeanResultantLenH',[],... % the mean resultant length of the mean phase direction (hilbert)
                              ...
                              'minPhaseFil',[],... % the phase which fires the least number of spikes 
                              'maxPhaseFil',[],... % the phase which fires the largest number of spikes 
                              'thetaModHist',[],... % theta modulation calculated based on theta phase histogram
                              'phaseMeanDire',[],... % the mean phase direction
                              'phaseMeanResultantLen',[],... % the mean resultant length of the mean phase direction
                              ...
                              'meanCorrDist',[],... % mean correlation between trials
                              'meanCorrDistNZ',[],... % mean correlation between trials with spikes
                              'adaptSpatialInfo',[],... % spatial information
                              'sparsity',[]); % sparsity
                          
    for i = 1:numRec
        disp(filenames(i,:));
        if(i == 7)
            a = 1;
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
        
        fullPathCorr = [filenames(i,:) '_spikesCorrDist_GoodTr_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fullPathCorr];
        if(exist(fullPath) == 0)
            disp('_spikesCorrDist_GoodTr_Run.mat file does not exist.');
            return;
        end
        if(goodTr == 1)
            load(fullPath,'indLapsGoodTr'); 
            indLaps = indLapsGoodTr;
        else
            load(fullPath,'indLapsBadTr'); 
            indLaps = indLapsBadTr;
        end
        if(mazeSess(i) == 0)
            indLaps = indLaps{1}; 
        else
            indLaps = indLaps{mazeSess(i)};
        end
                
        fullPathFRAll = [filenames(i,:) '_FR_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fullPathFRAll];
        if(exist(fullPath) == 0)
            disp('_FR_Run.mat file does not exist.');
            return;
        end
        load(fullPath,'mFRStruct','mFRStructSess'); 
        if(length(beh.mazeSessAll) > 1)
            mFRAll = mFRStructSess{mazeSess(i)};
        else
            mFRAll = mFRStruct;
        end
        
        fullPathFR = [filenames(i,:) '_FR_GoodTr_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fullPathFR];
        if(exist(fullPath) == 0)
            disp('_FR_GoodTr_Run.mat file does not exist.');
            return;
        end
        if(goodTr == 1)
            load(fullPath,'mFRStructSessGoodTr'); 
            mFR = mFRStructSessGoodTr;
        else
            load(fullPath,'mFRStructSessBadTr'); 
            mFR = mFRStructSessBadTr;
        end
        if(mazeSess(i) == 0)
            mFR = mFR{1}; 
        else
            mFR = mFR{mazeSess(i)};
        end
        
        fileNamePeakFR = [filenames(i,:) '_PeakFR_GoodTr' ...
                        num2str(spaceBin) 'mm_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNamePeakFR];
        if(exist(fullPath) == 0)
            disp(['The peak firing rate file does not exist. Please call ',...
                    'function "PeakFiringRate_smTr_GoodTr" first.']);
            return;
        end
        if(goodTr == 1)
            load(fullPath,'pFRStructSessGoodTr');
            pFR = pFRStructSessGoodTr;
        else
            load(fullPath,'pFRStructSessBadTr');
            pFR = pFRStructSessBadTr;
        end
        if(mazeSess(i) == 0)
            pFR = pFR{1}; 
        else
            pFR = pFR{mazeSess(i)};
        end
        
        fileNameThetaMod = [filenames(i,:) '_ThetaMod_GoodTr_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameThetaMod];
        if(exist(fullPath) == 0)
            disp('_ThetaMod_GoodTr file does not exist.');
            return;
        end
        load(fullPath,'thetaModSessGoodTr','thetaModSessBadTr');
        if(goodTr == 1)
            load(fullPath,'thetaModSessGoodTr');
            thetaMod = thetaModSessGoodTr;
        else
            load(fullPath,'thetaModSessBadTr');
            thetaMod = thetaModSessBadTr;
        end
        if(mazeSess(i) == 0)
            thetaMod = thetaMod{1}; 
        else
            thetaMod = thetaMod{mazeSess(i)};
        end
         
        fileNameSpInfo = [filenames(i,:) '_SpInfo_GoodTr_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameSpInfo];
        if(exist(fullPath) == 0)
            disp('_SpInfo_GoodTr file does not exist.');
            return;
        end
        if(goodTr == 1)
            load(fullPath,'spatialInfoSessGoodTr');
            spatialInfo = spatialInfoSessGoodTr;
        else
            load(fullPath,'spatialInfoSessBadTr');
            spatialInfo = spatialInfoSessBadTr;
        end
        if(mazeSess(i) == 0)
            spatialInfo = spatialInfo{1}; 
        else
            spatialInfo = spatialInfo{mazeSess(i)};
        end
        
        fileNameCorr = [filenames(i,:) '_meanSpikesCorrDist_GoodTr_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameCorr];
        if(exist(fullPath) == 0)
            disp('_meanSpikesCorrDist_GoodTr file does not exist.');
            return;
        end
        if(goodTr == 1)
            load(fullPath,'meanCorrDistGoodTr');
            meanCorrDist = meanCorrDistGoodTr;
        else
            load(fullPath,'meanCorrDistBadTr');
            meanCorrDist = meanCorrDistBadTr;
        end
        if(mazeSess(i) == 0)
            meanCorrDist = meanCorrDist{1}; 
        else
            meanCorrDist = meanCorrDist{mazeSess(i)};
        end
        
        if(methodTheta == 0)
            th = 'H';
        else
            th = 'L';
        end               
        fileNameBurst = [filenames(i,:) '_burstAll_GoodTr_TH' th '_Run' num2str(onlyRun) ...
                     '.mat'];
        fullPath = [paths(i,:) fileNameBurst];
        if(exist(fullPath) == 0)
            disp('_bustAll_GoodTr file does not exist.');
            return;
        end
        if(goodTr == 1)
            load(fullPath,'burstIsiPerNeuronSessGoodTr');
            burstIsi = burstIsiPerNeuronSessGoodTr;
        else
            load(fullPath,'burstIsiPerNeuronSessBadTr');
            burstIsi = burstIsiPerNeuronSessBadTr;
        end
        if(mazeSess(i) == 0)
            burstIsi = burstIsi{1}; 
        else
            burstIsi = burstIsi{mazeSess(i)};
        end
             
        fileNameThetaPhase = [filenames(i,:) '_ThetaPhaseL_GoodTr_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameThetaPhase];
        if(exist(fullPath) == 0)
            disp('_ThetaPhaseL_GoodTr file does not exist.');
            return;
        end
        if(goodTr == 1)
            load(fullPath,'spikeThetaPhaseStructSessGoodTr');
            spikeThetaPhase = spikeThetaPhaseStructSessGoodTr;
        else
            load(fullPath,'spikeThetaPhaseStructSessBadTr');
            spikeThetaPhase = spikeThetaPhaseStructSessBadTr;
        end
        if(mazeSess(i) == 0)
            spikeThetaPhase = spikeThetaPhase{1}; 
        else
            spikeThetaPhase = spikeThetaPhase{mazeSess(i)};
        end
        
        fileNameThetaPhase = [filenames(i,:) '_ThetaPhaseH_GoodTr_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameThetaPhase];
        if(exist(fullPath) == 0)
            disp('_ThetaPhaseH_GoodTr file does not exist.');
            return;
        end
        if(goodTr == 1)
            load(fullPath,'spikeThetaPhaseStructSessGoodTr');
            spikeThetaPhaseH = spikeThetaPhaseStructSessGoodTr;
        else
            load(fullPath,'spikeThetaPhaseStructSessBadTr');
            spikeThetaPhaseH = spikeThetaPhaseStructSessBadTr;
        end
        if(mazeSess(i) == 0)
            spikeThetaPhaseH = spikeThetaPhaseH{1}; 
        else
            spikeThetaPhaseH = spikeThetaPhaseH{mazeSess(i)};
        end
        trialLenMean = mean(beh.lenTrials(indLaps));
        
        fileNameFW = [filenames(i,:) '_FieldSpCorr_GoodTr_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameFW];
        if(exist(fullPath) == 0)
            disp('_FieldSpCorr_GoodTr file does not exist.');
            return;
        end
        if(goodTr == 1)
            load(fullPath,'fieldSpCorrSessGoodTr','paramF');
            fieldSpCorr = fieldSpCorrSessGoodTr;
        else
            load(fullPath,'fieldSpCorrSessBadTr','paramF');
            fieldSpCorr = fieldSpCorrSessBadTr;
        end
        if(mazeSess(i) == 0)
            fieldSpCorr = fieldSpCorr{1}; 
        else
            fieldSpCorr = fieldSpCorr{mazeSess(i)};
        end
                
%         indNeu = cluList.firingRate > minFR & cluList.firingRate < maxFR &...
%                     autoCorr.isPyrneuron == 1;
        indNeu = mFRAll.mFR > minFR & mFRAll.mFR < maxFR &...
                    autoCorr.isPyrneuron == 1;
        modPyr.task = [modPyr.task task*ones(1,sum(indNeu))];
        modPyr.indRec = [modPyr.indRec i*ones(1,sum(indNeu))];
        modPyr.indNeu = [modPyr.indNeu find(indNeu == 1)]; 
        modPyr.trialLenMean = [modPyr.trialLenMean trialLenMean*ones(1,sum(indNeu))];
        modPyr.thetaFreqHMean = [modPyr.thetaFreqHMean mean(beh.thetaFreqHMean(indLaps))*ones(1,sum(indNeu))];
        
        numNeu = sum(indNeu);
        if(length(indLaps) >= paramF.minNumTr)
            modPyr.thetaMod = [modPyr.thetaMod thetaMod.thetaMod(indNeu)];
            modPyr.trough = [modPyr.trough thetaMod.trough(indNeu)];
            modPyr.peak = [modPyr.peak thetaMod.peak(indNeu)];
            modPyr.thetaModInd = [modPyr.thetaModInd thetaMod.thetaModInd(indNeu)];
            modPyr.peakT3 = [modPyr.peakT3 thetaMod.peakT3(indNeu)];
            modPyr.troughT3 = [modPyr.troughT3 thetaMod.troughT3(indNeu)];
            modPyr.thetaAsym3 = [modPyr.thetaAsym3 ...
                (abs(thetaMod.peakT3(indNeu))-abs(thetaMod.troughT3(indNeu)))./...
                (abs(thetaMod.peakT3(indNeu)))];
            modPyr.thetaModFreq3 = [modPyr.thetaModFreq3 ...
                1000./(abs(thetaMod.peakT3(indNeu)))];
            modPyr.thetaModInd3 = [modPyr.thetaModInd3 thetaMod.thetaModInd3(indNeu)];
            modPyr.mFR = [modPyr.mFR mFR.mFR(indNeu)];

            nNeurons = length(cluList.firingRate);
            nFieldArr = zeros(1,sum(indNeu));
            indNeuWithField = zeros(1,nNeurons);
            fieldWidth = zeros(1,nNeurons);
            indStartField = zeros(1,nNeurons);
            indPeakField = zeros(1,nNeurons);
            
            [indNeuF,ia] = unique(fieldSpCorr.indNeuron);
            [indNeuF,ic] = intersect(indNeuF,find(autoCorr.isPyrneuron == 1));
            if(length(ia) ~= length(ic))
               ia = ia(ic); 
            end
            indNeuWithField(indNeuF) = 1;
            fieldWidth(indNeuF) = fieldSpCorr.FW(ia);
            indStartField(indNeuF) = fieldSpCorr.indStartField(ia);
            indPeakField(indNeuF) = fieldSpCorr.indPeakField(ia);
            numField = length(indNeuF);
            if(numField > 0)
                nFieldArr = numField * ones(1,sum(indNeu));
            end
            modPyr.nNeuWithField = [modPyr.nNeuWithField nFieldArr];
            modPyr.isNeuWithField = [modPyr.isNeuWithField indNeuWithField(indNeu)];
            modPyr.fieldWidth = [modPyr.fieldWidth fieldWidth(indNeu)];
            modPyr.indStartField = [modPyr.indStartField indStartField(indNeu)];
            modPyr.indPeakField = [modPyr.indPeakField indPeakField(indNeu)];
            modPyr.meanInstFR = [modPyr.meanInstFR pFR.meanInstFR(indNeu)];

            numSpPerBurst = zeros(1,nNeurons);
            for m = 1:nNeurons
                if(~isempty(burstIsi.numSpPerBurst{m}))
                    numSpPerBurst(m) = mean(burstIsi.numSpPerBurst{m});
                end
            end
            modPyr.fractBurst = [modPyr.fractBurst burstIsi.fractBurst(indNeu)];
            modPyr.numSpPerBurstMean = [modPyr.numSpPerBurstMean numSpPerBurst(indNeu)]; 
            modPyr.burstMeanDire = [modPyr.burstMeanDire burstIsi.meanDire(indNeu)];        
            modPyr.burstMeanResultantLen = [modPyr.burstMeanResultantLen burstIsi.meanResultantLen(indNeu)];

            modPyr.nonBurstMeanDire = [modPyr.nonBurstMeanDire burstIsi.meanDireNonBurst(indNeu)];
            modPyr.nonBurstMeanResultantLen = [modPyr.nonBurstMeanResultantLen burstIsi.meanResultantLenNonBurst(indNeu)];

            modPyr.burstMeanDireStart = [modPyr.burstMeanDireStart burstIsi.meanDireStart(indNeu)];        
            modPyr.burstMeanResultantLenStart = [modPyr.burstMeanResultantLenStart burstIsi.meanResultantLenStart(indNeu)];

            modPyr.minPhaseFil = [modPyr.minPhaseFil spikeThetaPhase.minPhaseFilArr(indNeu)];
            modPyr.maxPhaseFil = [modPyr.maxPhaseFil spikeThetaPhase.maxPhaseFilArr(indNeu)];
            modPyr.phaseMeanDire = [modPyr.phaseMeanDire spikeThetaPhase.meanDire(indNeu)];
            modPyr.phaseMeanResultantLen = [modPyr.phaseMeanResultantLen spikeThetaPhase.meanResultantLen(indNeu)];
            modPyr.thetaModHist = [modPyr.thetaModHist spikeThetaPhase.thetaMod(indNeu)]; 

            modPyr.minPhaseFilH = [modPyr.minPhaseFilH spikeThetaPhaseH.minPhaseFilArr(indNeu)];
            modPyr.maxPhaseFilH = [modPyr.maxPhaseFilH spikeThetaPhaseH.maxPhaseFilArr(indNeu)];
            modPyr.phaseMeanDireH = [modPyr.phaseMeanDireH spikeThetaPhaseH.meanDire(indNeu)];
            modPyr.phaseMeanResultantLenH = [modPyr.phaseMeanResultantLenH spikeThetaPhaseH.meanResultantLen(indNeu)];
            modPyr.thetaModHistH = [modPyr.thetaModHistH spikeThetaPhaseH.thetaMod(indNeu)]; 
            
            modPyr.meanCorrDist = [modPyr.meanCorrDist meanCorrDist.mean(indNeu)];
            modPyr.meanCorrDistNZ = [modPyr.meanCorrDistNZ meanCorrDist.meanNZ(indNeu)];
            modPyr.adaptSpatialInfo = [modPyr.adaptSpatialInfo spatialInfo.adaptSpatialInfo(indNeu)];
            modPyr.sparsity = [modPyr.sparsity spatialInfo.sparsity(indNeu)];
        else
            modPyr.thetaMod = [modPyr.thetaMod NaN(1,numNeu)];
            modPyr.trough = [modPyr.trough NaN(1,numNeu)];
            modPyr.peak = [modPyr.peak NaN(1,numNeu)];
            modPyr.thetaModInd = [modPyr.thetaModInd NaN(1,numNeu)];
            modPyr.peakT3 = [modPyr.peakT3 NaN(1,numNeu)];
            modPyr.troughT3 = [modPyr.troughT3 NaN(1,numNeu)];
            modPyr.thetaAsym3 = [modPyr.thetaAsym3 NaN(1,numNeu)];
            modPyr.thetaModFreq3 = [modPyr.thetaModFreq3 NaN(1,numNeu)];
            modPyr.thetaModInd3 = [modPyr.thetaModInd3 NaN(1,numNeu)];
            modPyr.mFR = [modPyr.mFR NaN(1,numNeu)];
            
            modPyr.nNeuWithField = [modPyr.nNeuWithField NaN(1,numNeu)];
            modPyr.isNeuWithField = [modPyr.isNeuWithField NaN(1,numNeu)];
            modPyr.fieldWidth = [modPyr.fieldWidth NaN(1,numNeu)];
            modPyr.indStartField = [modPyr.indStartField NaN(1,numNeu)];
            modPyr.indPeakField = [modPyr.indPeakField NaN(1,numNeu)];
            modPyr.meanInstFR = [modPyr.meanInstFR NaN(1,numNeu)];

            modPyr.fractBurst = [modPyr.fractBurst NaN(1,numNeu)];
            modPyr.numSpPerBurstMean = [modPyr.numSpPerBurstMean NaN(1,numNeu)]; 
            modPyr.burstMeanDire = [modPyr.burstMeanDire NaN(1,numNeu)];        
            modPyr.burstMeanResultantLen = [modPyr.burstMeanResultantLen NaN(1,numNeu)];

            modPyr.nonBurstMeanDire = [modPyr.nonBurstMeanDire NaN(1,numNeu)];
            modPyr.nonBurstMeanResultantLen = [modPyr.nonBurstMeanResultantLen NaN(1,numNeu)];

            modPyr.burstMeanDireStart = [modPyr.burstMeanDireStart NaN(1,numNeu)];        
            modPyr.burstMeanResultantLenStart = [modPyr.burstMeanResultantLenStart NaN(1,numNeu)];

            modPyr.minPhaseFil = [modPyr.minPhaseFil NaN(1,numNeu)];
            modPyr.maxPhaseFil = [modPyr.maxPhaseFil NaN(1,numNeu)];
            modPyr.phaseMeanDire = [modPyr.phaseMeanDire NaN(1,numNeu)];
            modPyr.phaseMeanResultantLen = [modPyr.phaseMeanResultantLen NaN(1,numNeu)];
            modPyr.thetaModHist = [modPyr.thetaModHist NaN(1,numNeu)]; 

            modPyr.minPhaseFilH = [modPyr.minPhaseFilH NaN(1,numNeu)];
            modPyr.maxPhaseFilH = [modPyr.maxPhaseFilH NaN(1,numNeu)];
            modPyr.phaseMeanDireH = [modPyr.phaseMeanDireH NaN(1,numNeu)];
            modPyr.phaseMeanResultantLenH = [modPyr.phaseMeanResultantLenH NaN(1,numNeu)];
            modPyr.thetaModHistH = [modPyr.thetaModHistH NaN(1,numNeu)]; 
            
            modPyr.meanCorrDist = [modPyr.meanCorrDist NaN(1,numNeu)];
            modPyr.meanCorrDistNZ = [modPyr.meanCorrDistNZ NaN(1,numNeu)];
            modPyr.adaptSpatialInfo = [modPyr.adaptSpatialInfo NaN(1,numNeu)];
            modPyr.sparsity = [modPyr.sparsity NaN(1,numNeu)];
        end
    end
    
    indDire = find(modPyr.burstMeanDire < 0);
    modPyr.burstMeanDire(indDire) = modPyr.burstMeanDire(indDire) + 2*pi;
    
    indDire = find(modPyr.nonBurstMeanDire < 0);
    modPyr.nonBurstMeanDire(indDire) = modPyr.nonBurstMeanDire(indDire) + 2*pi;
    
    indDire = find(modPyr.burstMeanDireStart < 0);
    modPyr.burstMeanDireStart(indDire) = modPyr.burstMeanDireStart(indDire) + 2*pi;
    
    indDire = find(modPyr.phaseMeanDire < 0);
    modPyr.phaseMeanDire(indDire) = modPyr.phaseMeanDire(indDire) + 2*pi;
    
    indDire = find(modPyr.phaseMeanDireH < 0);
    modPyr.phaseMeanDireH(indDire) = modPyr.phaseMeanDireH(indDire) + 2*pi;
end

function mod = modPyrAllCondi(modPyrNoCue, modPyrPL, modPyrAL)

    mod.task = [modPyrNoCue.task modPyrAL.task modPyrPL.task];
    mod.indRec = [modPyrNoCue.indRec modPyrAL.indRec modPyrPL.indRec];
    mod.indNeu = [modPyrNoCue.indNeu modPyrAL.indNeu modPyrPL.indNeu];
    mod.nNeuWithField = [modPyrNoCue.nNeuWithField modPyrAL.nNeuWithField modPyrPL.nNeuWithField];
    mod.isNeuWithField = [modPyrNoCue.isNeuWithField modPyrAL.isNeuWithField modPyrPL.isNeuWithField];
    mod.fieldWidth = [modPyrNoCue.fieldWidth modPyrAL.fieldWidth modPyrPL.fieldWidth];
    mod.indStartField = [modPyrNoCue.indStartField modPyrAL.indStartField modPyrPL.indStartField];
    mod.indPeakField = [modPyrNoCue.indPeakField modPyrAL.indPeakField modPyrPL.indPeakField];
    mod.percTrackStartField = [modPyrNoCue.indStartField./modPyrNoCue.trialLenMean...
        modPyrAL.indStartField./modPyrAL.trialLenMean modPyrPL.indStartField./modPyrPL.trialLenMean];
    mod.percTrackPeakField = [modPyrNoCue.indPeakField./modPyrNoCue.trialLenMean...
        modPyrAL.indPeakField./modPyrAL.trialLenMean modPyrPL.indPeakField./modPyrPL.trialLenMean];
       
    mod.mFR = [modPyrNoCue.mFR modPyrAL.mFR modPyrPL.mFR];
    mod.meanInstFR = [modPyrNoCue.meanInstFR modPyrAL.meanInstFR modPyrPL.meanInstFR];
     
    mod.burstMeanResultantLen = [modPyrNoCue.burstMeanResultantLen modPyrAL.burstMeanResultantLen ...
                modPyrPL.burstMeanResultantLen];
    mod.burstMeanDire = [modPyrNoCue.burstMeanDire modPyrAL.burstMeanDire modPyrPL.burstMeanDire];
    mod.nonBurstMeanDire = [modPyrNoCue.nonBurstMeanDire modPyrAL.nonBurstMeanDire modPyrPL.nonBurstMeanDire];
    mod.burstMeanDireStart = [modPyrNoCue.burstMeanDireStart modPyrAL.burstMeanDireStart modPyrPL.burstMeanDireStart];
    mod.numSpPerBurstMean = [modPyrNoCue.numSpPerBurstMean modPyrAL.numSpPerBurstMean modPyrPL.numSpPerBurstMean];
    mod.fractBurst = [modPyrNoCue.fractBurst modPyrAL.fractBurst modPyrPL.fractBurst];
        
    mod.thetaModHist = [modPyrNoCue.thetaModHist modPyrAL.thetaModHist modPyrPL.thetaModHist];
    mod.thetaModHistH = [modPyrNoCue.thetaModHistH modPyrAL.thetaModHistH modPyrPL.thetaModHistH];
    mod.phaseMeanDire = [modPyrNoCue.phaseMeanDire modPyrAL.phaseMeanDire modPyrPL.phaseMeanDire];
    mod.phaseMeanDireH = [modPyrNoCue.phaseMeanDireH modPyrAL.phaseMeanDireH modPyrPL.phaseMeanDireH];
    mod.maxPhaseArr = [modPyrNoCue.maxPhaseFil modPyrAL.maxPhaseFil modPyrPL.maxPhaseFil];
    mod.maxPhaseArrH = [modPyrNoCue.maxPhaseFilH modPyrAL.maxPhaseFilH modPyrPL.maxPhaseFilH];
    mod.minPhaseArr = [modPyrNoCue.minPhaseFil modPyrAL.minPhaseFil modPyrPL.minPhaseFil];
    mod.minPhaseArrH = [modPyrNoCue.minPhaseFilH modPyrAL.minPhaseFilH modPyrPL.minPhaseFilH];
    mod.phaseMeanResultantLen = [modPyrNoCue.phaseMeanResultantLen modPyrAL.phaseMeanResultantLen modPyrPL.phaseMeanResultantLen];
    
    phaseDiff = mod.maxPhaseArr - mod.minPhaseArr;
    phaseDiff(phaseDiff < 0) = phaseDiff(phaseDiff < 0) + 360;
    mod.phaseDiff = phaseDiff;
    phaseDiffH = mod.maxPhaseArrH - mod.minPhaseArrH;
    phaseDiffH(phaseDiffH < 0) = phaseDiffH(phaseDiffH < 0) + 360;
    mod.phaseDiffH = phaseDiffH;
    
    mod.diffNeuronLFPFreq = [modPyrNoCue.thetaModFreq3-modPyrNoCue.thetaFreqHMean...
            modPyrAL.thetaModFreq3-modPyrAL.thetaFreqHMean...
            modPyrPL.thetaModFreq3-modPyrPL.thetaFreqHMean];
    mod.thetaModFreq3 = [modPyrNoCue.thetaModFreq3 modPyrAL.thetaModFreq3 modPyrPL.thetaModFreq3];
    mod.thetaModInd3 = [modPyrNoCue.thetaModInd3 modPyrAL.thetaModInd3 modPyrPL.thetaModInd3];
    mod.thetaModInd = [modPyrNoCue.thetaModInd modPyrAL.thetaModInd modPyrPL.thetaModInd];   
    mod.thetaAsym3 = [modPyrNoCue.thetaAsym3 modPyrAL.thetaAsym3 modPyrPL.thetaAsym3];
    
    mod.meanCorrDist = [modPyrNoCue.meanCorrDist modPyrAL.meanCorrDist modPyrPL.meanCorrDist];
    mod.meanCorrDistNZ = [modPyrNoCue.meanCorrDistNZ modPyrAL.meanCorrDistNZ modPyrPL.meanCorrDistNZ];
    mod.adaptSpatialInfo = [modPyrNoCue.adaptSpatialInfo modPyrAL.adaptSpatialInfo modPyrPL.adaptSpatialInfo];
    mod.sparsity = [modPyrNoCue.sparsity modPyrAL.sparsity modPyrPL.sparsity];
    
    x = mod.burstMeanDire - mod.phaseMeanDire;    
    x(x < -pi) = x(x < -pi) + 2*pi;
    x(x > pi) = x(x > pi) - 2*pi;
    x((mod.fractBurst == 0)) = -100;
    mod.burstThetaDiff = x;
end

function modPyrStatsGoodVsBad = modPyrStats_GoodVsBad(modGood,modBad)
    modPyrStatsGoodVsBad = [];
    indNeuGood = ~isnan(modGood.phaseMeanDire);
    indNeuBad = ~isnan(modBad.phaseMeanDire);
    
    modPyrStatsGoodVsBad.isNeuWithField = [sum(modGood.isNeuWithField(indNeuGood)),sum(modBad.isNeuWithField(indNeuBad))];

    modPyrStatsGoodVsBad.mFR = [mean(modGood.mFR(indNeuGood)),mean(modBad.mFR(indNeuBad))];
    modPyrStatsGoodVsBad.meanInstFR = [mean(modGood.meanInstFR(indNeuGood)),mean(modBad.meanInstFR(indNeuBad))];

    modPyrStatsGoodVsBad.meanDiffNeuronLFPFreq = [mean(modGood.diffNeuronLFPFreq(indNeuGood)),mean(modBad.diffNeuronLFPFreq(indNeuBad))];

    indGoodBurst = modGood.fractBurst > 0;
    indBadBurst = modBad.fractBurst > 0;
    modPyrStatsGoodVsBad.meanBurstMeanDire = [circ_mean(modGood.burstMeanDire(indGoodBurst)'),circ_mean(modBad.burstMeanDire(indBadBurst)')];
    modPyrStatsGoodVsBad.meanNonBurstMeanDire = [circ_mean(modGood.nonBurstMeanDire(indNeuGood)'),circ_mean(modBad.nonBurstMeanDire(indNeuBad)')];     
    modPyrStatsGoodVsBad.meanBurstMeanDireStart = [circ_mean(modGood.burstMeanDireStart(indGoodBurst)'),circ_mean(modBad.burstMeanDireStart(indBadBurst)')];
    modPyrStatsGoodVsBad.meanFractBurst = [mean(modGood.fractBurst(indNeuGood)),mean(modBad.fractBurst(indNeuBad))];
    modPyrStatsGoodVsBad.meanNumSpPerBurstMean = [mean(modGood.numSpPerBurstMean(indGoodBurst)),...
        mean(modBad.numSpPerBurstMean(indBadBurst))];

    modPyrStatsGoodVsBad.meanPhaseMeanDire = [circ_mean(modGood.phaseMeanDire(indNeuGood)'),circ_mean(modBad.phaseMeanDire(indNeuBad)')];
    modPyrStatsGoodVsBad.meanPhaseMeanDireH = [circ_mean(modGood.phaseMeanDireH(indNeuGood)'),circ_mean(modBad.phaseMeanDireH(indNeuBad)')];
    modPyrStatsGoodVsBad.meanMaxPhase = [circ_mean(modGood.maxPhaseArr(indNeuGood)'/180*pi),circ_mean(modBad.maxPhaseArr(indNeuBad)'/180*pi)];
    modPyrStatsGoodVsBad.meanMaxPhaseH = [circ_mean(modGood.maxPhaseArrH(indNeuGood)'/180*pi),circ_mean(modBad.maxPhaseArrH(indNeuBad)'/180*pi)];
    modPyrStatsGoodVsBad.meanminPhase = [circ_mean(modGood.minPhaseArr(indNeuGood)'/180*pi),circ_mean(modBad.minPhaseArr(indNeuBad)'/180*pi)];
    modPyrStatsGoodVsBad.meanminPhaseH = [circ_mean(modGood.minPhaseArrH(indNeuGood)'/180*pi),circ_mean(modBad.minPhaseArrH(indNeuBad)'/180*pi)];
    modPyrStatsGoodVsBad.meanPhaseDiff = [mean(modGood.phaseDiff(indNeuGood)),mean(modBad.phaseDiff(indNeuBad))];
    modPyrStatsGoodVsBad.meanPhaseDiffH = [mean(modGood.phaseDiffH(indNeuGood)),mean(modBad.phaseDiffH(indNeuBad))];
    modPyrStatsGoodVsBad.meanThetaModHist = [mean(modGood.thetaModHist(indNeuGood)),mean(modBad.thetaModHist(indNeuBad))];
    modPyrStatsGoodVsBad.meanThetaModHistH = [mean(modGood.thetaModHistH(indNeuGood)),mean(modBad.thetaModHistH(indNeuBad))];

    modPyrStatsGoodVsBad.meanThetaModFreq3 = [mean(modGood.thetaModFreq3(indNeuGood)),mean(modBad.thetaModFreq3(indNeuBad))];
    modPyrStatsGoodVsBad.meanThetaAsym3 = [mean(modGood.thetaAsym3(indNeuGood)),mean(modBad.thetaAsym3(indNeuBad))];
    modPyrStatsGoodVsBad.meanThetaModInd3 = [mean(modGood.thetaModInd3(~isnan(modGood.thetaModInd3))),mean(modBad.thetaModInd3(~isnan(modBad.thetaModInd3)))];
    modPyrStatsGoodVsBad.meanThetaModInd = [mean(modGood.thetaModInd(~isnan(modGood.thetaModInd))),mean(modBad.thetaModInd(~isnan(modBad.thetaModInd)))];
    
    modPyrStatsGoodVsBad.meanCorrDist = [mean(modGood.meanCorrDist(indNeuGood)), mean(modBad.meanCorrDist(indNeuBad))];
    indGoodNZ = ~isnan(modGood.meanCorrDistNZ);
    indBadNZ = ~isnan(modBad.meanCorrDistNZ);
    modPyrStatsGoodVsBad.meanCorrDistNZ = [mean(modGood.meanCorrDistNZ(indGoodNZ)) mean(modBad.meanCorrDistNZ(indBadNZ))];
    modPyrStatsGoodVsBad.adaptSpatialInfo = [mean(modGood.adaptSpatialInfo(indGoodNZ)) mean(modBad.adaptSpatialInfo(indBadNZ))];
    modPyrStatsGoodVsBad.sparsity = [mean(modGood.sparsity(indGoodNZ)) mean(modBad.sparsity(indBadNZ))];

    modPyrStatsGoodVsBad.pRSDiffNeuronLFPFreq = ranksum(modGood.diffNeuronLFPFreq(indNeuGood),modBad.diffNeuronLFPFreq(indNeuBad));

    modPyrStatsGoodVsBad.pRSMFR = ranksum(modGood.mFR(indNeuGood),modBad.mFR(indNeuBad));
    modPyrStatsGoodVsBad.pRSMeanInstFR = ranksum(modGood.meanInstFR(indNeuGood),modBad.meanInstFR(indNeuBad));

    modPyrStatsGoodVsBad.pKBurstMeanDire = circ_ktest(modGood.burstMeanDire(indGoodBurst)',modBad.burstMeanDire(indBadBurst)');
    modPyrStatsGoodVsBad.pKNonBurstMeanDire = circ_ktest(modGood.nonBurstMeanDire(indNeuGood)',modBad.nonBurstMeanDire(indNeuBad)');
    modPyrStatsGoodVsBad.pKBurstMeanDireStart = circ_ktest(modGood.burstMeanDireStart(indGoodBurst)',modBad.burstMeanDireStart(indBadBurst)');

    modPyrStatsGoodVsBad.pWWBurstMeanDire = circ_wwtest(modGood.burstMeanDire(indGoodBurst)',modBad.burstMeanDire(indBadBurst)');
    modPyrStatsGoodVsBad.pWWNonBurstMeanDire = circ_wwtest(modGood.nonBurstMeanDire(indNeuGood)',modBad.nonBurstMeanDire(indNeuBad)');
    modPyrStatsGoodVsBad.pWWBurstMeanDireStart = circ_wwtest(modGood.burstMeanDireStart(indGoodBurst)',modBad.burstMeanDireStart(indBadBurst)');

    modPyrStatsGoodVsBad.pRSFractBurst = ranksum(modGood.fractBurst(indNeuGood),modBad.fractBurst(indNeuBad));        
    modPyrStatsGoodVsBad.pRSNumSpPerBurstMean = ranksum(modGood.numSpPerBurstMean(indGoodBurst),...
            modBad.numSpPerBurstMean(indBadBurst));

    modPyrStatsGoodVsBad.pWWBurstThetaDiff = circ_wwtest(modGood.burstThetaDiff(indGoodBurst),...
            modBad.burstThetaDiff(indBadBurst));
    modPyrStatsGoodVsBad.pKBurstThetaDiff = circ_ktest(modGood.burstThetaDiff(indGoodBurst),...
            modBad.burstThetaDiff(indBadBurst));

    modPyrStatsGoodVsBad.pKPhaseMeanDire = circ_ktest(modGood.phaseMeanDire(indNeuGood)',modBad.phaseMeanDire(indNeuBad)');
    modPyrStatsGoodVsBad.pKPhaseMeanDireH = circ_ktest(modGood.phaseMeanDireH(indNeuGood)',modBad.phaseMeanDireH(indNeuBad)');
    modPyrStatsGoodVsBad.pKMaxPhase = circ_ktest(modGood.maxPhaseArr(indNeuGood)'/180*pi,modBad.maxPhaseArr(indNeuBad)'/180*pi);
    modPyrStatsGoodVsBad.pKMaxPhaseH = circ_ktest(modGood.maxPhaseArrH(indNeuGood)'/180*pi,modBad.maxPhaseArrH(indNeuBad)'/180*pi);
    modPyrStatsGoodVsBad.pKMinPhase = circ_ktest(modGood.minPhaseArr(indNeuGood)'/180*pi,modBad.minPhaseArr(indNeuBad)'/180*pi);
    modPyrStatsGoodVsBad.pKMinPhaseH = circ_ktest(modGood.minPhaseArrH(indNeuGood)'/180*pi,modBad.minPhaseArrH(indNeuBad)'/180*pi);

    modPyrStatsGoodVsBad.pWWPhaseMeanDire = circ_wwtest(modGood.phaseMeanDire(indNeuGood)',modBad.phaseMeanDire(indNeuBad)');
    modPyrStatsGoodVsBad.pWWPhaseMeanDireH = circ_wwtest(modGood.phaseMeanDireH(indNeuGood)',modBad.phaseMeanDireH(indNeuBad)');
    modPyrStatsGoodVsBad.pWWMaxPhase = circ_wwtest(modGood.maxPhaseArr(indNeuGood)'/180*pi,modBad.maxPhaseArr(indNeuBad)'/180*pi);
    modPyrStatsGoodVsBad.pWWMaxPhaseH = circ_wwtest(modGood.maxPhaseArrH(indNeuGood)'/180*pi,modBad.maxPhaseArrH(indNeuBad)'/180*pi);
    modPyrStatsGoodVsBad.pWWMinPhase = circ_wwtest(modGood.minPhaseArr(indNeuGood)'/180*pi,modBad.minPhaseArr(indNeuBad)'/180*pi);
    modPyrStatsGoodVsBad.pWWMinPhaseH = circ_wwtest(modGood.minPhaseArrH(indNeuGood)'/180*pi,modBad.minPhaseArrH(indNeuBad)'/180*pi);

    modPyrStatsGoodVsBad.pRSPhaseDiff = ranksum(modGood.phaseDiff(indNeuGood),modBad.phaseDiff(indNeuBad));
    modPyrStatsGoodVsBad.pRSPhaseDiffH = ranksum(modGood.phaseDiffH(indNeuGood),modBad.phaseDiffH(indNeuBad));
    modPyrStatsGoodVsBad.pRSThetaModHist = ranksum(modGood.thetaModHist(indNeuGood),modBad.thetaModHist(indNeuBad));
    modPyrStatsGoodVsBad.pRSThetaModHistH = ranksum(modGood.thetaModHistH(indNeuGood),modBad.thetaModHistH(indNeuBad));

    modPyrStatsGoodVsBad.pRSThetaModFreq3 = ranksum(modGood.thetaModFreq3(indNeuGood),modBad.thetaModFreq3(indNeuBad));
    modPyrStatsGoodVsBad.pRSThetaAsym3 = ranksum(modGood.thetaAsym3(indNeuGood),modBad.thetaAsym3(indNeuBad));  
    modPyrStatsGoodVsBad.pRSThetaModInd3 = ranksum(modGood.thetaModInd3(indNeuGood),modBad.thetaModInd3(indNeuBad));
    modPyrStatsGoodVsBad.pRSThetaModInd = ranksum(modGood.thetaModInd(indNeuGood),modBad.thetaModInd(indNeuBad));
    
    modPyrStatsGoodVsBad.pRSMeanCorrDist = ranksum(modGood.meanCorrDist(indNeuGood), modBad.meanCorrDist(indNeuBad));
    modPyrStatsGoodVsBad.pRSMeanCorrDistNZ = ranksum(modGood.meanCorrDistNZ(indGoodNZ),modBad.meanCorrDistNZ(indBadNZ));
    modPyrStatsGoodVsBad.pRSAdaptSpatialInfo = ranksum(modGood.adaptSpatialInfo(indGoodNZ),modBad.adaptSpatialInfo(indBadNZ));
    modPyrStatsGoodVsBad.pRSSparsity = ranksum(modGood.sparsity(indGoodNZ),modBad.sparsity(indBadNZ));

end

function plotBurstVsTheta(burstMeanDire,phaseMeanDire,fractBurst,pathAnal,tit)

    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    x = burstMeanDire(fractBurst > 0) - phaseMeanDire(fractBurst > 0);
    x(x < -pi) = x(x < -pi) + 2*pi;
    x(x > pi) = x(x > pi) - 2*pi;
    X = [x' phaseMeanDire(fractBurst > 0)'];
    N = hist3(X,'CdataMode','auto','Ctrs',{-1.2*pi:pi/18:1.2*pi 0:pi/36:2*pi});
    imagesc(-1.2*pi:pi/36:1.2*pi,0:pi/36:2*pi,N');
    hold on;
    h = plot([0 0],[0 2*pi],'r-');
    xlabel('Burst mean phase - theta mean phase');
    ylabel('Theta mean phase');
    fileName1 = [pathAnal 'Pyr-Burst-ThetaMeanPhaseVsThetaPhase' tit];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
        
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
    x = [x1';x2'];
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

function plotBarPlot(x1,x2,yl,fn,pathAnal,ylimit,p,colorSel)
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
    meanx1 = mean(x1);
    meanx2 = mean(x2);
    sem1 = std(x1)/sqrt(size(x1,1));
    sem2 = std(x2)/sqrt(size(x2,1));
    hold on;
    h = bar(1,meanx1,0.6,'FaceColor',colorArr(2,:));
    set(h,'FaceAlpha',0.5);
    h = bar(2,meanx2,0.6,'FaceColor',colorArr(1,:));
    set(h,'FaceAlpha',0.5);
    h = errorbar(1,meanx1,sem1,'LineWidth',1,'Marker','.','Color','k');    
    h = errorbar(2,meanx2,sem2,'LineWidth',1,'Marker','.','Color','k');
    if(~isempty(ylimit))
        set(gca,'YLim',ylimit);
    end
    ylabel(yl);
    title(['p = ' num2str(p)]);
    
    fileName1 = [pathAnal fn];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotPolarPlot(x1,x2,ti,fn,pathAnal,p)
    [figNew,pos] = CreateFig();
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    colorArr = [163 207 98;...
                234 131 114]/255;
    
    polarhistogram(x1,0:pi/18:2*pi,'Normalization','probability','DisplayStyle','bar',....
        'FaceAlpha',0.5,'FaceColor',colorArr(2,:),'EdgeColor',[0.5 0.5 0.5]);
    hold on;
    if(~isempty(x2))
        polarhistogram(x2,0:pi/18:2*pi,'Normalization','probability','DisplayStyle','bar',....
            'FaceAlpha',0.5,'FaceColor',colorArr(1,:),'EdgeColor',[0.5 0.5 0.5]);
        title([ti ' p = ' num2str(p)])  
    else
        title(ti)
    end
    
    fileName1 = [pathAnal fn];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotDistri(x,y,histBins,xl,leg,pathAnal,fileN,xlimit)
    
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    colorArr = [234 131 114;...
                163 207 98]/255;
    hold on;
    hx = hist(x,histBins);
    hy = hist(y,histBins);
    h = plot(histBins,hx/sum(hx),'-');
    set(h,'LineWidth',2,'Color',colorArr(1,:));
    h = plot(histBins,hy/sum(hy),'-');
    set(h,'LineWidth',2,'Color',colorArr(2,:));
    if(~isempty(xlimit))
        set(gca,'xLim',xlimit);
    end
    xlabel(xl);
    ylabel('Probability');
    legend(leg);
    
    fileName1 = [pathAnal fileN];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotStackedBar(x,y,yl,fileN,pathAnal)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])

    % only for plotting two clusters with fields
    colorArr = [...
                234 131 114;...
                163 207 98]/255;
    a = [1 nan];
    b = [x,y;nan(1,2)];
    h = bar(a,b,'stacked');
    h(1).FaceColor = colorArr(1,:);
    h(1).BarWidth = 0.5;
    h(2).FaceColor = colorArr(2,:);
    h(2).BarWidth = 0.5;
    ylabel(yl);
    set(gca,'XLim',[0 2])
    
    fileName1 = [pathAnal fileN];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotCompXY(x,y,fileN,pathAnal,xl,yl,ti)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
        
    h = plot(x,y,'k+');
    set(h,'MarkerSize',9);
    hold on;
    minXY = min([x,y]);
    maxXY = max([x,y]);
    plot([minXY maxXY],[minXY maxXY],'r-');
    
    xlabel(xl);
    ylabel(yl);
    title(ti);
    
    fileName1 = [pathAnal fileN 'Field'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end
