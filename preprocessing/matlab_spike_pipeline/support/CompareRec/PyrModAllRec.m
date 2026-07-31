function PyrModAllRec(onlyRun,taskSel,methodKMean)
% e.g. PyrModAllRec(1,1,2)

    methodTheta = 1;
    minFR = 0.15;
    maxFR = 7;
    if(nargin == 1)
        methodKMean = 2; % which kmean method is used
    end
    
    RecordingList;
    pathAnal0 = 'Z:\Yingxue\DataAnalysisRaphi\Pyramidal\';
    if(taskSel == 1) % including all the neurons
        pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\Pyramidal\';
    elseif(taskSel == 2) % including AL and PL neurons
        pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\PyramidalALPL\';
    else % AL neurons only
        pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\PyramidalAL\';
    end
    
    if(exist([pathAnal0 'autoCorrPyrAllRec.mat']))
        load([pathAnal0 'autoCorrPyrAllRec.mat']);
    end
    
%     %% pyramidal neurons in no cue passive task
%     disp('No cue')
%     modPyrNoCue = accumPyrNeurons1(listRecordingsNoCuePath,...
%         listRecordingsNoCueFileName,mazeSessionNoCue,minFR,maxFR,1,methodTheta,onlyRun);
%     
%     disp('Active licking')
%     modPyrAL = accumPyrNeurons1(listRecordingsActiveLickPath,...
%         listRecordingsActiveLickFileName,mazeSessionActiveLick,minFR,maxFR,2,methodTheta,onlyRun);
%     
%     disp('Passive licking')
%     modPyrPL = accumPyrNeurons1(listRecordingsPassiveLickPath,...
%         listRecordingsPassiveLickFileName,mazeSessionPassiveLick,minFR,maxFR,3,methodTheta,onlyRun);
%     
%     save([pathAnal 'autoCorrPyrAllRec.mat'],'modPyrNoCue','modPyrAL','modPyrPL','-append'); 
    
%     plotConditions(modPyrNoCue.thetaModHist,modPyrNoCue.phaseMeanDire,...
%         modPyrAL.thetaModHist,modPyrAL.phaseMeanDire,...
%         modPyrPL.thetaModHist,modPyrPL.phaseMeanDire,...
%         'Theta modulation (hist)','Theta phase mean direction',[]);
%     
%     plotConditions(modPyrNoCue.thetaModHist,modPyrNoCue.thetaModInd3,...
%         modPyrAL.thetaModHist,modPyrAL.thetaModInd3,...
%         modPyrPL.thetaModHist,modPyrPL.thetaModInd3,...
%         'Theta modulation (hist)','Theta modulation 3',[]);
    
%     %% plot each Pyr clusters based on the task type
%     plotClusters(modPyrNoCue.burstMeanDire,modPyrNoCue.phaseMeanDire,...
%         autoCorrPyrAll.idxC2(autoCorrPyrAll.task == autoCorrPyrNoCue.task(1)),...
%         'Burst mean direction','Theta phase mean direction','No cue task')
%     plot([0 2*pi],[0 2*pi],'k-')
%     
%     plotClusters(modPyrAL.burstMeanDire,modPyrAL.phaseMeanDire,...
%         autoCorrPyrAll.idxC2(autoCorrPyrAll.task == autoCorrPyrAL.task(1)),...
%         'Burst mean direction','Theta phase mean direction','AL task')
%     plot([0 2*pi],[0 2*pi],'k-')
%     
%     plotClusters(modPyrPL.burstMeanDire,modPyrPL.phaseMeanDire,...
%         autoCorrPyrAll.idxC2(autoCorrPyrAll.task == autoCorrPyrPL.task(1)),...
%         'Burst mean direction','Theta phase mean direction','PL task')
%     plot([0 2*pi],[0 2*pi],'k-')
    
%     %% number of fields
%     plotClusters(modPyrAL.burstMeanDire,modPyrAL.nNeuWithField,...
%         autoCorrPyrAll.idxC2(autoCorrPyrAll.task == autoCorrPyrAL.task(1)),...
%         'Burst mean direction','Num. fields','AL task')
%     
%     plotClusters(modPyrPL.burstMeanDire,modPyrPL.nNeuWithField,...
%         autoCorrPyrAll.idxC2(autoCorrPyrAll.task == autoCorrPyrPL.task(1)),...
%         'Burst mean direction','Num. fields','PL task')
    if(taskSel == 1)
        mod.task = [autoCorrPyrNoCue.task autoCorrPyrAL.task autoCorrPyrPL.task];
        mod.indRec = [autoCorrPyrNoCue.indRec autoCorrPyrAL.indRec autoCorrPyrPL.indRec];
        mod.indNeu = [autoCorrPyrNoCue.indNeu autoCorrPyrAL.indNeu autoCorrPyrPL.indNeu];
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

    %     idxC = [autoCorrPyrAll.idxC(autoCorrPyrAll.task == autoCorrPyrNoCue.task(1))' ...
    %         autoCorrPyrAll.idxC(autoCorrPyrAll.task == autoCorrPyrAL.task(1))' ...
    %         autoCorrPyrAll.idxC(autoCorrPyrAll.task == autoCorrPyrPL.task(1))'];
        if(methodKMean == 1)
            mod.idxC = [autoCorrPyrAll.idxC1(autoCorrPyrAll.task == autoCorrPyrNoCue.task(1))' ...
                autoCorrPyrAll.idxC1(autoCorrPyrAll.task == autoCorrPyrAL.task(1))' ...
                autoCorrPyrAll.idxC1(autoCorrPyrAll.task == autoCorrPyrPL.task(1))'];
        elseif(methodKMean == 2)
            mod.idxC = [autoCorrPyrAll.idxC2(autoCorrPyrAll.task == autoCorrPyrNoCue.task(1))' ...
                autoCorrPyrAll.idxC2(autoCorrPyrAll.task == autoCorrPyrAL.task(1))' ...
                autoCorrPyrAll.idxC2(autoCorrPyrAll.task == autoCorrPyrPL.task(1))'];
        elseif(methodKMean == 3)
            mod.idxC = [autoCorrPyrAll.idxC3(autoCorrPyrAll.task == autoCorrPyrNoCue.task(1)) ...
                autoCorrPyrAll.idxC3(autoCorrPyrAll.task == autoCorrPyrAL.task(1)) ...
                autoCorrPyrAll.idxC3(autoCorrPyrAll.task == autoCorrPyrPL.task(1))];
        end
    elseif(taskSel == 2)
        mod.task = [autoCorrPyrAL.task autoCorrPyrPL.task];
        mod.indRec = [autoCorrPyrAL.indRec autoCorrPyrPL.indRec];
        mod.indNeu = [autoCorrPyrAL.indNeu autoCorrPyrPL.indNeu];
        mod.nNeuWithField = [modPyrAL.nNeuWithField modPyrPL.nNeuWithField];
        mod.isNeuWithField = [modPyrAL.isNeuWithField modPyrPL.isNeuWithField];
        mod.fieldWidth = [modPyrAL.fieldWidth modPyrPL.fieldWidth];
        mod.indStartField = [modPyrAL.indStartField modPyrPL.indStartField];
        mod.indPeakField = [modPyrAL.indPeakField modPyrPL.indPeakField];
        mod.percTrackStartField = [...
            modPyrAL.indStartField./modPyrAL.trialLenMean modPyrPL.indStartField./modPyrPL.trialLenMean];
        mod.percTrackPeakField = [...
            modPyrAL.indPeakField./modPyrAL.trialLenMean modPyrPL.indPeakField./modPyrPL.trialLenMean];

        mod.mFR = [modPyrAL.mFR modPyrPL.mFR];
        mod.meanInstFR = [modPyrAL.meanInstFR modPyrPL.meanInstFR];

        mod.burstMeanResultantLen = [modPyrAL.burstMeanResultantLen ...
                    modPyrPL.burstMeanResultantLen];
        mod.burstMeanDire = [modPyrAL.burstMeanDire modPyrPL.burstMeanDire];
        mod.nonBurstMeanDire = [modPyrAL.nonBurstMeanDire modPyrPL.nonBurstMeanDire];
        mod.burstMeanDireStart = [modPyrAL.burstMeanDireStart modPyrPL.burstMeanDireStart];
        mod.numSpPerBurstMean = [modPyrAL.numSpPerBurstMean modPyrPL.numSpPerBurstMean];
        mod.fractBurst = [modPyrAL.fractBurst modPyrPL.fractBurst];

        mod.thetaModHist = [modPyrAL.thetaModHist modPyrPL.thetaModHist];
        mod.thetaModHistH = [modPyrAL.thetaModHistH modPyrPL.thetaModHistH];
        mod.phaseMeanDire = [modPyrAL.phaseMeanDire modPyrPL.phaseMeanDire];
        mod.phaseMeanDireH = [modPyrAL.phaseMeanDireH modPyrPL.phaseMeanDireH];
        mod.maxPhaseArr = [modPyrAL.maxPhaseFil modPyrPL.maxPhaseFil];
        mod.maxPhaseArrH = [modPyrAL.maxPhaseFilH modPyrPL.maxPhaseFilH];
        mod.minPhaseArr = [modPyrAL.minPhaseFil modPyrPL.minPhaseFil];
        mod.minPhaseArrH = [modPyrAL.minPhaseFilH modPyrPL.minPhaseFilH];
        mod.phaseMeanResultantLen = [modPyrAL.phaseMeanResultantLen modPyrPL.phaseMeanResultantLen];

        phaseDiff = mod.maxPhaseArr - mod.minPhaseArr;
        phaseDiff(phaseDiff < 0) = phaseDiff(phaseDiff < 0) + 360;
        mod.phaseDiff = phaseDiff;
        phaseDiffH = mod.maxPhaseArrH - mod.minPhaseArrH;
        phaseDiffH(phaseDiffH < 0) = phaseDiffH(phaseDiffH < 0) + 360;
        mod.phaseDiffH = phaseDiffH;

        mod.diffNeuronLFPFreq = [...
                modPyrAL.thetaModFreq3-modPyrAL.thetaFreqHMean...
                modPyrPL.thetaModFreq3-modPyrPL.thetaFreqHMean];
        mod.thetaModFreq3 = [modPyrAL.thetaModFreq3 modPyrPL.thetaModFreq3];
        mod.thetaModInd3 = [modPyrAL.thetaModInd3 modPyrPL.thetaModInd3];
        mod.thetaModInd = [modPyrAL.thetaModInd modPyrPL.thetaModInd];   
        mod.thetaAsym3 = [modPyrAL.thetaAsym3 modPyrPL.thetaAsym3];

    %     idxC = [autoCorrPyrAll.idxC(autoCorrPyrAll.task == autoCorrPyrNoCue.task(1))' ...
    %         autoCorrPyrAll.idxC(autoCorrPyrAll.task == autoCorrPyrAL.task(1))' ...
    %         autoCorrPyrAll.idxC(autoCorrPyrAll.task == autoCorrPyrPL.task(1))'];
        if(methodKMean == 1)
            mod.idxC = [...
                autoCorrPyrAll.idxC1(autoCorrPyrAll.task == autoCorrPyrAL.task(1))' ...
                autoCorrPyrAll.idxC1(autoCorrPyrAll.task == autoCorrPyrPL.task(1))'];
        elseif(methodKMean == 2)
            mod.idxC = [...
                autoCorrPyrAll.idxC2(autoCorrPyrAll.task == autoCorrPyrAL.task(1))' ...
                autoCorrPyrAll.idxC2(autoCorrPyrAll.task == autoCorrPyrPL.task(1))'];
        elseif(methodKMean == 3)
            mod.idxC = [...
                autoCorrPyrAll.idxC3(autoCorrPyrAll.task == autoCorrPyrAL.task(1)) ...
                autoCorrPyrAll.idxC3(autoCorrPyrAll.task == autoCorrPyrPL.task(1))];
        end
    else
        mod.task = autoCorrPyrAL.task;
        mod.indRec = autoCorrPyrAL.indRec;
        mod.indNeu = autoCorrPyrAL.indNeu;
        mod.nNeuWithField = modPyrAL.nNeuWithField;
        mod.isNeuWithField = modPyrAL.isNeuWithField;
        mod.fieldWidth = modPyrAL.fieldWidth;
        mod.indStartField = modPyrAL.indStartField;
        mod.indPeakField = modPyrAL.indPeakField;
        mod.percTrackStartField = modPyrAL.indStartField./modPyrAL.trialLenMean;
        mod.percTrackPeakField = modPyrAL.indPeakField./modPyrAL.trialLenMean;

        mod.mFR = modPyrAL.mFR;
        mod.meanInstFR = modPyrAL.meanInstFR;

        mod.burstMeanResultantLen = modPyrAL.burstMeanResultantLen;
        mod.burstMeanDire = modPyrAL.burstMeanDire;
        mod.nonBurstMeanDire = modPyrAL.nonBurstMeanDire;
        mod.burstMeanDireStart = modPyrAL.burstMeanDireStart;
        mod.numSpPerBurstMean = modPyrAL.numSpPerBurstMean;
        mod.fractBurst = modPyrAL.fractBurst;

        mod.thetaModHist = modPyrAL.thetaModHist;
        mod.thetaModHistH = modPyrAL.thetaModHistH;
        mod.phaseMeanDire = modPyrAL.phaseMeanDire;
        mod.phaseMeanDireH = modPyrAL.phaseMeanDireH;
        mod.maxPhaseArr = modPyrAL.maxPhaseFil;
        mod.maxPhaseArrH = modPyrAL.maxPhaseFilH;
        mod.minPhaseArr = modPyrAL.minPhaseFil;
        mod.minPhaseArrH = modPyrAL.minPhaseFilH;
        mod.phaseMeanResultantLen = modPyrAL.phaseMeanResultantLen;

        phaseDiff = mod.maxPhaseArr - mod.minPhaseArr;
        phaseDiff(phaseDiff < 0) = phaseDiff(phaseDiff < 0) + 360;
        mod.phaseDiff = phaseDiff;
        phaseDiffH = mod.maxPhaseArrH - mod.minPhaseArrH;
        phaseDiffH(phaseDiffH < 0) = phaseDiffH(phaseDiffH < 0) + 360;
        mod.phaseDiffH = phaseDiffH;

        mod.diffNeuronLFPFreq = modPyrAL.thetaModFreq3-modPyrAL.thetaFreqHMean;
        mod.thetaModFreq3 = modPyrAL.thetaModFreq3;
        mod.thetaModInd3 = modPyrAL.thetaModInd3;
        mod.thetaModInd = modPyrAL.thetaModInd;   
        mod.thetaAsym3 = modPyrAL.thetaAsym3;

    %     idxC = [autoCorrPyrAll.idxC(autoCorrPyrAll.task == autoCorrPyrNoCue.task(1))' ...
    %         autoCorrPyrAll.idxC(autoCorrPyrAll.task == autoCorrPyrAL.task(1))' ...
    %         autoCorrPyrAll.idxC(autoCorrPyrAll.task == autoCorrPyrPL.task(1))'];
        if(methodKMean == 1)
            mod.idxC = autoCorrPyrAll.idxC1(autoCorrPyrAll.task == autoCorrPyrAL.task(1))';
        elseif(methodKMean == 2)
            mod.idxC = autoCorrPyrAll.idxC2(autoCorrPyrAll.task == autoCorrPyrAL.task(1))';
        elseif(methodKMean == 3)
            mod.idxC = autoCorrPyrAll.idxC3(autoCorrPyrAll.task == autoCorrPyrAL.task(1));
        end
    end
    
    x = mod.burstMeanDire - mod.phaseMeanDire;    
    x(x < -pi) = x(x < -pi) + 2*pi;
    x(x > pi) = x(x > pi) - 2*pi;
    x((mod.fractBurst == 0)) = -100;
    mod.burstThetaDiff = x;
    
%     
%     figure
%     for i = 1:max(idxC)
%         indCurC = idxC == i;
%         h = plot(phaseMeanDire(indCurC),phaseDiff(indCurC),'o');
%         set(gca,'XLim',[0 2*pi],'YLim',[0 300])
%         pause;
%     end

    % for each cluster, compare neurons in the recordings with fields vs.
    % recordings without field
    modPyrStatsField = modPyrStatsFieldPerC(mod,autoCorrPyrAll);
    
    % compare neurons between the two clusters
    modPyrStatsC = modPyrStatsClass(mod);
    
    % for each cluster, compare neurons with fields and without field
    modAlignedPyrStatsFNeuVsNoFNeu = modPyrStatsFieldNeu(mod,autoCorrPyrAll);
    
%     modAlignedPyrStatsFRecVsNoFRec = modPyrStatsFieldNeu(mod,autoCorrPyrAll);
    if(taskSel == 1)
        save([pathAnal 'autoCorrPyrAllRec_km' num2str(methodKMean) '.mat'],'modPyrStatsField','modPyrStatsC',...
            'modAlignedPyrStatsFNeuVsNoFNeu'); 
    elseif(taskSel == 2)
        save([pathAnal 'autoCorrPyrALPLRec_km' num2str(methodKMean) '.mat'],'modPyrStatsField','modPyrStatsC',...
            'modAlignedPyrStatsFNeuVsNoFNeu'); 
    else
        save([pathAnal 'autoCorrPyrALRec_km' num2str(methodKMean) '.mat'],'modPyrStatsField','modPyrStatsC',...
            'modAlignedPyrStatsFNeuVsNoFNeu'); 
    end
    
    if(taskSel == 1)
        diffFreqNoField = mod.diffNeuronLFPFreq(mod.isNeuWithField == 0 & (mod.task == 2 | mod.task == 3));
        diffFreqField = mod.diffNeuronLFPFreq(mod.isNeuWithField == 1 & (mod.task == 2 | mod.task == 3));
        pRSDiffFreq = ranksum(diffFreqNoField,diffFreqField);
        plotBarOnly([mean(diffFreqNoField),mean(diffFreqField)],...
            [std(diffFreqNoField)/sqrt(sum(mod.isNeuWithField == 0)),...
            std(diffFreqField)/sqrt(sum(mod.isNeuWithField))],...
            '','Neuron mod. freq. - LFP freq.', ['p=' num2str(pRSDiffFreq)],pathAnal,'ModFreqMinorsLFPFreq')
    
        diffFreqNoField = mod.diffNeuronLFPFreq(mod.isNeuWithField == 0 & mod.task == 1); % PL no cue no field neurons
        diffFreqField = mod.diffNeuronLFPFreq(mod.isNeuWithField == 1 & (mod.task == 2 | mod.task == 3));
                % PL and AL field neurons
        pRSDiffFreq = ranksum(diffFreqNoField,diffFreqField);
        plotBarOnly([mean(diffFreqNoField),mean(diffFreqField)],...
            [std(diffFreqNoField)/sqrt(sum(mod.isNeuWithField == 0)),...
            std(diffFreqField)/sqrt(sum(mod.isNeuWithField))],...
            '','Neuron mod. freq. - LFP freq.', ['p=' num2str(pRSDiffFreq)],pathAnal,'ModFreqMinorsLFPFreqPLALFieldVsPLNoCueNoField')
    end
    
    
    %% theta phase mean direction vs. 
    % diff between burst mean direction and theta mean direction
    plotBurstVsTheta(mod.burstMeanDire,mod.phaseMeanDire,mod.fractBurst,pathAnal);
   
    %% C1 and C2 contain what percentage of neurons with field 
    clu1 = 1; % deep
    clu2 = 2; % superficial
%     if(methodKMean == 3)
%         clu1 = 1; % deep
%         clu2 = 2; % superficial
%     end
    numNeuField = sum(mod.isNeuWithField);
    percFieldNeuC1 = sum(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu1})/numNeuField;
    percFieldNeuC2 = sum(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu2})/numNeuField;
    plotStackedBar(percFieldNeuC1,percFieldNeuC2,'Perc. neurons with field',...
        'Pyr_PercNeuWithFieldC1C2',pathAnal);
    
    %% plot polar plot of theta phase direction for the whole population
    plotPolarPlot(mod.phaseMeanDire,[],...
        'Mean theta phase direction',...
        'Pyr_ThetaMeanAllPolar',pathAnal,[]);
    
    plotPolarPlot(mod.burstMeanDire(mod.fractBurst > 0),[],...
        'Burst phase',...
        'Pyr_BurstMeanDireAllPolar',pathAnal,[]);
    
    plotPolarPlot(mod.nonBurstMeanDire,[],...
        'Non-burst phase',...
        'Pyr_NonBurstMeanDireAllPolar',pathAnal,[]);
    
    plotPolarPlot(mod.burstMeanDireStart(mod.fractBurst > 0),[],...
        'Burst start phase direction',...
        'Pyr_BurstMeanDireStartAllPolar',pathAnal,[]);
    
    %% compare two clusters
    colorSel = 0;
    idxC = mod.idxC;
        
    if(methodKMean == 1)
        plotBoxPlot(autoCorrPyrAll.relDepthNeuHDefC1{clu1},...
            autoCorrPyrAll.relDepthNeuHDefC1{clu2},'Depth','relDepthNeuHDefCmpC1C2Box',...
            pathAnal,[-6 6],autoCorrPyrAll.pRSRelDepthNeuHDefC1,colorSel);
    elseif(methodKMean == 2)
        plotBoxPlot(autoCorrPyrAll.relDepthNeuHDefC2{clu1},...
            autoCorrPyrAll.relDepthNeuHDefC2{clu2},'Depth','relDepthNeuHDefCmpC1C2Box',...
            pathAnal,[-6 6],autoCorrPyrAll.pRSRelDepthNeuHDefC2,colorSel);
    elseif(methodKMean == 3)
        clu1 = 1; % deep
        clu2 = 2; % superficial
        plotBoxPlot(autoCorrPyrAll.relDepthNeuHDefC3{clu1},...
            autoCorrPyrAll.relDepthNeuHDefC3{clu2},'Depth','relDepthNeuHDefCmpC1C2Box',...
            pathAnal,[-6 6],autoCorrPyrAll.pRSRelDepthNeuHDefC3,colorSel);
    end
    
    plotBoxPlot(mod.phaseMeanDire(idxC == clu1)/pi*180,...
        mod.phaseMeanDire(idxC == clu2)/pi*180,'Mean theta phase (deg.)','Pyr_ThetaMeanC1C2Box',...
        pathAnal,[],modPyrStatsC.pWWPhaseMeanDireC,colorSel);
    
    plotBoxPlot(mod.phaseMeanResultantLen(idxC == clu1),...
        mod.phaseMeanResultantLen(idxC == clu2),'Theta phase resultant length','Pyr_PhaseMeanResultantLenC1C2Box',...
        pathAnal,[-0.1 0.6],modPyrStatsC.pRSPhaseMeanResultantLenC,colorSel);
    
    plotBoxPlot(mod.burstMeanResultantLen(idxC == clu1),...
        mod.burstMeanResultantLen(idxC == clu2),'Burst phase resultant length','Pyr_BurstMeanResultantLenC1C2Box',...
        pathAnal,[-0.1 0.8],modPyrStatsC.pRSBurstMeanResultantLenC,colorSel);
    
    plotBoxPlot(mod.thetaModHist(idxC == clu1),...
        mod.thetaModHist(idxC == clu2),'Theta modulation','Pyr_ThetaModHistC1C2Box',...
        pathAnal,[-0.1 1.2],modPyrStatsC.pRSThetaModHistC,colorSel);
    
    plotBoxPlot(mod.burstThetaDiff(idxC == clu1 & mod.fractBurst > 0)/pi*180,...
        mod.burstThetaDiff(idxC == clu2 & mod.fractBurst > 0)/pi*180,'Burst phase - theta phase (deg.)','Pyr_BurstPhase-ThetaPhaseC1C2Box',...
        pathAnal,[],modPyrStatsC.pWWBurstThetaDiffC,colorSel);
    
    plotBoxPlot(mod.fractBurst(idxC == clu1),...
        mod.fractBurst(idxC == clu2),'Fract burst','Pyr_FractBurstC1C2Box',...
        pathAnal,[],modPyrStatsC.pRSFractBurstC,colorSel);
    
    plotBoxPlot(mod.numSpPerBurstMean(idxC == clu1),...
        mod.numSpPerBurstMean(idxC == clu2),'Num. spikes per burst','Pyr_NumSpPerBurstMeanC1C2Box',...
        pathAnal,[1.8 3.2],modPyrStatsC.pRSNumSpPerBurstMeanC,colorSel);
    
    plotBoxPlot(mod.thetaModFreq3(idxC == clu1),...
        mod.thetaModFreq3(idxC == clu2),'Theta modulation frequency (Hz)','Pyr_ThetaModFreq3C1C2Box',...
        pathAnal,[],modPyrStatsC.pRSThetaModFreq3C,colorSel);
    
    plotBoxPlot(mod.diffNeuronLFPFreq(idxC == clu1),...
        mod.diffNeuronLFPFreq(idxC == clu2),'Neuron mod freq. - LFP freq. (Hz)',...
        'Pyr_DiffNeuronLFPFreqC1C2Box',...
        pathAnal,[],modPyrStatsC.pRSDiffNeuronLFPFreq,colorSel);
    
    plotBoxPlot(mod.thetaAsym3(idxC == clu1),...
        mod.thetaAsym3(idxC == clu2),'Theta asymmetry','Pyr_ThetaAsym3C1C2Box',...
        pathAnal,[0.15 0.85],modPyrStatsC.pRSThetaAsym3C,colorSel);
        
    if(methodKMean == 1)
        plotDistri(autoCorrPyrAll.relDepthNeuHDefC1{clu1},...
            autoCorrPyrAll.relDepthNeuHDefC1{clu2},...
            [min(autoCorrPyrAll.relDepthNeuHDef):max(autoCorrPyrAll.relDepthNeuHDef)],...
            'Depth',['C1' 'C2'],pathAnal,'relDepthNeuHDefCmpC1C2',[-5 5]);
    elseif(methodKMean == 2)
        plotDistri(autoCorrPyrAll.relDepthNeuHDefC2{clu1},...
            autoCorrPyrAll.relDepthNeuHDefC2{clu2},...
            [min(autoCorrPyrAll.relDepthNeuHDef):max(autoCorrPyrAll.relDepthNeuHDef)],...
            'Depth',['C1' 'C2'],pathAnal,'relDepthNeuHDefCmpC1C2',[-5 5]);
    else
        plotDistri(autoCorrPyrAll.relDepthNeuHDefC3{clu1},...
            autoCorrPyrAll.relDepthNeuHDefC3{clu2},...
            [min(autoCorrPyrAll.relDepthNeuHDef):max(autoCorrPyrAll.relDepthNeuHDef)],...
            'Depth',['C1' 'C2'],pathAnal,'relDepthNeuHDefCmpC1C2',[-5 5]);
    end
    
    plotPolarPlot(mod.phaseMeanDire(idxC == clu1),...
        mod.phaseMeanDire(idxC == clu2),'Mean theta phase direction',...
        'Pyr_ThetaMeanC1C2Polar',pathAnal,modPyrStatsC.pWWPhaseMeanDireC);
    
    plotPolarPlot(mod.maxPhaseArr(idxC == clu1)/180*pi,...
        mod.maxPhaseArr(idxC == clu2)/180*pi,'Max theta phase',...
        'Pyr_MaxPhaseArrC1C2Polar',pathAnal,modPyrStatsC.pWWMaxPhaseC);
    
    plotPolarPlot(mod.minPhaseArr(idxC == clu1)/180*pi,...
        mod.minPhaseArr(idxC == clu2)/180*pi,'Min theta phase',...
        'Pyr_MinPhaseArrC1C2Polar',pathAnal,modPyrStatsC.pWWMinPhaseC);
    
    plotPolarPlot(mod.burstMeanDire(idxC == clu1 & mod.fractBurst > 0),...
        mod.burstMeanDire(idxC == clu2 & mod.fractBurst > 0),'Burst mean phase',...
        'Pyr_BurstMeanDireC1C2Polar',pathAnal,modPyrStatsC.pWWburstMeanDireC);
    
%     plotPolarPlot(mod.nonBurstMeanDire(idxC == clu1),...
%         mod.nonBurstMeanDire(idxC == clu2),'Non-burst mean phase',...
%         'Pyr_NonBurstMeanDireC1C2Polar',pathAnal,modPyrStatsC.pWWnonburstMeanDireC);
    
    plotPolarPlot(mod.burstMeanDireStart(idxC == clu1 & mod.fractBurst > 0),...
        mod.burstMeanDireStart(idxC == clu2 & mod.fractBurst > 0),'Burst start phase direction',...
        'Pyr_BurstMeanDireStartC1C2Polar',pathAnal,modPyrStatsC.pWWburstMeanDireStartC);
    
    %% compare field neurons vs. no field neuron in cluster 1
    colorSel = 1;   
    plotBarPlot(autoCorrPyrAll.relDepthNeuHDef(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu1}),...
        autoCorrPyrAll.relDepthNeuHDef(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu1}),'Depth',...
        'Pyr_RelDepthNeuHDefCmpC1FieldBar',pathAnal,[-5 6],modAlignedPyrStatsFNeuVsNoFNeu.pKWRelDepthNeuHDef(clu1),colorSel);
    
    plotBoxPlot(mod.fractBurst(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu1}),...
        mod.fractBurst(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu1}),'Fract burst',...
        'Pyr_FractBurstC1FieldBox',pathAnal,[-0.1 0.85],modAlignedPyrStatsFNeuVsNoFNeu.pRSMeanFractBurst(clu1),colorSel);
    
    plotBoxPlot(mod.numSpPerBurstMean(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu1}),...
        mod.numSpPerBurstMean(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu1}),'Burst length',...
        'Pyr_NumSpPerBurstMeanC1FieldBox',pathAnal,[1.9 3.5],modAlignedPyrStatsFNeuVsNoFNeu.pRSNumSpPerBurstMean(clu1),colorSel);
    
    plotBoxPlot(mod.thetaModHist(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu1}),...
        mod.thetaModHist(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu1}),'Theta modulation',...
        'Pyr_ThetaModHistC1FieldBox',pathAnal,[0 1.1],modAlignedPyrStatsFNeuVsNoFNeu.pRSThetaModHist(clu1),colorSel);
    
    x = mod.phaseMeanDire;
    x(x > pi) = x(x > pi) - 2*pi;
    x = x/pi*180;
    plotBoxPlot(x(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu1}),...
        x(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu1}),'Mean theta phase (deg.)',...
        'Pyr_ThetaMeanC1FieldBox',pathAnal,[-180 180],modAlignedPyrStatsFNeuVsNoFNeu.pWWPhaseMeanDire(clu1),colorSel);
    
    x = mod.burstMeanDire;
    x(x > pi) = x(x > pi) - 2*pi;
    x = x/pi*180;
    plotBoxPlot(x(modAlignedPyrStatsFNeuVsNoFNeu.indCurCFieldBurst{clu1}),...
        x(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoFieldBurst{clu1}),'Burst mean phase (deg.)',...
        'Pyr_BurstMeanDirC1FieldBox',pathAnal,[-180 180],modAlignedPyrStatsFNeuVsNoFNeu.pWWBurstMeanDire(clu1),colorSel);
    
    plotBoxPlot(mod.burstThetaDiff(modAlignedPyrStatsFNeuVsNoFNeu.indCurCFieldBurst{clu1})/pi*180,...
        mod.burstThetaDiff(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoFieldBurst{clu1})/pi*180,...
        'Burst phase - theta phase (deg.)','Pyr_BurstPhase-ThetaPhaseC1FieldBox',...
        pathAnal,[],modAlignedPyrStatsFNeuVsNoFNeu.pWWBurstThetaDiff(clu1),colorSel);
    
    x = mod.burstMeanDireStart;
    x(x > pi) = x(x > pi) - 2*pi;
    x = x/pi*180;
    plotBoxPlot(x(modAlignedPyrStatsFNeuVsNoFNeu.indCurCFieldBurst{clu1}),...
        x(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoFieldBurst{clu1}),'Burst start phase (deg.)',...
        'Pyr_BurstMeanDirStartC1FieldBox',pathAnal,[-150 180],modAlignedPyrStatsFNeuVsNoFNeu.pWWBurstMeanDireStart(clu1),colorSel);
    
    plotBoxPlot(mod.thetaModFreq3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu1}),...
        mod.thetaModFreq3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu1}),'Theta modulation frequency (Hz)',...
        'Pyr_ThetaModFreqC1FieldBox',pathAnal,[4 10.5],modAlignedPyrStatsFNeuVsNoFNeu.pRSThetaModFreq3(clu1),colorSel);
    
    plotBoxPlot(mod.diffNeuronLFPFreq(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu1}),...
        mod.diffNeuronLFPFreq(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu1}),'Neuron mod freq. - LFP freq. (Hz)',...
        'Pyr_DiffNeuronLFPFreqC1FieldBox',...
        pathAnal,[],modAlignedPyrStatsFNeuVsNoFNeu.pRSDiffNeuronLFPFreq(clu1),colorSel);
    
    plotBoxPlot(mod.thetaAsym3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu1}),...
        mod.thetaAsym3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu1}),'Theta asymmetry',...
        'Pyr_ThetaAsymC1FieldBox',pathAnal,[0.1 0.85],modAlignedPyrStatsFNeuVsNoFNeu.pRSThetaAsym3(clu1),colorSel);
    
    
    %% compare field neurons vs. no field neurons in cluster 2
    colorSel = 2;   
    plotBarPlot(autoCorrPyrAll.relDepthNeuHDef(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu2}),...
        autoCorrPyrAll.relDepthNeuHDef(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu2}),'Depth',...
        'Pyr_RelDepthNeuHDefCmpC2FieldBar',pathAnal,[-5 6],modAlignedPyrStatsFNeuVsNoFNeu.pRSRelDepthNeuHDef(clu2),colorSel);
    
    plotBoxPlot(mod.fractBurst(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu2}),...
        mod.fractBurst(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu2}),'Fract burst',...
        'Pyr_FractBurstC2FieldBox',pathAnal,[-0.1 0.85],modAlignedPyrStatsFNeuVsNoFNeu.pRSMeanFractBurst(clu2),colorSel);
    
    plotBoxPlot(mod.numSpPerBurstMean(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu2}),...
        mod.numSpPerBurstMean(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu2}),'Burst length',...
        'Pyr_NumSpPerBurstMeanC2FieldBox',pathAnal,[1.9 3.5],modAlignedPyrStatsFNeuVsNoFNeu.pRSNumSpPerBurstMean(clu2),colorSel);
    
    plotBoxPlot(mod.thetaModHist(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu2}),...
        mod.thetaModHist(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu2}),'Theta modulation',...
        'Pyr_ThetaModHistC2FieldBox',pathAnal,[0 1.1],modAlignedPyrStatsFNeuVsNoFNeu.pRSThetaModHist(clu2),colorSel);
    
    plotBoxPlot(mod.phaseMeanDire(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu2})/pi*180,...
        mod.phaseMeanDire(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu2})/pi*180,'Mean theta phase (deg.)',...
        'Pyr_ThetaMeanC2FieldBox',pathAnal,[-45 360],modAlignedPyrStatsFNeuVsNoFNeu.pWWPhaseMeanDire(clu2),colorSel);
    
    plotBoxPlot(mod.burstMeanDire(modAlignedPyrStatsFNeuVsNoFNeu.indCurCFieldBurst{clu2})/pi*180,...
        mod.burstMeanDire(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoFieldBurst{clu2})/pi*180,'Burst mean phase (deg.)',...
        'Pyr_BurstMeanDirC2FieldBox',pathAnal,[],modAlignedPyrStatsFNeuVsNoFNeu.pWWBurstMeanDire(clu2),colorSel);
    
    plotBoxPlot(mod.burstThetaDiff(modAlignedPyrStatsFNeuVsNoFNeu.indCurCFieldBurst{clu2})/pi*180,...
        mod.burstThetaDiff(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoFieldBurst{clu2})/pi*180,...
        'Burst phase - theta phase (deg.)','Pyr_BurstPhase-ThetaPhaseC2FieldBox',...
        pathAnal,[],modAlignedPyrStatsFNeuVsNoFNeu.pWWBurstThetaDiff(clu2),colorSel);
    
    plotBoxPlot(mod.burstMeanDireStart(modAlignedPyrStatsFNeuVsNoFNeu.indCurCFieldBurst{clu2})/pi*180,...
        mod.burstMeanDireStart(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoFieldBurst{clu2})/pi*180,'Burst start phase (deg.)',...
        'Pyr_BurstMeanDirStartC2FieldBox',pathAnal,[],modAlignedPyrStatsFNeuVsNoFNeu.pWWBurstMeanDireStart(clu2),colorSel);
    
    plotBoxPlot(mod.thetaModFreq3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu2}),...
        mod.thetaModFreq3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu2}),'Theta modulation frequency (Hz)',...
        'Pyr_ThetaModFreqC2FieldBox',pathAnal,[4.5 10],modAlignedPyrStatsFNeuVsNoFNeu.pRSThetaModFreq3(clu2),colorSel);
    
    plotBoxPlot(mod.diffNeuronLFPFreq(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu2}),...
        mod.diffNeuronLFPFreq(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu2}),'Neuron mod freq. - LFP freq. (Hz)',...
        'Pyr_DiffNeuronLFPFreqC2FieldBox',...
        pathAnal,[],modAlignedPyrStatsFNeuVsNoFNeu.pRSDiffNeuronLFPFreq(clu2),colorSel);
    
    plotBoxPlot(mod.thetaAsym3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{clu2}),...
        mod.thetaAsym3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{clu2}),'Theta asymmetry',...
        'Pyr_ThetaAsymC2FieldBox',pathAnal,[0.1 0.85],modAlignedPyrStatsFNeuVsNoFNeu.pRSThetaAsym3(clu2),colorSel);
    
    %% plot clusters with fields
    indCurCField1 = mod.nNeuWithField > 1 & mod.isNeuWithField == 1;
    indCurCNoField1 = mod.nNeuWithField <= 1 & mod.isNeuWithField == 1;
    indField = mod.nNeuWithField > 1;

    indCurCField2 = mod.isNeuWithField == 1;
    burstDiffDouble = [mod.burstThetaDiff/pi*180, mod.burstThetaDiff/pi*180];
    thetaDouble = [mod.phaseMeanDire/pi*180, mod.phaseMeanDire/pi*180+360];
    plotClustersFieldsPhaseLabelF(pathAnal,burstDiffDouble,thetaDouble,...
        [idxC idxC],[],[],...
        [],[],1:2*length(indCurCField2),...
        'Burst mean phase-Mean theta phase','Mean theta phase','All tasks',...
        'Pyr_BurstThetaDiffVsThetaMeanAll');
    set(gca,'Xlim',[-180 180],'XTick',[-180 0 180],'YLim',[0 720],'YTick',[0 360 720])
    plot([0 0],[0 720],'k-');
    fileName1 = [pathAnal 'Pyr_BurstThetaDiffVsThetaMeanAll'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    burstDiffDoubleField = [mod.burstThetaDiff(indCurCField2)/pi*180,...
        mod.burstThetaDiff(indCurCField2)/pi*180];
    thetaDoubleField = [mod.phaseMeanDire(indCurCField2)/pi*180,...
        mod.phaseMeanDire(indCurCField2)/pi*180+360];
    plotClustersFieldsPhaseLabelF(pathAnal,burstDiffDouble,thetaDouble,...
        [idxC idxC],burstDiffDoubleField,thetaDoubleField,...
        [],[],1:2*length(indCurCField2),...
        'Burst mean phase-Mean theta phase','Mean theta phase','All tasks',...
        'Pyr_BurstThetaDiffVsThetaMeanAll');
    set(gca,'Xlim',[-180 180],'XTick',[-180 0 180],'YLim',[0 720],'YTick',[0 360 720])
    plot([0 0],[0 720],'k-');
    fileName1 = [pathAnal 'Pyr_BurstThetaDiffVsThetaMeanAllWithField'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotClustersFieldsPhaseLabelF(pathAnal,mod.phaseMeanDire/pi*180,mod.thetaModFreq3,...
        idxC,mod.phaseMeanDire(indCurCField2)/pi*180,mod.thetaModFreq3(indCurCField2),...
        [],[],1:length(indCurCField2),...
        'Mean theta phase','Theta modulation freq. (Hz)','All tasks label Fields',...
        'Pyr_ThetaMeanVsModFreqAllVsField');

    %% compare neurons in the recordings with fields vs. in the recordings with no field in cluster 1
    colorSel = 1;   
    
    plotBarPlot(autoCorrPyrAll.relDepthNeuHDef(modPyrStatsField.indCurCField{clu1}),...
        autoCorrPyrAll.relDepthNeuHDef(modPyrStatsField.indCurCNoField{clu1}),'Depth',...
        'Pyr_RelDepthNeuHDefCmpC1FieldRecBar',pathAnal,[-5 6],modPyrStatsField.pKWRelDepthNeuHDef(clu1),colorSel);
        
    plotBoxPlot(mod.fractBurst(modPyrStatsField.indCurCField{clu1}),...
        mod.fractBurst(modPyrStatsField.indCurCNoField{clu1}),'Fract burst',...
        'Pyr_FractBurstC1FieldRecBox',pathAnal,[-0.1 0.85],modPyrStatsField.pRSFractBurst(clu1),colorSel);
    
    plotBoxPlot(mod.numSpPerBurstMean(modPyrStatsField.indCurCField{clu1}),...
        mod.numSpPerBurstMean(modPyrStatsField.indCurCNoField{clu1}),'Burst length',...
        'Pyr_NumSpPerBurstMeanC1FieldRecBox',pathAnal,[1.9 3.5],modPyrStatsField.pRSNumSpPerBurstMean(clu1),colorSel);
    
    plotBoxPlot(mod.thetaModHist(modPyrStatsField.indCurCField{clu1}),...
        mod.thetaModHist(modPyrStatsField.indCurCNoField{clu1}),'Theta modulation',...
        'Pyr_ThetaModHistC1FieldRecBox',pathAnal,[0 1.2],modPyrStatsField.pRSThetaModHist(clu1),colorSel);
    
    x = mod.phaseMeanDire;
    x(x > pi) = x(x > pi) - 2*pi;
    x = x/pi*180;
    plotBoxPlot(x(modPyrStatsField.indCurCField{clu1}),...
        x(modPyrStatsField.indCurCNoField{clu1}),'Mean theta phase (deg.)',...
        'Pyr_ThetaMeanC1FieldRecBox',pathAnal,[-120 180],modPyrStatsField.pWWPhaseMeanDire(clu1),colorSel);
    
    x = mod.burstMeanDire;
    x(x > pi) = x(x > pi) - 2*pi;
    x = x/pi*180;
    plotBoxPlot(x(modPyrStatsField.indCurCFieldBurst{clu1}),...
        x(modPyrStatsField.indCurCNoFieldBurst{clu1}),'Burst mean phase (deg.)',...
        'Pyr_BurstMeanDirC1FieldRecBox',pathAnal,[-180 180],modPyrStatsField.pWWBurstMeanDire(clu1),colorSel);
    
    x = mod.burstMeanDireStart;
    x(x > pi) = x(x > pi) - 2*pi;
    x = x/pi*180;
    plotBoxPlot(x(modPyrStatsField.indCurCFieldBurst{clu1}),...
        x(modPyrStatsField.indCurCNoFieldBurst{clu1}),'Burst start phase (deg.)',...
        'Pyr_BurstMeanDirStartC1FieldRecBox',pathAnal,[-180 180],modPyrStatsField.pWWBurstMeanDireStart(clu1),colorSel);
    
    plotBoxPlot(mod.burstThetaDiff(modPyrStatsField.indCurCFieldBurst{clu1})/pi*180,...
        mod.burstThetaDiff(modPyrStatsField.indCurCNoFieldBurst{clu1})/pi*180,...
        'Burst phase - theta phase (deg.)','Pyr_BurstPhase-ThetaPhaseC1FieldRecBox',...
        pathAnal,[],modPyrStatsField.pWWBurstThetaDiff(clu1),colorSel);
    
    plotBoxPlot(mod.thetaModFreq3(modPyrStatsField.indCurCField{clu1}),...
        mod.thetaModFreq3(modPyrStatsField.indCurCNoField{clu1}),'Theta modulation frequency (Hz)',...
        'Pyr_ThetaModFreqC1FieldRecBox',pathAnal,[4 10.5],modPyrStatsField.pRSThetaModFreq3(clu1),colorSel);
    
    plotBoxPlot(mod.diffNeuronLFPFreq(modPyrStatsField.indCurCField{clu1}),...
        mod.diffNeuronLFPFreq(modPyrStatsField.indCurCNoField{clu1}),'Neuron mod freq. - LFP freq. (Hz)',...
        'Pyr_DiffNeuronLFPFreqC1FieldRecBox',...
        pathAnal,[],modPyrStatsField.pRSDiffNeuronLFPFreq(clu1),colorSel);
    
    plotBoxPlot(mod.thetaAsym3(modPyrStatsField.indCurCField{clu1}),...
        mod.thetaAsym3(modPyrStatsField.indCurCNoField{clu1}),'Theta asymmetry',...
        'Pyr_ThetaAsymC1FieldRecBox',pathAnal,[0.1 0.85],modPyrStatsField.pRSThetaAsym3(clu1),colorSel);
    
    %% compare neurons in the recordings with fields vs. in the recordings with no field in cluster 2
    colorSel = 2;  
    plotBarPlot(autoCorrPyrAll.relDepthNeuHDef(modPyrStatsField.indCurCField{clu2}),...
        autoCorrPyrAll.relDepthNeuHDef(modPyrStatsField.indCurCNoField{clu2}),'Depth',...
        'Pyr_RelDepthNeuHDefCmpC2FieldRecBar',pathAnal,[-5 6],modPyrStatsField.pKWRelDepthNeuHDef(clu2),colorSel);
        
    plotBoxPlot(mod.fractBurst(modPyrStatsField.indCurCField{clu2}),...
        mod.fractBurst(modPyrStatsField.indCurCNoField{clu2}),'Fract burst',...
        'Pyr_FractBurstC2FieldRecBox',pathAnal,[-0.1 0.85],modPyrStatsField.pRSFractBurst(clu2),colorSel);
    
    plotBoxPlot(mod.numSpPerBurstMean(modPyrStatsField.indCurCField{clu2}),...
        mod.numSpPerBurstMean(modPyrStatsField.indCurCNoField{clu2}),'Burst length',...
        'Pyr_NumSpPerBurstMeanC2FieldRecBox',pathAnal,[1.9 3.5],modPyrStatsField.pRSNumSpPerBurstMean(clu2),colorSel);
    
    plotBoxPlot(mod.thetaModHist(modPyrStatsField.indCurCField{clu2}),...
        mod.thetaModHist(modPyrStatsField.indCurCNoField{clu2}),'Theta modulation',...
        'Pyr_ThetaModHistC2FieldRecBox',pathAnal,[0 1.1],modPyrStatsField.pRSThetaModHist(clu2),colorSel);
    
    plotBoxPlot(mod.phaseMeanDire(modPyrStatsField.indCurCField{clu2})/pi*180,...
        mod.phaseMeanDire(modPyrStatsField.indCurCNoField{clu2})/pi*180,'Mean theta phase (deg.)',...
        'Pyr_ThetaMeanC2FieldRecBox',pathAnal,[-45 360],modPyrStatsField.pWWPhaseMeanDire(clu2),colorSel);
    
    plotBoxPlot(mod.burstMeanDire(modPyrStatsField.indCurCFieldBurst{clu2})/pi*180,...
        mod.burstMeanDire(modPyrStatsField.indCurCNoFieldBurst{clu2})/pi*180,'Burst mean phase (deg.)',...
        'Pyr_BurstMeanDirC2FieldRecBox',pathAnal,[],modPyrStatsField.pWWBurstMeanDire(clu2),colorSel);
    
    plotBoxPlot(mod.burstMeanDireStart(modPyrStatsField.indCurCFieldBurst{clu2})/pi*180,...
        mod.burstMeanDireStart(modPyrStatsField.indCurCNoFieldBurst{clu2})/pi*180,'Burst start phase (deg.)',...
        'Pyr_BurstMeanDirStartC2FieldRecBox',pathAnal,[],modPyrStatsField.pWWBurstMeanDireStart(clu2),colorSel);
    
    plotBoxPlot(mod.burstThetaDiff(modPyrStatsField.indCurCFieldBurst{clu2})/pi*180,...
        mod.burstThetaDiff(modPyrStatsField.indCurCNoFieldBurst{clu2})/pi*180,...
        'Burst phase - theta phase (deg.)','Pyr_BurstPhase-ThetaPhaseC2FieldRecBox',...
        pathAnal,[],modPyrStatsField.pWWBurstThetaDiff(clu2),colorSel);
    
    plotBoxPlot(mod.thetaModFreq3(modPyrStatsField.indCurCField{clu2}),...
        mod.thetaModFreq3(modPyrStatsField.indCurCNoField{clu2}),'Theta modulation frequency (Hz)',...
        'Pyr_ThetaModFreqC2FieldRecBox',pathAnal,[4.5 10.5],modPyrStatsField.pRSThetaModFreq3(clu2),colorSel);
    
    plotBoxPlot(mod.diffNeuronLFPFreq(modPyrStatsField.indCurCField{clu2}),...
        mod.diffNeuronLFPFreq(modPyrStatsField.indCurCNoField{clu2}),'Neuron mod freq. - LFP freq. (Hz)',...
        'Pyr_DiffNeuronLFPFreqC2FieldRecBox',...
        pathAnal,[],modPyrStatsField.pRSDiffNeuronLFPFreq(clu2),colorSel);
    
    plotBoxPlot(mod.thetaAsym3(modPyrStatsField.indCurCField{clu2}),...
        mod.thetaAsym3(modPyrStatsField.indCurCNoField{clu2}),'Theta asymmetry',...
        'Pyr_ThetaAsymC2FieldRecBox',pathAnal,[0.1 0.85],modPyrStatsField.pRSThetaAsym3(clu2),colorSel);
    
    %%
    plotClustersFields(pathAnal,mod.thetaModHist,mod.phaseMeanDire,...
        idxC,indField,...
        'Theta modulation (hist)','Theta phase mean direction','All tasks',...
        'Pyr_ThetaModHistVsThetaMean');
    
end

function modPyr = accumPyrNeurons1(paths,filenames,mazeSess,minFR,maxFR,task,methodTheta,onlyRun)

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
                              'nNeuWithFieldAligned',[],... % number of neurons with fields in the recording
                              'isNeuWithFieldAligned',[],... % does the neuron have field(s)
                              'fieldWidthAligned',[],... % field width
                              'indStartFieldAligned',[],... % start index of a field
                              'indPeakFieldAligned',[],... % peak index of a field
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
                              'phaseMeanResultantLen',[]); % the mean resultant length of the mean phase direction
                          
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
        
        fileNameFW = [filenames(i,:) '_FieldSpCorr_GoodTr_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameFW];
        if(exist(fullPath) == 0)
            disp(['The field detection file does not exist. Please call ',...
                    'function "FieldDetection_GoodTr" first.']);
            return;
        end
        load(fullPath,'fieldSpCorrSessGoodTr','paramF'); 
        if(mazeSess(i) == 0)
            mazeSessTmp = 1;
        else
            mazeSessTmp = mazeSess(i);
        end
        fieldSpCorrSessGoodTr = fieldSpCorrSessGoodTr{mazeSessTmp};
        
        fileNameFW = [filenames(i,:) '_FieldSpCorrAligned_Run' num2str(mazeSess(i)) ...
                        '_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameFW];
        if(exist(fullPath) == 0)
            disp(['The field detection file does not exist. Please call ',...
                    'function "FieldDetectionAligned" first.']);
            return;
        end
        load(fullPath,'fieldSpCorrSessNonStimGood');   
                
        fullPathFR = [filenames(i,:) '_FR_Run' num2str(onlyRun) '.mat'];
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
        
        fileNamePeakFR = [filenames(i,:) '_PeakFRAligned_msess' num2str(mazeSess(i)) ...
                        '_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNamePeakFR];
        if(exist(fullPath) == 0)
            disp(['The peak firing rate file does not exist. Please call ',...
                    'function "PeakFiringRate_Aligned" first.']);
            return;
        end
        load(fullPath,'pFRNonStimGoodStruct');
        
        fileNameThetaMod = [filenames(i,:) '_ThetaMod_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameThetaMod];
        if(exist(fullPath) == 0)
            disp('_ThetaMod file does not exist.');
            return;
        end
        load(fullPath,'thetaModSess');
        if(length(beh.mazeSessAll) == 1)
           thetaModSessTmp = thetaModSess{1};
        else
           thetaModSessTmp = thetaModSess{mazeSess(i)}; 
        end
        
        if(methodTheta == 0)
            th = 'H';
        else
            th = 'L';
        end
        fileNameBurst = [filenames(i,:) '_burstAll_TH' th '_Run' num2str(onlyRun) ...
                     '.mat'];
        fullPath = [paths(i,:) fileNameBurst];
        if(exist(fullPath) == 0)
            disp('_bustAll file does not exist.');
            return;
        end
        load(fullPath,'burstIsiPerNeuron','burstIsiPerNeuronSess');
        if(length(beh.mazeSessAll) > 1)
            burstIsi = burstIsiPerNeuronSess{mazeSess(i)};
        else
            burstIsi = burstIsiPerNeuron;
        end
        
        fileNameThetaPhase = [filenames(i,:) '_ThetaPhaseL_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameThetaPhase];
        if(exist(fullPath) == 0)
            disp('_ThetaPhaseL file does not exist.');
            return;
        end
        load(fullPath,'spikeThetaPhaseStruct','spikeThetaPhaseStructSess');
        if(length(beh.mazeSessAll) > 1)
            spikeThetaPhase = spikeThetaPhaseStructSess{mazeSess(i)};
        else
            spikeThetaPhase = spikeThetaPhaseStruct;
        end
        
        fileNameThetaPhase = [filenames(i,:) '_ThetaPhaseH_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameThetaPhase];
        if(exist(fullPath) == 0)
            disp('_ThetaPhaseH file does not exist.');
            return;
        end
        load(fullPath,'spikeThetaPhaseStruct','spikeThetaPhaseStructSess');
        if(length(beh.mazeSessAll) > 1)
            spikeThetaPhaseH = spikeThetaPhaseStructSess{mazeSess(i)};
        else
            spikeThetaPhaseH = spikeThetaPhaseStruct;
        end
        trialLenMean = mean(beh.lenTrials(spikeThetaPhaseH.indLapList));
        
%         indNeu = cluList.firingRate > minFR & cluList.firingRate < maxFR &...
%                     autoCorr.isPyrneuron == 1;
        indNeu = mFR.mFR > minFR & mFR.mFR < maxFR &...
                    autoCorr.isPyrneuron == 1;
        modPyr.task = [modPyr.task task*ones(1,sum(indNeu))];
        modPyr.indRec = [modPyr.indRec i*ones(1,sum(indNeu))];
        modPyr.indNeu = [modPyr.indNeu find(indNeu == 1)]; 
        modPyr.trialLenMean = [modPyr.trialLenMean trialLenMean*ones(1,sum(indNeu))];
        modPyr.thetaFreqHMean = [modPyr.thetaFreqHMean mean(beh.thetaFreqHMean(spikeThetaPhaseH.indLapList))*ones(1,sum(indNeu))];
                
        modPyr.thetaMod = [modPyr.thetaMod thetaModSessTmp.thetaMod(indNeu)];
        modPyr.trough = [modPyr.trough thetaModSessTmp.trough(indNeu)];
        modPyr.peak = [modPyr.peak thetaModSessTmp.peak(indNeu)];
        modPyr.thetaModInd = [modPyr.thetaModInd thetaModSessTmp.thetaModInd(indNeu)];
        modPyr.peakT3 = [modPyr.peakT3 thetaModSessTmp.peakT3(indNeu)];
        modPyr.troughT3 = [modPyr.troughT3 thetaModSessTmp.troughT3(indNeu)];
        modPyr.thetaAsym3 = [modPyr.thetaAsym3 ...
            (abs(thetaModSessTmp.peakT3(indNeu))-abs(thetaModSessTmp.troughT3(indNeu)))./...
            (abs(thetaModSessTmp.peakT3(indNeu)))];
        modPyr.thetaModFreq3 = [modPyr.thetaModFreq3 ...
            1000./(abs(thetaModSessTmp.peakT3(indNeu)))];
        modPyr.thetaModInd3 = [modPyr.thetaModInd3 thetaModSessTmp.thetaModInd3(indNeu)];
        modPyr.mFR = [modPyr.mFR mFR.mFR(indNeu)];
        
        %% neurons with fields
        nNeurons = length(cluList.firingRate);
        nFieldArr = zeros(1,sum(indNeu));
        indNeuWithField = zeros(1,nNeurons);
        fieldWidth = zeros(1,nNeurons);
        indStartField = zeros(1,nNeurons);
        indPeakField = zeros(1,nNeurons);
        if(length(pFRNonStimGoodStruct.indLapList) > paramF.minNumTr) % more than 15 trirals
            if(~isempty(fieldSpCorrSessGoodTr))
                [indNeuF,ia] = unique(fieldSpCorrSessGoodTr.indNeuron); 
                indNeuWithField(indNeuF) = 1;
                fieldWidth(indNeuF) = fieldSpCorrSessGoodTr.FW(ia);
                indStartField(indNeuF) = fieldSpCorrSessGoodTr.indStartField(ia);
                indPeakField(indNeuF) = fieldSpCorrSessGoodTr.indPeakField(ia);
                numField = length(indNeuF);
                if(numField > 0)
                    nFieldArr = numField * ones(1,sum(indNeu));
                end
            end
        end
        modPyr.nNeuWithField = [modPyr.nNeuWithField nFieldArr];
        modPyr.isNeuWithField = [modPyr.isNeuWithField indNeuWithField(indNeu)]; % the total number of neurons with field for each recording here might not be the same as in nNeuWithField, because of the subselection of neurons
        modPyr.fieldWidth = [modPyr.fieldWidth fieldWidth(indNeu)];
        modPyr.indStartField = [modPyr.indStartField indStartField(indNeu)];
        modPyr.indPeakField = [modPyr.indPeakField indPeakField(indNeu)];
        modPyr.meanInstFR = [modPyr.meanInstFR pFRNonStimGoodStruct.meanInstFR(indNeu)];
        
        %% neurons with fields after alignment
        nNeurons = length(cluList.firingRate);
        nFieldArr = zeros(1,sum(indNeu));
        indNeuWithField = zeros(1,nNeurons);
        fieldWidth = zeros(1,nNeurons);
        indStartField = zeros(1,nNeurons);
        indPeakField = zeros(1,nNeurons);
        if(length(pFRNonStimGoodStruct.indLapList) > paramF.minNumTr) % more than 15 trirals
            if(~isempty(fieldSpCorrSessNonStimGood))
                [indNeuF,ia] = unique(fieldSpCorrSessNonStimGood.indNeuron); 
                indNeuWithField(indNeuF) = 1;
                fieldWidth(indNeuF) = fieldSpCorrSessNonStimGood.FW(ia);
                indStartField(indNeuF) = fieldSpCorrSessNonStimGood.indStartField(ia);
                indPeakField(indNeuF) = fieldSpCorrSessNonStimGood.indPeakField(ia);
                numField = length(indNeuF);
                if(numField > 0)
                    nFieldArr = numField * ones(1,sum(indNeu));
                end
            end
        end
        modPyr.nNeuWithFieldAligned = [modPyr.nNeuWithFieldAligned nFieldArr];
        modPyr.isNeuWithFieldAligned = [modPyr.isNeuWithFieldAligned indNeuWithField(indNeu)];
        modPyr.fieldWidthAligned = [modPyr.fieldWidthAligned fieldWidth(indNeu)];
        modPyr.indStartFieldAligned = [modPyr.indStartFieldAligned indStartField(indNeu)];
        modPyr.indPeakFieldAligned = [modPyr.indPeakFieldAligned indPeakField(indNeu)];
        
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

function modPyrStatsField = modPyrStatsFieldPerC(mod,autoCorrPyrAll)
    modPyrStatsField = [];
    idxC = mod.idxC;
    for i = 1:max(idxC)
        indCurCField = idxC == i & mod.nNeuWithField > 1;
        indCurCNoField = idxC == i & mod.nNeuWithField < 1;
        modPyrStatsField.indCurCField{i} = indCurCField;
        modPyrStatsField.indCurCNoField{i} = indCurCNoField;
        modPyrStatsField.isNeuWithField(i,:) = [sum(mod.isNeuWithField(indCurCField)),sum(mod.isNeuWithField(indCurCNoField))];
              
        modPyrStatsField.meanPeakTo40ms(i,:) = [mean(autoCorrPyrAll.peakTo40ms(indCurCField)),mean(autoCorrPyrAll.peakTo40ms(indCurCNoField))];
        modPyrStatsField.meanPeakTime(i,:) = [mean(autoCorrPyrAll.peakTime(indCurCField)),mean(autoCorrPyrAll.peakTime(indCurCNoField))];
        
        modPyrStatsField.mFR(i,:) = [mean(mod.mFR(indCurCField)),mean(mod.mFR(indCurCNoField))];
        modPyrStatsField.meanInstFR(i,:) = [mean(mod.meanInstFR(indCurCField)),mean(mod.meanInstFR(indCurCNoField))];
        
        modPyrStatsField.meanDiffNeuronLFPFreq(i,:) = [mean(mod.diffNeuronLFPFreq(indCurCField)),mean(mod.diffNeuronLFPFreq(indCurCNoField))];
                
        indCurCFieldBurst = idxC == i & mod.nNeuWithField > 1 & mod.fractBurst > 0;
        indCurCNoFieldBurst = idxC == i & mod.nNeuWithField <= 1 & mod.fractBurst > 0;
        modPyrStatsField.indCurCFieldBurst{i} = indCurCFieldBurst;
        modPyrStatsField.indCurCNoFieldBurst{i} = indCurCNoFieldBurst;
        modPyrStatsField.meanBurstMeanDire(i,:) = [circ_mean(mod.burstMeanDire(indCurCFieldBurst)'),circ_mean(mod.burstMeanDire(indCurCNoFieldBurst)')];
        modPyrStatsField.meanNonBurstMeanDire(i,:) = [circ_mean(mod.nonBurstMeanDire(indCurCField)'),circ_mean(mod.nonBurstMeanDire(indCurCNoField)')];     
        modPyrStatsField.meanBurstMeanDireStart(i,:) = [circ_mean(mod.burstMeanDireStart(indCurCFieldBurst)'),circ_mean(mod.burstMeanDireStart(indCurCNoFieldBurst)')];
        modPyrStatsField.meanFractBurst(i,:) = [mean(mod.fractBurst(indCurCField)),mean(mod.fractBurst(indCurCNoField))];
        nonZeroF = mod.numSpPerBurstMean(indCurCField);
        nonZeroF = nonZeroF(nonZeroF > 0);
        nonZeroNoF = mod.numSpPerBurstMean(indCurCNoField);
        nonZeroNoF = nonZeroNoF(nonZeroNoF > 0);
        modPyrStatsField.meanNumSpPerBurstMean(i,:) = [mean(nonZeroF),mean(nonZeroNoF)];
        
        modPyrStatsField.meanPhaseMeanDire(i,:) = [circ_mean(mod.phaseMeanDire(indCurCField)'),circ_mean(mod.phaseMeanDire(indCurCNoField)')];
        modPyrStatsField.meanPhaseMeanDireH(i,:) = [circ_mean(mod.phaseMeanDireH(indCurCField)'),circ_mean(mod.phaseMeanDireH(indCurCNoField)')];
        modPyrStatsField.meanMaxPhase(i,:) = [circ_mean(mod.maxPhaseArr(indCurCField)'/180*pi),circ_mean(mod.maxPhaseArr(indCurCNoField)'/180*pi)];
        modPyrStatsField.meanMaxPhaseH(i,:) = [circ_mean(mod.maxPhaseArrH(indCurCField)'/180*pi),circ_mean(mod.maxPhaseArrH(indCurCNoField)'/180*pi)];
        modPyrStatsField.meanminPhase(i,:) = [circ_mean(mod.minPhaseArr(indCurCField)'/180*pi),circ_mean(mod.minPhaseArr(indCurCNoField)'/180*pi)];
        modPyrStatsField.meanminPhaseH(i,:) = [circ_mean(mod.minPhaseArrH(indCurCField)'/180*pi),circ_mean(mod.minPhaseArrH(indCurCNoField)'/180*pi)];
        modPyrStatsField.meanPhaseDiff(i,:) = [mean(mod.phaseDiff(indCurCField)),mean(mod.phaseDiff(indCurCNoField))];
        modPyrStatsField.meanPhaseDiffH(i,:) = [mean(mod.phaseDiffH(indCurCField)),mean(mod.phaseDiffH(indCurCNoField))];
        modPyrStatsField.meanThetaModHist(i,:) = [mean(mod.thetaModHist(indCurCField)),mean(mod.thetaModHist(indCurCNoField))];
        modPyrStatsField.meanThetaModHistH(i,:) = [mean(mod.thetaModHistH(indCurCField)),mean(mod.thetaModHistH(indCurCNoField))];
        
        modPyrStatsField.meanThetaModFreq3(i,:) = [mean(mod.thetaModFreq3(indCurCField)),mean(mod.thetaModFreq3(indCurCNoField))];
        modPyrStatsField.meanThetaAsym3(i,:) = [mean(mod.thetaAsym3(indCurCField)),mean(mod.thetaAsym3(indCurCNoField))];
        modPyrStatsField.meanThetaModInd3(i,:) = [mean(mod.thetaModInd3(indCurCField)),mean(mod.thetaModInd3(indCurCNoField))];
        modPyrStatsField.meanThetaModInd(i,:) = [mean(mod.thetaModInd(indCurCField)),mean(mod.thetaModInd(indCurCNoField))];
        
        modPyrStatsField.pRSDiffNeuronLFPFreq(i) = ranksum(mod.diffNeuronLFPFreq(indCurCField),mod.diffNeuronLFPFreq(indCurCNoField));
        
        modPyrStatsField.pRSPeakTo40ms(i) = ranksum(autoCorrPyrAll.peakTo40ms(indCurCField),autoCorrPyrAll.peakTo40ms(indCurCNoField));
        modPyrStatsField.pRSPeakTime(i) = ranksum(autoCorrPyrAll.peakTime(indCurCField),autoCorrPyrAll.peakTime(indCurCNoField));
        modPyrStatsField.pRSRelDepthNeuHDef(i) = ranksum(autoCorrPyrAll.relDepthNeuHDef(indCurCField),...
                            autoCorrPyrAll.relDepthNeuHDef(indCurCNoField));
        modPyrStatsField.pKWRelDepthNeuHDef(i) = kruskalwallis([autoCorrPyrAll.relDepthNeuHDef(indCurCField),...
                            autoCorrPyrAll.relDepthNeuHDef(indCurCNoField)],...
                            [ones(1,sum(indCurCField)),2*ones(1,sum(indCurCNoField))]);
                        
        modPyrStatsField.pRSMFR(i) = ranksum(mod.mFR(indCurCField),mod.mFR(indCurCNoField));
        modPyrStatsField.pRSMeanInstFR(i) = ranksum(mod.meanInstFR(indCurCField),mod.meanInstFR(indCurCNoField));
        
        modPyrStatsField.pKBurstMeanDire(i) = circ_ktest(mod.burstMeanDire(indCurCFieldBurst)',mod.burstMeanDire(indCurCNoFieldBurst)');
        modPyrStatsField.pKNonBurstMeanDire(i) = circ_ktest(mod.nonBurstMeanDire(indCurCField)',mod.nonBurstMeanDire(indCurCNoField)');
        modPyrStatsField.pKBurstMeanDireStart(i) = circ_ktest(mod.burstMeanDireStart(indCurCFieldBurst)',mod.burstMeanDireStart(indCurCNoFieldBurst)');
        
        modPyrStatsField.pWWBurstMeanDire(i) = circ_wwtest(mod.burstMeanDire(indCurCFieldBurst)',mod.burstMeanDire(indCurCNoFieldBurst)');
        modPyrStatsField.pWWNonBurstMeanDire(i) = circ_wwtest(mod.nonBurstMeanDire(indCurCField)',mod.nonBurstMeanDire(indCurCNoField)');
        modPyrStatsField.pWWBurstMeanDireStart(i) = circ_wwtest(mod.burstMeanDireStart(indCurCFieldBurst)',mod.burstMeanDireStart(indCurCNoFieldBurst)');
        
        modPyrStatsField.pRSFractBurst(i) = ranksum(mod.fractBurst(indCurCField),mod.fractBurst(indCurCNoField));        
        modPyrStatsField.pRSNumSpPerBurstMean(i) = ranksum(nonZeroF,nonZeroNoF);
        
        modPyrStatsField.pWWBurstThetaDiff(i) = circ_wwtest(mod.burstThetaDiff(indCurCFieldBurst),...
                mod.burstThetaDiff(indCurCNoFieldBurst));
        modPyrStatsField.pKBurstThetaDiff(i) = circ_ktest(mod.burstThetaDiff(indCurCFieldBurst),...
                mod.burstThetaDiff(indCurCNoFieldBurst));
        
        modPyrStatsField.pKPhaseMeanDire(i) = circ_ktest(mod.phaseMeanDire(indCurCField)',mod.phaseMeanDire(indCurCNoField)');
        modPyrStatsField.pKPhaseMeanDireH(i) = circ_ktest(mod.phaseMeanDireH(indCurCField)',mod.phaseMeanDireH(indCurCNoField)');
        modPyrStatsField.pKMaxPhase(i) = circ_ktest(mod.maxPhaseArr(indCurCField)'/180*pi,mod.maxPhaseArr(indCurCNoField)'/180*pi);
        modPyrStatsField.pKMaxPhaseH(i) = circ_ktest(mod.maxPhaseArrH(indCurCField)'/180*pi,mod.maxPhaseArrH(indCurCNoField)'/180*pi);
        modPyrStatsField.pKMinPhase(i) = circ_ktest(mod.minPhaseArr(indCurCField)'/180*pi,mod.minPhaseArr(indCurCNoField)'/180*pi);
        modPyrStatsField.pKMinPhaseH(i) = circ_ktest(mod.minPhaseArrH(indCurCField)'/180*pi,mod.minPhaseArrH(indCurCNoField)'/180*pi);
        
        modPyrStatsField.pWWPhaseMeanDire(i) = circ_wwtest(mod.phaseMeanDire(indCurCField)',mod.phaseMeanDire(indCurCNoField)');
        modPyrStatsField.pWWPhaseMeanDireH(i) = circ_wwtest(mod.phaseMeanDireH(indCurCField)',mod.phaseMeanDireH(indCurCNoField)');
        modPyrStatsField.pWWMaxPhase(i) = circ_wwtest(mod.maxPhaseArr(indCurCField)'/180*pi,mod.maxPhaseArr(indCurCNoField)'/180*pi);
        modPyrStatsField.pWWMaxPhaseH(i) = circ_wwtest(mod.maxPhaseArrH(indCurCField)'/180*pi,mod.maxPhaseArrH(indCurCNoField)'/180*pi);
        modPyrStatsField.pWWMinPhase(i) = circ_wwtest(mod.minPhaseArr(indCurCField)'/180*pi,mod.minPhaseArr(indCurCNoField)'/180*pi);
        modPyrStatsField.pWWMinPhaseH(i) = circ_wwtest(mod.minPhaseArrH(indCurCField)'/180*pi,mod.minPhaseArrH(indCurCNoField)'/180*pi);
        
        modPyrStatsField.pRSPhaseDiff(i) = ranksum(mod.phaseDiff(indCurCField),mod.phaseDiff(indCurCNoField));
        modPyrStatsField.pRSPhaseDiffH(i) = ranksum(mod.phaseDiffH(indCurCField),mod.phaseDiffH(indCurCNoField));
        modPyrStatsField.pRSThetaModHist(i) = ranksum(mod.thetaModHist(indCurCField),mod.thetaModHist(indCurCNoField));
        modPyrStatsField.pRSThetaModHistH(i) = ranksum(mod.thetaModHistH(indCurCField),mod.thetaModHistH(indCurCNoField));
        
        modPyrStatsField.pRSThetaModFreq3(i) = ranksum(mod.thetaModFreq3(indCurCField),mod.thetaModFreq3(indCurCNoField));
        modPyrStatsField.pRSThetaAsym3(i) = ranksum(mod.thetaAsym3(indCurCField),mod.thetaAsym3(indCurCNoField));  
        modPyrStatsField.pRSThetaModInd3(i) = ranksum(mod.thetaModInd3(indCurCField),mod.thetaModInd3(indCurCNoField));
        modPyrStatsField.pRSThetaModInd(i) = ranksum(mod.thetaModInd(indCurCField),mod.thetaModInd(indCurCNoField));
        
        modPyrStatsField.pKBurstVsNonBurstMeanDire(i) = circ_ktest(mod.burstMeanDire(indCurCFieldBurst)',mod.nonBurstMeanDire(indCurCField)');
        modPyrStatsField.pKBurstVsThetaMeanDire(i) = circ_ktest(mod.burstMeanDire(indCurCFieldBurst)',mod.phaseMeanDire(indCurCField)');
        modPyrStatsField.pKNonBurstVsThetaMeanDire(i) = circ_ktest(mod.nonBurstMeanDire(indCurCField)',mod.phaseMeanDire(indCurCField)');
        modPyrStatsField.pKBurstStartVsThetaMeanDire(i) = circ_ktest(mod.burstMeanDireStart(indCurCFieldBurst)',mod.phaseMeanDire(indCurCField)');
        
        modPyrStatsField.pWWBurstVsNonBurstMeanDire(i) = circ_wwtest(mod.burstMeanDire(indCurCFieldBurst)',mod.nonBurstMeanDire(indCurCField)');
        modPyrStatsField.pWWBurstVsThetaMeanDire(i) = circ_wwtest(mod.burstMeanDire(indCurCFieldBurst)',mod.phaseMeanDire(indCurCField)');
        modPyrStatsField.pWWNonBurstVsThetaMeanDire(i) = circ_wwtest(mod.nonBurstMeanDire(indCurCField)',mod.phaseMeanDire(indCurCField)');
        modPyrStatsField.pWWBurstStartVsThetaMeanDire(i) = circ_wwtest(mod.burstMeanDireStart(indCurCFieldBurst)',mod.phaseMeanDire(indCurCField)');
    end
end

function modPyrStatsField = modPyrStatsClass(mod)
    idxC = mod.idxC;
    modPyrStatsField.pRSDiffNeuronLFPFreq = ranksum(mod.diffNeuronLFPFreq(idxC == 1),mod.diffNeuronLFPFreq(idxC == 2));
    nonZeroF1 = mod.numSpPerBurstMean(idxC == 1);
    nonZeroF1 = nonZeroF1(nonZeroF1 > 0);
    nonZeroF2 = mod.numSpPerBurstMean(idxC == 2);
    nonZeroF2 = nonZeroF2(nonZeroF2 > 0);
    modPyrStatsField.pRSNumSpPerBurstMeanC = ranksum(nonZeroF1,nonZeroF2);
    modPyrStatsField.pRSFractBurstC = ranksum(mod.fractBurst(idxC == 1),mod.fractBurst(idxC == 2));    
    
    modPyrStatsField.pWWBurstThetaDiffC = circ_wwtest(mod.burstThetaDiff(idxC == 1 & mod.fractBurst > 0),...
        mod.burstThetaDiff(idxC == 2 & mod.fractBurst > 0));
    modPyrStatsField.pKBurstThetaDiffC = circ_ktest(mod.burstThetaDiff(idxC == 1 & mod.fractBurst > 0),...
        mod.burstThetaDiff(idxC == 2 & mod.fractBurst > 0));
    
    modPyrStatsField.pWWPhaseMeanDireC = circ_wwtest(mod.phaseMeanDire(idxC == 1)',mod.phaseMeanDire(idxC == 2)');
    modPyrStatsField.pKPhaseMeanDireC = circ_ktest(mod.phaseMeanDire(idxC == 1)',mod.phaseMeanDire(idxC == 2)');
    modPyrStatsField.pWWPhaseMeanDireHC = circ_wwtest(mod.phaseMeanDireH(idxC == 1)',mod.phaseMeanDireH(idxC == 2)');
    modPyrStatsField.pKPhaseMeanDireHC = circ_ktest(mod.phaseMeanDireH(idxC == 1)',mod.phaseMeanDireH(idxC == 2)');
    modPyrStatsField.pWWMaxPhaseC = circ_wwtest(mod.maxPhaseArr(idxC == 1)'/180*pi,mod.maxPhaseArr(idxC == 2)'/180*pi);
    modPyrStatsField.pKMaxPhaseC = circ_ktest(mod.maxPhaseArr(idxC == 1)'/180*pi,mod.maxPhaseArr(idxC == 2)'/180*pi);
    modPyrStatsField.pWWMinPhaseC = circ_wwtest(mod.minPhaseArr(idxC == 1)'/180*pi,mod.minPhaseArr(idxC == 2)'/180*pi);
    modPyrStatsField.pKMinPhaseC = circ_ktest(mod.minPhaseArr(idxC == 1)'/180*pi,mod.minPhaseArr(idxC == 2)'/180*pi);
    modPyrStatsField.pWWburstMeanDireC = circ_wwtest(mod.burstMeanDire(idxC == 1 & mod.fractBurst > 0)',...
        mod.burstMeanDire(idxC == 2 & mod.fractBurst > 0)');
    modPyrStatsField.pKburstMeanDireC = circ_ktest(mod.burstMeanDire(idxC == 1 & mod.fractBurst > 0)',...
        mod.burstMeanDire(idxC == 2 & mod.fractBurst > 0)');
    modPyrStatsField.pWWburstMeanDireStartC = circ_wwtest(mod.burstMeanDireStart(idxC == 1 & mod.fractBurst > 0)',...
        mod.burstMeanDireStart(idxC == 2 & mod.fractBurst > 0)');
    modPyrStatsField.pKburstMeanDirStartC = circ_ktest(mod.burstMeanDireStart(idxC == 1 & mod.fractBurst > 0)',...
        mod.burstMeanDireStart(idxC == 2 & mod.fractBurst > 0)');
    modPyrStatsField.pRSBurstMeanResultantLenC = ranksum(mod.burstMeanResultantLen(idxC == 1),mod.burstMeanResultantLen(idxC == 2));
    modPyrStatsField.pRSPhaseMeanResultantLenC = ranksum(mod.phaseMeanResultantLen(idxC == 1),mod.phaseMeanResultantLen(idxC == 2));
    
    modPyrStatsField.pRSThetaModHistC = ranksum(mod.thetaModHist(idxC == 1),mod.thetaModHist(idxC == 2));
    modPyrStatsField.pRSThetaModHistHC = ranksum(mod.thetaModHistH(idxC == 1),mod.thetaModHistH(idxC == 2));
    modPyrStatsField.pRSThetaModFreq3C = ranksum(mod.thetaModFreq3(idxC == 1),mod.thetaModFreq3(idxC == 2));
    modPyrStatsField.pRSThetaAsym3C = ranksum(mod.thetaAsym3(idxC == 1),mod.thetaAsym3(idxC == 2));
    
    modPyrStatsField.pRSMFRC = ranksum(mod.mFR(idxC == 1),mod.mFR(idxC == 2));
    modPyrStatsField.pRSMeanInstFR = ranksum(mod.meanInstFR(idxC == 1),mod.meanInstFR(idxC == 2));
end

function modPyrStatsFieldF = modPyrStatsFieldNeu(mod,autoCorrPyrAll)
    %% field vs no field neurons for each cluster
    modPyrStatsFieldF = [];
    idxC = mod.idxC;
    for i = 1:max(idxC)        
        indCurCField = idxC == i & mod.isNeuWithField == 1;
        indCurCNoField = idxC == i & mod.isNeuWithField == 0;
        indCurCFieldBurst = idxC == i & mod.isNeuWithField == 1 & mod.fractBurst > 0;
        indCurCNoFieldBurst = idxC == i & mod.isNeuWithField == 0 & mod.fractBurst > 0;
        modPyrStatsFieldF.indCurCField{i} = indCurCField;
        modPyrStatsFieldF.indCurCNoField{i} = indCurCNoField;
        modPyrStatsFieldF.indCurCFieldBurst{i} = indCurCFieldBurst;
        modPyrStatsFieldF.indCurCNoFieldBurst{i} = indCurCNoFieldBurst;
        modPyrStatsFieldF.fieldWidthCFRec{i} = mod.fieldWidth(indCurCField);
        modPyrStatsFieldF.indPeakFieldFRec{i} = mod.indPeakField(indCurCField);
        modPyrStatsFieldF.indStartFieldFRec{i} = mod.indStartField(indCurCField);
        modPyrStatsFieldF.percTrackStartField{i} = mod.percTrackStartField(indCurCField);
        modPyrStatsFieldF.percTrackPeakField{i} = mod.percTrackPeakField(indCurCField);
        modPyrStatsFieldF.skewness{i} = (mod.indPeakField(indCurCField) - mod.indStartField(indCurCField) + 1)./mod.fieldWidth(indCurCField);
      
        modPyrStatsFieldF.meanDiffNeuronLFPFreq(i,:) = [mean(mod.diffNeuronLFPFreq(indCurCField)),mean(mod.diffNeuronLFPFreq(indCurCNoField))];
        modPyrStatsFieldF.meanRelDepthNeuHDef(i,:) = [mean(autoCorrPyrAll.relDepthNeuHDef(indCurCField)),...
                            mean(autoCorrPyrAll.relDepthNeuHDef(indCurCNoField))];
          
        modPyrStatsFieldF.meanPhaseMeanDire(i,:) = [circ_mean(mod.phaseMeanDire(indCurCField)'),circ_mean(mod.phaseMeanDire(indCurCNoField)')];
        modPyrStatsFieldF.meanBurstMeanDire(i,:) = [circ_mean(mod.burstMeanDire(indCurCFieldBurst)'),circ_mean(mod.burstMeanDire(indCurCNoFieldBurst)')];  
        modPyrStatsFieldF.meanNonBurstMeanDire(i,:) = [circ_mean(mod.nonBurstMeanDire(indCurCField)'),circ_mean(mod.nonBurstMeanDire(indCurCNoField)')]; 
        modPyrStatsFieldF.meanBurstMeanDireStart(i,:) = [circ_mean(mod.burstMeanDireStart(indCurCFieldBurst)'),circ_mean(mod.burstMeanDireStart(indCurCNoFieldBurst)')];  
        modPyrStatsFieldF.meanFractBurst(i,:) = [mean(mod.fractBurst(indCurCField)),mean(mod.fractBurst(indCurCNoField))];
        modPyrStatsFieldF.meanNumSpPerBurstMean(i,:) = [mean(mod.numSpPerBurstMean(indCurCField)),mean(mod.numSpPerBurstMean(indCurCNoField))];
        modPyrStatsFieldF.meanBurstMeanResultantLen(i,:) = [mean(mod.burstMeanResultantLen(indCurCFieldBurst)),mean(mod.burstMeanResultantLen(indCurCNoFieldBurst))];
        modPyrStatsFieldF.meanPhaseMeanResultantLen(i,:) = [mean(mod.phaseMeanResultantLen(indCurCField)),mean(mod.phaseMeanResultantLen(indCurCNoField))];
            
        modPyrStatsFieldF.meanThetaModHist(i,:) = [mean(mod.thetaModHist(indCurCField)),mean(mod.thetaModHist(indCurCNoField))];
        modPyrStatsFieldF.meanThetaModHistH(i,:) = [mean(mod.thetaModHistH(indCurCField)),mean(mod.thetaModHistH(indCurCNoField))];
        modPyrStatsFieldF.meanThetaModFreq3(i,:) = [mean(mod.thetaModFreq3(indCurCField)),mean(mod.thetaModFreq3(indCurCNoField))];
        modPyrStatsFieldF.meanThetaAsym3(i,:) = [mean(mod.thetaAsym3(indCurCField)),mean(mod.thetaAsym3(indCurCNoField))];
        
        modPyrStatsFieldF.pRSDiffNeuronLFPFreq(i,:) = ranksum(mod.diffNeuronLFPFreq(indCurCField),mod.diffNeuronLFPFreq(indCurCNoField));
        modPyrStatsFieldF.pRSRelDepthNeuHDef(i) = ranksum(autoCorrPyrAll.relDepthNeuHDef(indCurCField),...
                            autoCorrPyrAll.relDepthNeuHDef(indCurCNoField));
        modPyrStatsFieldF.pKWRelDepthNeuHDef(i) = kruskalwallis([autoCorrPyrAll.relDepthNeuHDef(indCurCField),...
                            autoCorrPyrAll.relDepthNeuHDef(indCurCNoField)],...
                            [ones(1,sum(indCurCField)),2*ones(1,sum(indCurCNoField))]);
        
        modPyrStatsFieldF.pRSThetaModHist(i) = ranksum(mod.thetaModHist(indCurCField),mod.thetaModHist(indCurCNoField));
        modPyrStatsFieldF.pRSThetaModHistH(i) = ranksum(mod.thetaModHistH(indCurCField),mod.thetaModHistH(indCurCNoField));
        modPyrStatsFieldF.pRSThetaModFreq3(i) = ranksum(mod.thetaModFreq3(indCurCField),mod.thetaModFreq3(indCurCNoField));
        modPyrStatsFieldF.pRSThetaAsym3(i) = ranksum(mod.thetaAsym3(indCurCField),mod.thetaAsym3(indCurCNoField));
    
        modPyrStatsFieldF.pKPhaseMeanDire(i) = circ_ktest(mod.phaseMeanDire(indCurCField)',mod.phaseMeanDire(indCurCNoField)');
        modPyrStatsFieldF.pKBurstMeanDire(i) = circ_ktest(mod.burstMeanDire(indCurCFieldBurst)',mod.burstMeanDire(indCurCNoFieldBurst)'); 
        modPyrStatsFieldF.pKNonBurstMeanDire(i) = circ_ktest(mod.nonBurstMeanDire(indCurCField)',mod.nonBurstMeanDire(indCurCNoField)'); 
        modPyrStatsFieldF.pKBurstMeanDireStart(i) = circ_ktest(mod.burstMeanDireStart(indCurCFieldBurst)',mod.burstMeanDireStart(indCurCNoFieldBurst)'); 
        modPyrStatsFieldF.pWWPhaseMeanDire(i) = circ_wwtest(mod.phaseMeanDire(indCurCField)',mod.phaseMeanDire(indCurCNoField)');
        modPyrStatsFieldF.pWWBurstMeanDire(i) = circ_wwtest(mod.burstMeanDire(indCurCFieldBurst)',mod.burstMeanDire(indCurCNoFieldBurst)'); 
        modPyrStatsFieldF.pWWNonBurstMeanDire(i) = circ_wwtest(mod.nonBurstMeanDire(indCurCField)',mod.nonBurstMeanDire(indCurCNoField)');
        modPyrStatsFieldF.pWWBurstMeanDireStart(i) = circ_wwtest(mod.burstMeanDireStart(indCurCFieldBurst)',mod.burstMeanDireStart(indCurCNoFieldBurst)'); 
        modPyrStatsFieldF.pMEDBurstVsThetaMeanDire(i) = circ_medtest(mod.burstMeanDire(indCurCFieldBurst)'- mod.phaseMeanDire(indCurCFieldBurst)',0);  
        modPyrStatsFieldF.pMEDNonBurstVsThetaMeanDire(i) = circ_medtest(mod.nonBurstMeanDire(indCurCField)'- mod.phaseMeanDire(indCurCField)',0);    
        modPyrStatsFieldF.pRSBurstMeanResultantLen(i) = ranksum(mod.burstMeanResultantLen(indCurCFieldBurst),mod.burstMeanResultantLen(indCurCNoFieldBurst));
        modPyrStatsFieldF.pRSPhaseMeanResultantLen(i) = ranksum(mod.phaseMeanResultantLen(indCurCField),mod.phaseMeanResultantLen(indCurCNoField));
        
        modPyrStatsFieldF.pWWBurstThetaDiff(i) = circ_wwtest(mod.burstThetaDiff(indCurCFieldBurst),...
                mod.burstThetaDiff(indCurCNoFieldBurst));
        modPyrStatsFieldF.pKBurstThetaDiff(i) = circ_ktest(mod.burstThetaDiff(indCurCFieldBurst),...
                mod.burstThetaDiff(indCurCNoFieldBurst));
    
        modPyrStatsFieldF.pRSMeanFractBurst(i) = ranksum(mod.fractBurst(indCurCField),mod.fractBurst(indCurCNoField));
        modPyrStatsFieldF.pRSNumSpPerBurstMean(i) = ranksum(mod.numSpPerBurstMean(indCurCField),mod.numSpPerBurstMean(indCurCNoField));
    end
    modPyrStatsFieldF.pRSFieldWidthC = ranksum(modPyrStatsFieldF.fieldWidthCFRec{1},modPyrStatsFieldF.fieldWidthCFRec{2});
    modPyrStatsFieldF.pRSIndPeakFieldFRec = ranksum(modPyrStatsFieldF.indPeakFieldFRec{1},modPyrStatsFieldF.indPeakFieldFRec{2});
    modPyrStatsFieldF.pRSSkewnessC = ranksum(modPyrStatsFieldF.skewness{1},modPyrStatsFieldF.skewness{2});
    modPyrStatsFieldF.pRSThetaModFreq3C = ranksum(mod.thetaModFreq3(modPyrStatsFieldF.indCurCField{1}),...
                mod.thetaModFreq3(modPyrStatsFieldF.indCurCField{2}));
    modPyrStatsFieldF.pRSPercTrackStartFieldC = ranksum(modPyrStatsFieldF.percTrackStartField{1},...
                modPyrStatsFieldF.percTrackStartField{2});
    modPyrStatsFieldF.pRSPercTrackPeakFieldC = ranksum(modPyrStatsFieldF.percTrackPeakField{1},...
                modPyrStatsFieldF.percTrackPeakField{2});
end

function plotBurstVsTheta(burstMeanDire,phaseMeanDire,fractBurst,pathAnal)

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
    fileName1 = [pathAnal 'Pyr-Burst-ThetaMeanPhaseVsThetaPhase'];
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

function plotConditions(x1,y1,x2,y2,x3,y3,xl,yl,ti)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    colorArr = [0.5 0.5 0.9;...
                0.9 0.5 0.5;...
                0.2 0.8 0.5];
    h = plot(x1,y1,'o');
    set(h,'MarkerSize',6,'Color',colorArr(1,:));
    hold on;
    h = plot(x2,y2,'o');
    set(h,'MarkerSize',6,'Color',colorArr(2,:));
    h = plot(x3,y3,'o');
    set(h,'MarkerSize',6,'Color',colorArr(3,:));
    maxX = max([x1,x2,x3]);
    maxY = max([y1,y2,y3]);
    minX = min([x1,x2,x3]);
    minY = min([y1,y2,y3]);
    set(gca,'XLim',[minX maxX],'YLim',[minY maxY]);
    xlabel(xl)
    ylabel(yl)
    title(ti)
    legend('NoCue','AL','PL')
end

function plotClusters(x,y,idx,xl,yl,ti)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])

    % only for plotting two clusters with fields
    colorArr = [...
                163 207 98;...
                234 131 114;...
                ]/255;
            
%     colorArr = [0.5 0.5 0.9;...
%         0.9 0.5 0.5;...
%         0.3 0.3 0.7;...
%         0.7 0.3 0.3;...
%         0.5 0.9 0.5;...
%         0.2 0.5 0.8;...
%         0.2 0.8 0.5;...
%         0.8 0.5 0.2;...
%         0.3 0.7 0.3];
    hold on;
    for i = 1:max(idx)
        indTmp = idx == i;
        h = plot(x(indTmp),y(indTmp),'.');
        set(h,'MarkerSize',8,'Color',colorArr(mod(i,2)+1,:));
    end
%     h = plot(xt,yt,'k+');
%     set(h,'MarkerSize',7);
    maxX = max(x);
    maxY = max(y);
    minX = min(x);
    minY = min(y);
    set(gca,'XLim',[minX maxX],'YLim',[minY maxY]);
    xlabel(xl)
    ylabel(yl)
    title(ti);
end

function plotClustersFieldsPhase(pathAnal,x,y,idx,indField,xl,yl,ti,fileN)
    plotClusters(x(indField),y(indField),...
        idx(indField),xl,yl,[ti ' - field']);
    plot([0 2*pi],[0 2*pi],'k-')
    
    fileName1 = [pathAnal fileN 'Field'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    if(sum(~indField) > 0)
        plotClusters(x(~indField),y(~indField),...
            idx(~indField),xl,yl,[ti ' - no field']);
        plot([0 2*pi],[0 2*pi],'k-')

        fileName1 = [pathAnal fileN 'NoField'];
        saveas(gcf,fileName1);
        print('-painters', '-dpdf', fileName1, '-r600')
    end
end

function plotClustersFieldsPhaseLabelF(pathAnal,x,y,idx,xfF,yfF,xfNoF,yfNoF,indField,xl,yl,ti,fileN)
    plotClusters(x(indField),y(indField),...
        idx(indField),xl,yl,[ti ' - field']);
%     plot([0 2*pi],[0 2*pi],'k-')
    h = plot(xfF,yfF,'k+');
    set(h,'MarkerSize',9);
    
    fileName1 = [pathAnal fileN 'Field'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    if(sum(~indField) > 0)
        plotClusters(x(~indField),y(~indField),...
            idx(~indField),xl,yl,[ti ' - no field']);
    %     plot([0 2*pi],[0 2*pi],'k-')
        h = plot(xfNoF,yfNoF,'k+');
        set(h,'MarkerSize',9);

        fileName1 = [pathAnal fileN 'NoField'];
        saveas(gcf,fileName1);
        print('-painters', '-dpdf', fileName1, '-r600')
    end
end

function plotClustersFields(pathAnal,x,y,idx,indField,xl,yl,ti,fileN)
    plotClusters(x(indField),y(indField),...
        idx(indField),xl,yl,[ti ' - field']);
    
    fileName1 = [pathAnal fileN 'Field'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotClusters(x(~indField),y(~indField),...
        idx(~indField),xl,yl,[ti ' - no field']);
        
    fileName1 = [pathAnal fileN 'NoField'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotClustersFieldsLabelF(pathAnal,x,y,idx,xfF,yfF,xfNoF,yfNoF,indField,xl,yl,ti,fileN)
    plotClusters(x(indField),y(indField),...
        idx(indField),xl,yl,[ti ' - field']);
    h = plot(xfF,yfF,'k+');
    set(h,'MarkerSize',9);
    
    fileName1 = [pathAnal fileN 'Field'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotClusters(x(~indField),y(~indField),...
        idx(~indField),xl,yl,[ti ' - no field']);
    h = plot(xfNoF,yfNoF,'k+');
    set(h,'MarkerSize',9);
    
    fileName1 = [pathAnal fileN 'NoField'];
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
