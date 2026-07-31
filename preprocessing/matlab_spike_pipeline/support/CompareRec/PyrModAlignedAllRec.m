function PyrModAlignedAllRec(onlyRun,methodKMean)

    methodTheta = 1;
    minFR = 0.15;
    maxFR = 7;
    if(nargin == 1)
        methodKMean = 2;
    end
    
    RecordingList;
    pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\PyramidalAligned\';
    pathAnal1 = 'Z:\Yingxue\DataAnalysisRaphi\Pyramidal\';
    
    load([pathAnal1 'autoCorrPyrAllRec.mat']);
    
    %% pyramidal neurons in no cue passive task
%     disp('No cue')
%     modAlignedPyrNoCue = accumPyrNeurons1(listRecordingsNoCuePath,...
%         listRecordingsNoCueFileName,mazeSessionNoCue,minFR,maxFR,1,methodTheta,onlyRun);
%     
%     disp('Active licking')
%     modAlignedPyrAL = accumPyrNeurons1(listRecordingsActiveLickPath,...
%         listRecordingsActiveLickFileName,mazeSessionActiveLick,minFR,maxFR,2,methodTheta,onlyRun);
%     
%     disp('Passive licking')
%     modAlignedPyrPL = accumPyrNeurons1(listRecordingsPassiveLickPath,...
%         listRecordingsPassiveLickFileName,mazeSessionPassiveLick,minFR,maxFR,3,methodTheta,onlyRun);
%     
%     save([pathAnal1 'autoCorrPyrAllRec.mat'],'modAlignedPyrNoCue','modAlignedPyrAL','modAlignedPyrPL','-append'); 

    mod.task = [autoCorrPyrNoCue.task autoCorrPyrAL.task autoCorrPyrPL.task];
    mod.indRec = [autoCorrPyrNoCue.indRec autoCorrPyrAL.indRec autoCorrPyrPL.indRec];
    mod.nNeuWithField = [modPyrNoCue.nNeuWithField modPyrAL.nNeuWithField modPyrPL.nNeuWithField];
    mod.isNeuWithField = [modPyrNoCue.isNeuWithField modPyrAL.isNeuWithField modPyrPL.isNeuWithField];
    mod.fieldWidth = [modPyrNoCue.fieldWidth modPyrAL.fieldWidth modPyrPL.fieldWidth];
    mod.indStartField = [modPyrNoCue.indStartField modPyrAL.indStartField modPyrPL.indStartField];
    mod.indPeakField = [modPyrNoCue.indPeakField modPyrAL.indPeakField modPyrPL.indPeakField];
    mod.percTrackStartField = [modPyrNoCue.indStartField./modPyrNoCue.trialLenMean...
        modPyrAL.indStartField./modPyrAL.trialLenMean modPyrPL.indStartField./modPyrPL.trialLenMean];
    mod.percTrackPeakField = [modPyrNoCue.indPeakField./modPyrNoCue.trialLenMean...
        modPyrAL.indPeakField./modPyrAL.trialLenMean modPyrPL.indPeakField./modPyrPL.trialLenMean];
       
    mod.mFR = [modAlignedPyrNoCue.mFR modAlignedPyrAL.mFR modAlignedPyrPL.mFR];
     
    mod.burstMeanResultantLen = [modAlignedPyrNoCue.burstMeanResultantLen modAlignedPyrAL.burstMeanResultantLen ...
                modAlignedPyrPL.burstMeanResultantLen];
    mod.burstMeanDire = [modAlignedPyrNoCue.burstMeanDire modAlignedPyrAL.burstMeanDire modAlignedPyrPL.burstMeanDire];
    mod.nonBurstMeanDire = [modAlignedPyrNoCue.nonBurstMeanDire modAlignedPyrAL.nonBurstMeanDire modAlignedPyrPL.nonBurstMeanDire];
    mod.burstMeanDireStart = [modAlignedPyrNoCue.burstMeanDireStart modAlignedPyrAL.burstMeanDireStart modAlignedPyrPL.burstMeanDireStart];
    mod.numSpPerBurstMean = [modAlignedPyrNoCue.numSpPerBurstMean modAlignedPyrAL.numSpPerBurstMean modAlignedPyrPL.numSpPerBurstMean];
    mod.fractBurst = [modAlignedPyrNoCue.fractBurst modAlignedPyrAL.fractBurst modAlignedPyrPL.fractBurst];
            
    mod.thetaModHist = [modAlignedPyrNoCue.thetaModHist modAlignedPyrAL.thetaModHist modAlignedPyrPL.thetaModHist];
    mod.thetaModHistH = [modAlignedPyrNoCue.thetaModHistH modAlignedPyrAL.thetaModHistH modAlignedPyrPL.thetaModHistH];
    mod.phaseMeanDire = [modAlignedPyrNoCue.phaseMeanDire modAlignedPyrAL.phaseMeanDire modAlignedPyrPL.phaseMeanDire];
    mod.phaseMeanDireH = [modAlignedPyrNoCue.phaseMeanDireH modAlignedPyrAL.phaseMeanDireH modAlignedPyrPL.phaseMeanDireH];
    mod.maxPhaseArr = [modAlignedPyrNoCue.maxPhaseFil modAlignedPyrAL.maxPhaseFil modAlignedPyrPL.maxPhaseFil];
    mod.maxPhaseArrH = [modAlignedPyrNoCue.maxPhaseFilH modAlignedPyrAL.maxPhaseFilH modAlignedPyrPL.maxPhaseFilH];
    mod.minPhaseArr = [modAlignedPyrNoCue.minPhaseFil modAlignedPyrAL.minPhaseFil modAlignedPyrPL.minPhaseFil];
    mod.minPhaseArrH = [modAlignedPyrNoCue.minPhaseFilH modAlignedPyrAL.minPhaseFilH modAlignedPyrPL.minPhaseFilH];
    mod.phaseMeanResultantLen = [modAlignedPyrNoCue.phaseMeanResultantLen modAlignedPyrAL.phaseMeanResultantLen modAlignedPyrPL.phaseMeanResultantLen];
        
    phaseDiff = mod.maxPhaseArr - mod.minPhaseArr;
    phaseDiff(phaseDiff < 0) = phaseDiff(phaseDiff < 0) + 360;
    mod.phaseDiff = phaseDiff;
    phaseDiffH = mod.maxPhaseArrH - mod.minPhaseArrH;
    phaseDiffH(phaseDiffH < 0) = phaseDiffH(phaseDiffH < 0) + 360;
    mod.phaseDiffH = phaseDiffH;
    
    mod.diffNeuronLFPFreq = [modAlignedPyrNoCue.thetaModFreq3-modAlignedPyrNoCue.thetaFreqHMean...
            modAlignedPyrAL.thetaModFreq3-modAlignedPyrAL.thetaFreqHMean...
            modAlignedPyrPL.thetaModFreq3-modAlignedPyrPL.thetaFreqHMean];
    mod.thetaModFreq3 = [modAlignedPyrNoCue.thetaModFreq3 modAlignedPyrAL.thetaModFreq3 modAlignedPyrPL.thetaModFreq3];
    mod.thetaModInd3 = [modAlignedPyrNoCue.thetaModInd3 modAlignedPyrAL.thetaModInd3 modAlignedPyrPL.thetaModInd3];
    mod.thetaModInd = [modAlignedPyrNoCue.thetaModInd modAlignedPyrAL.thetaModInd modAlignedPyrPL.thetaModInd];   
    mod.thetaAsym3 = [modAlignedPyrNoCue.thetaAsym3 modAlignedPyrAL.thetaAsym3 modAlignedPyrPL.thetaAsym3];
    
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

    x = mod.burstMeanDire(mod.fractBurst > 0) - mod.phaseMeanDire(mod.fractBurst > 0);
    mod.idxCBurst = mod.idxC(mod.fractBurst > 0);
    x(x < -pi) = x(x < -pi) + 2*pi;
    x(x > pi) = x(x > pi) - 2*pi;
    mod.burstThetaDiff = x;
 
    % for each cluster, compare neurons in the recordings with fields vs.
    % recordings without field
    modAlignedPyrStatsField = modPyrStatsFieldPerC(mod,autoCorrPyrAll);
    
    % compare neurons between the two clusters
    modAlignedPyrStatsC = modPyrStatsC(mod);
    
    % for each cluster, compare neurons with fields and without field
    modAlignedPyrStatsFNeuVsNoFNeu = modPyrStatsFieldNeu(mod,autoCorrPyrAll);
    
    save([pathAnal1 'autoCorrPyrAllRec_km' num2str(methodKMean) '.mat'],'modAlignedPyrStatsField','modAlignedPyrStatsFNeuVsNoFNeu','modAlignedPyrStatsC','-append'); 
    
    %% theta phase mean direction vs. 
    % diff between burst mean direction and theta mean direction
    plotBurstVsTheta(mod.burstMeanDire,mod.phaseMeanDire,mod.fractBurst,pathAnal);
   
    %% compare two clusters
    colorSel = 0;
    idxC = mod.idxC;
    plotBoxPlot(autoCorrPyrAll.relDepthNeuHDefC1{1},...
        autoCorrPyrAll.relDepthNeuHDefC1{2},'Depth','relDepthNeuHDefCmpC1C2Box',...
        pathAnal,[-6 6],autoCorrPyrAll.pRSRelDepthNeuHDefC1,colorSel);

    plotBoxPlot(mod.phaseMeanDire(idxC == 1),...
        mod.phaseMeanDire(idxC == 2),'Mean theta phase direction','Pyr_ThetaMeanC1C2Box',...
        pathAnal,[],modAlignedPyrStatsC.pWWPhaseMeanDireC,colorSel);
    
    plotBoxPlot(mod.phaseMeanResultantLen(idxC == 1),...
        mod.phaseMeanResultantLen(idxC == 2),'Theta phase resultant length','Pyr_PhaseMeanResultantLenC1C2Box',...
        pathAnal,[-0.1 0.6],modAlignedPyrStatsC.pRSPhaseMeanResultantLenC,colorSel);
    
    plotBoxPlot(mod.burstMeanResultantLen(idxC == 1),...
        mod.burstMeanResultantLen(idxC == 2),'Burst phase resultant length','Pyr_BurstMeanResultantLenC1C2Box',...
        pathAnal,[-0.1 1],modAlignedPyrStatsC.pRSBurstMeanResultantLenC,colorSel);
    
    plotBoxPlot(mod.thetaModHist(idxC == 1),...
        mod.thetaModHist(idxC == 2),'Theta modulation','Pyr_ThetaModHistC1C2Box',...
        pathAnal,[-0.1 1.2],modAlignedPyrStatsC.pRSThetaModHistC,colorSel);
    
    plotBoxPlot(mod.burstThetaDiff(mod.idxCBurst == 1),...
        mod.burstThetaDiff(mod.idxCBurst == 2),'Burst phase - theta phase','Pyr_BurstPhase-ThetaPhaseC1C2Box',...
        pathAnal,[],modAlignedPyrStatsC.pWWBurstThetaDiffC,colorSel);
    
    plotBoxPlot(mod.fractBurst(idxC == 1),...
        mod.fractBurst(idxC == 2),'Fract burst','Pyr_FractBurstC1C2Box',...
        pathAnal,[],modAlignedPyrStatsC.pRSFractBurstC,colorSel);
    
    plotBoxPlot(mod.numSpPerBurstMean(idxC == 1),...
        mod.numSpPerBurstMean(idxC == 2),'Num. spikes per burst','Pyr_NumSpPerBurstMeanC1C2Box',...
        pathAnal,[1.8 3.2],modAlignedPyrStatsC.pRSNumSpPerBurstMeanC,colorSel);
    
    plotBoxPlot(mod.thetaModFreq3(idxC == 1),...
        mod.thetaModFreq3(idxC == 2),'Theta modulation frequency (Hz)','Pyr_ThetaModFreq3C1C2Box',...
        pathAnal,[],modAlignedPyrStatsC.pRSThetaModFreq3C,colorSel);
    
    plotBoxPlot(mod.diffNeuronLFPFreq(idxC == 1),...
        mod.diffNeuronLFPFreq(idxC == 2),'Neuron mod freq. - LFP freq. (Hz)','Pyr_DiffNeuronLFPFreqC1C2Box',...
        pathAnal,[],modAlignedPyrStatsC.pRSDiffNeuronLFPFreq,colorSel);
    
    plotBoxPlot(mod.thetaAsym3(idxC == 1),...
        mod.thetaAsym3(idxC == 2),'Theta asymmetry','Pyr_ThetaAsym3C1C2Box',...
        pathAnal,[0.15 0.85],modAlignedPyrStatsC.pRSThetaAsym3C,colorSel);
        
    plotDistri(autoCorrPyrAll.relDepthNeuHDefC1{1},...
        autoCorrPyrAll.relDepthNeuHDefC1{2},...
        [min(autoCorrPyrAll.relDepthNeuHDef):max(autoCorrPyrAll.relDepthNeuHDef)],...
        'Depth',['C1' 'C2'],pathAnal,'relDepthNeuHDefCmpC1C2',[-5 5]);
    
    plotPolarPlot(mod.phaseMeanDire(idxC == 1),...
        mod.phaseMeanDire(idxC == 2),'Mean theta phase direction',...
        'Pyr_ThetaMeanC1C2Polar',pathAnal,modAlignedPyrStatsC.pWWPhaseMeanDireC);
    
    plotPolarPlot(mod.maxPhaseArr(idxC == 1)/180*pi,...
        mod.maxPhaseArr(idxC == 2)/180*pi,'Max theta phase',...
        'Pyr_MaxPhaseArrC1C2Polar',pathAnal,modAlignedPyrStatsC.pWWMaxPhaseC);
    
    plotPolarPlot(mod.minPhaseArr(idxC == 1)/180*pi,...
        mod.minPhaseArr(idxC == 2)/180*pi,'Min theta phase',...
        'Pyr_MinPhaseArrC1C2Polar',pathAnal,modAlignedPyrStatsC.pWWMinPhaseC);
    
    x = mod.burstMeanDire(mod.fractBurst > 0);
    plotPolarPlot(x(mod.idxCBurst == 1),...
        x(mod.idxCBurst == 2),'Burst phase',...
        'Pyr_BurstMeanDireC1C2Polar',pathAnal,modAlignedPyrStatsC.pWWburstMeanDireC);
    
    x = mod.burstMeanDireStart(mod.fractBurst > 0);
    plotPolarPlot(x(mod.idxCBurst == 1),...
        x(mod.idxCBurst == 2),'Burst start phase direction',...
        'Pyr_BurstMeanDireStartC1C2Polar',pathAnal,modAlignedPyrStatsC.pWWburstMeanDireStartC);
    
    %% compare field vs. no field cluster 1 (deep cells)
    colorSel = 1;   
    cluster = 1;
    plotBarPlot(autoCorrPyrAll.relDepthNeuHDef(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        autoCorrPyrAll.relDepthNeuHDef(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Depth',...
        'Pyr_RelDepthNeuHDefCmpC1FieldBar',pathAnal,[-5 6],modAlignedPyrStatsFNeuVsNoFNeu.pKWRelDepthNeuHDef(cluster),colorSel);
    
    plotBoxPlot(mod.fractBurst(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        mod.fractBurst(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Fract burst',...
        'Pyr_FractBurstC1FieldBox',pathAnal,[-0.1 0.85],modAlignedPyrStatsFNeuVsNoFNeu.pRSMeanFractBurst(cluster),colorSel);
    
    plotBoxPlot(mod.numSpPerBurstMean(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        mod.numSpPerBurstMean(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Burst length',...
        'Pyr_NumSpPerBurstMeanC1FieldBox',pathAnal,[1.9 3.5],modAlignedPyrStatsFNeuVsNoFNeu.pRSNumSpPerBurstMean(cluster),colorSel);
    
    plotBoxPlot(mod.thetaModHist(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        mod.thetaModHist(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Theta modulation',...
        'Pyr_ThetaModHistC1FieldBox',pathAnal,[0 1.1],modAlignedPyrStatsFNeuVsNoFNeu.pRSThetaModHist(cluster),colorSel);
    
    x = mod.phaseMeanDire;
    x(x > pi) = x(x > pi) - 2*pi;
    plotBoxPlot(x(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        x(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Mean theta phase',...
        'Pyr_ThetaMeanC1FieldBox',pathAnal,[-2 3],modAlignedPyrStatsFNeuVsNoFNeu.pWWPhaseMeanDire(cluster),colorSel);
    
    x = mod.burstMeanDire;
    x(x > pi) = x(x > pi) - 2*pi;
    plotBoxPlot(x(modAlignedPyrStatsFNeuVsNoFNeu.indCurCFieldBurst{cluster}),...
        x(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoFieldBurst{cluster}),'Burst phase',...
        'Pyr_BurstMeanDirC1FieldBox',pathAnal,[-3 pi],modAlignedPyrStatsFNeuVsNoFNeu.pWWBurstMeanDire(cluster),colorSel);
    
    x = mod.burstMeanDireStart;
    x(x > pi) = x(x > pi) - 2*pi;
    plotBoxPlot(x(modAlignedPyrStatsFNeuVsNoFNeu.indCurCFieldBurst{cluster}),...
        x(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoFieldBurst{cluster}),'Burst start phase',...
        'Pyr_BurstMeanDirStartC1FieldBox',pathAnal,[-3 pi],modAlignedPyrStatsFNeuVsNoFNeu.pWWBurstMeanDireStart(cluster),colorSel);
    
    plotBoxPlot(mod.thetaModFreq3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        mod.thetaModFreq3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Theta modulation frequency (Hz)',...
        'Pyr_ThetaModFreqC1FieldBox',pathAnal,[4 10.5],modAlignedPyrStatsFNeuVsNoFNeu.pRSThetaModFreq3(cluster),colorSel);
    
    plotBoxPlot(mod.diffNeuronLFPFreq(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        mod.diffNeuronLFPFreq(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Neuron mod freq. - LFP freq. (Hz)',...
        'Pyr_DiffNeuronLFPFreqC1C2Box',...
        pathAnal,[],modAlignedPyrStatsFNeuVsNoFNeu.pRSDiffNeuronLFPFreq(cluster),colorSel);
    
    plotBoxPlot(mod.thetaAsym3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        mod.thetaAsym3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Theta asymmetry',...
        'Pyr_ThetaAsymC1FieldBox',pathAnal,[0.1 0.85],modAlignedPyrStatsFNeuVsNoFNeu.pRSThetaAsym3(cluster),colorSel);
    
    
    %% compare field vs. no field cluster 2
    colorSel = 2;   
    cluster = 2;
    plotBarPlot(autoCorrPyrAll.relDepthNeuHDef(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        autoCorrPyrAll.relDepthNeuHDef(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Depth',...
        'Pyr_RelDepthNeuHDefCmpC2FieldBar',pathAnal,[-5 6],modAlignedPyrStatsFNeuVsNoFNeu.pRSRelDepthNeuHDef(cluster),colorSel);
    
    plotBoxPlot(mod.fractBurst(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        mod.fractBurst(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Fract burst',...
        'Pyr_FractBurstC2FieldBox',pathAnal,[-0.1 0.85],modAlignedPyrStatsFNeuVsNoFNeu.pRSMeanFractBurst(cluster),colorSel);
    
    plotBoxPlot(mod.numSpPerBurstMean(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        mod.numSpPerBurstMean(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Burst length',...
        'Pyr_NumSpPerBurstMeanC2FieldBox',pathAnal,[1.9 3.5],modAlignedPyrStatsFNeuVsNoFNeu.pRSNumSpPerBurstMean(cluster),colorSel);
    
    plotBoxPlot(mod.thetaModHist(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        mod.thetaModHist(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Theta modulation',...
        'Pyr_ThetaModHistC2FieldBox',pathAnal,[0 1.1],modAlignedPyrStatsFNeuVsNoFNeu.pRSThetaModHist(cluster),colorSel);
    
    plotBoxPlot(mod.phaseMeanDire(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        mod.phaseMeanDire(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Mean theta phase',...
        'Pyr_ThetaMeanC2FieldBox',pathAnal,[-0.5 6],modAlignedPyrStatsFNeuVsNoFNeu.pWWPhaseMeanDire(cluster),colorSel);
    
    plotBoxPlot(mod.burstMeanDire(modAlignedPyrStatsFNeuVsNoFNeu.indCurCFieldBurst{cluster}),...
        mod.burstMeanDire(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoFieldBurst{cluster}),'Burst phase',...
        'Pyr_BurstMeanDirC2FieldBox',pathAnal,[],modAlignedPyrStatsFNeuVsNoFNeu.pWWBurstMeanDire(cluster),colorSel);
    
    plotBoxPlot(mod.burstMeanDireStart(modAlignedPyrStatsFNeuVsNoFNeu.indCurCFieldBurst{cluster}),...
        mod.burstMeanDireStart(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoFieldBurst{cluster}),'Burst start phase',...
        'Pyr_BurstMeanDirStartC2FieldBox',pathAnal,[],modAlignedPyrStatsFNeuVsNoFNeu.pWWBurstMeanDireStart(cluster),colorSel);
    
    plotBoxPlot(mod.thetaModFreq3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        mod.thetaModFreq3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Theta modulation frequency (Hz)',...
        'Pyr_ThetaModFreqC2FieldBox',pathAnal,[4.5 10],modAlignedPyrStatsFNeuVsNoFNeu.pRSThetaModFreq3(cluster),colorSel);
    
    plotBoxPlot(mod.diffNeuronLFPFreq(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        mod.diffNeuronLFPFreq(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Neuron mod freq. - LFP freq. (Hz)',...
        'Pyr_DiffNeuronLFPFreqC1C2Box',...
        pathAnal,[],modAlignedPyrStatsFNeuVsNoFNeu.pRSDiffNeuronLFPFreq(cluster),colorSel);
    
    plotBoxPlot(mod.thetaAsym3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCField{cluster}),...
        mod.thetaAsym3(modAlignedPyrStatsFNeuVsNoFNeu.indCurCNoField{cluster}),'Theta asymmetry',...
        'Pyr_ThetaAsymC2FieldBox',pathAnal,[0.1 0.85],modAlignedPyrStatsFNeuVsNoFNeu.pRSThetaAsym3(cluster),colorSel);
    
    %% plot clusters with fields
    indCurCField1 = mod.nNeuWithField > 1 & mod.isNeuWithField == 1;
    indCurCNoField1 = mod.nNeuWithField <= 1 & mod.isNeuWithField == 1;
    indField = mod.nNeuWithField > 1;

    indCurCField2 = mod.isNeuWithField == 1;
    plotClustersFieldsPhaseLabelF(pathAnal,mod.phaseMeanDire,mod.burstMeanDire,...
        idxC,mod.phaseMeanDire(indCurCField2),mod.burstMeanDire(indCurCField2),...
        [],[],1:length(indCurCField2),...
        'Mean theta phase','Burst mean phase','All tasks label Fields',...
        'Pyr_BurstMeanVsThetaMeanAllVsField');
    plot([0 2*pi],[0 2*pi],'k-');
    
    plotClustersFieldsPhaseLabelF(pathAnal,mod.phaseMeanDire,mod.thetaModFreq3,...
        idxC,mod.phaseMeanDire(indCurCField2),mod.thetaModFreq3(indCurCField2),...
        [],[],1:length(indCurCField2),...
        'Mean theta phase','Theta modulation freq. (Hz)','All tasks label Fields',...
        'Pyr_ThetaMeanVsModFreqAllVsField');

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
        load(fullPath,'autoCorr','beh','fieldStructSess'); 
        
        fullPath = [paths(i,:) filenames(i,:) '_alignRun_msess' num2str(mazeSess(i)) '.mat']; 
        if(exist(fullPath) == 0)
            disp('The aligned to run file does not exist');
            return;
        end
        load(fullPath,'trialsRun');
                
        fullPathFR = [filenames(i,:) '_FRAlignedRun_msess' num2str(mazeSess(i)) '_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fullPathFR];
        if(exist(fullPath) == 0)
            disp('_FRAlignedRun.mat file does not exist.');
            return;
        end
        load(fullPath,'mFRStructNonStimGood'); 
        mFR = mFRStructNonStimGood;
        
        fileNamePeakFR = [filenames(i,:) '_PeakFRAligned_msess' num2str(mazeSess(i)) ...
                        '_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNamePeakFR];
        if(exist(fullPath) == 0)
            disp(['The peak firing rate file does not exist. Please call ',...
                    'function "PeakFiringRate_Aligned" first.']);
            return;
        end
        load(fullPath,'pFRNonStimGoodStruct','trialNoNonStimGood');
        
        fileNameThetaMod = [filenames(i,:) '_ThetaModAlignedRun_msess' num2str(mazeSess(i)) ...
                        '_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameThetaMod];
        if(exist(fullPath) == 0)
            disp('_ThetaModAlignedRun file does not exist.');
            return;
        end
        load(fullPath,'thetaModNonStimGood');
        thetaModSessTmp = thetaModNonStimGood;
                
        if(methodTheta == 0)
            th = 'H';
        else
            th = 'L';
        end
        fileNameBurst = [filenames(i,:) '_burstAllAlignedRun_TH' th '_msess' num2str(mazeSess(i)) ...
                     '_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameBurst];
        if(exist(fullPath) == 0)
            disp('_burstAllAlignedRun file does not exist.');
            return;
        end
        load(fullPath,'burstIsiPerNeuronNonStimGood');
        burstIsi = burstIsiPerNeuronNonStimGood;
        
        fileNameThetaPhase = [filenames(i,:) '_ThetaPhaseLAligned_msess' num2str(mazeSess(i)) ...
                    '_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameThetaPhase];
        if(exist(fullPath) == 0)
            disp('_ThetaPhaseLAligned file does not exist.');
            return;
        end
        load(fullPath,'spikeThetaPhaseRunNoStimGood');
        spikeThetaPhase = spikeThetaPhaseRunNoStimGood;
                
        fileNameThetaPhase = [filenames(i,:) '_ThetaPhaseHAligned_msess' num2str(mazeSess(i)) ...
                    '_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameThetaPhase];
        if(exist(fullPath) == 0)
            disp('_ThetaPhaseHAligned file does not exist.');
            return;
        end
        load(fullPath,'spikeThetaPhaseRunNoStimGood');
        spikeThetaPhaseH = spikeThetaPhaseRunNoStimGood;
        
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
        
        trialLenMean = mean(trialsRun.numSamples(trialNoNonStimGood));
        
%         indNeu = cluList.firingRate > minFR & cluList.firingRate < maxFR &...
%                     autoCorr.isPyrneuron == 1;
        indNeu = mFR.mFR > minFR & mFR.mFR < maxFR &...
                    autoCorr.isPyrneuron == 1;
        modPyr.task = [modPyr.task task*ones(1,sum(indNeu))];
        modPyr.indRec = [modPyr.indRec i*ones(1,sum(indNeu))];
        modPyr.indNeu = [modPyr.indNeu find(indNeu == 1)]; 
        modPyr.trialLenMean = [modPyr.trialLenMean trialLenMean*ones(1,sum(indNeu))];
        modPyr.thetaFreqHMean = [modPyr.thetaFreqHMean mean(beh.thetaFreqHMean(trialNoNonStimGood))*ones(1,sum(indNeu))];
                
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
        modPyr.mFR = [modPyr.mFR mFR.mFR(indNeu)]; % ??
        
        nNeurons = length(cluList.firingRate);
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

function modAlignedPyrStatsField = modPyrStatsFieldPerC(mod,autoCorrPyrAll)
    modAlignedPyrStatsField = [];
    for i = 1:max(mod.idxC)
        indCurCField = mod.idxC == i & mod.nNeuWithField > 1;
        indCurCNoField = mod.idxC == i & mod.nNeuWithField < 1;
        modAlignedPyrStatsField.indCurCField{i} = indCurCField;
        modAlignedPyrStatsField.indCurCNoField{i} = indCurCNoField;
        modAlignedPyrStatsField.isNeuWithField(i,:) = [sum(mod.isNeuWithField(indCurCField)),sum(mod.isNeuWithField(indCurCNoField))];
              
        modAlignedPyrStatsField.meanPeakTo40ms(i,:) = [mean(autoCorrPyrAll.peakTo40ms(indCurCField)),mean(autoCorrPyrAll.peakTo40ms(indCurCNoField))];
        modAlignedPyrStatsField.meanPeakTime(i,:) = [mean(autoCorrPyrAll.peakTime(indCurCField)),mean(autoCorrPyrAll.peakTime(indCurCNoField))];
        
        modAlignedPyrStatsField.mFR(i,:) = [mean(mod.mFR(indCurCField)),mean(mod.mFR(indCurCNoField))];
        
        modAlignedPyrStatsField.meanDiffNeuronLFPFreq(i,:) = [mean(mod.diffNeuronLFPFreq(indCurCField)),mean(mod.diffNeuronLFPFreq(indCurCNoField))];
        
        indCurCFieldBurst = mod.idxC == i & mod.nNeuWithField > 1 & mod.fractBurst > 0;
        indCurCNoFieldBurst = mod.idxC == i & mod.nNeuWithField <= 1 & mod.fractBurst > 0;
        modAlignedPyrStatsField.meanBurstMeanDire(i,:) = [circ_mean(mod.burstMeanDire(indCurCFieldBurst)'),circ_mean(mod.burstMeanDire(indCurCNoFieldBurst)')];
        modAlignedPyrStatsField.meanNonBurstMeanDire(i,:) = [circ_mean(mod.nonBurstMeanDire(indCurCField)'),circ_mean(mod.nonBurstMeanDire(indCurCNoField)')];     
        modAlignedPyrStatsField.meanBurstMeanDireStart(i,:) = [circ_mean(mod.burstMeanDireStart(indCurCFieldBurst)'),circ_mean(mod.burstMeanDireStart(indCurCNoFieldBurst)')];
        modAlignedPyrStatsField.meanFractBurst(i,:) = [mean(mod.fractBurst(indCurCField)),mean(mod.fractBurst(indCurCNoField))];
        nonZeroF = mod.numSpPerBurstMean(indCurCField);
        nonZeroF = nonZeroF(nonZeroF > 0);
        nonZeroNoF = mod.numSpPerBurstMean(indCurCNoField);
        nonZeroNoF = nonZeroNoF(nonZeroNoF > 0);
        modAlignedPyrStatsField.meanNumSpPerBurstMean(i,:) = [mean(nonZeroF),mean(nonZeroNoF)];
        
        modAlignedPyrStatsField.meanPhaseMeanDire(i,:) = [circ_mean(mod.phaseMeanDire(indCurCField)'),circ_mean(mod.phaseMeanDire(indCurCNoField)')];
        modAlignedPyrStatsField.meanPhaseMeanDireH(i,:) = [circ_mean(mod.phaseMeanDireH(indCurCField)'),circ_mean(mod.phaseMeanDireH(indCurCNoField)')];
        modAlignedPyrStatsField.meanMaxPhase(i,:) = [circ_mean(mod.maxPhaseArr(indCurCField)'/180*pi),circ_mean(mod.maxPhaseArr(indCurCNoField)'/180*pi)];
        modAlignedPyrStatsField.meanMaxPhaseH(i,:) = [circ_mean(mod.maxPhaseArrH(indCurCField)'/180*pi),circ_mean(mod.maxPhaseArrH(indCurCNoField)'/180*pi)];
        modAlignedPyrStatsField.meanminPhase(i,:) = [circ_mean(mod.minPhaseArr(indCurCField)'/180*pi),circ_mean(mod.minPhaseArr(indCurCNoField)'/180*pi)];
        modAlignedPyrStatsField.meanminPhaseH(i,:) = [circ_mean(mod.minPhaseArrH(indCurCField)'/180*pi),circ_mean(mod.minPhaseArrH(indCurCNoField)'/180*pi)];
        modAlignedPyrStatsField.meanPhaseDiff(i,:) = [mean(mod.phaseDiff(indCurCField)),mean(mod.phaseDiff(indCurCNoField))];
        modAlignedPyrStatsField.meanPhaseDiffH(i,:) = [mean(mod.phaseDiffH(indCurCField)),mean(mod.phaseDiffH(indCurCNoField))];
        modAlignedPyrStatsField.meanThetaModHist(i,:) = [mean(mod.thetaModHist(indCurCField)),mean(mod.thetaModHist(indCurCNoField))];
        modAlignedPyrStatsField.meanThetaModHistH(i,:) = [mean(mod.thetaModHistH(indCurCField)),mean(mod.thetaModHistH(indCurCNoField))];
        
        modAlignedPyrStatsField.meanThetaModFreq3(i,:) = [mean(mod.thetaModFreq3(indCurCField)),mean(mod.thetaModFreq3(indCurCNoField))];
        modAlignedPyrStatsField.meanThetaAsym3(i,:) = [mean(mod.thetaAsym3(indCurCField)),mean(mod.thetaAsym3(indCurCNoField))];
        modAlignedPyrStatsField.meanThetaModInd3(i,:) = [mean(mod.thetaModInd3(indCurCField)),mean(mod.thetaModInd3(indCurCNoField))];
        modAlignedPyrStatsField.meanThetaModInd(i,:) = [mean(mod.thetaModInd(indCurCField)),mean(mod.thetaModInd(indCurCNoField))];
        
        modAlignedPyrStatsField.pRSDiffNeuronLFPFreq(i) = ranksum(mod.diffNeuronLFPFreq(indCurCField),mod.diffNeuronLFPFreq(indCurCNoField));     
        
        modAlignedPyrStatsField.pRSPeakTo40ms(i) = ranksum(autoCorrPyrAll.peakTo40ms(indCurCField),autoCorrPyrAll.peakTo40ms(indCurCNoField));
        modAlignedPyrStatsField.pRSPeakTime(i) = ranksum(autoCorrPyrAll.peakTime(indCurCField),autoCorrPyrAll.peakTime(indCurCNoField));
        
        modAlignedPyrStatsField.pRSMFR(i) = ranksum(mod.mFR(indCurCField),mod.mFR(indCurCNoField));
        
        modAlignedPyrStatsField.pKBurstMeanDire(i) = circ_ktest(mod.burstMeanDire(indCurCFieldBurst)',mod.burstMeanDire(indCurCNoFieldBurst)');
        modAlignedPyrStatsField.pKNonBurstMeanDire(i) = circ_ktest(mod.nonBurstMeanDire(indCurCField)',mod.nonBurstMeanDire(indCurCNoField)');
        modAlignedPyrStatsField.pKBurstMeanDireStart(i) = circ_ktest(mod.burstMeanDireStart(indCurCFieldBurst)',mod.burstMeanDireStart(indCurCNoFieldBurst)');
        
        modAlignedPyrStatsField.pWWBurstMeanDire(i) = circ_wwtest(mod.burstMeanDire(indCurCFieldBurst)',mod.burstMeanDire(indCurCNoFieldBurst)');
        modAlignedPyrStatsField.pWWNonBurstMeanDire(i) = circ_wwtest(mod.nonBurstMeanDire(indCurCField)',mod.nonBurstMeanDire(indCurCNoField)');
        modAlignedPyrStatsField.pWWBurstMeanDireStart(i) = circ_wwtest(mod.burstMeanDireStart(indCurCFieldBurst)',mod.burstMeanDireStart(indCurCNoFieldBurst)');
        
        modAlignedPyrStatsField.pRSFractBurst(i) = ranksum(mod.fractBurst(indCurCField),mod.fractBurst(indCurCNoField));        
        modAlignedPyrStatsField.pRSNumSpPerBurstMean(i) = ranksum(nonZeroF,nonZeroNoF);
        
        modAlignedPyrStatsField.pKPhaseMeanDire(i) = circ_ktest(mod.phaseMeanDire(indCurCField)',mod.phaseMeanDire(indCurCNoField)');
        modAlignedPyrStatsField.pKPhaseMeanDireH(i) = circ_ktest(mod.phaseMeanDireH(indCurCField)',mod.phaseMeanDireH(indCurCNoField)');
        modAlignedPyrStatsField.pKMaxPhase(i) = circ_ktest(mod.maxPhaseArr(indCurCField)'/180*pi,mod.maxPhaseArr(indCurCNoField)'/180*pi);
        modAlignedPyrStatsField.pKMaxPhaseH(i) = circ_ktest(mod.maxPhaseArrH(indCurCField)'/180*pi,mod.maxPhaseArrH(indCurCNoField)'/180*pi);
        modAlignedPyrStatsField.pKMinPhase(i) = circ_ktest(mod.minPhaseArr(indCurCField)'/180*pi,mod.minPhaseArr(indCurCNoField)'/180*pi);
        modAlignedPyrStatsField.pKMinPhaseH(i) = circ_ktest(mod.minPhaseArrH(indCurCField)'/180*pi,mod.minPhaseArrH(indCurCNoField)'/180*pi);
        
        modAlignedPyrStatsField.pWWPhaseMeanDire(i) = circ_wwtest(mod.phaseMeanDire(indCurCField)',mod.phaseMeanDire(indCurCNoField)');
        modAlignedPyrStatsField.pWWPhaseMeanDireH(i) = circ_wwtest(mod.phaseMeanDireH(indCurCField)',mod.phaseMeanDireH(indCurCNoField)');
        modAlignedPyrStatsField.pWWMaxPhase(i) = circ_wwtest(mod.maxPhaseArr(indCurCField)'/180*pi,mod.maxPhaseArr(indCurCNoField)'/180*pi);
        modAlignedPyrStatsField.pWWMaxPhaseH(i) = circ_wwtest(mod.maxPhaseArrH(indCurCField)'/180*pi,mod.maxPhaseArrH(indCurCNoField)'/180*pi);
        modAlignedPyrStatsField.pWWMinPhase(i) = circ_wwtest(mod.minPhaseArr(indCurCField)'/180*pi,mod.minPhaseArr(indCurCNoField)'/180*pi);
        modAlignedPyrStatsField.pWWMinPhaseH(i) = circ_wwtest(mod.minPhaseArrH(indCurCField)'/180*pi,mod.minPhaseArrH(indCurCNoField)'/180*pi);
        
        modAlignedPyrStatsField.pRSPhaseDiff(i) = ranksum(mod.phaseDiff(indCurCField),mod.phaseDiff(indCurCNoField));
        modAlignedPyrStatsField.pRSPhaseDiffH(i) = ranksum(mod.phaseDiffH(indCurCField),mod.phaseDiffH(indCurCNoField));
        modAlignedPyrStatsField.pRSThetaModHist(i) = ranksum(mod.thetaModHist(indCurCField),mod.thetaModHist(indCurCNoField));
        modAlignedPyrStatsField.pRSThetaModHistH(i) = ranksum(mod.thetaModHistH(indCurCField),mod.thetaModHistH(indCurCNoField));
        
        modAlignedPyrStatsField.pRSThetaModFreq3(i) = ranksum(mod.thetaModFreq3(indCurCField),mod.thetaModFreq3(indCurCNoField));
        modAlignedPyrStatsField.pRSThetaAsym3(i) = ranksum(mod.thetaAsym3(indCurCField),mod.thetaAsym3(indCurCNoField));  
        modAlignedPyrStatsField.pRSThetaModInd3(i) = ranksum(mod.thetaModInd3(indCurCField),mod.thetaModInd3(indCurCNoField));
        modAlignedPyrStatsField.pRSThetaModInd(i) = ranksum(mod.thetaModInd(indCurCField),mod.thetaModInd(indCurCNoField));
        
        modAlignedPyrStatsField.pKBurstVsNonBurstMeanDire(i) = circ_ktest(mod.burstMeanDire(indCurCFieldBurst)',mod.nonBurstMeanDire(indCurCField)');
        modAlignedPyrStatsField.pKBurstVsThetaMeanDire(i) = circ_ktest(mod.burstMeanDire(indCurCFieldBurst)',mod.phaseMeanDire(indCurCField)');
        modAlignedPyrStatsField.pKNonBurstVsThetaMeanDire(i) = circ_ktest(mod.nonBurstMeanDire(indCurCField)',mod.phaseMeanDire(indCurCField)');
        modAlignedPyrStatsField.pKBurstStartVsThetaMeanDire(i) = circ_ktest(mod.burstMeanDireStart(indCurCFieldBurst)',mod.phaseMeanDire(indCurCField)');
        
        modAlignedPyrStatsField.pWWBurstVsNonBurstMeanDire(i) = circ_wwtest(mod.burstMeanDire(indCurCFieldBurst)',mod.nonBurstMeanDire(indCurCField)');
        modAlignedPyrStatsField.pWWBurstVsThetaMeanDire(i) = circ_wwtest(mod.burstMeanDire(indCurCFieldBurst)',mod.phaseMeanDire(indCurCField)');
        modAlignedPyrStatsField.pWWNonBurstVsThetaMeanDire(i) = circ_wwtest(mod.nonBurstMeanDire(indCurCField)',mod.phaseMeanDire(indCurCField)');
        modAlignedPyrStatsField.pWWBurstStartVsThetaMeanDire(i) = circ_wwtest(mod.burstMeanDireStart(indCurCFieldBurst)',mod.phaseMeanDire(indCurCField)');
    end
end

function modAlignedPyrStatsC = modPyrStatsC(mod)
    idxC = mod.idxC;
    modAlignedPyrStatsC.pRSDiffNeuronLFPFreq = ranksum(mod.diffNeuronLFPFreq(idxC == 1),mod.diffNeuronLFPFreq(idxC == 2));
    nonZeroF1 = mod.numSpPerBurstMean(idxC == 1);
    nonZeroF1 = nonZeroF1(nonZeroF1 > 0);
    nonZeroF2 = mod.numSpPerBurstMean(idxC == 2);
    nonZeroF2 = nonZeroF2(nonZeroF2 > 0);
    modAlignedPyrStatsC.pRSNumSpPerBurstMeanC = ranksum(nonZeroF1,nonZeroF2);
    modAlignedPyrStatsC.pRSFractBurstC = ranksum(mod.fractBurst(idxC == 1),mod.fractBurst(idxC == 2));    
    
    modAlignedPyrStatsC.pWWBurstThetaDiffC = circ_wwtest(mod.burstThetaDiff(mod.idxCBurst == 1),...
        mod.burstThetaDiff(mod.idxCBurst == 2));
    modAlignedPyrStatsC.pKBurstThetaDiffC = circ_ktest(mod.burstThetaDiff(mod.idxCBurst == 1),...
        mod.burstThetaDiff(mod.idxCBurst == 2));
    
    modAlignedPyrStatsC.pWWPhaseMeanDireC = circ_wwtest(mod.phaseMeanDire(idxC == 1)',mod.phaseMeanDire(idxC == 2)');
    modAlignedPyrStatsC.pKPhaseMeanDireC = circ_ktest(mod.phaseMeanDire(idxC == 1)',mod.phaseMeanDire(idxC == 2)');
    modAlignedPyrStatsC.pWWPhaseMeanDireHC = circ_wwtest(mod.phaseMeanDireH(idxC == 1)',mod.phaseMeanDireH(idxC == 2)');
    modAlignedPyrStatsC.pKPhaseMeanDireHC = circ_ktest(mod.phaseMeanDireH(idxC == 1)',mod.phaseMeanDireH(idxC == 2)');
    modAlignedPyrStatsC.pWWMaxPhaseC = circ_wwtest(mod.maxPhaseArr(idxC == 1)'/180*pi,mod.maxPhaseArr(idxC == 2)'/180*pi);
    modAlignedPyrStatsC.pKMaxPhaseC = circ_ktest(mod.maxPhaseArr(idxC == 1)'/180*pi,mod.maxPhaseArr(idxC == 2)'/180*pi);
    modAlignedPyrStatsC.pWWMinPhaseC = circ_wwtest(mod.minPhaseArr(idxC == 1)'/180*pi,mod.minPhaseArr(idxC == 2)'/180*pi);
    modAlignedPyrStatsC.pKMinPhaseC = circ_ktest(mod.minPhaseArr(idxC == 1)'/180*pi,mod.minPhaseArr(idxC == 2)'/180*pi);
    modAlignedPyrStatsC.pWWburstMeanDireC = circ_wwtest(mod.burstMeanDire(idxC == 1 & mod.fractBurst > 0)',...
        mod.burstMeanDire(idxC == 2 & mod.fractBurst > 0)');
    modAlignedPyrStatsC.pKburstMeanDireC = circ_ktest(mod.burstMeanDire(idxC == 1 & mod.fractBurst > 0)',...
        mod.burstMeanDire(idxC == 2 & mod.fractBurst > 0)');
    modAlignedPyrStatsC.pWWburstMeanDireStartC = circ_wwtest(mod.burstMeanDireStart(idxC == 1 & mod.fractBurst > 0)',...
        mod.burstMeanDireStart(idxC == 2 & mod.fractBurst > 0)');
    modAlignedPyrStatsC.pKburstMeanDirStartC = circ_ktest(mod.burstMeanDireStart(idxC == 1 & mod.fractBurst > 0)',...
        mod.burstMeanDireStart(idxC == 2 & mod.fractBurst > 0)');
    modAlignedPyrStatsC.pRSBurstMeanResultantLenC = ranksum(mod.burstMeanResultantLen(idxC == 1),mod.burstMeanResultantLen(idxC == 2));
    modAlignedPyrStatsC.pRSPhaseMeanResultantLenC = ranksum(mod.phaseMeanResultantLen(idxC == 1),mod.phaseMeanResultantLen(idxC == 2));
    
    modAlignedPyrStatsC.pRSThetaModHistC = ranksum(mod.thetaModHist(idxC == 1),mod.thetaModHist(idxC == 2));
    modAlignedPyrStatsC.pRSThetaModHistHC = ranksum(mod.thetaModHistH(idxC == 1),mod.thetaModHistH(idxC == 2));
    modAlignedPyrStatsC.pRSThetaModFreq3C = ranksum(mod.thetaModFreq3(idxC == 1),mod.thetaModFreq3(idxC == 2));
    modAlignedPyrStatsC.pRSThetaAsym3C = ranksum(mod.thetaAsym3(idxC == 1),mod.thetaAsym3(idxC == 2));
    
    modAlignedPyrStatsC.pRSMFRC = ranksum(mod.mFR(idxC == 1),mod.mFR(idxC == 2));
end

function modAlignedPyrStatsFieldF = modPyrStatsFieldNeu(mod,autoCorrPyrAll)
    %% field vs no field neurons for each cluster
    modAlignedPyrStatsFieldF = [];
    for i = 1:max(mod.idxC)        
        indCurCField = mod.idxC == i & mod.isNeuWithField == 1;
        indCurCNoField = mod.idxC == i & mod.isNeuWithField == 0;
        indCurCFieldBurst = mod.idxC == i & mod.isNeuWithField == 1 & mod.fractBurst > 0;
        indCurCNoFieldBurst = mod.idxC == i & mod.isNeuWithField == 0 & mod.fractBurst > 0;
        modAlignedPyrStatsFieldF.indCurCField{i} = indCurCField;
        modAlignedPyrStatsFieldF.indCurCNoField{i} = indCurCNoField;
        modAlignedPyrStatsFieldF.indCurCFieldBurst{i} = indCurCFieldBurst;
        modAlignedPyrStatsFieldF.indCurCNoFieldBurst{i} = indCurCNoFieldBurst;
        modAlignedPyrStatsFieldF.fieldWidthCFRec{i} = mod.fieldWidth(indCurCField);
        modAlignedPyrStatsFieldF.indPeakFieldFRec{i} = mod.indPeakField(indCurCField);
        modAlignedPyrStatsFieldF.indStartFieldFRec{i} = mod.indStartField(indCurCField);
        modAlignedPyrStatsFieldF.percTrackStartField{i} = mod.percTrackStartField(indCurCField);
        modAlignedPyrStatsFieldF.percTrackPeakField{i} = mod.percTrackPeakField(indCurCField);
        modAlignedPyrStatsFieldF.skewness{i} = (mod.indPeakField(indCurCField) - mod.indStartField(indCurCField) + 1)./mod.fieldWidth(indCurCField);
      
        modAlignedPyrStatsFieldF.meanDiffNeuronLFPFreq(i,:) = [mean(mod.diffNeuronLFPFreq(indCurCField)),mean(mod.diffNeuronLFPFreq(indCurCNoField))];
        modAlignedPyrStatsFieldF.meanRelDepthNeuHDef(i,:) = [mean(autoCorrPyrAll.relDepthNeuHDef(indCurCField)),...
                            mean(autoCorrPyrAll.relDepthNeuHDef(indCurCNoField))];
          
        modAlignedPyrStatsFieldF.meanPhaseMeanDire(i,:) = [circ_mean(mod.phaseMeanDire(indCurCField)'),circ_mean(mod.phaseMeanDire(indCurCNoField)')];
        modAlignedPyrStatsFieldF.meanBurstMeanDire(i,:) = [circ_mean(mod.burstMeanDire(indCurCFieldBurst)'),circ_mean(mod.burstMeanDire(indCurCNoFieldBurst)')];  
        modAlignedPyrStatsFieldF.meanNonBurstMeanDire(i,:) = [circ_mean(mod.nonBurstMeanDire(indCurCFieldBurst)'),circ_mean(mod.nonBurstMeanDire(indCurCNoFieldBurst)')]; 
        modAlignedPyrStatsFieldF.meanBurstMeanDireStart(i,:) = [circ_mean(mod.burstMeanDireStart(indCurCFieldBurst)'),circ_mean(mod.burstMeanDireStart(indCurCNoFieldBurst)')];  
        modAlignedPyrStatsFieldF.meanFractBurst(i,:) = [mean(mod.fractBurst(indCurCField)),mean(mod.fractBurst(indCurCNoField))];
        modAlignedPyrStatsFieldF.meanNumSpPerBurstMean(i,:) = [mean(mod.numSpPerBurstMean(indCurCField)),mean(mod.numSpPerBurstMean(indCurCNoField))];
        modAlignedPyrStatsFieldF.meanBurstMeanResultantLen(i,:) = [mean(mod.burstMeanResultantLen(indCurCFieldBurst)),mean(mod.burstMeanResultantLen(indCurCNoFieldBurst))];
        modAlignedPyrStatsFieldF.meanPhaseMeanResultantLen(i,:) = [mean(mod.phaseMeanResultantLen(indCurCField)),mean(mod.phaseMeanResultantLen(indCurCNoField))];
            
        modAlignedPyrStatsFieldF.meanThetaModHist(i,:) = [mean(mod.thetaModHist(indCurCField)),mean(mod.thetaModHist(indCurCNoField))];
        modAlignedPyrStatsFieldF.meanThetaModHistH(i,:) = [mean(mod.thetaModHistH(indCurCField)),mean(mod.thetaModHistH(indCurCNoField))];
        modAlignedPyrStatsFieldF.meanThetaModFreq3(i,:) = [mean(mod.thetaModFreq3(indCurCField)),mean(mod.thetaModFreq3(indCurCNoField))];
        modAlignedPyrStatsFieldF.meanThetaAsym3(i,:) = [mean(mod.thetaAsym3(indCurCField)),mean(mod.thetaAsym3(indCurCNoField))];
        
        modAlignedPyrStatsFieldF.pRSDiffNeuronLFPFreq(i,:) = ranksum(mod.diffNeuronLFPFreq(indCurCField),mod.diffNeuronLFPFreq(indCurCNoField));
        modAlignedPyrStatsFieldF.pRSRelDepthNeuHDef(i) = ranksum(autoCorrPyrAll.relDepthNeuHDef(indCurCField),...
                            autoCorrPyrAll.relDepthNeuHDef(indCurCNoField));
        modAlignedPyrStatsFieldF.pKWRelDepthNeuHDef(i) = kruskalwallis([autoCorrPyrAll.relDepthNeuHDef(indCurCField),...
                            autoCorrPyrAll.relDepthNeuHDef(indCurCNoField)],...
                            [ones(1,sum(indCurCField)),2*ones(1,sum(indCurCNoField))]);
        
        modAlignedPyrStatsFieldF.pRSThetaModHist(i) = ranksum(mod.thetaModHist(indCurCField),mod.thetaModHist(indCurCNoField));
        modAlignedPyrStatsFieldF.pRSThetaModHistH(i) = ranksum(mod.thetaModHistH(indCurCField),mod.thetaModHistH(indCurCNoField));
        modAlignedPyrStatsFieldF.pRSThetaModFreq3(i) = ranksum(mod.thetaModFreq3(indCurCField),mod.thetaModFreq3(indCurCNoField));
        modAlignedPyrStatsFieldF.pRSThetaAsym3(i) = ranksum(mod.thetaAsym3(indCurCField),mod.thetaAsym3(indCurCNoField));
    
        modAlignedPyrStatsFieldF.pKPhaseMeanDire(i) = circ_ktest(mod.phaseMeanDire(indCurCField)',mod.phaseMeanDire(indCurCNoField)');
        modAlignedPyrStatsFieldF.pKBurstMeanDire(i) = circ_ktest(mod.burstMeanDire(indCurCFieldBurst)',mod.burstMeanDire(indCurCNoFieldBurst)'); 
        modAlignedPyrStatsFieldF.pKNonBurstMeanDire(i) = circ_ktest(mod.nonBurstMeanDire(indCurCField)',mod.nonBurstMeanDire(indCurCNoField)'); 
        modAlignedPyrStatsFieldF.pKBurstMeanDireStart(i) = circ_ktest(mod.burstMeanDireStart(indCurCFieldBurst)',mod.burstMeanDireStart(indCurCNoFieldBurst)'); 
        modAlignedPyrStatsFieldF.pWWPhaseMeanDire(i) = circ_wwtest(mod.phaseMeanDire(indCurCField)',mod.phaseMeanDire(indCurCNoField)');
        modAlignedPyrStatsFieldF.pWWBurstMeanDire(i) = circ_wwtest(mod.burstMeanDire(indCurCFieldBurst)',mod.burstMeanDire(indCurCNoFieldBurst)'); 
        modAlignedPyrStatsFieldF.pWWNonBurstMeanDire(i) = circ_wwtest(mod.nonBurstMeanDire(indCurCField)',mod.nonBurstMeanDire(indCurCNoField)');
        modAlignedPyrStatsFieldF.pWWBurstMeanDireStart(i) = circ_wwtest(mod.burstMeanDireStart(indCurCFieldBurst)',mod.burstMeanDireStart(indCurCNoFieldBurst)'); 
        modAlignedPyrStatsFieldF.pMEDBurstVsThetaMeanDire(i) = circ_medtest(mod.burstMeanDire(indCurCFieldBurst)'- mod.phaseMeanDire(indCurCFieldBurst)',0);  
        modAlignedPyrStatsFieldF.pMEDNonBurstVsThetaMeanDire(i) = circ_medtest(mod.nonBurstMeanDire(indCurCField)'- mod.phaseMeanDire(indCurCField)',0);    
        modAlignedPyrStatsFieldF.pRSBurstMeanResultantLen(i) = ranksum(mod.burstMeanResultantLen(indCurCFieldBurst),mod.burstMeanResultantLen(indCurCNoFieldBurst));
        modAlignedPyrStatsFieldF.pRSPhaseMeanResultantLen(i) = ranksum(mod.phaseMeanResultantLen(indCurCField),mod.phaseMeanResultantLen(indCurCNoField));
        
        modAlignedPyrStatsFieldF.pWWBurstThetaDiff(i) = circ_wwtest(mod.burstThetaDiff(indCurCFieldBurst(mod.fractBurst > 0)),...
                mod.burstThetaDiff(indCurCNoFieldBurst(mod.fractBurst > 0)));
        modAlignedPyrStatsFieldF.pKBurstThetaDiff(i) = circ_ktest(mod.burstThetaDiff(indCurCFieldBurst(mod.fractBurst > 0)),...
                mod.burstThetaDiff(indCurCNoFieldBurst(mod.fractBurst > 0)));
    
        modAlignedPyrStatsFieldF.pRSMeanFractBurst(i) = ranksum(mod.fractBurst(indCurCField),mod.fractBurst(indCurCNoField));
        modAlignedPyrStatsFieldF.pRSNumSpPerBurstMean(i) = ranksum(mod.numSpPerBurstMean(indCurCField),mod.numSpPerBurstMean(indCurCNoField));
    end
    modAlignedPyrStatsFieldF.pRSFieldWidthC = ranksum(modAlignedPyrStatsFieldF.fieldWidthCFRec{1},modAlignedPyrStatsFieldF.fieldWidthCFRec{2});
    modAlignedPyrStatsFieldF.pRSIndPeakFieldFRec = ranksum(modAlignedPyrStatsFieldF.indPeakFieldFRec{1},modAlignedPyrStatsFieldF.indPeakFieldFRec{2});
    modAlignedPyrStatsFieldF.pRSSkewnessC = ranksum(modAlignedPyrStatsFieldF.skewness{1},modAlignedPyrStatsFieldF.skewness{2});
    modAlignedPyrStatsFieldF.pRSThetaModFreq3C = ranksum(mod.thetaModFreq3(modAlignedPyrStatsFieldF.indCurCField{1}),...
                mod.thetaModFreq3(modAlignedPyrStatsFieldF.indCurCField{2}));
    modAlignedPyrStatsFieldF.pRSPercTrackStartFieldC = ranksum(modAlignedPyrStatsFieldF.percTrackStartField{1},...
                modAlignedPyrStatsFieldF.percTrackStartField{2});
    modAlignedPyrStatsFieldF.pRSPercTrackPeakFieldC = ranksum(modAlignedPyrStatsFieldF.percTrackPeakField{1},...
                modAlignedPyrStatsFieldF.percTrackPeakField{2});
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
    polarhistogram(x2,0:pi/18:2*pi,'Normalization','probability','DisplayStyle','bar',....
        'FaceAlpha',0.5,'FaceColor',colorArr(1,:),'EdgeColor',[0.5 0.5 0.5]);
    title([ti ' p = ' num2str(p)])    
    
    fileName1 = [pathAnal fn];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotDistri(x,y,histBins,xl,leg,pathAnal,fileN,xlimit)
    
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    colorArr = [163 207 98;...
                234 131 114]/255;
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
    colorArr = [163 207 98;...
                234 131 114;...
                163 207 98]/255;
            
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
        set(h,'MarkerSize',8,'Color',colorArr(mod(i,6)+1,:));
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
