function InterneuronModAllRec(task)

    onlyRun = 1;
    methodTheta = 1;
    minFRInt = 3;
    sampleFq = 1250;
    
    RecordingList;
    pathAnal0 = 'Z:\Yingxue\DataAnalysisRaphi\Interneuron\';
    if(task == 1) % including all the neurons
        pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\interneuron\';
    elseif(task == 2) % including AL and PL neurons
        pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\interneuronALPL\';
    else % AL neurons only
        pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\interneuronAL\';
    end
    
    load([pathAnal0 'autoCorrIntAllRec.mat']);
    
%     %% interneurons in no cue passive task
%     disp('No cue')
%     modIntNoCue = accumInterneurons1(listRecordingsNoCuePath,...
%         listRecordingsNoCueFileName,mazeSessionNoCue,minFRInt,1,methodTheta,onlyRun);
%     
%     disp('Active licking')
%     modIntAL = accumInterneurons1(listRecordingsActiveLickPath,...
%         listRecordingsActiveLickFileName,mazeSessionActiveLick,minFRInt,2,methodTheta,onlyRun);
%     
%     disp('Passive licking')
%     modIntPL = accumInterneurons1(listRecordingsPassiveLickPath,...
%         listRecordingsPassiveLickFileName,mazeSessionPassiveLick,minFRInt,3,methodTheta,onlyRun);
%     
%     save([pathAnal 'autoCorrIntAllRec.mat'],'modIntNoCue','modIntAL','modIntPL','-append');     
%     
%     plotConditions(modIntNoCue.thetaModHist,modIntNoCue.phaseMeanDire,...
%         modIntAL.thetaModHist,modIntAL.phaseMeanDire,...
%         modIntPL.thetaModHist,modIntPL.phaseMeanDire,...
%         'Theta modulation (hist)','Theta phase mean direction',[]);
%     
%     plotConditions(modIntNoCue.thetaModHist,modIntNoCue.thetaModInd3,...
%         modIntAL.thetaModHist,modIntAL.thetaModInd3,...
%         modIntPL.thetaModHist,modIntPL.thetaModInd3,...
%         'Theta modulation (hist)','Theta modulation 3',[]);
%     
%     %% plot each interneuron clusters based on the task type
%     plotClusters(modIntNoCue.burstMeanDire,modIntNoCue.phaseMeanDire,...
%         autoCorrIntAll.idxC(autoCorrIntAll.task == autoCorrIntNoCue.task(1)),...
%         'Burst mean direction','Theta phase mean direction','No cue task')
%     plot([0 2*pi],[0 2*pi],'k-')
%     
%     plotClusters(modIntAL.burstMeanDire,modIntAL.phaseMeanDire,...
%         autoCorrIntAll.idxC(autoCorrIntAll.task == autoCorrIntAL.task(1)),...
%         'Burst mean direction','Theta phase mean direction','AL task')
%     plot([0 2*pi],[0 2*pi],'k-')
%     
%     plotClusters(modIntPL.burstMeanDire,modIntPL.phaseMeanDire,...
%         autoCorrIntAll.idxC(autoCorrIntAll.task == autoCorrIntPL.task(1)),...
%         'Burst mean direction','Theta phase mean direction','PL task')
%     plot([0 2*pi],[0 2*pi],'k-')
%     
%     %% number of fields
%     plotClusters(modIntAL.burstMeanDire,modIntAL.nNeuWithField,...
%         autoCorrIntAll.idxC(autoCorrIntAll.task == autoCorrIntAL.task(1)),...
%         'Burst mean direction','Num. fields','AL task')
%     
%     plotClusters(modIntPL.burstMeanDire,modIntPL.nNeuWithField,...
%         autoCorrIntAll.idxC(autoCorrIntAll.task == autoCorrIntPL.task(1)),...
%         'Burst mean direction','Num. fields','PL task')
    
    if(task == 1)
        nNeuWithField = [modIntNoCue.nNeuWithField modIntAL.nNeuWithField modIntPL.nNeuWithField];
        burstMeanDire = [modIntNoCue.burstMeanDire modIntAL.burstMeanDire modIntPL.burstMeanDire];
        nonBurstMeanDire = [modIntNoCue.nonBurstMeanDire modIntAL.nonBurstMeanDire modIntPL.nonBurstMeanDire];
        burstMeanDireStart = [modIntNoCue.burstMeanDireStart modIntAL.burstMeanDireStart modIntPL.burstMeanDireStart];
        numSpPerBurstMean = [modIntNoCue.numSpPerBurstMean modIntAL.numSpPerBurstMean modIntPL.numSpPerBurstMean];
        fractBurst = [modIntNoCue.fractBurst modIntAL.fractBurst modIntPL.fractBurst];

        thetaModHist = [modIntNoCue.thetaModHist modIntAL.thetaModHist modIntPL.thetaModHist];
        thetaModHistH = [modIntNoCue.thetaModHistH modIntAL.thetaModHistH modIntPL.thetaModHistH];
        phaseMeanDire = [modIntNoCue.phaseMeanDire modIntAL.phaseMeanDire modIntPL.phaseMeanDire];
        phaseMeanDireH = [modIntNoCue.phaseMeanDireH modIntAL.phaseMeanDireH modIntPL.phaseMeanDireH];
        maxPhaseArr = [modIntNoCue.maxPhaseFil modIntAL.maxPhaseFil modIntPL.maxPhaseFil];
        maxPhaseArrH = [modIntNoCue.maxPhaseFilH modIntAL.maxPhaseFilH modIntPL.maxPhaseFilH];
        minPhaseArr = [modIntNoCue.minPhaseFil modIntAL.minPhaseFil modIntPL.minPhaseFil];
        minPhaseArrH = [modIntNoCue.minPhaseFilH modIntAL.minPhaseFilH modIntPL.minPhaseFilH];
        phaseDiff = abs(maxPhaseArr - minPhaseArr);
        phaseDiff(phaseDiff < 0) = phaseDiff(phaseDiff < 0) + 360;
        phaseDiffH = abs(maxPhaseArrH - minPhaseArrH);
        phaseDiffH(phaseDiffH < 0) = phaseDiffH(phaseDiffH < 0) + 360;

        thetaModFreq3 = [modIntNoCue.thetaModFreq3 modIntAL.thetaModFreq3 modIntPL.thetaModFreq3];
        thetaModInd3 = [modIntNoCue.thetaModInd3 modIntAL.thetaModInd3 modIntPL.thetaModInd3];
        thetaModInd = [modIntNoCue.thetaModInd modIntAL.thetaModInd modIntPL.thetaModInd];   
        thetaAsym3 = [modIntNoCue.thetaAsym3 modIntAL.thetaAsym3 modIntPL.thetaAsym3];

    %     idxC = [autoCorrIntAll.idxC(autoCorrIntAll.task == autoCorrIntNoCue.task(1))' ...
    %         autoCorrIntAll.idxC(autoCorrIntAll.task == autoCorrIntAL.task(1))' ...
    %         autoCorrIntAll.idxC(autoCorrIntAll.task == autoCorrIntPL.task(1))'];
    %     idxC = [autoCorrIntAll.idxC1(autoCorrIntAll.task == autoCorrIntNoCue.task(1))' ...
    %         autoCorrIntAll.idxC1(autoCorrIntAll.task == autoCorrIntAL.task(1))' ...
    %         autoCorrIntAll.idxC1(autoCorrIntAll.task == autoCorrIntPL.task(1))'];
        idxC = [autoCorrIntAll.idxC2(autoCorrIntAll.task == autoCorrIntNoCue.task(1))' ...
            autoCorrIntAll.idxC2(autoCorrIntAll.task == autoCorrIntAL.task(1))' ...
            autoCorrIntAll.idxC2(autoCorrIntAll.task == autoCorrIntPL.task(1))'];
    elseif(task == 2)
        nNeuWithField = [modIntPL.nNeuWithField modIntAL.nNeuWithField];
        burstMeanDire = [modIntPL.burstMeanDire modIntAL.burstMeanDire];
        nonBurstMeanDire = [modIntPL.nonBurstMeanDire modIntAL.nonBurstMeanDire];
        burstMeanDireStart = [modIntPL.burstMeanDireStart modIntAL.burstMeanDireStart];
        numSpPerBurstMean = [modIntPL.numSpPerBurstMean modIntAL.numSpPerBurstMean];
        fractBurst = [modIntPL.fractBurst modIntAL.fractBurst];

        thetaModHist = [modIntPL.thetaModHist modIntAL.thetaModHist];
        thetaModHistH = [modIntPL.thetaModHistH modIntAL.thetaModHistH];
        phaseMeanDire = [modIntPL.phaseMeanDire modIntAL.phaseMeanDire];
        phaseMeanDireH = [modIntPL.phaseMeanDireH modIntAL.phaseMeanDireH];
        maxPhaseArr = [modIntPL.maxPhaseFil modIntAL.maxPhaseFil];
        maxPhaseArrH = [modIntPL.maxPhaseFilH modIntAL.maxPhaseFilH];
        minPhaseArr = [modIntPL.minPhaseFil modIntAL.minPhaseFil];
        minPhaseArrH = [modIntPL.minPhaseFilH modIntAL.minPhaseFilH];
        phaseDiff = abs(maxPhaseArr - minPhaseArr);
        phaseDiff(phaseDiff < 0) = phaseDiff(phaseDiff < 0) + 360;
        phaseDiffH = abs(maxPhaseArrH - minPhaseArrH);
        phaseDiffH(phaseDiffH < 0) = phaseDiffH(phaseDiffH < 0) + 360;

        thetaModFreq3 = [modIntPL.thetaModFreq3 modIntAL.thetaModFreq3];
        thetaModInd3 = [modIntPL.thetaModInd3 modIntAL.thetaModInd3];
        thetaModInd = [modIntPL.thetaModInd modIntAL.thetaModInd];   
        thetaAsym3 = [modIntPL.thetaAsym3 modIntAL.thetaAsym3];

        idxC = [autoCorrIntAll.idxC2(autoCorrIntAll.task == autoCorrIntPL.task(1))'...
            autoCorrIntAll.idxC2(autoCorrIntAll.task == autoCorrIntAL.task(1))'];
    else
        nNeuWithField = modIntAL.nNeuWithField;
        burstMeanDire = modIntAL.burstMeanDire;
        nonBurstMeanDire = modIntAL.nonBurstMeanDire;
        burstMeanDireStart = modIntAL.burstMeanDireStart;
        numSpPerBurstMean = modIntAL.numSpPerBurstMean;
        fractBurst = modIntAL.fractBurst;

        thetaModHist = modIntAL.thetaModHist;
        thetaModHistH = modIntAL.thetaModHistH;
        phaseMeanDire = modIntAL.phaseMeanDire;
        phaseMeanDireH = modIntAL.phaseMeanDireH;
        maxPhaseArr = modIntAL.maxPhaseFil;
        maxPhaseArrH = modIntAL.maxPhaseFilH;
        minPhaseArr = modIntAL.minPhaseFil;
        minPhaseArrH = modIntAL.minPhaseFilH;
        phaseDiff = abs(maxPhaseArr - minPhaseArr);
        phaseDiff(phaseDiff < 0) = phaseDiff(phaseDiff < 0) + 360;
        phaseDiffH = abs(maxPhaseArrH - minPhaseArrH);
        phaseDiffH(phaseDiffH < 0) = phaseDiffH(phaseDiffH < 0) + 360;

        thetaModFreq3 = modIntAL.thetaModFreq3;
        thetaModInd3 = modIntAL.thetaModInd3;
        thetaModInd = modIntAL.thetaModInd;   
        thetaAsym3 = modIntAL.thetaAsym3;

        idxC = autoCorrIntAll.idxC2(autoCorrIntAll.task == autoCorrIntAL.task(1))';
    end
        
%     figure
%     for i = 1:max(idxC)
%         indCurC = idxC == i;
%         h = plot(phaseMeanDire(indCurC),phaseDiff(indCurC),'o');
%         set(gca,'XLim',[0 2*pi],'YLim',[0 300])
%         pause;
%     end

    modIntStatsField = [];
    for i = 1:max(idxC)
        indCurCField = idxC == i & nNeuWithField >= 2;
        indCurCNoField = idxC == i & nNeuWithField < 1;
        modIntStatsField.meanPeakTo40ms(i,:) = [mean(autoCorrIntAll.peakTo40ms(indCurCField)),mean(autoCorrIntAll.peakTo40ms(indCurCNoField))];
        modIntStatsField.meanPeakTime(i,:) = [mean(autoCorrIntAll.peakTime(indCurCField)),mean(autoCorrIntAll.peakTime(indCurCNoField))];
        
        indCurCFieldBurst = idxC == i & nNeuWithField >= 2 & fractBurst > 0;
        indCurCNoFieldBurst = idxC == i & nNeuWithField < 1 & fractBurst > 0;
        modIntStatsField.meanBurstMeanDire(i,:) = [circ_mean(burstMeanDire(indCurCFieldBurst)'),circ_mean(burstMeanDire(indCurCNoFieldBurst)')];
        modIntStatsField.meanNonBurstMeanDire(i,:) = [circ_mean(nonBurstMeanDire(indCurCFieldBurst)'),circ_mean(nonBurstMeanDire(indCurCNoFieldBurst)')];        
        modIntStatsField.meanBurstMeanDireStart(i,:) = [circ_mean(burstMeanDireStart(indCurCFieldBurst)'),circ_mean(burstMeanDireStart(indCurCNoFieldBurst)')];
        modIntStatsField.meanFractBurst(i,:) = [mean(fractBurst(indCurCField)),mean(fractBurst(indCurCNoField))];
        nonZeroF = numSpPerBurstMean(indCurCField);
        nonZeroF = nonZeroF(nonZeroF > 0);
        nonZeroNoF = numSpPerBurstMean(indCurCNoField);
        nonZeroNoF = nonZeroNoF(nonZeroNoF > 0);
        modIntStatsField.meanNumSpPerBurstMean(i,:) = [mean(nonZeroF),mean(nonZeroNoF)];
        
        modIntStatsField.meanPhaseMeanDire(i,:) = [circ_mean(phaseMeanDire(indCurCField)'),circ_mean(phaseMeanDire(indCurCNoField)')];
        modIntStatsField.meanPhaseMeanDireH(i,:) = [circ_mean(phaseMeanDireH(indCurCField)'),circ_mean(phaseMeanDireH(indCurCNoField)')];
        modIntStatsField.meanMaxPhase(i,:) = [circ_mean(maxPhaseArr(indCurCField)'/180*pi),circ_mean(maxPhaseArr(indCurCNoField)'/180*pi)];
        modIntStatsField.meanMaxPhaseH(i,:) = [circ_mean(maxPhaseArrH(indCurCField)'/180*pi),circ_mean(maxPhaseArrH(indCurCNoField)'/180*pi)];
        modIntStatsField.meanminPhase(i,:) = [circ_mean(minPhaseArr(indCurCField)'/180*pi),circ_mean(minPhaseArr(indCurCNoField)'/180*pi)];
        modIntStatsField.meanminPhaseH(i,:) = [circ_mean(minPhaseArrH(indCurCField)'/180*pi),circ_mean(minPhaseArrH(indCurCNoField)'/180*pi)];
        modIntStatsField.meanPhaseDiff(i,:) = [mean(phaseDiff(indCurCField)),mean(phaseDiff(indCurCNoField))];
        modIntStatsField.meanPhaseDiffH(i,:) = [mean(phaseDiffH(indCurCField)),mean(phaseDiffH(indCurCNoField))];
        modIntStatsField.meanThetaModHist(i,:) = [mean(thetaModHist(indCurCField)),mean(thetaModHist(indCurCNoField))];
        modIntStatsField.meanThetaModHistH(i,:) = [mean(thetaModHistH(indCurCField)),mean(thetaModHistH(indCurCNoField))];
        
        modIntStatsField.meanThetaModFreq3(i,:) = [mean(thetaModFreq3(indCurCField)),mean(thetaModFreq3(indCurCNoField))];
        modIntStatsField.meanThetaAsym3(i,:) = [mean(thetaAsym3(indCurCField)),mean(thetaAsym3(indCurCNoField))];
        modIntStatsField.meanThetaModInd3(i,:) = [mean(thetaModInd3(indCurCField)),mean(thetaModInd3(indCurCNoField))];
        modIntStatsField.meanThetaModInd(i,:) = [mean(thetaModInd(indCurCField)),mean(thetaModInd(indCurCNoField))];
        
        modIntStatsField.pRSPeakTo40ms(i) = ranksum(autoCorrIntAll.peakTo40ms(indCurCField),autoCorrIntAll.peakTo40ms(indCurCNoField));
        modIntStatsField.pRSPeakTime(i) = ranksum(autoCorrIntAll.peakTime(indCurCField),autoCorrIntAll.peakTime(indCurCNoField));
        
        modIntStatsField.pWWBurstMeanDire(i) = circ_wwtest(burstMeanDire(indCurCField)',burstMeanDire(indCurCNoField)');
        modIntStatsField.pWWNonBurstMeanDire(i) = circ_wwtest(nonBurstMeanDire(indCurCField)',nonBurstMeanDire(indCurCNoField)');
        modIntStatsField.pWWBurstMeanDireStart(i) = circ_wwtest(burstMeanDireStart(indCurCFieldBurst)',burstMeanDireStart(indCurCNoFieldBurst)');
        modIntStatsField.pRSFractBurst(i) = ranksum(fractBurst(indCurCField),fractBurst(indCurCNoField));
        modIntStatsField.pRSNumSpPerBurstMean(i) = ranksum(nonZeroF,nonZeroNoF);
        
        modIntStatsField.pKBurstMeanDire(i) = circ_ktest(burstMeanDire(indCurCField)',burstMeanDire(indCurCNoField)');
        modIntStatsField.pKNonBurstMeanDire(i) = circ_ktest(nonBurstMeanDire(indCurCField)',nonBurstMeanDire(indCurCNoField)');
        modIntStatsField.pKBurstMeanDireStart(i) = circ_ktest(burstMeanDireStart(indCurCFieldBurst)',burstMeanDireStart(indCurCNoFieldBurst)');
        
        modIntStatsField.pWWPhaseMeanDire(i) = circ_wwtest(phaseMeanDire(indCurCField)',phaseMeanDire(indCurCNoField)');
        modIntStatsField.pWWPhaseMeanDireH(i) = circ_wwtest(phaseMeanDireH(indCurCField)',phaseMeanDireH(indCurCNoField)');
        modIntStatsField.pWWMaxPhase(i) = circ_wwtest(maxPhaseArr(indCurCField)'/180*pi,maxPhaseArr(indCurCNoField)'/180*pi);
        modIntStatsField.pWWMaxPhaseH(i) = circ_wwtest(maxPhaseArrH(indCurCField)'/180*pi,maxPhaseArrH(indCurCNoField)'/180*pi);
        modIntStatsField.pWWMinPhase(i) = circ_wwtest(minPhaseArr(indCurCField)'/180*pi,minPhaseArr(indCurCNoField)'/180*pi);
        modIntStatsField.pWWMinPhaseH(i) = circ_wwtest(minPhaseArrH(indCurCField)'/180*pi,minPhaseArrH(indCurCNoField)'/180*pi);
        
        modIntStatsField.pKPhaseMeanDire(i) = circ_ktest(phaseMeanDire(indCurCField)',phaseMeanDire(indCurCNoField)');
        modIntStatsField.pKPhaseMeanDireH(i) = circ_ktest(phaseMeanDireH(indCurCField)',phaseMeanDireH(indCurCNoField)');
        modIntStatsField.pKMaxPhase(i) = circ_ktest(maxPhaseArr(indCurCField)'/180*pi,maxPhaseArr(indCurCNoField)'/180*pi);
        modIntStatsField.pKMaxPhaseH(i) = circ_ktest(maxPhaseArrH(indCurCField)'/180*pi,maxPhaseArrH(indCurCNoField)'/180*pi);
        modIntStatsField.pKMinPhase(i) = circ_ktest(minPhaseArr(indCurCField)'/180*pi,minPhaseArr(indCurCNoField)'/180*pi);
        modIntStatsField.pKMinPhaseH(i) = circ_ktest(minPhaseArrH(indCurCField)'/180*pi,minPhaseArrH(indCurCNoField)'/180*pi);
        
        modIntStatsField.pRSPhaseDiff(i) = ranksum(phaseDiff(indCurCField),phaseDiff(indCurCNoField));
        modIntStatsField.pRSPhaseDiffH(i) = ranksum(phaseDiffH(indCurCField),phaseDiffH(indCurCNoField));
        modIntStatsField.pRSThetaModHist(i) = ranksum(thetaModHist(indCurCField),thetaModHist(indCurCNoField));
        modIntStatsField.pRSThetaModHistH(i) = ranksum(thetaModHistH(indCurCField),thetaModHistH(indCurCNoField));
        
        modIntStatsField.pRSThetaModFreq3(i) = ranksum(thetaModFreq3(indCurCField),thetaModFreq3(indCurCNoField));
        modIntStatsField.pRSThetaAsym3(i) = ranksum(thetaAsym3(indCurCField),thetaAsym3(indCurCNoField));  
        modIntStatsField.pRSThetaModInd3(i) = ranksum(thetaModInd3(indCurCField),thetaModInd3(indCurCNoField));
        modIntStatsField.pRSThetaModInd(i) = ranksum(thetaModInd(indCurCField),thetaModInd(indCurCNoField));
        
        modIntStatsField.pWWBurstVsNonBurstMeanDire(i) = circ_wwtest(burstMeanDire(indCurCField)',nonBurstMeanDire(indCurCField)');
        modIntStatsField.pWWBurstVsThetaMeanDire(i) = circ_wwtest(burstMeanDire(indCurCField)',phaseMeanDire(indCurCField)');
        modIntStatsField.pWWNonBurstVsThetaMeanDire(i) = circ_wwtest(nonBurstMeanDire(indCurCField)',phaseMeanDire(indCurCField)');
        modIntStatsField.pWWBurstStartVsThetaMeanDire(i) = circ_wwtest(burstMeanDireStart(indCurCFieldBurst)',phaseMeanDire(indCurCField)');
        
        modIntStatsField.pKBurstVsNonBurstMeanDire(i) = circ_ktest(burstMeanDire(indCurCField)',nonBurstMeanDire(indCurCField)');
        modIntStatsField.pKBurstVsThetaMeanDire(i) = circ_ktest(burstMeanDire(indCurCField)',phaseMeanDire(indCurCField)');
        modIntStatsField.pKNonBurstVsThetaMeanDire(i) = circ_ktest(nonBurstMeanDire(indCurCField)',phaseMeanDire(indCurCField)');
        modIntStatsField.pKBurstStartVsThetaMeanDire(i) = circ_ktest(burstMeanDireStart(indCurCFieldBurst)',phaseMeanDire(indCurCField)');
    end
    
    if(task == 1)
        save([pathAnal 'autoCorrIntAllRec.mat'],'modIntStatsField','-append'); 
    elseif(task == 2)
        save([pathAnal 'autoCorrIntALPLRec.mat'],'modIntStatsField'); 
    else
        save([pathAnal 'autoCorrIntALRec.mat'],'modIntStatsField'); 
    end
    
    %% compare cluster 1, recordings with fields vs. without
    colorSel = 1;
    cluster = 4;
    indCurCField = idxC == cluster & nNeuWithField >= 2;
    indCurCNoField = idxC == cluster & nNeuWithField < 1;
    plotBoxPlot(phaseMeanDire(indCurCField),...
        phaseMeanDire(indCurCNoField),'Mean theta phase',['Int_ThetaMeanC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[0 2*pi],modIntStatsField.pWWPhaseMeanDire(cluster),colorSel);
    
    plotBoxPlot(maxPhaseArr(indCurCField)/180*pi,...
        maxPhaseArr(indCurCNoField)/180*pi,'Max theta phase',['Int_ThetaMaxC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[0 2*pi],modIntStatsField.pWWMaxPhase(cluster),colorSel);
  
    minPhaseField = minPhaseArr(indCurCField)/180*pi;
    minPhaseNoField = minPhaseArr(indCurCNoField)/180*pi;
    minPhaseField(minPhaseField > pi) = minPhaseField(minPhaseField > pi) - 2*pi;
    minPhaseNoField(minPhaseNoField > pi) = minPhaseNoField(minPhaseNoField > pi) - 2*pi;
    plotBoxPlot(minPhaseField,...
        minPhaseNoField,'Min theta phase',['Int_ThetaMinC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[-pi pi],modIntStatsField.pWWMinPhase(cluster),colorSel);
    
    plotBoxPlot(fractBurst(indCurCField),...
        fractBurst(indCurCNoField),'Fract. burst',['Int_FractBurstC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[-0.1 0.4],modIntStatsField.pRSFractBurst(cluster),colorSel);
    
    plotBoxPlot(thetaModHist(indCurCField),...
        thetaModHist(indCurCNoField),'Theta modulation',['Int_ThetaModHistC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[0 1.2],modIntStatsField.pRSThetaModHist(cluster),colorSel);
    
    nonZeroF = numSpPerBurstMean(indCurCField);
    nonZeroF = nonZeroF(nonZeroF > 0);
    nonZeroNoF = numSpPerBurstMean(indCurCNoField);
    nonZeroNoF = nonZeroNoF(nonZeroNoF > 0);
    plotBoxPlot(nonZeroF,...
        nonZeroNoF,'Burst length',['Int_NumSpPerBurstMeanC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[1 4],modIntStatsField.pRSNumSpPerBurstMean(cluster),colorSel);
    
    plotBoxPlot(nonBurstMeanDire(indCurCField),...
        nonBurstMeanDire(indCurCNoField),'Non-burst phase',['Int_NonBurstMeanDireC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[1 5],modIntStatsField.pWWNonBurstMeanDire(cluster),colorSel);
    
    plotPolarPlot(phaseMeanDire(indCurCField),...
        phaseMeanDire(indCurCNoField),'Mean theta phase',...
        ['Int_ThetaMeanC1' num2str(cluster) 'Field-NoFieldPolar'],pathAnal,modIntStatsField.pWWPhaseMeanDire(cluster),colorSel);
    
    plotPolarPlot(maxPhaseArr(indCurCField)/180*pi,...
        maxPhaseArr(indCurCNoField)/180*pi,'Max theta phase',['Int_MaxPhaseC' num2str(cluster) 'Field-NoFieldPolar'],...
        pathAnal,modIntStatsField.pWWMaxPhase(cluster),colorSel);
    
    plotPolarPlot(minPhaseArr(indCurCField)/180*pi,...
        minPhaseArr(indCurCNoField)/180*pi,'Min theta phase',['Int_MinPhaseC' num2str(cluster) 'Field-NoFieldPolar'],...
        pathAnal,modIntStatsField.pWWMinPhase(cluster),colorSel);
    
    plotPolarPlot(burstMeanDire(indCurCField),...
        burstMeanDire(indCurCNoField),'Burst phase',['Int_BurstMeanDireC' num2str(cluster) 'Field-NoFieldPolar'],...
        pathAnal,modIntStatsField.pKBurstMeanDire(cluster),colorSel);
    
    %% compare cluster 2, recordings with fields vs. without
    colorSel = 6;
    cluster = 2;
    indCurCField = idxC == cluster & nNeuWithField >= 2;
    indCurCNoField = idxC == cluster & nNeuWithField < 1;
    plotBoxPlot(phaseMeanDire(indCurCField),...
        phaseMeanDire(indCurCNoField),'Mean theta phase',['Int_ThetaMeanC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[0 pi],modIntStatsField.pWWPhaseMeanDire(cluster),colorSel);
    
    plotBoxPlot(thetaModHist(indCurCField),...
        thetaModHist(indCurCNoField),'Theta modulation',['Int_ThetaModHistC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[0 1.2],modIntStatsField.pRSThetaModHist(cluster),colorSel);
    
    plotBoxPlot(fractBurst(indCurCField),...
        fractBurst(indCurCNoField),'Fract. burst',['Int_FractBurstC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[-0.1 0.4],modIntStatsField.pRSFractBurst(cluster),colorSel);
    
    nonZeroF = numSpPerBurstMean(indCurCField);
    nonZeroF = nonZeroF(nonZeroF > 0);
    nonZeroNoF = numSpPerBurstMean(indCurCNoField);
    nonZeroNoF = nonZeroNoF(nonZeroNoF > 0);
    plotBoxPlot(nonZeroF,...
        nonZeroNoF,'Burst length',['Int_NumSpPerBurstMeanC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[1 4],modIntStatsField.pRSNumSpPerBurstMean(cluster),colorSel);
    
    plotBoxPlot(nonBurstMeanDire(indCurCField),...
        nonBurstMeanDire(indCurCNoField),'Non-burst phase',['Int_NonBurstMeanDireC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[0 pi],modIntStatsField.pWWNonBurstMeanDire(cluster),colorSel);
    
    plotBoxPlot(burstMeanDire(indCurCField),...
        burstMeanDire(indCurCNoField),'Burst phase',['Int_BurstMeanDireC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[0 2*pi],modIntStatsField.pWWBurstMeanDire(cluster),colorSel);
    
    plotPolarPlot(phaseMeanDire(indCurCField),...
        phaseMeanDire(indCurCNoField),'Mean theta phase',...
        ['Int_ThetaMeanC' num2str(cluster) 'Field-NoFieldPolar'],pathAnal,modIntStatsField.pWWPhaseMeanDire(cluster),colorSel);
    
    plotPolarPlot(phaseMeanDire(idxC == cluster),...
        [],'Mean theta phase',...
        ['Int_ThetaMeanCPolar' num2str(cluster)],pathAnal,1,colorSel);
    
    plotPolarPlot(burstMeanDire(idxC == cluster),...
        [],'Burst phase',...
        ['Int_BurstMeanDireCPolar' num2str(cluster)],pathAnal,1,colorSel);
    
    plotPolarPlot(maxPhaseArr(indCurCField)/180*pi,...
        maxPhaseArr(indCurCNoField)/180*pi,'Max theta phase',['Int_MaxPhaseC' num2str(cluster) 'Field-NoFieldPolar'],...
        pathAnal,modIntStatsField.pWWMaxPhase(cluster),colorSel);
    
    plotPolarPlot(minPhaseArr(indCurCField)/180*pi,...
        minPhaseArr(indCurCNoField)/180*pi,'Min theta phase',['Int_MinPhaseC' num2str(cluster) 'Field-NoFieldPolar'],...
        pathAnal,modIntStatsField.pWWMinPhase(cluster),colorSel);
    
    plotPolarPlot(burstMeanDire(indCurCField),...
        burstMeanDire(indCurCNoField),'Burst phase',['Int_BurstMeanDireC' num2str(cluster) 'Field-NoFieldPolar'],...
        pathAnal,modIntStatsField.pWWBurstMeanDire(cluster),colorSel);
    
    %% compare cluster 4, recordings with fields vs. without
    colorSel = 4;
    cluster = 1;
    indCurCField = idxC == cluster & nNeuWithField >= 2;
    indCurCNoField = idxC == cluster & nNeuWithField < 1;
    plotBoxPlot(phaseMeanDire(indCurCField),...
        phaseMeanDire(indCurCNoField),'Mean theta phase',['Int_ThetaMeanC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[0 pi],modIntStatsField.pWWPhaseMeanDire(cluster),colorSel);
    
    plotBoxPlot(thetaModHist(indCurCField),...
        thetaModHist(indCurCNoField),'Theta modulation',['Int_ThetaModHistC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[0 1.2],modIntStatsField.pRSThetaModHist(cluster),colorSel);
    
    plotBoxPlot(fractBurst(indCurCField),...
        fractBurst(indCurCNoField),'Fract. burst',['Int_FractBurstC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[-0.1 0.2],modIntStatsField.pRSFractBurst(cluster),colorSel);
    
    nonZeroF = numSpPerBurstMean(indCurCField);
    nonZeroF = nonZeroF(nonZeroF > 0);
    nonZeroNoF = numSpPerBurstMean(indCurCNoField);
    nonZeroNoF = nonZeroNoF(nonZeroNoF > 0);
    plotBoxPlot(nonZeroF,...
        nonZeroNoF,'Burst length',['Int_NumSpPerBurstMeanC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[1 4],modIntStatsField.pRSNumSpPerBurstMean(cluster),colorSel);
    
    plotBoxPlot(nonBurstMeanDire(indCurCField),...
        nonBurstMeanDire(indCurCNoField),'Non-burst phase',['Int_NonBurstMeanDireC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[0 pi],modIntStatsField.pWWNonBurstMeanDire(cluster),colorSel);
    
    plotBoxPlot(burstMeanDire(indCurCField),...
        burstMeanDire(indCurCNoField),'Burst phase',['Int_BurstMeanDireC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[0 2*pi],modIntStatsField.pWWBurstMeanDire(cluster),colorSel);
    
    plotPolarPlot(phaseMeanDire(indCurCField),...
        phaseMeanDire(indCurCNoField),'Mean theta phase',...
        ['Int_ThetaMeanC' num2str(cluster) 'Field-NoFieldPolar'],pathAnal,modIntStatsField.pWWPhaseMeanDire(cluster),colorSel);
    
    plotPolarPlot(maxPhaseArr(indCurCField)/180*pi,...
        maxPhaseArr(indCurCNoField)/180*pi,'Max theta phase',['Int_MaxPhaseC' num2str(cluster) 'Field-NoFieldPolar'],...
        pathAnal,modIntStatsField.pWWMaxPhase(cluster),colorSel);
    
    plotPolarPlot(minPhaseArr(indCurCField)/180*pi,...
        minPhaseArr(indCurCNoField)/180*pi,'Min theta phase',['Int_MinPhaseC' num2str(cluster) 'Field-NoFieldPolar'],...
        pathAnal,modIntStatsField.pWWMinPhase(cluster),colorSel);
    
    plotPolarPlot(burstMeanDire(indCurCField),...
        burstMeanDire(indCurCNoField),'Burst phase',['Int_BurstMeanDireC' num2str(cluster) 'Field-NoFieldPolar'],...
        pathAnal,modIntStatsField.pWWBurstMeanDire(cluster),colorSel);
    
    %% compare cluster 6, recordings with fields vs. without
    colorSel = 2;
    cluster = 3;
    indCurCField = idxC == cluster & nNeuWithField >= 2;
    indCurCNoField = idxC == cluster & nNeuWithField < 1;
    plotBoxPlot(phaseMeanDire(indCurCField),...
        phaseMeanDire(indCurCNoField),'Mean theta phase',['Int_ThetaMeanC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[0 pi*2],modIntStatsField.pWWPhaseMeanDire(cluster),colorSel);
    
    plotBoxPlot(thetaModHist(indCurCField),...
        thetaModHist(indCurCNoField),'Theta modulation',['Int_ThetaModHistC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[0 1.2],modIntStatsField.pRSThetaModHist(cluster),colorSel);
    
    plotBoxPlot(fractBurst(indCurCField),...
        fractBurst(indCurCNoField),'Fract. burst',['Int_FractBurstC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[-0.1 0.2],modIntStatsField.pRSFractBurst(cluster),colorSel);
    
    nonZeroF = numSpPerBurstMean(indCurCField);
    nonZeroF = nonZeroF(nonZeroF > 0);
    nonZeroNoF = numSpPerBurstMean(indCurCNoField);
    nonZeroNoF = nonZeroNoF(nonZeroNoF > 0);
    plotBoxPlot(nonZeroF,...
        nonZeroNoF,'Burst length',['Int_NumSpPerBurstMeanC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[1 4],modIntStatsField.pRSNumSpPerBurstMean(cluster),colorSel);
    
    plotBoxPlot(nonBurstMeanDire(indCurCField),...
        nonBurstMeanDire(indCurCNoField),'Non-burst phase',['Int_NonBurstMeanDireC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[0 pi*2],modIntStatsField.pWWNonBurstMeanDire(cluster),colorSel);
    
    plotBoxPlot(burstMeanDire(indCurCField),...
        burstMeanDire(indCurCNoField),'Burst phase',['Int_BurstMeanDireC' num2str(cluster) 'Field-NoFieldBox'],...
        pathAnal,[0 2*pi],modIntStatsField.pWWBurstMeanDire(cluster),colorSel);
    
    plotPolarPlot(phaseMeanDire(indCurCField),...
        phaseMeanDire(indCurCNoField),'Mean theta phase',...
        ['Int_ThetaMeanC' num2str(cluster) 'Field-NoFieldPolar'],pathAnal,modIntStatsField.pWWPhaseMeanDire(cluster),colorSel);
    
    plotPolarPlot(maxPhaseArr(indCurCField)/180*pi,...
        maxPhaseArr(indCurCNoField)/180*pi,'Max theta phase',['Int_MaxPhaseC' num2str(cluster) 'Field-NoFieldPolar'],...
        pathAnal,modIntStatsField.pWWMaxPhase(cluster),colorSel);
    
    plotPolarPlot(minPhaseArr(indCurCField)/180*pi,...
        minPhaseArr(indCurCNoField)/180*pi,'Min theta phase',['Int_MinPhaseC' num2str(cluster) 'Field-NoFieldPolar'],...
        pathAnal,modIntStatsField.pWWMinPhase(cluster),colorSel);
    
    plotPolarPlot(burstMeanDire(indCurCField),...
        burstMeanDire(indCurCNoField),'Burst phase',['Int_BurstMeanDireC' num2str(cluster) 'Field-NoFieldPolar'],...
        pathAnal,modIntStatsField.pWWBurstMeanDire(cluster),colorSel);
    
    %%
%     plotClustersFieldsPhase(pathAnal,burstMeanDire,phaseMeanDire,...
%         idxC,indField,...
%         'Burst mean direction','Theta phase mean direction','All tasks',...
%         'Interneuron_BurstMeanVsThetaMean');
%     
%     plotClustersFields(pathAnal,thetaModHist,phaseMeanDire,...
%         idxC,indField,...
%         'Theta modulation (hist)','Theta phase mean direction','All tasks',...
%         'Interneuron_ThetaModHistVsThetaMean');
%     
%     plotClustersFields(pathAnal,thetaModHistH,phaseMeanDireH,...
%         idxC,indField,...
%         'Theta modulation (hist hilbert)','Theta phase mean direction (hilbert)','All tasks',...
%         'Interneuron_ThetaModHistHVsThetaMean');
%    
%     plotClustersFields(pathAnal,burstMeanDire,fractBurst,...
%         idxC,indField,...
%         'Burst mean direction','Fract burst','All tasks',...
%         'Interneuron_BurstMeanVsFractBurst');
%     
% %     plotClustersFields(pathAnal,phaseMeanDire,thetaModFreq3,...
% %         idxC,indField,...
% %         'Theta phase mean direction','Theta modulation freq. (Hz)','All tasks',...
% %         'Interneuron_ThetaMeanVsModFreq');
%     
%     plotClustersFields(pathAnal,phaseMeanDire,phaseDiff,...
%         idxC,indField,...
%         'Theta phase mean direction','Max-min theta phase','All tasks',...
%         'Interneuron_ThetaMeanVsPhaseDiff');
%     
%     plotClustersFields(pathAnal,phaseMeanDire,minPhaseArr,...
%         idxC,indField,...
%         'Theta phase mean direction','Min theta phase','All tasks',...
%         'Interneuron_ThetaMeanVsMinPhase');
%     
%     plotClustersFields(pathAnal,phaseMeanDire,thetaAsym3,...
%         idxC,indField,...
%         'Theta phase mean direction','Theta asymmetry','All tasks',...
%         'Interneuron_ThetaMeanVsThetaAsym');
%     
%     plotClustersFieldsPhase(pathAnal,burstMeanDire,nonBurstMeanDire,...
%         idxC,indField,...
%         'Burst mean direction','Nonburst mean direction','All tasks',...
%         'Interneuron_BurstMeanVsNonBurstMean');
%     
%     plotClustersFields(pathAnal,thetaModFreq3,thetaAsym3,...
%         idxC,indField,...
%         'Theta modulation freq. (Hz)','Theta asymmetry','All tasks',...
%         'Interneuron_ModFreqVsThetaAsym');

end

function modInt = accumInterneurons1(paths,filenames,mazeSess,minFRInt,task,methodTheta,onlyRun)

    %% interneurons in no cue passive task
    numRec = size(paths,1);
    modInt = struct('task',[],... % no cue - 1, AL - 2, PL - 3
                              'indRec',[],... % recording index
                              'indNeu',[],... % neuron indices
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
                              'nNeuWithField',[],... % number of neurons with fields in the recording
                              'nNeuWithFieldAligned',[],... % number of neurons with fields in the recording
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
        if(i == 3)
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
        
        fileNameFR = [filenames(i,:) '_FR_Run' num2str(onlyRun) '.mat'];
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
                                      
%         indNeu = autoCorr.isInterneuron == 1 & cluList.firingRate > minFRInt;
        indNeu = autoCorr.isInterneuron == 1 & mFR.mFR > minFRInt;         
        modInt.task = [modInt.task task*ones(1,sum(indNeu))];
        modInt.indRec = [modInt.indRec i*ones(1,sum(indNeu))];
        modInt.indNeu = [modInt.indNeu find(indNeu == 1)]; 
        
        modInt.thetaMod = [modInt.thetaMod thetaModSessTmp.thetaMod(indNeu)];
        modInt.trough = [modInt.trough thetaModSessTmp.trough(indNeu)];
        modInt.peak = [modInt.peak thetaModSessTmp.peak(indNeu)];
        modInt.thetaModInd = [modInt.thetaModInd thetaModSessTmp.thetaModInd(indNeu)];
        modInt.peakT3 = [modInt.peakT3 thetaModSessTmp.peakT3(indNeu)];
        modInt.troughT3 = [modInt.troughT3 thetaModSessTmp.troughT3(indNeu)];
        modInt.thetaAsym3 = [modInt.thetaAsym3 ...
            (abs(thetaModSessTmp.peakT3(indNeu))-abs(thetaModSessTmp.troughT3(indNeu)))./...
            (abs(thetaModSessTmp.peakT3(indNeu)))];
        modInt.thetaModFreq3 = [modInt.thetaModFreq3 ...
            1000./(abs(thetaModSessTmp.peakT3(indNeu)))];
        modInt.thetaModInd3 = [modInt.thetaModInd3 thetaModSessTmp.thetaModInd3(indNeu)];
                
        nFieldArr = zeros(1,sum(indNeu));
        if(length(pFRNonStimGoodStruct.indLapList) >= paramF.minNumTr) % more than 15 trirals
            if(~isempty(fieldSpCorrSessGoodTr))
                numField = length(unique(fieldSpCorrSessGoodTr.indNeuron));
            else
                numField = [];
            end
            if(~isempty(numField))
                nFieldArr = numField * ones(1,sum(indNeu));
            end
        end
        modInt.nNeuWithField = [modInt.nNeuWithField nFieldArr];
        
        nFieldArr = zeros(1,sum(indNeu));
        if(length(pFRNonStimGoodStruct.indLapList) >= paramF.minNumTr) % more than 15 trirals
            if(~isempty(fieldSpCorrSessNonStimGood))
                numField = length(unique(fieldSpCorrSessNonStimGood.indNeuron));
            else
                numField = [];
            end                
            if(~isempty(numField))
                nFieldArr = numField * ones(1,sum(indNeu));
            end
        end
        modInt.nNeuWithFieldAligned = [modInt.nNeuWithFieldAligned nFieldArr];
        
        nNeurons = length(cluList.firingRate);
        numSpPerBurst = zeros(1,nNeurons);
        for m = 1:nNeurons
            if(~isempty(burstIsi.numSpPerBurst{m}))
                numSpPerBurst(m) = mean(burstIsi.numSpPerBurst{m});
            end
        end
        modInt.fractBurst = [modInt.fractBurst burstIsi.fractBurst(indNeu)];
        modInt.numSpPerBurstMean = [modInt.numSpPerBurstMean numSpPerBurst(indNeu)]; 
        modInt.burstMeanDire = [modInt.burstMeanDire burstIsi.meanDire(indNeu)];        
        modInt.burstMeanResultantLen = [modInt.burstMeanResultantLen burstIsi.meanResultantLen(indNeu)];
        
        modInt.nonBurstMeanDire = [modInt.nonBurstMeanDire burstIsi.meanDireNonBurst(indNeu)];
        modInt.nonBurstMeanResultantLen = [modInt.nonBurstMeanResultantLen burstIsi.meanResultantLenNonBurst(indNeu)];
        
        modInt.burstMeanDireStart = [modInt.burstMeanDireStart burstIsi.meanDireStart(indNeu)];        
        modInt.burstMeanResultantLenStart = [modInt.burstMeanResultantLenStart burstIsi.meanResultantLenStart(indNeu)];
        
        modInt.minPhaseFil = [modInt.minPhaseFil spikeThetaPhase.minPhaseFilArr(indNeu)];
        modInt.maxPhaseFil = [modInt.maxPhaseFil spikeThetaPhase.maxPhaseFilArr(indNeu)];
        modInt.phaseMeanDire = [modInt.phaseMeanDire spikeThetaPhase.meanDire(indNeu)];
        modInt.phaseMeanResultantLen = [modInt.phaseMeanResultantLen spikeThetaPhase.meanResultantLen(indNeu)];
        modInt.thetaModHist = [modInt.thetaModHist spikeThetaPhase.thetaMod(indNeu)]; 
        
        modInt.minPhaseFilH = [modInt.minPhaseFilH spikeThetaPhaseH.minPhaseFilArr(indNeu)];
        modInt.maxPhaseFilH = [modInt.maxPhaseFilH spikeThetaPhaseH.maxPhaseFilArr(indNeu)];
        modInt.phaseMeanDireH = [modInt.phaseMeanDireH spikeThetaPhaseH.meanDire(indNeu)];
        modInt.phaseMeanResultantLenH = [modInt.phaseMeanResultantLenH spikeThetaPhaseH.meanResultantLen(indNeu)];
        modInt.thetaModHistH = [modInt.thetaModHistH spikeThetaPhaseH.thetaMod(indNeu)]; 
    end
    
    indDire = find(modInt.burstMeanDire < 0);
    modInt.burstMeanDire(indDire) = modInt.burstMeanDire(indDire) + 2*pi;
    
    indDire = find(modInt.nonBurstMeanDire < 0);
    modInt.nonBurstMeanDire(indDire) = modInt.nonBurstMeanDire(indDire) + 2*pi;
    
    indDire = find(modInt.burstMeanDireStart < 0);
    modInt.burstMeanDireStart(indDire) = modInt.burstMeanDireStart(indDire) + 2*pi;
    
    indDire = find(modInt.phaseMeanDire < 0);
    modInt.phaseMeanDire(indDire) = modInt.phaseMeanDire(indDire) + 2*pi;
    
    indDire = find(modInt.phaseMeanDireH < 0);
    modInt.phaseMeanDireH(indDire) = modInt.phaseMeanDireH(indDire) + 2*pi;
end

function plotBoxPlot(x1,x2,yl,fn,pathAnal,ylimit,p,colorSel)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 200 400])
    if(colorSel == 6)
        colorArr = [...
                127 127 229;...
                125 125 168]/255;
    elseif(colorSel == 1)            
        colorArr = [...
                230 127 128;...
                153 85 86]/255;
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

function plotPolarPlot(x1,x2,ti,fn,pathAnal,p,colorSel)
    [figNew,pos] = CreateFig();
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    if(colorSel == 6)
        colorArr = [...
                121 181 226;...
                0 0 200]/255;
    elseif(colorSel == 1)            
        colorArr = [...
                230 127 128;...
                153 85 86]/255;
    else        
        colorArr = [163 207 98;... 
            63 79 37]/255;
    end
    
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
%     colorArr = [0.2 0.2 0.9;...
%                 0.9 0.2 0.2;...
%                 0.2 0.9 0.2;...
%                 0.6 0.6 0.6;...
%                 0.3 0.3 0.7;...
%                 0.3 0.7 0.3;...
%                 0.2 0.5 0.8;...
%                 0.2 0.8 0.5;...
%                 0.8 0.5 0.2];
            
    colorArr = [0.5 0.5 0.9;...
        0.9 0.5 0.5;...
        0.3 0.3 0.7;...
        0.7 0.3 0.3;...
        0.5 0.9 0.5;...
        0.2 0.5 0.8;...
        0.2 0.8 0.5;...
        0.8 0.5 0.2;...
        0.3 0.7 0.3];
    hold on;
    for i = 1:max(idx)
        indTmp = idx == i;
        h = plot(x(indTmp),y(indTmp),'.');
        set(h,'MarkerSize',20,'Color',colorArr(mod(i,6)+1,:));
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
    
    plotClusters(x(~indField),y(~indField),...
        idx(~indField),xl,yl,[ti ' - no field']);
    plot([0 2*pi],[0 2*pi],'k-')
    
    fileName1 = [pathAnal fileN 'NoField'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
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
