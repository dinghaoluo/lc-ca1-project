function plotPVSSTPutativeCells()

    pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\interneuron\';
    if(exist([pathAnal 'autoCorrIntAllRec.mat']))
        load([pathAnal 'autoCorrIntAllRec.mat']);
        load([pathAnal 'initPeakIntAllRec.mat'],'modInt1NoCue','modInt1AL','modInt1PL');
    end
    
    posPhase = 2.5:5:1080;
    indPhase = find(posPhase>= 360); 
    
    avgFRProfile = [modInt1NoCue.avgFRProfile; modInt1AL.avgFRProfile; modInt1PL.avgFRProfile];
    
    % SST
    cluSST1 = [2];
    peakTime = autoCorrIntAll.peakTime(autoCorrIntAll.idxC2 == cluSST1);
    [~,indPeakTime] = sort(peakTime);
    SSTCCGVal = autoCorrIntAll.ccgVal(autoCorrIntAll.idxC2 == cluSST1,:);
    SSTCCGVal = SSTCCGVal(indPeakTime,:);
    for i = 1:size(SSTCCGVal,1)
        SSTCCGVal(i,:) = SSTCCGVal(i,:)/max(SSTCCGVal(i,:));
    end
    SSTHistPhaseFil = autoCorrIntAll.histPhaseFil(autoCorrIntAll.idxC2 == cluSST1,indPhase);
    SSTHistPhaseFil = SSTHistPhaseFil(indPeakTime,:);
    for i = 1:size(SSTHistPhaseFil,1)
        SSTHistPhaseFil(i,:) = SSTHistPhaseFil(i,:)/max(SSTHistPhaseFil(i,:));
    end
    
    idxTagSST1 = autoCorrIntTag.cellType == 2 & autoCorrIntTag.idxC2 == cluSST1;
    peakTime = autoCorrIntTag.peakTime(idxTagSST1);
    indRec = autoCorrIntTag.indRec(idxTagSST1);
    indNeu = autoCorrIntTag.indNeu(idxTagSST1);
    indSST = [];
    for i = 1:length(indRec)
        ind = find(autoCorrIntAll.indNeu == indNeu(i) & autoCorrIntAll.indRec == indRec(i) ...
            & autoCorrIntAll.task == 2);
        indSST = [indSST ind];
    end
    [~,indPeakTime] = sort(peakTime);
    SSTCCGValTag = autoCorrIntTag.ccgVal(idxTagSST1,:);
    SSTCCGValTag = SSTCCGValTag(indPeakTime,:);
    for i = 1:size(SSTCCGValTag,1)
        SSTCCGValTag(i,:) = SSTCCGValTag(i,:)/max(SSTCCGValTag(i,:));
    end
    SSTHistPhaseFilTag = autoCorrIntAll.histPhaseFil(indSST,indPhase);
    SSTHistPhaseFilTag = SSTHistPhaseFilTag(indPeakTime,:);
    for i = 1:size(SSTHistPhaseFilTag,1)
        SSTHistPhaseFilTag(i,:) = SSTHistPhaseFilTag(i,:)/max(SSTHistPhaseFilTag(i,:));
    end
    
    p = randperm(size(SSTCCGVal,1));
    p = sort(p(1:15));
    plotCCG(-50:50,SSTCCGVal(1:20,:),'Putative SST',pathAnal,['putativeSSTCCG' num2str(cluSST1)]);
    plotCCG(-50:50,SSTCCGValTag(1:5,:),'Tagged SST',pathAnal,['taggedSSTCCG'  num2str(cluSST1)]);
    plotHistPhase(posPhase(indPhase)-360,SSTHistPhaseFil(1:20,:),'Putative SST',pathAnal,['putativeSSTHistPhase' num2str(cluSST1)]);
    plotHistPhase(posPhase(indPhase)-360,SSTHistPhaseFilTag(1:5,:),'Tagged SST',pathAnal,['taggedSSTHistPhase'  num2str(cluSST1)]);
    
    ind = find(modInt1NoCue.timeStepRun == -1);
    indSSTTmp = find(autoCorrIntAll.idxC2 == cluSST1);
    for i = 1:size(avgFRProfile,1)
        avgFRProfileNorm(i,:) = avgFRProfile(i,:)/max(avgFRProfile(i,:));
    end
    plotAvgFRProfile(modInt1NoCue.timeStepRun(ind:end),...
        avgFRProfileNorm(indSSTTmp,ind:end),['C' num2str(cluSST1) ' FR (Hz)'],...
        ['Int_FRProfilePutativeSSTC' num2str(cluSST1)],...
        pathAnal,[0.4 0.9])
    plotAvgFRProfile(modInt1NoCue.timeStepRun(ind:end),...
        avgFRProfileNorm(indSST,ind:end),['Tagged C' num2str(cluSST1) ' FR (Hz)'],...
        ['Int_FRProfileTaggedSSTC' num2str(cluSST1)],...
        pathAnal,[0.4 1])
    plotAvgFRProfileCmp(modInt1NoCue.timeStepRun(ind:end),...
        avgFRProfileNorm(indSSTTmp,ind:end),avgFRProfileNorm(indSST,ind:end),...
        ['C' num2str(cluSST1) ' FR (Hz)'],...
        ['Int_FRProfilePutativeTaggedSSTC' num2str(cluSST1)],...
        pathAnal,[])
    
    %% SST
    cluSST2 = [3];
    peakTime = autoCorrIntAll.peakTime(autoCorrIntAll.idxC2 == cluSST2);
    [~,indPeakTime] = sort(peakTime);
    SSTCCGVal = autoCorrIntAll.ccgVal(autoCorrIntAll.idxC2 == cluSST2,:);
    SSTCCGVal = SSTCCGVal(indPeakTime,:);
    for i = 1:size(SSTCCGVal,1)
        SSTCCGVal(i,:) = SSTCCGVal(i,:)/max(SSTCCGVal(i,:));
    end
    SSTHistPhaseFil = autoCorrIntAll.histPhaseFil(autoCorrIntAll.idxC2 == cluSST2,indPhase);
    SSTHistPhaseFil = SSTHistPhaseFil(indPeakTime,:);
    for i = 1:size(SSTHistPhaseFil,1)
        SSTHistPhaseFil(i,:) = SSTHistPhaseFil(i,:)/max(SSTHistPhaseFil(i,:));
    end
    
    idxTagSST1 = autoCorrIntTag.cellType == 2 & autoCorrIntTag.idxC2 == cluSST2;
    peakTime = autoCorrIntTag.peakTime(idxTagSST1);
    indRec = autoCorrIntTag.indRec(idxTagSST1);
    indNeu = autoCorrIntTag.indNeu(idxTagSST1);
    indSST = [];
    for i = 1:length(indRec)
        ind = find(autoCorrIntAll.indNeu == indNeu(i) & autoCorrIntAll.indRec == indRec(i) ...
            & autoCorrIntAll.task == 2);
        indSST = [indSST ind];
    end
    [~,indPeakTime] = sort(peakTime);
    SSTCCGValTag = autoCorrIntTag.ccgVal(idxTagSST1,:);
    SSTCCGValTag = SSTCCGValTag(indPeakTime,:);
    for i = 1:size(SSTCCGValTag,1)
        SSTCCGValTag(i,:) = SSTCCGValTag(i,:)/max(SSTCCGValTag(i,:));
    end
    SSTHistPhaseFilTag = autoCorrIntAll.histPhaseFil(indSST,indPhase);
    SSTHistPhaseFilTag = SSTHistPhaseFilTag(indPeakTime,:);
    for i = 1:size(SSTHistPhaseFilTag,1)
        SSTHistPhaseFilTag(i,:) = SSTHistPhaseFilTag(i,:)/max(SSTHistPhaseFilTag(i,:));
    end
    
    p = randperm(size(SSTCCGVal,1));
    p = sort(p(1:15));
    plotCCG(-50:50,SSTCCGVal(:,:),'Putative SST',pathAnal,['putativeSSTCCG' num2str(cluSST2)])
    plotCCG(-50:50,SSTCCGValTag(:,:),'Tagged SST',pathAnal,['taggedSSTCCG' num2str(cluSST2)])
    plotHistPhase(posPhase(indPhase)-360,SSTHistPhaseFil(:,:),'Putative SST',pathAnal,['putativeSSTHistPhase' num2str(cluSST2)])
    plotHistPhase(posPhase(indPhase)-360,SSTHistPhaseFilTag(:,:),'Tagged SST',pathAnal,['taggedSSTHistPhase' num2str(cluSST2)])
    
    ind = find(modInt1NoCue.timeStepRun == -2);
    indSSTTmp = find(autoCorrIntAll.idxC2 == cluSST2);
    for i = 1:size(avgFRProfile,1)
        avgFRProfileNorm(i,:) = avgFRProfile(i,:)/max(avgFRProfile(i,:));
    end
    plotAvgFRProfile(modInt1NoCue.timeStepRun(ind:end),...
        avgFRProfileNorm(indSSTTmp(1:end),ind:end),['C' num2str(cluSST2) ' FR (Hz)'],...
        ['Int_FRProfilePutativeSSTC' num2str(cluSST2)],...
        pathAnal,[])
    plotAvgFRProfile(modInt1NoCue.timeStepRun(ind:end),...
        avgFRProfileNorm(indSST,ind:end),['Tagged C' num2str(cluSST2) ' FR (Hz)'],...
        ['Int_FRProfileTaggedSSTC' num2str(cluSST2)],...
        pathAnal,[])
    
        
    %% PV
    cluPV1 = [4];
    peakTime = autoCorrIntAll.peakTime(autoCorrIntAll.idxC2 == cluPV1);
    [~,indPeakTime] = sort(peakTime);
    PVCCGVal = autoCorrIntAll.ccgVal(autoCorrIntAll.idxC2 == cluPV1,:);
    PVCCGVal = PVCCGVal(indPeakTime,:);
    for i = 1:size(PVCCGVal,1)
        PVCCGVal(i,:) = PVCCGVal(i,:)/max(PVCCGVal(i,:));
    end
    PVHistPhaseFil = autoCorrIntAll.histPhaseFil(autoCorrIntAll.idxC2 == cluPV1,indPhase);
    PVHistPhaseFil = PVHistPhaseFil(indPeakTime,:);
    for i = 1:size(PVHistPhaseFil,1)
        PVHistPhaseFil(i,:) = PVHistPhaseFil(i,:)/max(PVHistPhaseFil(i,:));
    end
    
    idxTagPV1 = autoCorrIntTag.cellType == 1 & autoCorrIntTag.idxC2 == cluPV1;
    peakTime = autoCorrIntTag.peakTime(idxTagPV1);
    indRec = autoCorrIntTag.indRec(idxTagPV1);
    indNeu = autoCorrIntTag.indNeu(idxTagPV1);
    indPV = [];
    for i = 1:length(indRec)
        ind = find(autoCorrIntAll.indNeu == indNeu(i) & autoCorrIntAll.indRec == indRec(i) ...
            & autoCorrIntAll.task == 2);
        indPV = [indPV ind];
    end
    [~,indPeakTime] = sort(peakTime);
    PVCCGValTag = autoCorrIntTag.ccgVal(idxTagPV1,:);
    PVCCGValTag = PVCCGValTag(indPeakTime,:);
    for i = 1:size(PVCCGValTag,1)
        PVCCGValTag(i,:) = PVCCGValTag(i,:)/max(PVCCGValTag(i,:));
    end
    PVHistPhaseFilTag = autoCorrIntAll.histPhaseFil(indPV,indPhase);
    PVHistPhaseFilTag = PVHistPhaseFilTag(indPeakTime,:);
    for i = 1:size(PVHistPhaseFilTag,1)
        PVHistPhaseFilTag(i,:) = PVHistPhaseFilTag(i,:)/max(PVHistPhaseFilTag(i,:));
    end
    p = randperm(size(PVCCGVal,1));
    p = sort(p(1:15));
    plotCCG(-50:50,PVCCGVal(1:20,:),'Putative PV',pathAnal,['putativePVCCG' num2str(cluPV1)]);
    plotCCG(-50:50,PVCCGValTag(1:5,:),'Tagged PV',pathAnal,['taggedPVCCG' num2str(cluPV1)]);
    plotHistPhase(posPhase(indPhase)-360,PVHistPhaseFil(1:20,:),'Putative PV',pathAnal,['putativePVHistPhase' num2str(cluPV1)]);
    plotHistPhase(posPhase(indPhase)-360,PVHistPhaseFilTag(1:5,:),'Tagged PV',pathAnal,['taggedPVHistPhase' num2str(cluPV1)]);
    
    ind = find(modInt1NoCue.timeStepRun == -1);
    indPVTmp = find(autoCorrIntAll.idxC2 == cluPV1);
    plotAvgFRProfile(modInt1NoCue.timeStepRun(ind:end),...
        avgFRProfileNorm(indPVTmp(1:end),ind:end),['C' num2str(cluPV1) ' FR (Hz)'],...
        ['Int_FRProfilePutativePVC' num2str(cluPV1)],...
        pathAnal,[0.4 0.9])
    plotAvgFRProfile(modInt1NoCue.timeStepRun(ind:end),...
        avgFRProfileNorm(indPV,ind:end),['Tagged C' num2str(cluPV1) ' FR (Hz)'],...
        ['Int_FRProfileTaggedPVC' num2str(cluPV1)],...
        pathAnal,[0.4 1])
    plotAvgFRProfileCmp(modInt1NoCue.timeStepRun(ind:end),...
        avgFRProfileNorm(indPVTmp,ind:end),avgFRProfileNorm(indPV,ind:end),...
        ['C' num2str(cluPV1) ' FR (Hz)'],...
        ['Int_FRProfilePutativeTaggedPVC' num2str(cluPV1)],...
        pathAnal,[])
    
    %% PV
    cluPV2 = [6];
    peakTime = autoCorrIntAll.peakTime(autoCorrIntAll.idxC2 == cluPV2);
    [~,indPeakTime] = sort(peakTime);
    PVCCGVal = autoCorrIntAll.ccgVal(autoCorrIntAll.idxC2 == cluPV2,:);
    PVCCGVal = PVCCGVal(indPeakTime,:);
    for i = 1:size(PVCCGVal,1)
        PVCCGVal(i,:) = PVCCGVal(i,:)/max(PVCCGVal(i,:));
    end
    PVHistPhaseFil = autoCorrIntAll.histPhaseFil(autoCorrIntAll.idxC2 == cluPV2,indPhase);
    PVHistPhaseFil = PVHistPhaseFil(indPeakTime,:);
    for i = 1:size(PVHistPhaseFil,1)
        PVHistPhaseFil(i,:) = PVHistPhaseFil(i,:)/max(PVHistPhaseFil(i,:));
    end
    
    idxTagPV1 = autoCorrIntTag.cellType == 1 & autoCorrIntTag.idxC2 == cluPV2;
    peakTime = autoCorrIntTag.peakTime(idxTagPV1);
    indRec = autoCorrIntTag.indRec(idxTagPV1);
    indNeu = autoCorrIntTag.indNeu(idxTagPV1);
    indPV = [];
    for i = 1:length(indRec)
        ind = find(autoCorrIntAll.indNeu == indNeu(i) & autoCorrIntAll.indRec == indRec(i) ...
            & autoCorrIntAll.task == 2);
        indPV = [indPV ind];
    end
    [~,indPeakTime] = sort(peakTime);
    PVCCGValTag = autoCorrIntTag.ccgVal(idxTagPV1,:);
    PVCCGValTag = PVCCGValTag(indPeakTime,:);
    for i = 1:size(PVCCGValTag,1)
        PVCCGValTag(i,:) = PVCCGValTag(i,:)/max(PVCCGValTag(i,:));
    end
    PVHistPhaseFilTag = autoCorrIntAll.histPhaseFil(indPV,indPhase);
    PVHistPhaseFilTag = PVHistPhaseFilTag(indPeakTime,:);
    for i = 1:size(PVHistPhaseFilTag,1)
        PVHistPhaseFilTag(i,:) = PVHistPhaseFilTag(i,:)/max(PVHistPhaseFilTag(i,:));
    end
    p = randperm(size(PVCCGVal,1));
    p = sort(p(1:15));
    plotCCG(-50:50,PVCCGVal(:,:),'Putative PV',pathAnal,['putativePVCCG' num2str(cluPV2)]);
    plotCCG(-50:50,PVCCGValTag(:,:),'Tagged PV',pathAnal,['taggedPVCCG' num2str(cluPV2)]);
    plotHistPhase(posPhase(indPhase)-360,PVHistPhaseFil(:,:),'Putative PV',pathAnal,['putativePVHistPhase' num2str(cluPV2)]);
    plotHistPhase(posPhase(indPhase)-360,PVHistPhaseFilTag(:,:),'Tagged PV',pathAnal,['taggedPVHistPhase' num2str(cluPV2)]);
    
    ind = find(modInt1NoCue.timeStepRun == -2);
    indPVTmp = find(autoCorrIntAll.idxC2 == cluPV2);
    plotAvgFRProfile(modInt1NoCue.timeStepRun(ind:end),...
        avgFRProfileNorm(indPVTmp(1:end),ind:end),['C' num2str(cluPV2) ' FR (Hz)'],...
        ['Int_FRProfilePutativePVC' num2str(cluPV2)],...
        pathAnal,[])
    plotAvgFRProfile(modInt1NoCue.timeStepRun(ind:end),...
        avgFRProfileNorm(indPV,ind:end),['Tagged PV C' num2str(cluPV2) ' FR (Hz)'],...
        ['Int_FRProfileTaggedPVC' num2str(cluPV2)],...
        pathAnal,[])
    
    idxTagSST1 = autoCorrIntTag.cellType == 2 & autoCorrIntTag.idxC2 == cluPV2;
    peakTime = autoCorrIntTag.peakTime(idxTagSST1);
    indRec = autoCorrIntTag.indRec(idxTagSST1);
    indNeu = autoCorrIntTag.indNeu(idxTagSST1);
    indSST = [];
    for i = 1:length(indRec)
        ind = find(autoCorrIntAll.indNeu == indNeu(i) & autoCorrIntAll.indRec == indRec(i) ...
            & autoCorrIntAll.task == 2);
        indSST = [indSST ind];
    end
    ind = find(modInt1NoCue.timeStepRun == -2);
    plotAvgFRProfile(modInt1NoCue.timeStepRun(ind:end),...
        avgFRProfileNorm(indSST,ind:end),['Tagged SST C' num2str(cluPV2) ' FR (Hz)'],...
        ['Int_FRProfileTaggedSSTC' num2str(cluPV2)],...
        pathAnal,[])
    
    %% Unknown
    cluPV2 = [3];
    peakTime = autoCorrIntAll.peakTime(autoCorrIntAll.idxC2 == cluPV2);
    [~,indPeakTime] = sort(peakTime);
    PVCCGVal = autoCorrIntAll.ccgVal(autoCorrIntAll.idxC2 == cluPV2,:);
    PVCCGVal = PVCCGVal(indPeakTime,:);
    for i = 1:size(PVCCGVal,1)
        PVCCGVal(i,:) = PVCCGVal(i,:)/max(PVCCGVal(i,:));
    end
    PVHistPhaseFil = autoCorrIntAll.histPhaseFil(autoCorrIntAll.idxC2 == cluPV2,indPhase);
    PVHistPhaseFil = PVHistPhaseFil(indPeakTime,:);
    for i = 1:size(PVHistPhaseFil,1)
        PVHistPhaseFil(i,:) = PVHistPhaseFil(i,:)/max(PVHistPhaseFil(i,:));
    end
    
    idxTagPV1 = autoCorrIntTag.cellType == 1 & autoCorrIntTag.idxC2 == cluPV2;
    peakTime = autoCorrIntTag.peakTime(idxTagPV1);
    indRec = autoCorrIntTag.indRec(idxTagPV1);
    indNeu = autoCorrIntTag.indNeu(idxTagPV1);
    indPV = [];
    for i = 1:length(indRec)
        ind = find(autoCorrIntAll.indNeu == indNeu(i) & autoCorrIntAll.indRec == indRec(i) ...
            & autoCorrIntAll.task == 2);
        indPV = [indPV ind];
    end
    [~,indPeakTime] = sort(peakTime);
    PVCCGValTag = autoCorrIntTag.ccgVal(idxTagPV1,:);
    PVCCGValTag = PVCCGValTag(indPeakTime,:);
    for i = 1:size(PVCCGValTag,1)
        PVCCGValTag(i,:) = PVCCGValTag(i,:)/max(PVCCGValTag(i,:));
    end
    PVHistPhaseFilTag = autoCorrIntAll.histPhaseFil(indPV,indPhase);
    PVHistPhaseFilTag = PVHistPhaseFilTag(indPeakTime,:);
    for i = 1:size(PVHistPhaseFilTag,1)
        PVHistPhaseFilTag(i,:) = PVHistPhaseFilTag(i,:)/max(PVHistPhaseFilTag(i,:));
    end
    p = randperm(size(PVCCGVal,1));
    p = sort(p(1:15));
    plotCCG(-50:50,PVCCGVal(:,:),'Putative PV',pathAnal,['putativePVCCG' num2str(cluPV2)]);
    plotCCG(-50:50,PVCCGValTag(:,:),'Tagged PV',pathAnal,['taggedPVCCG' num2str(cluPV2)]);
    plotHistPhase(posPhase(indPhase)-360,PVHistPhaseFil(:,:),'Putative PV',pathAnal,['putativePVHistPhase' num2str(cluPV2)]);
    plotHistPhase(posPhase(indPhase)-360,PVHistPhaseFilTag(:,:),'Tagged PV',pathAnal,['taggedPVHistPhase' num2str(cluPV2)]);
    
    ind = find(modInt1NoCue.timeStepRun == -2);
    indPVTmp = find(autoCorrIntAll.idxC2 == cluPV2);
    plotAvgFRProfile(modInt1NoCue.timeStepRun(ind:end),...
        avgFRProfileNorm(indPVTmp(1:end),ind:end),['C' num2str(cluPV2) ' FR (Hz)'],...
        ['Int_FRProfilePutativeC' num2str(cluPV2)],...
        pathAnal,[])

    
    %% Unknown
    cluPV2 = [5];
    peakTime = autoCorrIntAll.peakTime(autoCorrIntAll.idxC2 == cluPV2);
    [~,indPeakTime] = sort(peakTime);
    PVCCGVal = autoCorrIntAll.ccgVal(autoCorrIntAll.idxC2 == cluPV2,:);
    PVCCGVal = PVCCGVal(indPeakTime,:);
    for i = 1:size(PVCCGVal,1)
        PVCCGVal(i,:) = PVCCGVal(i,:)/max(PVCCGVal(i,:));
    end
    PVHistPhaseFil = autoCorrIntAll.histPhaseFil(autoCorrIntAll.idxC2 == cluPV2,indPhase);
    PVHistPhaseFil = PVHistPhaseFil(indPeakTime,:);
    for i = 1:size(PVHistPhaseFil,1)
        PVHistPhaseFil(i,:) = PVHistPhaseFil(i,:)/max(PVHistPhaseFil(i,:));
    end
    
    idxTagPV1 = autoCorrIntTag.cellType == 1 & autoCorrIntTag.idxC2 == cluPV2;
    peakTime = autoCorrIntTag.peakTime(idxTagPV1);
    indRec = autoCorrIntTag.indRec(idxTagPV1);
    indNeu = autoCorrIntTag.indNeu(idxTagPV1);
    indPV = [];
    for i = 1:length(indRec)
        ind = find(autoCorrIntAll.indNeu == indNeu(i) & autoCorrIntAll.indRec == indRec(i) ...
            & autoCorrIntAll.task == 2);
        indPV = [indPV ind];
    end
    [~,indPeakTime] = sort(peakTime);
    PVCCGValTag = autoCorrIntTag.ccgVal(idxTagPV1,:);
    PVCCGValTag = PVCCGValTag(indPeakTime,:);
    for i = 1:size(PVCCGValTag,1)
        PVCCGValTag(i,:) = PVCCGValTag(i,:)/max(PVCCGValTag(i,:));
    end
    PVHistPhaseFilTag = autoCorrIntAll.histPhaseFil(indPV,indPhase);
    PVHistPhaseFilTag = PVHistPhaseFilTag(indPeakTime,:);
    for i = 1:size(PVHistPhaseFilTag,1)
        PVHistPhaseFilTag(i,:) = PVHistPhaseFilTag(i,:)/max(PVHistPhaseFilTag(i,:));
    end
    p = randperm(size(PVCCGVal,1));
    p = sort(p(1:15));
    plotCCG(-50:50,PVCCGVal(:,:),'Putative PV',pathAnal,['putativePVCCG' num2str(cluPV2)]);
    plotCCG(-50:50,PVCCGValTag(:,:),'Tagged PV',pathAnal,['taggedPVCCG' num2str(cluPV2)]);
    plotHistPhase(posPhase(indPhase)-360,PVHistPhaseFil(:,:),'Putative PV',pathAnal,['putativePVHistPhase' num2str(cluPV2)]);
    plotHistPhase(posPhase(indPhase)-360,PVHistPhaseFilTag(:,:),'Tagged PV',pathAnal,['taggedPVHistPhase' num2str(cluPV2)]);
    
    ind = find(modInt1NoCue.timeStepRun == -2);
    indPVTmp = find(autoCorrIntAll.idxC2 == cluPV2);
    plotAvgFRProfile(modInt1NoCue.timeStepRun(ind:end),...
        avgFRProfileNorm(indPVTmp(1:end),ind:end),['C' num2str(cluPV2) ' FR (Hz)'],...
        ['Int_FRProfilePutativeC' num2str(cluPV2)],...
        pathAnal,[])
    
end

function plotCCG(time,ccg,ylab,path,fileName)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
        
    imagesc(time,1:size(ccg,1),ccg);
%     colormap jet
    xlabel('Time (ms)')
    ylabel(ylab)
    
    fileName1 = [path fileName];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotHistPhase(phase,ccg,ylab,path,fileName)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
        
    imagesc(phase,1:size(ccg,1),ccg);
%     colormap jet
    set(gca,'XTick',[0 180 360 540 720]);
    xlabel('Theta phase (deg.)')
    ylabel(ylab)
    
    fileName1 = [path fileName];
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
    h = plot([0 0],[0 1],'r-');
    set(h,'LineWidth',1)
    set(gca,'XLim',[timeStepRun(1) 3],'XTick',timeStepRun(1):1:4);
    if(~isempty(ylimit))
        set(gca,'YLim',ylimit);
    else
        set(gca,'YLim',[0 1]);
%         set(gca,'YLim',[min(mean(avgFRProfile)-std(avgFRProfile)/sqrt(size(avgFRProfile,1)))*0.95 ...
%         max(mean(avgFRProfile)+std(avgFRProfile)/sqrt(size(avgFRProfile,1)))*1.05]);
    end
    xlabel('Time (s)')
    ylabel(yl)
    
    fileName1 = [pathAnal fileName];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
        
end

function plotAvgFRProfileCmp(timeStepRun,avgFRProfilex,avgFRProfiley,yl,fileName,pathAnal,ylimit)
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
    h = plot([0 0],[0 1],'r-');
    set(h,'LineWidth',1)
    set(gca,'XLim',[timeStepRun(1) 3],'XTick',timeStepRun(1):1:4,'YLim',[0 1]);
    if(~isempty(ylimit))
        set(gca,'YLim',ylimit);
    else
        set(gca,'YLim',[0.4 1]);
%         set(gca,'YLim',[min([minX minY])*0.95 ...
%         max([maxX maxY])*1.05]);
    end
    xlabel('Time (s)')
    ylabel(yl)
    
    fileName1 = [pathAnal fileName];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
        
end
