function InterneuronInitPeakAllRec(taskSel)
    
    methodTheta = 1;
    minFRInt = 3;
    
    RecordingList;
    pathAnal0 = 'Z:\Yingxue\DataAnalysisRaphi\interneuron\';
    if(taskSel == 1) % including all the neurons
        pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\interneuron\';
    elseif(taskSel == 2) % including AL and PL neurons
        pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\interneuronALPL\';
    else
        pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\interneuronAL\';
    end
    
    sampleFq = 1250;
    
    load([pathAnal0 'autoCorrIntAllRec.mat']);
    if(exist([pathAnal0 'initPeakIntAllRec.mat']))
        load([pathAnal0 'initPeakIntAllRec.mat']);
    end
    
%     nNeuWithField = [modIntNoCue.nNeuWithField modIntAL.nNeuWithField modIntPL.nNeuWithField];
%     nNeuWithFieldAligned = [modIntNoCue.nNeuWithFieldAligned modIntAL.nNeuWithFieldAligned ...
%             modIntPL.nNeuWithFieldAligned];
%         
%     disp('No cue - peak firing rate')
%     modInt1NoCue = accumInterneurons2(listRecordingsNoCuePath,...
%         listRecordingsNoCueFileName,mazeSessionNoCue,autoCorrIntAll,...
%         nNeuWithField,nNeuWithFieldAligned,minFRInt,1,sampleFq);
%     
%     disp('Active licking - peak firing rate')
%     modInt1AL = accumInterneurons2(listRecordingsActiveLickPath,...
%         listRecordingsActiveLickFileName,mazeSessionActiveLick,autoCorrIntAll,...
%         nNeuWithField,nNeuWithFieldAligned,minFRInt,2,sampleFq);
%     
%     disp('Passive licking - peak firing rate')
%     modInt1PL = accumInterneurons2(listRecordingsPassiveLickPath,...
%         listRecordingsPassiveLickFileName,mazeSessionPassiveLick,autoCorrIntAll,...
%         nNeuWithField,nNeuWithFieldAligned,minFRInt,3,sampleFq);
%     
%     save([pathAnal 'initPeakIntAllRec.mat'],'modInt1NoCue','modInt1AL','modInt1PL');
%     
    if(taskSel == 1) 
        task = [modInt1NoCue.task modInt1AL.task modInt1PL.task];
        taskBad = [modInt1NoCue.taskBad modInt1AL.taskBad modInt1PL.taskBad];
        indRec = [modInt1NoCue.indRec modInt1AL.indRec modInt1PL.indRec];
        indRecBad = [modInt1NoCue.indRecBad modInt1AL.indRecBad modInt1PL.indRecBad];
        indNeu = [modInt1NoCue.indNeu modInt1AL.indNeu modInt1PL.indNeu];
        indNeuBad = [modInt1NoCue.indNeuBad modInt1AL.indNeuBad modInt1PL.indNeuBad];
        idxC = [modInt1NoCue.idxC2 modInt1AL.idxC2 modInt1PL.idxC2];
        idxCBad = [modInt1NoCue.idxC2Bad modInt1AL.idxC2Bad modInt1PL.idxC2Bad];
        nNeuWithField = [modInt1NoCue.nNeuWithField modInt1AL.nNeuWithField modInt1PL.nNeuWithField];
        nNeuWithFieldAligned = [modInt1NoCue.nNeuWithFieldAligned modInt1AL.nNeuWithFieldAligned modInt1PL.nNeuWithFieldAligned];
        avgFRProfile = [modInt1NoCue.avgFRProfile; modInt1AL.avgFRProfile; modInt1PL.avgFRProfile];
        avgFRProfileBad = [modInt1NoCue.avgFRProfileBad; modInt1AL.avgFRProfileBad; modInt1PL.avgFRProfileBad];
    elseif(taskSel == 2)
        task = [modInt1AL.task modInt1PL.task];
        taskBad = [modInt1AL.taskBad modInt1PL.taskBad];
        indRec = [modInt1AL.indRec modInt1PL.indRec];
        indRecBad = [modInt1AL.indRecBad modInt1PL.indRecBad];
        indNeu = [modInt1AL.indNeu modInt1PL.indNeu];
        indNeuBad = [modInt1AL.indNeuBad modInt1PL.indNeuBad];
        idxC = [modInt1AL.idxC2 modInt1PL.idxC2];
        idxCBad = [modInt1AL.idxC2Bad modInt1PL.idxC2Bad];
        nNeuWithField = [modInt1AL.nNeuWithField modInt1PL.nNeuWithField];
        nNeuWithFieldAligned = [modInt1AL.nNeuWithFieldAligned modInt1PL.nNeuWithFieldAligned];
        avgFRProfile =[modInt1AL.avgFRProfile;  modInt1PL.avgFRProfile];
        avgFRProfileBad = [modInt1AL.avgFRProfileBad; modInt1PL.avgFRProfileBad];
    else
        task = modInt1AL.task;
        taskBad = modInt1AL.taskBad;
        indRec = modInt1AL.indRec;
        indRecBad = modInt1AL.indRecBad;
        indNeu = modInt1AL.indNeu;
        indNeuBad = modInt1AL.indNeuBad;
        idxC = modInt1AL.idxC2;
        idxCBad = modInt1AL.idxC2Bad;
        nNeuWithField = modInt1AL.nNeuWithField;
        nNeuWithFieldAligned = modInt1AL.nNeuWithFieldAligned;
        avgFRProfile = modInt1AL.avgFRProfile;
        avgFRProfileBad = modInt1AL.avgFRProfileBad;
    end
    avgFRProfileNorm = zeros(size(avgFRProfile,1),size(avgFRProfile,2));
    avgFRProfileNormBad = zeros(size(avgFRProfileBad,1),size(avgFRProfileBad,2));
    for i = 1:size(avgFRProfile,1)
        if(max(avgFRProfile(i,:)) ~= 0)
            avgFRProfileNorm(i,:) = avgFRProfile(i,:)/max(avgFRProfile(i,:));
        end
    end
    for i = 1:size(avgFRProfileBad,1)
        if(max(avgFRProfileBad(i,:)) ~= 0)
            avgFRProfileNormBad(i,:) = avgFRProfileBad(i,:)/max(avgFRProfileBad(i,:));
        end
    end

    FRProfileMean = accumMean(avgFRProfile,modInt1AL.timeStepRun);
    
    FRProfileMeanBad = accumMean(avgFRProfileBad,modInt1AL.timeStepRun);
    
    FRProfileMeanStatC = accumMeanStatC(FRProfileMean,idxC,nNeuWithField);
    
    FRProfileMeanStatCGoodBad = accumMeanStatCGoodBad(FRProfileMean,FRProfileMeanBad,idxC,idxCBad);
    
    if(taskSel == 1)
        save([pathAnal 'initPeakIntAllRec.mat'],'FRProfileMean','FRProfileMeanBad',...
            'FRProfileMeanStatC','FRProfileMeanStatCGoodBad','-append'); 
    elseif(taskSel == 2)
        save([pathAnal 'initPeakIntALPLRec.mat'],'FRProfileMean','FRProfileMeanBad',...
            'FRProfileMeanStatC','FRProfileMeanStatCGoodBad'); 
    else
        save([pathAnal 'initPeakIntALRec.mat'],'FRProfileMean','FRProfileMeanBad',...
            'FRProfileMeanStatC','FRProfileMeanStatCGoodBad'); 
    end
    
    colorSel = 0;
    
    ind = find(modInt1NoCue.timeStepRun == -1.5);
    
    % firing profile field vs. non-field
    for i = 1:max(idxC)
        plotAvgFRProfile(modInt1NoCue.timeStepRun(ind:end),...
            avgFRProfile(idxC == i,ind:end),['C' num2str(i) ' FR (Hz)'],...
            ['Int_FRProfileSecC' num2str(i)],...
            pathAnal,[])
        
        indCurCField = idxC == i & nNeuWithField > 1;
        indCurCNoField = idxC == i & nNeuWithField < 1;
        plotAvgFRProfileCmp(modInt1NoCue.timeStepRun(ind:end),...
            avgFRProfile(indCurCField,ind:end),avgFRProfile(indCurCNoField,ind:end),...
            ['C' num2str(i) ' FR (Hz) F vs. NoF'],...
            ['Int_FRProfileCFieldNoField' num2str(i)],...
            pathAnal,[])
        
        plotAvgFRProfileCmp(modInt1NoCue.timeStepRun(ind:end),...
            avgFRProfileNorm(indCurCField,ind:end),avgFRProfileNorm(indCurCNoField,ind:end),...
            ['C' num2str(i) ' Norm FR (Hz) F vs. NoF'],...
            ['Int_FRProfileNormCFieldNoField' num2str(i)],...
            pathAnal,[])
    end
    
    % firing profile field vs. non-field after aligning to run onset
    for i = 1:max(idxC)        
        indCurCField = idxC == i & nNeuWithFieldAligned > 1;
        indCurCNoField = idxC == i & nNeuWithFieldAligned < 1;
        plotAvgFRProfileCmp(modInt1NoCue.timeStepRun(ind:end),...
            avgFRProfile(indCurCField,ind:end),avgFRProfile(indCurCNoField,ind:end),...
            ['C' num2str(i) ' FR (Hz) F vs. NoF'],...
            ['Int_FRProfileAlignedCFieldNoField' num2str(i)],...
            pathAnal,[])
        
        plotAvgFRProfileCmp(modInt1NoCue.timeStepRun(ind:end),...
            avgFRProfileNorm(indCurCField,ind:end),avgFRProfileNorm(indCurCNoField,ind:end),...
            ['C' num2str(i) ' Norm FR (Hz) F vs. NoF'],...
            ['Int_FRProfileAlignedNormCFieldNoField' num2str(i)],...
            pathAnal,[])
    end
        
    % firing profile good vs. bad trials
    for i = 1:max(idxC)        
        plotAvgFRProfileCmp(modInt1NoCue.timeStepRun(ind:end),...
            avgFRProfile(idxC == i,ind:end),avgFRProfileBad(idxCBad == i,ind:end),...
            ['C' num2str(i) ' FR (Hz) Good vs. Bad'],...
            ['Int_FRProfileCGoodVsBad' num2str(i)],...
            pathAnal,[])
        
        plotAvgFRProfileCmp(modInt1NoCue.timeStepRun(ind:end),...
            avgFRProfileNorm(idxC == i,ind:end),avgFRProfileNormBad(idxCBad == i,ind:end),...
            ['C' num2str(i) ' Norm FR Good vs. Bad'],...
            ['Int_FRProfileNormCGoodVsBad' num2str(i)],...
            pathAnal,[])
    end
    
    % firing profile good vs. bad trials for tagged cells
    for i = 1:max(idxC)        
        idxTmp = find(autoCorrIntTag.idxC2 == i & autoCorrIntTag.cellType == 1);
        taskTag = autoCorrIntTag.task(idxTmp);
        recTag = autoCorrIntTag.indRec(idxTmp);
        neuTag = autoCorrIntTag.indNeu(idxTmp);
        idxTag = [];        
        for j = 1:length(idxTmp)
            idxTagTmp = find(task == taskTag(j) & indRec == recTag(j) & indNeu == neuTag(j));
            idxTag = [idxTag idxTagTmp];
        end
        idxTagBad = []; 
        for j = 1:length(idxTmp)
            idxTagTmp = find(taskBad == taskTag(j) & indRecBad == recTag(j) & indNeuBad == neuTag(j));
            idxTagBad = [idxTagBad idxTagTmp];
        end
        plotAvgFRProfileCmp(modInt1NoCue.timeStepRun(ind:end),...
            avgFRProfile(idxTag,ind:end),...
            avgFRProfileBad(idxTagBad,ind:end),...
            ['C' num2str(i) ' FR (Hz) Good vs. Bad Tagged PV'],...
            ['Int_FRProfileCGoodVsBadTaggedPV' num2str(i)],...
            pathAnal,[])
        
        plotAvgFRProfileCmp(modInt1NoCue.timeStepRun(ind:end),...
            avgFRProfileNorm(idxTag,ind:end),avgFRProfileNormBad(idxTagBad,ind:end),...
            ['C' num2str(i) ' Norm FR Good vs. Bad Tagged PV'],...
            ['Int_FRProfileNormCGoodVsBadTaggedPV' num2str(i)],...
            pathAnal,[])
        
        idxTmp = find(autoCorrIntTag.idxC2 == i & autoCorrIntTag.cellType == 2);
        taskTag = autoCorrIntTag.task(idxTmp);
        recTag = autoCorrIntTag.indRec(idxTmp);
        neuTag = autoCorrIntTag.indNeu(idxTmp);
        idxTag = [];        
        for j = 1:length(idxTmp)
            idxTagTmp = find(task == taskTag(j) & indRec == recTag(j) & indNeu == neuTag(j));
            idxTag = [idxTag idxTagTmp];
        end
        idxTagBad = []; 
        for j = 1:length(idxTmp)
            idxTagTmp = find(taskBad == taskTag(j) & indRecBad == recTag(j) & indNeuBad == neuTag(j));
            idxTagBad = [idxTagBad idxTagTmp];
        end
        plotAvgFRProfileCmp(modInt1NoCue.timeStepRun(ind:end),...
            avgFRProfile(idxTag,ind:end),...
            avgFRProfileBad(idxTagBad,ind:end),...
            ['C' num2str(i) ' FR (Hz) Good vs. Bad Tagged SST'],...
            ['Int_FRProfileCGoodVsBadTaggedSST' num2str(i)],...
            pathAnal,[])
        
        plotAvgFRProfileCmp(modInt1NoCue.timeStepRun(ind:end),...
            avgFRProfileNorm(idxTag,ind:end),avgFRProfileNormBad(idxTagBad,ind:end),...
            ['C' num2str(i) ' Norm FR Good vs. Bad Tagged SST'],...
            ['Int_FRProfileNormCGoodVsBadTaggedSST' num2str(i)],...
            pathAnal,[])
    end
    
    % firing profile of individual neurons ordered
    for i = 1:max(idxC)
        plotIndFRProfile(modInt1NoCue.timeStepRun(ind:end),...
            avgFRProfileNorm(idxC == i,ind:end),['C' num2str(i) ' FR (Hz)'],...
            ['Int_IndFRProfileNormCPeak' num2str(i)],...
            pathAnal,[],1,[],[]) % ordered based on peak time
        
%         plotIndFRProfile(modInt1NoCue.timeStepRun,...
%             avgFRProfileNorm(idxC == i,:),['C' num2str(i) ' FR (Hz)'],...
%             ['Int_IndFRProfileNormCFR0to1' num2str(i)],...
%             pathAnal,[],2,FRProfileMean.indFR0to1,[]) % ordered based on 0to1 mean
%         
%         plotIndFRProfile(modInt1NoCue.timeStepRun,...
%             avgFRProfileNorm(idxC == i,:),['C' num2str(i) ' FR (Hz)'],...
%             ['Int_IndFRProfileNormCFRBefRun' num2str(i)],...
%             pathAnal,[],3,FRProfileMean.indFRBefRun,[]) % ordered based on -1to0 mean
%         
        plotIndFRProfile(modInt1NoCue.timeStepRun,...
            avgFRProfileNorm(idxC == i,:),['C' num2str(i) ' FR (Hz)'],...
            ['Int_IndFRProfileNormCFR0to1VsBefRun' num2str(i)],...
            pathAnal,[],4,FRProfileMean.indFRBefRun,FRProfileMean.indFR0to1) % ordered based on -1to0 to 0to1 mean ratio
        
        idxTmp = find(autoCorrIntTag.idxC2 == i & autoCorrIntTag.cellType == 1);
        taskTag = autoCorrIntTag.task(idxTmp);
        recTag = autoCorrIntTag.indRec(idxTmp);
        neuTag = autoCorrIntTag.indNeu(idxTmp);
        idxTag = [];        
        for j = 1:length(idxTmp)
            idxTagTmp = find(task == taskTag(j) & indRec == recTag(j) & indNeu == neuTag(j));
            idxTag = [idxTag idxTagTmp];
        end
        disp(['No tagged PV neurons in C' num2str(i) ' is: ' num2str(length(idxTmp)) ...
            '; No tagged neurons included here is :' num2str(length(idxTag))]);
        plotIndFRProfile(modInt1NoCue.timeStepRun,...
            avgFRProfileNorm(idxTag,:),['C' num2str(i) ' FR (Hz) Tagged PV'],...
            ['Int_IndFRProfileNormCFR0to1VsBefRunTagPV' num2str(i)],...
            pathAnal,[],4,FRProfileMean.indFRBefRun,FRProfileMean.indFR0to1) % tagged cells ordered based on -1to0 to 0to1 mean ratio
         
        idxTmp = find(autoCorrIntTag.idxC2 == i & autoCorrIntTag.cellType == 2);
        taskTag = autoCorrIntTag.task(idxTmp);
        recTag = autoCorrIntTag.indRec(idxTmp);
        neuTag = autoCorrIntTag.indNeu(idxTmp);
        idxTag = [];
        for j = 1:length(idxTmp)
            idxTagTmp = find(task == taskTag(j) & indRec == recTag(j) & indNeu == neuTag(j));
            idxTag = [idxTag idxTagTmp];
        end
        disp(['No tagged SST neurons in C' num2str(i) ' is: ' num2str(length(idxTmp)) ...
            '; No tagged neurons included here is :' num2str(length(idxTag))]);
        plotIndFRProfile(modInt1NoCue.timeStepRun,...
            avgFRProfileNorm(idxTag,:),['C' num2str(i) ' FR (Hz)  Tagged SST'],...
            ['Int_IndFRProfileNormCFR0to1VsBefRunTagSST' num2str(i)],...
            pathAnal,[],4,FRProfileMean.indFRBefRun,FRProfileMean.indFR0to1) % tagged cells ordered based on -1to0 to 0to1 mean ratio
        
%         plotAvgFRProfile(modInt1NoCue.timeStepRun,...
%             avgFRProfile(idxC == i,:),['C' num2str(i) ' FR (Hz)'],...
%             ['Int_FRProfileC' num2str(i)],...
%             pathAnal,[])
        
    end
    
    %% firing profile of individual neurons with field and without field ordered
    for i = 1:max(idxC)
        indCurCField = idxC == i & nNeuWithField > 1;
        plotIndFRProfile(modInt1NoCue.timeStepRun,...
            avgFRProfileNorm(indCurCField,:),['C' num2str(i) ' FR (Hz) F'],...
            ['Int_IndFRProfileNormCField0to1VsBefRun' num2str(i)],...
            pathAnal,[],4,FRProfileMean.indFRBefRun,FRProfileMean.indFR0to1) % ordered based on peak time
        
        indCurCNoField = idxC == i & nNeuWithField < 1;
        plotIndFRProfile(modInt1NoCue.timeStepRun,...
            avgFRProfileNorm(indCurCNoField,:),['C' num2str(i) ' FR (Hz) NoF'],...
            ['Int_IndFRProfileNormCNoField0to1VsBefRun' num2str(i)],...
            pathAnal,[],4,FRProfileMean.indFRBefRun,FRProfileMean.indFR0to1) % ordered based on peak time
        
    end
    
    % compare the mean of each segment of avgFRProfile
    for i = 1:max(idxC)        
        plotBoxPlot(FRProfileMean.meanAvgFRProfileBaseline(idxC == i),...
            FRProfileMean.meanAvgFRProfile0to1(idxC == i),...
            ['C' num2str(i) ' average FR (Hz) BL vs 0to1s'],...
            ['Int_FRMeanC' num2str(i) 'BL-0to1Box'],...
            pathAnal,[],FRProfileMeanStatC.pRS0to1VsBL(i),colorSel);
        
        plotBoxPlot(FRProfileMean.meanAvgFRProfileBaseline(idxC == i),...
            FRProfileMean.meanAvgFRProfileBefRun(idxC == i),...
            ['C' num2str(i) ' average FR (Hz) BL vs BefRun'],...
            ['Int_FRMeanC' num2str(i) 'BL-BefRunBox'],...
            pathAnal,[],FRProfileMeanStatC.pRSBefRunVsBL(i),colorSel);
%         
%         plotBoxPlot(FRProfileMean.meanAvgFRProfileBaseline(idxC == i),...
%             FRProfileMean.meanAvgFRProfile1to5(idxC == i),...
%             ['C' num2str(i) ' average FR (Hz) BL vs 1to5s'],...
%             ['Int_FRMeanC' num2str(i) 'BL-1to5Box'],...
%             pathAnal,[],FRProfileMeanStatC.pRS1to5VsBL(i),colorSel);
%         
%         plotBoxPlot(FRProfileMean.meanAvgFRProfileBefRun(idxC == i),...
%             FRProfileMean.meanAvgFRProfile1to5(idxC == i),...
%             ['C' num2str(i) ' average FR (Hz) BefRun vs 1to5s'],...
%             ['Int_FRMeanC' num2str(i) 'BefRun-1to5Box'],...
%             pathAnal,[],FRProfileMeanStatC.pRS1to5VsBefRun(i),colorSel);
%         
        plotBoxPlot(FRProfileMean.meanAvgFRProfileBefRun(idxC == i),...
            FRProfileMean.meanAvgFRProfile0to1(idxC == i),...
            ['C' num2str(i) ' average FR (Hz) BefRun vs 0to1s'],...
            ['Int_FRMeanC' num2str(i) 'BefRun-0to1Box'],...
            pathAnal,[],FRProfileMeanStatC.pRS0to1VsBefRun(i),colorSel);
%         
%         plotBoxPlot(FRProfileMean.meanAvgFRProfile0to1(idxC == i),...
%             FRProfileMean.meanAvgFRProfile1to5(idxC == i),...
%             ['C' num2str(i) ' average FR (Hz) 0to1s vs 1to5s'],...
%             ['Int_FRMeanC' num2str(i) '0to1-1to5Box'],...
%             pathAnal,[],FRProfileMeanStatC.pRS1to5Vs0to1(i),colorSel);
%         
    end
    
    % compare each segment of avgFRProfile
    for i = 1:max(idxC)    
        FRBLArr = avgFRProfile(idxC == i,FRProfileMean.indFRBaseline);
        FR0to1Arr = avgFRProfile(idxC == i,FRProfileMean.indFR0to1);
        plotBoxPlot(FR0to1Arr(:),FRBLArr(:),...
            ['C' num2str(i) ' average FR (Hz) All BL vs 0to1s All'],...
            ['Int_FRMeanC' num2str(i) 'BL-0to1BoxAll'],...
            pathAnal,[],FRProfileMeanStatC.pRS0to1VsBL(i),colorSel);
        
        FRBefRunArr = avgFRProfile(idxC == i,FRProfileMean.indFRBefRun);
        plotBoxPlot(FRBLArr(:),FRBefRunArr(:),...
            ['C' num2str(i) ' average FR (Hz) All BL vs BefRun All'],...
            ['Int_FRMeanC' num2str(i) 'BL-BefRunBoxAll'],...
            pathAnal,[],FRProfileMeanStatC.pRSBefRunVsBL(i),colorSel);
         
        plotBoxPlot(FRBefRunArr(:),FR0to1Arr(:),...
            ['C' num2str(i) ' average FR (Hz) All BefRun vs 0to1s All'],...
            ['Int_FRMeanC' num2str(i) 'BefRun-0to1BoxAll'],...
            pathAnal,[],FRProfileMeanStatC.pRS0to1VsBefRun(i),colorSel);
         
    end
        
    % compare each segment of avgFRProfile (field vs. no field)
    for i = 1:max(idxC) 
        indCurCField = idxC == i & nNeuWithField > 1;
        indCurCNoField = idxC == i & nNeuWithField < 1;
        FRBLArrF = avgFRProfile(indCurCField,FRProfileMean.indFRBaseline);
        FRBLArrNoF = avgFRProfile(indCurCNoField,FRProfileMean.indFRBaseline);
        FR0to1ArrF = avgFRProfile(indCurCField,FRProfileMean.indFR0to1);
        FR0to1ArrNoF = avgFRProfile(indCurCNoField,FRProfileMean.indFR0to1);
        plotBoxPlot(FR0to1ArrF(:)./FRBLArrF(:),...
            FR0to1ArrNoF(:)./FRBLArrNoF(:),...
            ['C' num2str(i) ' average FR change 0-1s/BL F/NoF'],...
            ['Int_FRChangeC' num2str(i) '0to1-BLFieldNoFieldBox'],...
            pathAnal,[0 3],FRProfileMeanStatC.pRSPercChange0to1VsBLFieldVsNoField(i),colorSel);
        
        FRBefRunArrF = avgFRProfile(indCurCField,FRProfileMean.indFRBefRun);
        FRBefRunArrNoF = avgFRProfile(indCurCNoField,FRProfileMean.indFRBefRun);
        plotBoxPlot(FRBefRunArrF(:)./FRBLArrF(:),...
            FRBefRunArrNoF(:)./FRBLArrNoF(:),...
            ['C' num2str(i) ' average FR change BefRun/BL F/NoF'],...
            ['Int_FRChangeC' num2str(i) 'BefRun-BLFieldNoFieldBox'],...
            pathAnal,[0 2],FRProfileMeanStatC.pRSPercChangeBefRunVsBLFieldVsNoField(i),colorSel);
        
        plotBoxPlot(FR0to1ArrF(:)./FRBefRunArrF(:),...
            FR0to1ArrNoF(:)./FRBefRunArrNoF(:),...
            ['C' num2str(i) ' average FR change 0-1s/BefRun F/NoF'],...
            ['Int_FRChangeC' num2str(i) '0to1-BefRunFieldNoFieldBox'],...
            pathAnal,[0 3.5],FRProfileMeanStatC.pRSPercChange0to1VsBefRunFieldVsNoField(i),colorSel);

    end
        
    % FR good vs bad trials
    for i = 1:max(idxC) 
        FRBLArrG = avgFRProfile(idxC == i,FRProfileMean.indFRBaseline);
        FRBLArrB = avgFRProfileBad(idxCBad == i,FRProfileMean.indFRBaseline); 
        plotBoxPlot(FRBLArrG(:),FRBLArrB(:),...
            ['C' num2str(i) ' average FR BL Good/Bad'],...
            ['Int_FRC' num2str(i) 'BLGoodBadTrBox'],...
            pathAnal,[],FRProfileMeanStatCGoodBad.pRSBL(i),colorSel);
        
        FRBefRunArrG = avgFRProfile(idxC == i,FRProfileMean.indFRBefRun);
        FRBefRunArrB = avgFRProfileBad(idxCBad == i,FRProfileMean.indFRBefRun);
        plotBoxPlot(FRBefRunArrG(:),FRBefRunArrB(:),...
            ['C' num2str(i) ' average FR BefRun Good/Bad'],...
            ['Int_FRC' num2str(i) 'BefRunGoodBadTrdBox'],...
            pathAnal,[],FRProfileMeanStatCGoodBad.pRSBefRun(i),colorSel);
        
        FR0to1ArrG = avgFRProfile(idxC == i,FRProfileMean.indFR0to1);
        FR0to1ArrB = avgFRProfileBad(idxCBad == i,FRProfileMean.indFR0to1);
        plotBoxPlot(FR0to1ArrG(:),FR0to1ArrB(:),...
            ['C' num2str(i) ' average FR 0-1s Good/Bad'],...
            ['Int_FRC' num2str(i) '0to1GoodBadTrBox'],...
            pathAnal,[],FRProfileMeanStatCGoodBad.pRS0to1(i),colorSel);
    end
    
    % FR change good vs bad trials
    for i = 1:max(idxC) 
        FRBLArrG = avgFRProfile(idxC == i,FRProfileMean.indFRBaseline);
        FRBLArrB = avgFRProfileBad(idxCBad == i,FRProfileMean.indFRBaseline); 
        FR0to1ArrG = avgFRProfile(idxC == i,FRProfileMean.indFR0to1);
        FR0to1ArrB = avgFRProfileBad(idxCBad == i,FRProfileMean.indFR0to1);
        plotBoxPlot(FR0to1ArrG(:)./FRBLArrG(:),...
            FR0to1ArrB(:)./FRBLArrB(:),...
            ['C' num2str(i) ' average FR change 0-1s/BL Good/Bad'],...
            ['Int_FRChangeC' num2str(i) '0to1-BLGoodBadTrBox'],...
            pathAnal,[0 3],FRProfileMeanStatCGoodBad.pRSPercChange0to1VsBL(i),colorSel);
        
        FRBefRunArrG = avgFRProfile(idxC == i,FRProfileMean.indFRBefRun);
        FRBefRunArrB = avgFRProfileBad(idxCBad == i,FRProfileMean.indFRBefRun);
        plotBoxPlot(FRBefRunArrG(:)./FRBLArrG(:),...
            FRBefRunArrB(:)./FRBLArrB(:),...
            ['C' num2str(i) ' average FR change BefRun/BL Good/Bad'],...
            ['Int_FRChangeC' num2str(i) 'BefRun-BLGoodBadTrBox'],...
            pathAnal,[0 2],FRProfileMeanStatCGoodBad.pRSPercChangeBefRunVsBL(i),colorSel);
        
        plotBoxPlot(FR0to1ArrG(:)./FRBefRunArrG(:),...
            FR0to1ArrB(:)./FRBefRunArrB(:),...
            ['C' num2str(i) ' average FR change 0-1s/BefRun Good/Bad'],...
            ['Int_FRChangeC' num2str(i) '0to1-BefRunGoodBadTrBox'],...
            pathAnal,[0 3.5],FRProfileMeanStatCGoodBad.pRSPercChange0to1VsBefRun(i),colorSel);
    end
    
end

function modInt1 = accumInterneurons2(paths,filenames,mazeSess,autoCorrIntAll,nNeuWithField,...
            nNeuWithFieldAligned,minFRInt,task,sampleFq)
    numRec = size(paths,1);
    modInt1 = struct('task',[],... % no cue - 1, AL - 2, PL - 3 Good trials
                              'taskBad',[],... % no cue - 1, AL - 2, PL - 3 Bad trials
                              'indRec',[],... % recording index  Good trials
                              'indRecBad',[],... % recording index Bad trials
                              'indNeu',[],... % neuron indices   Good trials
                              'indNeuBad',[],... % neuron indices  Bad trials
                              'idxC2',[],...  % cluster no. Good trials
                              'idxC2Bad',[],... % cluster no. for bad trials
                              'nNeuWithField',[],... % number of neurons with fields
                              'nNeuWithFieldAligned',[],... % number of neurons with fields after aligning to run onset
                              ...
                              'timeStepRun',[],...
                              'avgFRProfile',[],...% average firing rate profile good trials
                              'avgFRProfileBad',[]); % average firing rate profile bad trials
    
    totExcNeu = 0;
    for i = 1:numRec
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
                    'function "PeakFiringRate_AlignedRunOnset" first.']);
            return;
        end
        load(fullPath,'pFRNonStimGoodStruct','pFRNonStimBadStruct');
        
        fileNameConv = [filenames(i,:) '_convSpikesAligned_msess' num2str(mazeSess(i)) '_BefRun0.mat'];
        fullPath = [paths(i,:) fileNameConv];
        if(exist(fullPath) == 0)
            disp(['The convSpikesAligned file does not exist. Please call ',...
                    'function "ConvSpikeTrain_AlignedRunOnset" first.']);
            return;
        end
        load(fullPath,'timeStepRun');
        
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
        
        indTmp = find(autoCorrIntAll.task == task & autoCorrIntAll.indRec == i);
        if(length(indTmp) ~= sum(indNeu))
            disp(['the number of animals in recording task = ' num2str(task) ' rec. no. = ' num2str(indRec)...
                    'does not match that in the autoCorrIntAll struct.']);
        end
        
        if(length(pFRNonStimGoodStruct.indLapList) >= paramF.minNumTr)
            modInt1.avgFRProfile = [modInt1.avgFRProfile; pFRNonStimGoodStruct.avgFRProfile(indNeu,:)];
            modInt1.task = [modInt1.task task*ones(1,sum(indNeu))];
            modInt1.indRec = [modInt1.indRec i*ones(1,sum(indNeu))];
            modInt1.indNeu = [modInt1.indNeu find(indNeu == 1)];            
            if(length(indTmp) == sum(indNeu))
                modInt1.idxC2 = [modInt1.idxC2 autoCorrIntAll.idxC2(indTmp)'];
                modInt1.nNeuWithField = [modInt1.nNeuWithField nNeuWithField(indTmp)];
                modInt1.nNeuWithFieldAligned = [modInt1.nNeuWithFieldAligned nNeuWithFieldAligned(indTmp)];
            end
        else
            disp([filenames(i,:) ' only has ' num2str(length(pFRNonStimGoodStruct.indLapList)) ...
                ' good trials.']);
            disp(['No. interneurons in this recording is ' num2str(sum(indNeu))]);
        end
        if(~isempty(pFRNonStimBadStruct) && length(pFRNonStimBadStruct.indLapList) >= paramF.minNumTr)
            modInt1.avgFRProfileBad = [modInt1.avgFRProfileBad; pFRNonStimBadStruct.avgFRProfile(indNeu,:)];
            modInt1.taskBad = [modInt1.taskBad task*ones(1,sum(indNeu))];
            modInt1.indRecBad = [modInt1.indRecBad i*ones(1,sum(indNeu))];
            modInt1.indNeuBad = [modInt1.indNeuBad find(indNeu == 1)];        
            if(length(indTmp) == sum(indNeu))
                modInt1.idxC2Bad = [modInt1.idxC2Bad autoCorrIntAll.idxC2(indTmp)'];
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
        modInt1.timeStepRun = timeStepRun/sampleFq;
        
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
    indFR1to5 = find(timeStep >= 3 & timeStep < 5);  
    FRProfileMean.indFR1to5 = indFR1to5;
    FRProfileMean.meanAvgFRProfile1to5 = mean(avgFRProfile(:,indFR1to5),2);
    
    % perc change from 0.5-1.5 s to baseline
    FRProfileMean.percChange0to1VsBL = FRProfileMean.meanAvgFRProfile0to1...
        ./FRProfileMean.meanAvgFRProfileBaseline;
    
    % perc change -1.5- -0.5 s to baseline
    FRProfileMean.percChangeBefRunVsBL = FRProfileMean.meanAvgFRProfileBefRun...
        ./FRProfileMean.meanAvgFRProfileBaseline;
    
    % perc change 0.5 to 1.5 s to -1.5- -0.5 s 
    FRProfileMean.percChange0to1VsBefRun = FRProfileMean.meanAvgFRProfile0to1...
        ./FRProfileMean.meanAvgFRProfileBefRun;
    
    % perc change from 0.5-1.5 s to 3-5s
    FRProfileMean.percChange0to1Vs1to5 = FRProfileMean.meanAvgFRProfile0to1...
        ./FRProfileMean.meanAvgFRProfile1to5;
    
    % perc change -1.5- -0.5 s to 3-5s
    FRProfileMean.percChangeBefRunVs1to5 = FRProfileMean.meanAvgFRProfileBefRun...
        ./FRProfileMean.meanAvgFRProfile1to5;
end

function FRProfileMeanStatC = accumMeanStatC(FRProfileMean,idxC,nNeuWithField)
    
    numC = max(idxC);
    for i = 1:numC
        idxCI = idxC == i;
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
                
        FRProfileMeanStatC.pRS1to5VsBL(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(idxCI),...
                    FRProfileMean.meanAvgFRProfile1to5(idxCI));
        FRProfileMeanStatC.pRS1to5VsBLField(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(indCurCField),...
                    FRProfileMean.meanAvgFRProfile1to5(indCurCField));
        FRProfileMeanStatC.pRS1to5VsBLNoField(i) = ranksum(FRProfileMean.meanAvgFRProfileBaseline(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfile1to5(indCurCNoField));
                
        FRProfileMeanStatC.pRS0to1VsBefRun(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(idxCI),...
                    FRProfileMean.meanAvgFRProfileBefRun(idxCI));
        FRProfileMeanStatC.pRS0to1VsBefRunField(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(indCurCField),...
                    FRProfileMean.meanAvgFRProfileBefRun(indCurCField));
        FRProfileMeanStatC.pRS0to1VsBefRunNoField(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfileBefRun(indCurCNoField));
                
        FRProfileMeanStatC.pRS1to5Vs0to1(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(idxCI),...
                    FRProfileMean.meanAvgFRProfile1to5(idxCI));
        FRProfileMeanStatC.pRS1to5Vs0to1Field(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(indCurCField),...
                    FRProfileMean.meanAvgFRProfile1to5(indCurCField));
        FRProfileMeanStatC.pRS1to5Vs0to1NoField(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfile1to5(indCurCNoField));
                
        FRProfileMeanStatC.pRS1to5VsBefRun(i) = ranksum(FRProfileMean.meanAvgFRProfileBefRun(idxCI),...
                    FRProfileMean.meanAvgFRProfile1to5(idxCI));
        FRProfileMeanStatC.pRS1to5VsBefRunField(i) = ranksum(FRProfileMean.meanAvgFRProfileBefRun(indCurCField),...
                    FRProfileMean.meanAvgFRProfile1to5(indCurCField));
        FRProfileMeanStatC.pRS1to5VsBefRunNoField(i) = ranksum(FRProfileMean.meanAvgFRProfileBefRun(indCurCNoField),...
                    FRProfileMean.meanAvgFRProfile1to5(indCurCNoField));
                
        FRProfileMeanStatC.pTTPercChange0to1VsBL(i) = ttest(FRProfileMean.percChange0to1VsBL(idxCI));
        FRProfileMeanStatC.pRSPercChange0to1VsBLFieldVsNoField(i) = ...
            ranksum(FRProfileMean.percChange0to1VsBL(indCurCField),...
            FRProfileMean.percChange0to1VsBL(indCurCNoField));
        
        FRProfileMeanStatC.pTTPercChangeBefRunVsBL(i) = ttest(FRProfileMean.percChangeBefRunVsBL(idxCI));
        FRProfileMeanStatC.pRSPercChangeBefRunVsBLFieldVsNoField(i) = ...
            ranksum(FRProfileMean.percChangeBefRunVsBL(indCurCField),...
            FRProfileMean.percChangeBefRunVsBL(indCurCNoField));
        
        FRProfileMeanStatC.pTTPercChange0to1VsBefRun(i) = ttest(FRProfileMean.percChange0to1VsBefRun(idxCI));
        FRProfileMeanStatC.pRSPercChange0to1VsBefRunFieldVsNoField(i) = ...
            ranksum(FRProfileMean.percChange0to1VsBefRun(indCurCField),...
            FRProfileMean.percChange0to1VsBefRun(indCurCNoField));
        
        FRProfileMeanStatC.pTTPercChange0to1Vs1to5(i) = ttest(FRProfileMean.percChange0to1Vs1to5(idxCI));
        FRProfileMeanStatC.pRSPercChange0to1Vs1to5FieldVsNoField(i) = ...
            ranksum(FRProfileMean.percChange0to1Vs1to5(indCurCField),...
            FRProfileMean.percChange0to1Vs1to5(indCurCNoField));
        
        FRProfileMeanStatC.pTTPercChangeBefRunVs1to5(i) = ttest(FRProfileMean.percChangeBefRunVs1to5(idxCI));
        FRProfileMeanStatC.pRSPercChangeBefRunVs1to5FieldVsNoField(i) = ...
            ranksum(FRProfileMean.percChangeBefRunVs1to5(indCurCField),...
            FRProfileMean.percChangeBefRunVs1to5(indCurCNoField));
    end
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
               
        FRProfileMeanStatC.pRS1to5(i) = ranksum(FRProfileMean.meanAvgFRProfile1to5(idxCI),...
                    FRProfileMeanBad.meanAvgFRProfile1to5(idxCIBad));
               
        FRProfileMeanStatC.pRS0to1(i) = ranksum(FRProfileMean.meanAvgFRProfile0to1(idxCI),...
                    FRProfileMeanBad.meanAvgFRProfile0to1(idxCIBad));
                
        % perc change from 0-1 s to baseline
        FRProfileMeanStatC.pRSPercChange0to1VsBL(i) = ranksum(FRProfileMean.percChange0to1VsBL(idxCI),...
                    FRProfileMeanBad.percChange0to1VsBL(idxCIBad));

        % perc change -1- -0.5 s to baseline
        FRProfileMeanStatC.pRSPercChangeBefRunVsBL(i) = ranksum(FRProfileMean.percChangeBefRunVsBL(idxCI),...
                    FRProfileMeanBad.percChangeBefRunVsBL(idxCIBad));

        % perc change 0-1 s to -1- -0.5 s 
        FRProfileMeanStatC.pRSPercChange0to1VsBefRun(i) = ranksum(FRProfileMean.percChange0to1VsBefRun(idxCI),...
                    FRProfileMeanBad.percChange0to1VsBefRun(idxCIBad));

        % perc change from 0-1 s to 1-5s
        FRProfileMeanStatC.pRSPercChange0to1Vs1to5(i) = ranksum(FRProfileMean.percChange0to1Vs1to5(idxCI),...
                    FRProfileMeanBad.percChange0to1Vs1to5(idxCIBad));

        % perc change -1- -0.5 s to 1-5s
        FRProfileMeanStatC.pRSPercChangeBefRunVs1to5(i) = ranksum(FRProfileMean.percChangeBefRunVs1to5(idxCI),...
                    FRProfileMeanBad.percChangeBefRunVs1to5(idxCIBad));
        
    end
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
    set(gca,'XLim',[timeStepRun(1) 4]);
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

function plotIndFRProfile(timeStepRun,avgFRProfile,yl,fileName,pathAnal,ylimit,ordMethod,indT,indT1)
    if(isempty(avgFRProfile))
        return;
    end
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 200 400]);
    numNeurons = size(avgFRProfile,1);
    if(ordMethod == 1)
        [~,indMax] = max(avgFRProfile');
    elseif(ordMethod == 2)
        indMax = mean(avgFRProfile(:,indT)');
    elseif(ordMethod == 3)
        indMax = mean(avgFRProfile(:,indT)');
    elseif(ordMethod == 4)
        indMax1 = mean(avgFRProfile(:,indT1)');
        indMax2 = mean(avgFRProfile(:,indT)');
        indMax = indMax1./indMax2;
    end
    if(ordMethod == 4)
        [~,indOrd] = sort(indMax,'descend');
    else
        [~,indOrd] = sort(indMax);
    end
    h = imagesc(timeStepRun,1:numNeurons,avgFRProfile(indOrd,:));
%     set(h,'LineWidth',0.1)
    set(gca,'XLim',[-1 4]);
%     if(~isempty(ylimit))
%         set(gca,'YLim',ylimit);
%     else
%         set(gca,'YLim',[min(avgFRProfile(:)) max(avgFRProfile(:))]);
%     end
    xlabel('Time (s)')
    ylabel(yl)
    
    fileName1 = [pathAnal fileName];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
        
end

function plotAvgFRProfileCmp(timeStepRun,avgFRProfilex,avgFRProfiley, yl,fileName,pathAnal,ylimit)
    if(isempty(avgFRProfilex))
        return;
    end
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
    h = plot([0 0],[min([minX minY])*0.95 ...
        max([maxX maxY])*1.05],'r-');
    set(h,'LineWidth',1)
    set(gca,'XLim',[timeStepRun(1) 4]);
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