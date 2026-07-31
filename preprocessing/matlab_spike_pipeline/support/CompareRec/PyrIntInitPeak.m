function PyrIntInitPeak(taskSel)
% compare Pyramidal neurons and PV interneurons on their initial peak
    
    pathAnal0 = 'Z:\Yingxue\DataAnalysisRaphi\Pyramidal\';
    pathAnalInt0 = 'Z:\Yingxue\DataAnalysisRaphi\Interneuron\';
    
    if(exist([pathAnal0 'initPeakPyrAllRec.mat']))
        load([pathAnal0 'initPeakPyrAllRec.mat']);
    end
    if(exist([pathAnalInt0 'initPeakIntAllRec.mat']))
        load([pathAnalInt0 'initPeakIntAllRec.mat']);
    end
    
    if(taskSel == 1) % including all the neurons
        pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\Pyramidal\';
    elseif(taskSel == 2) % including AL and PL neurons
        pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\PyramidalALPL\';
    else
        pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\PyramidalAL\';
    end
    if(exist([pathAnal 'initPeakPyrIntAllRec.mat']))
        load([pathAnal 'initPeakPyrIntAllRec.mat']);
    end
    
    if(taskSel == 1)   
        task = [modPyr1NoCue.task modPyr1AL.task modPyr1PL.task];
        taskBad = [modPyr1NoCue.taskBad modPyr1AL.taskBad modPyr1PL.taskBad];
        indRec = [modPyr1NoCue.indRec modPyr1AL.indRec modPyr1PL.indRec];
        indRecBad = [modPyr1NoCue.indRecBad modPyr1AL.indRecBad modPyr1PL.indRecBad];
        indNeu = [modPyr1NoCue.indNeu modPyr1AL.indNeu modPyr1PL.indNeu];
        indNeuBad = [modPyr1NoCue.indNeuBad modPyr1AL.indNeuBad modPyr1PL.indNeuBad];
        idxC2 = [modPyr1NoCue.idxC2 modPyr1AL.idxC2 modPyr1PL.idxC2];
        avgFRProfile = [modPyr1NoCue.avgFRProfile; modPyr1AL.avgFRProfile; modPyr1PL.avgFRProfile];
        avgFRProfileBad = [modPyr1NoCue.avgFRProfileBad; modPyr1AL.avgFRProfileBad; modPyr1PL.avgFRProfileBad];
        relDepthNeuHDef = [modPyr1NoCue.relDepthNeuHDef modPyr1AL.relDepthNeuHDef modPyr1PL.relDepthNeuHDef];
        isNeuWithFieldAligned = [modPyr1NoCue.isNeuWithFieldAligned modPyr1AL.isNeuWithFieldAligned...
            modPyr1PL.isNeuWithFieldAligned];
        taskInt = [modInt1NoCue.task modInt1AL.task modInt1PL.task];
        taskIntBad = [modInt1NoCue.taskBad modInt1AL.taskBad modInt1PL.taskBad];
        indRecInt = [modInt1NoCue.indRec modInt1AL.indRec modInt1PL.indRec];
        indRecIntBad = [modInt1NoCue.indRecBad modInt1AL.indRecBad modInt1PL.indRecBad];
        indNeuInt = [modInt1NoCue.indNeu modInt1AL.indNeu modInt1PL.indNeu];
        indNeuIntBad = [modInt1NoCue.indNeuBad modInt1AL.indNeuBad modInt1PL.indNeuBad];
        idxC2Int = [modInt1NoCue.idxC2 modInt1AL.idxC2 modInt1PL.idxC2];
        avgFRProfileInt = [modInt1NoCue.avgFRProfile; modInt1AL.avgFRProfile; modInt1PL.avgFRProfile];
        avgFRProfileBadInt = [modInt1NoCue.avgFRProfileBad; modInt1AL.avgFRProfileBad; modInt1PL.avgFRProfileBad];
    elseif(taskSel == 2)
        task = [modPyr1AL.task modPyr1PL.task];
        taskBad = [modPyr1AL.taskBad modPyr1PL.taskBad];
        indRec = [modPyr1AL.indRec modPyr1PL.indRec];
        indRecBad = [modPyr1AL.indRecBad modPyr1PL.indRecBad];
        indNeu = [modPyr1AL.indNeu modPyr1PL.indNeu];
        indNeuBad = [modPyr1AL.indNeuBad modPyr1PL.indNeuBad];
        idxC2 = [modPyr1AL.idxC2 modPyr1PL.idxC2];
        avgFRProfile = [modPyr1AL.avgFRProfile; modPyr1PL.avgFRProfile];
        avgFRProfileBad = [modPyr1AL.avgFRProfileBad; modPyr1PL.avgFRProfileBad];
        relDepthNeuHDef = [modPyr1AL.relDepthNeuHDef modPyr1PL.relDepthNeuHDef];
        isNeuWithFieldAligned = [modPyr1AL.isNeuWithFieldAligned...
            modPyr1PL.isNeuWithFieldAligned];
        taskInt = [modInt1AL.task modInt1PL.task];
        taskIntBad = [modInt1AL.taskBad modInt1PL.taskBad];
        indRecInt = [modInt1AL.indRec modInt1PL.indRec];
        indRecIntBad = [modInt1AL.indRecBad modInt1PL.indRecBad];
        indNeuInt = [modInt1AL.indNeu modInt1PL.indNeu];
        indNeuIntBad = [modInt1AL.indNeuBad modInt1PL.indNeuBad];
        idxC2Int = [modInt1AL.idxC2 modInt1PL.idxC2];
        avgFRProfileInt = [modInt1AL.avgFRProfile; modInt1PL.avgFRProfile];
        avgFRProfileBadInt = [modInt1AL.avgFRProfileBad; modInt1PL.avgFRProfileBad];
    else
        task = modPyr1AL.task;
        taskBad = modPyr1AL.taskBad;
        indRec = modPyr1AL.indRec;
        indRecBad = modPyr1AL.indRecBad;
        indNeu = modPyr1AL.indNeu;
        indNeuBad = modPyr1AL.indNeuBad;
        idxC2 = modPyr1AL.idxC2;
        avgFRProfile = modPyr1AL.avgFRProfile;
        avgFRProfileBad = modPyr1AL.avgFRProfileBad;
        relDepthNeuHDef = modPyr1AL.relDepthNeuHDef;
        isNeuWithFieldAligned = modPyr1AL.isNeuWithFieldAligned;
        taskInt = modInt1AL.task;
        taskIntBad = modInt1AL.taskBad;
        indRecInt = modInt1AL.indRec;
        indRecIntBad = modInt1AL.indRecBad;
        indNeuInt = modInt1AL.indNeu;
        indNeuIntBad = modInt1AL.indNeuBad;
        idxC2Int = modInt1AL.idxC2;
        avgFRProfileInt = modInt1AL.avgFRProfile;
        avgFRProfileBadInt = modInt1AL.avgFRProfileBad;        
    end
    
    indTmp = modPyr1NoCue.timeStepRun >=-1 & modPyr1NoCue.timeStepRun <= 4 ;
    avgFRProfileNorm = zeros(size(avgFRProfile,1),size(avgFRProfile,2));
    for i = 1:size(avgFRProfile,1)
        if(max(avgFRProfile(i,:)) ~= 0)
            avgFRProfileNorm(i,:) = avgFRProfile(i,:)/max(avgFRProfile(i,:));
        end
    end
    n = 0;
    avgFRProfileNorm1 = zeros(size(avgFRProfile,1),size(avgFRProfile,2));
    for i = 1:size(avgFRProfile,1)
        if(max(avgFRProfile(i,:)) ~= 0)
            rangeFR = max(avgFRProfile(i,indTmp))-min(avgFRProfile(i,indTmp));
            if(rangeFR ~= 0)
                avgFRProfileNorm1(i,:) = (avgFRProfile(i,:)-min(avgFRProfile(i,indTmp)))...
                    /rangeFR;
            else
                n = n+1;
            end
        end
    end
    avgFRProfileNormBad = zeros(size(avgFRProfileBad,1),size(avgFRProfileBad,2));
    for i = 1:size(avgFRProfileBad,1)
        if(max(avgFRProfileBad(i,:)) ~= 0)
            avgFRProfileNormBad(i,:) = avgFRProfileBad(i,:)/max(avgFRProfileBad(i,:));
        end
    end
    
    avgFRProfileNormInt = zeros(size(avgFRProfileInt,1),size(avgFRProfileInt,2));
    for i = 1:size(avgFRProfileInt,1)
        if(max(avgFRProfileInt(i,:)) ~= 0)
            avgFRProfileNormInt(i,:) = avgFRProfileInt(i,:)/max(avgFRProfileInt(i,:));
        end
    end
    avgFRProfileNormInt1 = zeros(size(avgFRProfileInt,1),size(avgFRProfileInt,2));
    m = 0;
    for i = 1:size(avgFRProfileInt,1)
        if(max(avgFRProfileInt(i,:)) ~= 0)
            rangeFR = max(avgFRProfileInt(i,indTmp))-min(avgFRProfileInt(i,indTmp));
            if(rangeFR ~= 0)
                avgFRProfileNormInt1(i,:) = (avgFRProfileInt(i,:)-min(avgFRProfileInt(i,indTmp)))...
                    /rangeFR;
            else
                m = m+1;
            end
        end
    end
    avgFRProfileNormBadInt = zeros(size(avgFRProfileBadInt,1),size(avgFRProfileBadInt,2));
    for i = 1:size(avgFRProfileBadInt,1)
        if(max(avgFRProfileBadInt(i,:)) ~= 0)
            avgFRProfileNormBadInt(i,:) = avgFRProfileBadInt(i,:)/max(avgFRProfileBadInt(i,:));
        end
    end
        
    mean0to1 = mean(avgFRProfile(:,FRProfileMean.indFR0to1),2);
    meanBefRun = mean(avgFRProfile(:,FRProfileMean.indFRBefRun),2);
    ratio0to1BefRun = mean0to1./meanBefRun;
    [ratio0to1BefRunOrd,indOrd] = sort(ratio0to1BefRun,'descend');
    idxNan = isnan(ratio0to1BefRunOrd);
    idxInf = isinf(ratio0to1BefRunOrd);
%     idx = find(ratio0to1BefRunOrd < 1.25,1);
%     idx1 = find(ratio0to1BefRunOrd >= 0.8,1,'last');
    
%     idx = find(ratio0to1BefRunOrd < 2,1);
%     idx1 = find(ratio0to1BefRunOrd >= 0.5,1,'last');
    
    idx = find(ratio0to1BefRunOrd < 1.5,1);
    idx1 = find(ratio0to1BefRunOrd >= 2/3,1,'last');
    
    %% neurons with FR increase around 0
    indOrdTmp = indOrd(1:idx);
    idxNanTmp = idxNan(1:idx);
    indOrdTmp = indOrdTmp(idxNanTmp == 0);
    PyrRise.idxRise = indOrdTmp;
    PyrRise.task = task(indOrdTmp);
    PyrRise.indRec = indRec(indOrdTmp);
    PyrRise.indNeu = indNeu(indOrdTmp);
    PyrRise.idxRiseBad = [];
    PyrRise.taskBad = [];
    PyrRise.indRecBad = [];
    PyrRise.indNeuBad = [];
    for i = 1:length(indOrdTmp)
        idxBad = find(taskBad == PyrRise.task(i) & indRecBad == PyrRise.indRec(i) ...
            & indNeuBad == PyrRise.indNeu(i));
        if(sum(idxBad) > 0)
            PyrRise.idxRiseBad = [PyrRise.idxRiseBad idxBad];
            PyrRise.taskBad = [PyrRise.taskBad taskBad(idxBad)];
            PyrRise.indRecBad = [PyrRise.indRecBad indRecBad(idxBad)];
            PyrRise.indNeuBad = [PyrRise.indNeuBad indNeuBad(idxBad)];
        end
    end
    
    %% neurons with FR decrease around 0
    indOrdTmp = indOrd(idx1:end);
    idxNanTmp = idxNan(idx1:end);
    indOrdTmp = indOrdTmp(idxNanTmp == 0);
    PyrDown.idxDown = indOrdTmp;
    PyrDown.task = task(indOrdTmp);
    PyrDown.indRec = indRec(indOrdTmp);
    PyrDown.indNeu = indNeu(indOrdTmp);
    PyrDown.idxDownBad = [];
    PyrDown.taskBad = [];
    PyrDown.indRecBad = [];
    PyrDown.indNeuBad = [];
    for i = 1:length(indOrdTmp)
        idxBad = find(taskBad == PyrDown.task(i) & indRecBad == PyrDown.indRec(i) ...
            & indNeuBad == PyrDown.indNeu(i));
        if(sum(idxBad) > 0)
            PyrDown.idxDownBad = [PyrDown.idxDownBad idxBad];
            PyrDown.taskBad = [PyrDown.taskBad taskBad(idxBad)];
            PyrDown.indRecBad = [PyrDown.indRecBad indRecBad(idxBad)];
            PyrDown.indNeuBad = [PyrDown.indNeuBad indNeuBad(idxBad)];
        end
    end
    
    for i = 1:max(idxC2Int)
        idxCTmp = find(idxC2Int == i);
        mean0to1Int = mean(avgFRProfileInt(idxC2Int == i,FRProfileMean.indFR0to1),2);
        meanBefRunInt = mean(avgFRProfileInt(idxC2Int == i,FRProfileMean.indFRBefRun),2);
        ratio0to1BefRunInt = mean0to1Int./meanBefRunInt;
        [ratio0to1BefRunIntOrd,indOrdInt] = sort(ratio0to1BefRunInt,'descend');
    %     idxInt = find(ratio0to1BefRunIntOrd < 1.25,1);
    %     idxInt1 = find(ratio0to1BefRunIntOrd >= 0.8,1,'last');

    %     idxInt = find(ratio0to1BefRunIntOrd < 2,1);
    %     idxInt1 = find(ratio0to1BefRunIntOrd >= 0.5,1,'last');

        idxNan = isnan(ratio0to1BefRunIntOrd);
        idxInt = find(ratio0to1BefRunIntOrd < 1.5,1);
        idxInt1 = find(ratio0to1BefRunIntOrd >= 2/3,1,'last');
        
        %% neurons with FR increase around 0
        indOrdTmp = idxCTmp(indOrdInt(1:idxInt));
        idxNanTmp = idxNan(1:idxInt);
        indOrdTmp = indOrdTmp(idxNanTmp == 0);
        IntRise.idxRise{i} = indOrdTmp;
        IntRise.task{i} = taskInt(indOrdTmp);
        IntRise.indRec{i} = indRecInt(indOrdTmp);
        IntRise.indNeu{i} = indNeuInt(indOrdTmp);
        IntRise.idxRiseBad{i} = [];
        IntRise.taskBad{i} = [];
        IntRise.indRecBad{i} = [];
        IntRise.indNeuBad{i} = [];
        for j = 1:length(indOrdTmp)
            idxBad = find(taskIntBad == IntRise.task{i}(j) & indRecIntBad == IntRise.indRec{i}(j) ...
                & indNeuIntBad == IntRise.indNeu{i}(j));
            if(sum(idxBad) > 0)
                IntRise.idxRiseBad{i} = [IntRise.idxRiseBad{i} idxBad];
                IntRise.taskBad{i} = [IntRise.taskBad{i} taskIntBad(idxBad)];
                IntRise.indRecBad{i} = [IntRise.indRecBad{i} indRecIntBad(idxBad)];
                IntRise.indNeuBad{i} = [IntRise.indNeuBad{i} indNeuIntBad(idxBad)];
            end
        end

        %% neurons with FR decrease around 0
        indOrdTmp = idxCTmp(indOrdInt(idxInt1:end));
        idxNanTmp = idxNan(idxInt1:end);
        indOrdTmp = indOrdTmp(idxNanTmp == 0);
        IntDown.idxDown{i} = indOrdTmp;
        IntDown.task{i} = taskInt(indOrdTmp);
        IntDown.indRec{i} = indRecInt(indOrdTmp);
        IntDown.indNeu{i} = indNeuInt(indOrdTmp);
        IntDown.idxDownBad{i} = [];
        IntDown.taskBad{i} = [];
        IntDown.indRecBad{i} = [];
        IntDown.indNeuBad{i} = [];
        for j = 1:length(indOrdTmp)
            idxBad = find(taskIntBad == IntDown.task{i}(j) & indRecIntBad == IntDown.indRec{i}(j) ...
                & indNeuIntBad == IntDown.indNeu{i}(j));
            if(sum(idxBad) > 0)
                IntDown.idxDownBad{i} = [IntDown.idxDownBad{i} idxBad];
                IntDown.taskBad{i} = [IntDown.taskBad{i} taskIntBad(idxBad)];
                IntDown.indRecBad{i} = [IntDown.indRecBad{i} indRecIntBad(idxBad)];
                IntDown.indNeuBad{i} = [IntDown.indNeuBad{i} indNeuIntBad(idxBad)];
            end
        end
    end
    
    %% pyramidal neurons
    FRProfileMeanRise = accumMean(avgFRProfile(PyrRise.idxRise,:),modPyr1NoCue.timeStepRun);
    FRProfileMeanDown = accumMean(avgFRProfile(PyrDown.idxDown,:),modPyr1NoCue.timeStepRun);
    
    relDepthNeuHDefRise.depth = relDepthNeuHDef(PyrRise.idxRise);
    relDepthNeuHDefDown.depth = relDepthNeuHDef(PyrDown.idxDown);
    relDepthNeuHDefRise.depthMean = mean(relDepthNeuHDefRise.depth);
    relDepthNeuHDefDown.depthMean = mean(relDepthNeuHDefDown.depth);
    relDepthNeuHDefRise.depthSem = std(relDepthNeuHDefRise.depth)/...
            sqrt(length(relDepthNeuHDefRise.depth));
    relDepthNeuHDefDown.depthSem = std(relDepthNeuHDefDown.depth)/...
            sqrt(length(relDepthNeuHDefDown.depth));  
    relDepthNeuHDefRise.pRSRelDepthNeuHDef = ranksum(relDepthNeuHDefRise.depth,relDepthNeuHDefDown.depth);
    
    %% for field
    isNeuWithFieldAlignedRise.isField = isNeuWithFieldAligned(PyrRise.idxRise);
    isNeuWithFieldAlignedDown.isField = isNeuWithFieldAligned(PyrDown.idxDown);
    isNeuWithFieldAlignedRise.isFieldBad = isNeuWithFieldAligned(PyrRise.idxRiseBad);
    isNeuWithFieldAlignedDown.isFieldBad = isNeuWithFieldAligned(PyrDown.idxDownBad);
    isNeuWithFieldAlignedRise.numField = sum(isNeuWithFieldAlignedRise.isField);
    isNeuWithFieldAlignedDown.numField = sum(isNeuWithFieldAlignedDown.isField);
    isNeuWithFieldAlignedRise.idxRise = PyrRise.idxRise(isNeuWithFieldAlignedRise.isField == 1);
    isNeuWithFieldAlignedDown.idxDown = PyrDown.idxDown(isNeuWithFieldAlignedDown.isField == 1);
    isNeuWithFieldAlignedRise.idxRiseBad = PyrRise.idxRiseBad(isNeuWithFieldAlignedRise.isFieldBad == 1);
    isNeuWithFieldAlignedDown.idxDownBad = PyrDown.idxDownBad(isNeuWithFieldAlignedDown.isFieldBad == 1);
    
    FRProfileMeanRiseField = accumMean(avgFRProfile(isNeuWithFieldAlignedRise.idxRise,:),...
        modPyr1NoCue.timeStepRun);
    FRProfileMeanDownField = accumMean(avgFRProfile(isNeuWithFieldAlignedDown.idxDown,:),...
        modPyr1NoCue.timeStepRun);
    FRProfileMeanRiseFieldBad = accumMean(avgFRProfileBad(isNeuWithFieldAlignedRise.idxRiseBad,:),...
        modPyr1NoCue.timeStepRun);
    FRProfileMeanDownFieldBad = accumMean(avgFRProfileBad(isNeuWithFieldAlignedDown.idxDownBad,:),...
        modPyr1NoCue.timeStepRun);
    
    FRProfileMeanStatRiseField = accumMeanStatC(FRProfileMeanRiseField);
    FRProfileMeanStatDownField = accumMeanStatC(FRProfileMeanDownField);
    
    FRProfileMeanStatRiseFieldBad = accumMeanStatC(FRProfileMeanRiseFieldBad);
    FRProfileMeanStatDownFieldBad = accumMeanStatC(FRProfileMeanDownFieldBad);
    
    % compare good and bad trials
    FRProfileMeanStatRiseFieldGoodBad = accumMeanStatCGoodBad(FRProfileMeanRiseField,FRProfileMeanRiseFieldBad);
    FRProfileMeanStatDownFieldGoodBad = accumMeanStatCGoodBad(FRProfileMeanDownField,FRProfileMeanDownFieldBad);
    
    %% for all the pyramidal neurons
    FRProfileMeanRiseBad = accumMean(avgFRProfileBad(PyrRise.idxRiseBad,:),modPyr1NoCue.timeStepRun);
    FRProfileMeanDownBad = accumMean(avgFRProfileBad(PyrDown.idxDownBad,:),modPyr1NoCue.timeStepRun);
    
    FRProfileMeanStatRise = accumMeanStatC(FRProfileMeanRise);
    FRProfileMeanStatDown = accumMeanStatC(FRProfileMeanDown);
    
    FRProfileMeanStatRiseBad = accumMeanStatC(FRProfileMeanRiseBad);
    FRProfileMeanStatDownBad = accumMeanStatC(FRProfileMeanDownBad);
    
    % compare good and bad trials
    FRProfileMeanStatRiseGoodBad = accumMeanStatCGoodBad(FRProfileMeanRise,FRProfileMeanRiseBad);
    FRProfileMeanStatDownGoodBad = accumMeanStatCGoodBad(FRProfileMeanDown,FRProfileMeanDownBad);
    
    for i = 1:max(idxC2Int)
        FRProfileMeanRiseInt{i} = accumMean(avgFRProfileInt(IntRise.idxRise{i},:),modInt1NoCue.timeStepRun);
        FRProfileMeanDownInt{i} = accumMean(avgFRProfileInt(IntDown.idxDown{i},:),modInt1NoCue.timeStepRun);

        FRProfileMeanRiseIntBad{i} = accumMean(avgFRProfileBadInt(IntRise.idxRiseBad{i},:),modInt1NoCue.timeStepRun);
        FRProfileMeanDownIntBad{i} = accumMean(avgFRProfileBadInt(IntDown.idxDownBad{i},:),modInt1NoCue.timeStepRun);

        FRProfileMeanStatRiseInt{i} = accumMeanStatC(FRProfileMeanRiseInt{i});
        FRProfileMeanStatDownInt{i} = accumMeanStatC(FRProfileMeanDownInt{i});

        FRProfileMeanStatRiseIntBad{i} = accumMeanStatC(FRProfileMeanRiseIntBad{i});
        FRProfileMeanStatDownIntBad{i} = accumMeanStatC(FRProfileMeanDownIntBad{i});

        % compare good and bad trials
        FRProfileMeanStatRiseIntGoodBad{i} = accumMeanStatCGoodBad(FRProfileMeanRiseInt{i},FRProfileMeanRiseIntBad{i});
        FRProfileMeanStatDownIntGoodBad{i} = accumMeanStatCGoodBad(FRProfileMeanDownInt{i},FRProfileMeanDownIntBad{i});
    end
    
    save([pathAnal 'initPeakPyrIntAllRec.mat'],'PyrRise','PyrDown',...
        'FRProfileMeanRise','FRProfileMeanRiseBad',...
        'FRProfileMeanDown','FRProfileMeanDownBad',...
        'FRProfileMeanStatRise','FRProfileMeanStatDown',...
        'FRProfileMeanStatRiseBad','FRProfileMeanStatDownBad',...
        'FRProfileMeanStatRiseGoodBad','FRProfileMeanStatDownGoodBad',...
        'IntRise','IntDown','FRProfileMeanRiseInt','FRProfileMeanRiseIntBad',...
        'FRProfileMeanDownInt','FRProfileMeanDownIntBad',...
        'FRProfileMeanStatRiseInt','FRProfileMeanStatDownInt',...
        'FRProfileMeanStatRiseIntBad','FRProfileMeanStatDownIntBad',...
        'FRProfileMeanStatRiseIntGoodBad','FRProfileMeanStatDownIntGoodBad',...
        'relDepthNeuHDefRise','relDepthNeuHDefDown',...
        'isNeuWithFieldAlignedRise','isNeuWithFieldAlignedDown',...
        'FRProfileMeanRiseField','FRProfileMeanDownField','FRProfileMeanRiseFieldBad','FRProfileMeanDownFieldBad',...
        'FRProfileMeanStatRiseField','FRProfileMeanStatDownField','FRProfileMeanStatRiseFieldBad',...
        'FRProfileMeanStatDownFieldBad','FRProfileMeanStatRiseFieldGoodBad','FRProfileMeanStatDownFieldGoodBad'); 
    
    %% plot depth of rise and down pyramidal neurons
    colorSel = 0;
    plotBoxPlot(relDepthNeuHDefRise.depth',...
        relDepthNeuHDefDown.depth','Depth',...
        'Pyr_relDepthNeuHDefRiseVsDown',pathAnal,[],relDepthNeuHDefRise.pRSRelDepthNeuHDef,colorSel);
    
    %% order pyramidal neurons with field based on the peak firing rate after 0
    plotIndFRProfile(modInt1NoCue.timeStepRun,...
            avgFRProfileNorm(isNeuWithFieldAligned == 1,:),['FR (Hz)'],...
            ['Pyr_IndFRProfileNormFRPeakAftRunNeuField'],...
            pathAnal,[],5,[],[])
        
    %% order pyramidal neurons with field based on before and after run FR ratio  
    plotIndFRProfile(modInt1NoCue.timeStepRun,...
            avgFRProfileNorm(isNeuWithFieldAligned == 1,:),['FR (Hz)'],...
            ['Pyr_IndFRProfileNormFR0to1VsBefRunNeuField'],...
            pathAnal,[],4,FRProfileMean.indFRBefRun,...
            FRProfileMean.indFR0to1) % ordered based on -1to0 to 0to1 mean ratio
    
    plotAvgFRProfileCmp(modPyr1NoCue.timeStepRun,...
            avgFRProfile(isNeuWithFieldAlignedRise.idxRise,:),...
            avgFRProfileBad(isNeuWithFieldAlignedRise.idxRiseBad,:),...
            ['FR PyrRise Good/Bad Field'],...
            ['Pyr_FRProfilePyrFieldRiseGoodBad'],...
            pathAnal,[])
        
    plotAvgFRProfileCmp(modPyr1NoCue.timeStepRun,...
            avgFRProfile(isNeuWithFieldAlignedDown.idxDown,:),...
            avgFRProfileBad(isNeuWithFieldAlignedDown.idxDownBad,:),...
            ['FR PyrDown Good/Bad Field'],...
            ['Pyr_FRProfilePyrFieldDownGoodBad'],...
            pathAnal,[])
        
    plotBoxPlot(FRProfileMeanRiseField.percChangeBefRunVs0to1,...
        FRProfileMeanRiseFieldBad.percChangeBefRunVs0to1,...
        ['average FR change BefRun/0to1s G/B Rise Field'],...
        ['Pyr_FRChangeBefRun-0to1RiseGoodBadFieldBox'],...
        pathAnal,[-1 14],FRProfileMeanStatRiseFieldGoodBad.pRSPercChangeBefRunVs0to1All,colorSel);
    
    plotBoxPlot(FRProfileMeanDownField.percChangeBefRunVs0to1,...
        FRProfileMeanDownFieldBad.percChangeBefRunVs0to1,...
        ['average FR change BefRun/0to1s G/B Down Field'],...
        ['Pyr_FRChangeBefRun-0to1DownFieldGoodBadBox'],...
        pathAnal,[-1 3],FRProfileMeanStatDownFieldGoodBad.pRSPercChangeBefRunVs0to1All,colorSel);
        
    %% order pyramidal cells based on before and after run FR ratio
    plotIndFRProfile(modInt1NoCue.timeStepRun,...
            avgFRProfileNorm,['FR (Hz)'],...
            ['Pyr_IndFRProfileNormFR0to1VsBefRunAll'],...
            pathAnal,[],4,FRProfileMean.indFRBefRun,FRProfileMean.indFR0to1) % ordered based on -1to0 to 0to1 mean ratio
    
    plotAvgFRProfileCmp(modPyr1NoCue.timeStepRun,...
            avgFRProfile(PyrRise.idxRise,:),avgFRProfileBad(PyrRise.idxRiseBad,:),...
            ['FR PyrRise Good/Bad'],...
            ['Pyr_FRProfilePyrRiseGoodBad'],...
            pathAnal,[])
        
    plotAvgFRProfileCmp(modPyr1NoCue.timeStepRun,...
            avgFRProfileNorm(PyrRise.idxRise,:),avgFRProfileNormBad(PyrRise.idxRiseBad,:),...
            ['Norm FR PyrRise Good/Bad'],...
            ['Pyr_FRProfileNormPyrRiseGoodBad'],...
            pathAnal,[])
        
    plotAvgFRProfileCmp(modPyr1NoCue.timeStepRun,...
            avgFRProfile(PyrDown.idxDown,:),avgFRProfileBad(PyrDown.idxDownBad,:),...
            ['FR PyrDown Good/Bad'],...
            ['Pyr_FRProfilePyrDownGoodBad'],...
            pathAnal,[])
        
    plotAvgFRProfileCmp(modPyr1NoCue.timeStepRun,...
            avgFRProfileNorm(PyrDown.idxDown,:),avgFRProfileNormBad(PyrDown.idxDownBad,:),...
            ['Norm FR PyrDown Good/Bad'],...
            ['Pyr_FRProfileNormPyrDownGoodBad'],...
            pathAnal,[])
    
    for i = 1:max(idxC2Int)
        plotAvgFRProfileCmp(modInt1NoCue.timeStepRun,...
            avgFRProfileInt(IntRise.idxRise{i},:),avgFRProfileBadInt(IntRise.idxRiseBad{i},:),...
            ['FR IntC' num2str(i) ' Rise Good/Bad'],...
            ['Int_FRProfileIntRiseGoodBadC' num2str(i)],...
            pathAnal,[])
        
        plotAvgFRProfileCmp(modInt1NoCue.timeStepRun,...
                avgFRProfileNormInt(IntRise.idxRise{i},:),avgFRProfileNormBadInt(IntRise.idxRiseBad{i},:),...
                ['Norm FR IntC' num2str(i) ' Rise Good/Bad'],...
                ['Int_FRProfileNormIntRiseGoodBadC' num2str(i)],...
                pathAnal,[])

%         plotAvgFRProfileCmp(modInt1NoCue.timeStepRun,...
%                 avgFRProfileInt(IntDown.idxDown{i},:),avgFRProfileIntBad(IntDown.idxDownBad{i},:),...
%                 ['FR IntC' num2str(i) ' Down Good/Bad'],...
%                 ['Int_FRProfileIntDownGoodBadC' num2str(i)],...
%                 pathAnal,[])
% 
%         plotAvgFRProfileCmp(modInt1NoCue.timeStepRun,...
%                 avgFRProfileNormInt(IntDown.idxDown{i},:),avgFRProfileNormBadInt(IntDown.idxDownBad{i},:),...
%                 ['Norm FR IntC' num2str(i) ' Down Good/Bad'],...
%                 ['Int_FRProfileNormIntDownGoodBadC' num2str(i)],...
%                 pathAnal,[])
    end
        
    %% compare pyramidal cells profile with each type of interneuron
    for i = 1:max(idxC2Int)
        mean0to1Int = mean(avgFRProfileInt(idxC2Int == i,FRProfileMean.indFR0to1),2);
        meanBefRunInt = mean(avgFRProfileInt(idxC2Int == i,FRProfileMean.indFRBefRun),2);
        ratio0to1BefRunInt = mean0to1Int./meanBefRunInt;
        [ratio0to1BefRunIntOrd,indOrdInt] = sort(ratio0to1BefRunInt,'descend');
    %     idxInt = find(ratio0to1BefRunIntOrd < 1.25,1);
    %     idxInt1 = find(ratio0to1BefRunIntOrd >= 0.8,1,'last');

    %     idxInt = find(ratio0to1BefRunIntOrd < 2,1);
    %     idxInt1 = find(ratio0to1BefRunIntOrd >= 0.5,1,'last');

        idxInt = find(ratio0to1BefRunIntOrd < 1.5,1);
        idxInt1 = find(ratio0to1BefRunIntOrd >= 2/3,1,'last');
       
        meanAvgFRProfile = mean(avgFRProfileNorm(indOrd(1:idx),:));
        meanAvgFRProfile1 = mean(avgFRProfileNorm(indOrd(idx1:end),:));
        avgFRProfileIntC = avgFRProfileNormInt(idxC2Int == i,:);
        meanAvgFRProfileInt = mean(avgFRProfileIntC(indOrdInt(1:idxInt),:));
        
        plotMeanFRProfilePyrVsInt(modPyr1NoCue.timeStepRun,...
            (meanAvgFRProfile-min(meanAvgFRProfile(indTmp)))...
            /(max(meanAvgFRProfile(indTmp))-min(meanAvgFRProfile(indTmp))),...
            (meanAvgFRProfileInt-min(meanAvgFRProfileInt(indTmp)))...
            /(max(meanAvgFRProfileInt(indTmp))-min(meanAvgFRProfileInt(indTmp))),...
            ['Norm. FR PyrInc Vs Int C' num2str(i)],...
            ['FRProfileNormPyrIncVsIntC' num2str(i)],pathAnal,[]);
        
        plotMeanFRProfilePyrVsInt(modPyr1NoCue.timeStepRun,...
            (meanAvgFRProfile1-min(meanAvgFRProfile1(indTmp)))...
            /(max(meanAvgFRProfile1(indTmp))-min(meanAvgFRProfile1(indTmp))),...
            (meanAvgFRProfileInt-min(meanAvgFRProfileInt(indTmp)))...
            /(max(meanAvgFRProfileInt(indTmp))-min(meanAvgFRProfileInt(indTmp))),...
            ['Norm. FR PyrDec Vs Int C' num2str(i)],...
            ['FRProfileNormPyrDecVsIntC' num2str(i)],pathAnal,[]);
        
        plotAvgFRProfileCmp(modPyr1NoCue.timeStepRun,...
            avgFRProfileNorm(indOrd(1:idx),:),avgFRProfileIntC(indOrdInt(1:idxInt),:),...
            ['Norm FR Pyr vs. IntC' num2str(i)],...
            ['FRProfileNormPyrVsIntC' num2str(i)],...
            pathAnal,[0 1])
    end
    
    %% FR change good vs. bad trials
    colorSel = 0;
    plotBoxPlot(FRProfileMeanRise.percChange0to1VsBL,...
        FRProfileMeanRiseBad.percChange0to1VsBL,...
        ['average FR change 0-1s/BL G/B Rise'],...
        ['Pyr_FRChange0to1-BLRiseGoodBadBox'],...
        pathAnal,[-1 5],FRProfileMeanStatRiseGoodBad.pRSPercChange0to1VsBLAll,colorSel);

    plotBoxPlot(FRProfileMeanRise.percChangeBefRunVsBL,...
        FRProfileMeanRiseBad.percChangeBefRunVsBL,...
        ['average FR change BefRun/BL G/B Rise'],...
        ['Pyr_FRChangeBefRun-BLRiseGoodBadBox'],...
        pathAnal,[-1 3],FRProfileMeanStatRiseGoodBad.pRSPercChangeBefRunVsBLAll,colorSel);

    plotBoxPlot(FRProfileMeanRise.percChangeBefRunVs0to1,...
        FRProfileMeanRiseBad.percChangeBefRunVs0to1,...
        ['average FR change BefRun/0to1s G/B Rise'],...
        ['Pyr_FRChangeBefRun-0to1RiseGoodBadBox'],...
        pathAnal,[-1 7],FRProfileMeanStatRiseGoodBad.pRSPercChangeBefRunVs0to1All,colorSel);
    
    plotBoxPlot(FRProfileMeanRise.percChange0to1Vs3to5,...
        FRProfileMeanRiseBad.percChange0to1Vs3to5,...
        ['average FR change 0-1s/3-5s G/B Rise'],...
        ['Pyr_FRChange0to1-3to5RiseGoodBadBox'],...
        pathAnal,[-1 7],FRProfileMeanStatRiseGoodBad.pRSPercChange0to1Vs3to5All,colorSel);

    plotBoxPlot(FRProfileMeanRise.percChangeBefRunVs3to5,...
        FRProfileMeanRiseBad.percChangeBefRunVs3to5,...
        ['average FR change BefRun/3-5s G/B Rise'],...
        ['Pyr_FRChangeBefRun-3to5RiseGoodBadBox'],...
        pathAnal,[-1 3],FRProfileMeanStatRiseGoodBad.pRSPercChangeBefRunVs3to5All,colorSel);
    
    plotBoxPlot(FRProfileMeanDown.percChange0to1VsBL,...
        FRProfileMeanDownBad.percChange0to1VsBL,...
        ['average FR change 0-1s/BL G/B Down'],...
        ['Pyr_FRChange0to1-BLDownGoodBadBox'],...
        pathAnal,[-1 3],FRProfileMeanStatDownGoodBad.pRSPercChange0to1VsBLAll,colorSel);

    plotBoxPlot(FRProfileMeanDown.percChangeBefRunVsBL,...
        FRProfileMeanDownBad.percChangeBefRunVsBL,...
        ['average FR change BefRun/BL G/B Down'],...
        ['Pyr_FRChangeBefRun-BLDownGoodBadBox'],...
        pathAnal,[-1 4],FRProfileMeanStatDownGoodBad.pRSPercChangeBefRunVsBLAll,colorSel);

    plotBoxPlot(FRProfileMeanDown.percChangeBefRunVs0to1,...
        FRProfileMeanDownBad.percChangeBefRunVs0to1,...
        ['average FR change BefRun/0to1s G/B Down'],...
        ['Pyr_FRChangeBefRun-0to1DownGoodBadBox'],...
        pathAnal,[-1 3],FRProfileMeanStatDownGoodBad.pRSPercChangeBefRunVs0to1All,colorSel);
    
    plotBoxPlot(FRProfileMeanDown.percChange0to1Vs3to5,...
        FRProfileMeanDownBad.percChange0to1Vs3to5,...
        ['average FR change 0-1s/3-5s G/B Down'],...
        ['Pyr_FRChange0to1-3to5DownGoodBadBox'],...
        pathAnal,[-1 3],FRProfileMeanStatDownGoodBad.pRSPercChange0to1Vs3to5All,colorSel);

    plotBoxPlot(FRProfileMeanDown.percChangeBefRunVs3to5,...
        FRProfileMeanDownBad.percChangeBefRunVs3to5,...
        ['average FR change BefRun/3-5s G/B Down'],...
        ['Pyr_FRChangeBefRun-3to5DownGoodBadBox'],...
        pathAnal,[-1 4],FRProfileMeanStatDownGoodBad.pRSPercChangeBefRunVs3to5All,colorSel);
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

function FRProfileMeanStat = accumMeanStatC(FRProfileMean)
      
    if(isempty(FRProfileMean.meanAvgFRProfileBaseline))
        FRProfileMeanStat = [];
        return;
    end
    FRProfileMeanStat.pRS0to1VsBL = ranksum(FRProfileMean.meanAvgFRProfileBaseline,...
                FRProfileMean.meanAvgFRProfile0to1);    
    FRProfileMeanStat.pRSBefRunVsBL = ranksum(FRProfileMean.meanAvgFRProfileBaseline,...
                FRProfileMean.meanAvgFRProfileBefRun);    
    FRProfileMeanStat.pRS3to5VsBL = ranksum(FRProfileMean.meanAvgFRProfileBaseline,...
                FRProfileMean.meanAvgFRProfile3to5);
    FRProfileMeanStat.pRSBefRunVs0to1 = ranksum(FRProfileMean.meanAvgFRProfile0to1,...
                FRProfileMean.meanAvgFRProfileBefRun);
    FRProfileMeanStat.pRS3to5Vs0to1 = ranksum(FRProfileMean.meanAvgFRProfile0to1,...
                FRProfileMean.meanAvgFRProfile3to5);
    FRProfileMeanStat.pRS3to5VsBefRun = ranksum(FRProfileMean.meanAvgFRProfileBefRun,...
                FRProfileMean.meanAvgFRProfile3to5);

    FRProfileMeanStat.pTTPercChange0to1VsBL = ttest(FRProfileMean.percChange0to1VsBL);
    FRProfileMeanStat.pTTPercChangeBefRunVsBL = ttest(FRProfileMean.percChangeBefRunVsBL);
    FRProfileMeanStat.pTTPercChangeBefRunVs0to1 = ttest(FRProfileMean.percChangeBefRunVs0to1);
    FRProfileMeanStat.pTTPercChange0to1Vs3to5 = ttest(FRProfileMean.percChange0to1Vs3to5);
    FRProfileMeanStat.pTTPercChangeBefRunVs3to5 = ttest(FRProfileMean.percChangeBefRunVs3to5);
    
end

function FRProfileMeanStatC = accumMeanStatCGoodBad(FRProfileMean,FRProfileMeanBad)
        
    if(isempty(FRProfileMeanBad.meanAvgFRProfileBaseline))
        FRProfileMeanStatC = [];
        return;
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
%     set(gca,'XLim',[timeStepRun(1) 7]);
    set(gca,'XLim',[-1 4]);
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
    elseif(ordMethod == 5)
        indTmp = timeStepRun > 0;
        [~,indMax] = max(avgFRProfile(:,indTmp)');
    end
    if(ordMethod == 4 | ordMethod == 5)
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

function plotMeanFRProfilePyrVsInt(timeStepRun,avgFRProfilex,avgFRProfiley,yl,fileName,pathAnal,ylimit)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 200]);
    
    h = plot(timeStepRun,avgFRProfilex);
    set(h,'LineWidth',1, 'Color',[ 39 169 225]./255);
    hold on;
    h = plot(timeStepRun,avgFRProfiley);
    set(h,'LineWidth',1, 'Color',[167 169  171]./255);
    
    if(~isempty(ylimit))
        h = plot([0 0],ylimit,'r-');
    else
        h = plot([0 0],[0 1],'r-');
    end
    set(h,'LineWidth',1)
    set(gca,'XLim',[-1 4]);
    
    xlabel('Time (s)')
    ylabel(yl)
    
    fileName1 = [pathAnal fileName];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end