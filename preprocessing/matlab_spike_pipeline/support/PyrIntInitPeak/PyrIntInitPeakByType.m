function InitAll = PyrIntInitPeakByType(taskSel,modPyr1NoCue,modPyr1AL,modPyr1PL,...
                                modInt1NoCue,modInt1AL,modInt1PL)

    if(taskSel == 1)   
        InitAll.task = [modPyr1NoCue.taskGood modPyr1AL.taskGood modPyr1PL.taskGood];
        InitAll.taskBad = [modPyr1NoCue.taskBad modPyr1AL.taskBad modPyr1PL.taskBad];
        InitAll.indRec = [modPyr1NoCue.indRecGood modPyr1AL.indRecGood modPyr1PL.indRecGood];
        InitAll.indRecBad = [modPyr1NoCue.indRecBad modPyr1AL.indRecBad modPyr1PL.indRecBad];
        InitAll.indNeu = [modPyr1NoCue.indNeuGood modPyr1AL.indNeuGood modPyr1PL.indNeuGood];
        InitAll.indNeuBad = [modPyr1NoCue.indNeuBad modPyr1AL.indNeuBad modPyr1PL.indNeuBad];
        InitAll.idxC2 = [modPyr1NoCue.idxC2Good modPyr1AL.idxC2Good modPyr1PL.idxC2Good];
        InitAll.idxC2Bad = [modPyr1NoCue.idxC2Bad modPyr1AL.idxC2Bad modPyr1PL.idxC2Bad];
        InitAll.avgFRProfile = [modPyr1NoCue.avgFRProfileGood; modPyr1AL.avgFRProfileGood; modPyr1PL.avgFRProfileGood];
        InitAll.avgFRProfileBad = [modPyr1NoCue.avgFRProfileBad; modPyr1AL.avgFRProfileBad; modPyr1PL.avgFRProfileBad];
        InitAll.relDepthNeuHDef = [modPyr1NoCue.relDepthNeuHDefGood modPyr1AL.relDepthNeuHDefGood modPyr1PL.relDepthNeuHDefGood];
        InitAll.isNeuWithFieldAligned = [modPyr1NoCue.isNeuWithFieldAlignedGood modPyr1AL.isNeuWithFieldAlignedGood...
            modPyr1PL.isNeuWithFieldAlignedGood];
        InitAll.taskInt = [modInt1NoCue.taskGood modInt1AL.taskGood modInt1PL.taskGood];
        InitAll.taskIntBad = [modInt1NoCue.taskBad modInt1AL.taskBad modInt1PL.taskBad];
        InitAll.indRecInt = [modInt1NoCue.indRecGood modInt1AL.indRecGood modInt1PL.indRecGood];
        InitAll.indRecIntBad = [modInt1NoCue.indRecBad modInt1AL.indRecBad modInt1PL.indRecBad];
        InitAll.indNeuInt = [modInt1NoCue.indNeuGood modInt1AL.indNeuGood modInt1PL.indNeuGood];
        InitAll.indNeuIntBad = [modInt1NoCue.indNeuBad modInt1AL.indNeuBad modInt1PL.indNeuBad];
        InitAll.idxC2Int = [modInt1NoCue.idxC2Good modInt1AL.idxC2Good modInt1PL.idxC2Good];
        InitAll.idxC2BadInt = [modInt1NoCue.idxC2Bad modInt1AL.idxC2Bad modInt1PL.idxC2Bad];
        InitAll.avgFRProfileInt = [modInt1NoCue.avgFRProfileGood; modInt1AL.avgFRProfileGood; modInt1PL.avgFRProfileGood];
        InitAll.avgFRProfileBadInt = [modInt1NoCue.avgFRProfileBad; modInt1AL.avgFRProfileBad; modInt1PL.avgFRProfileBad];
    elseif(taskSel == 2)
        InitAll.task = [modPyr1AL.taskGood modPyr1PL.taskGood];
        InitAll.taskBad = [modPyr1AL.taskBad modPyr1PL.taskBad];
        InitAll.indRec = [modPyr1AL.indRecGood modPyr1PL.indRecGood];
        InitAll.indRecBad = [modPyr1AL.indRecBad modPyr1PL.indRecBad];
        InitAll.indNeu = [modPyr1AL.indNeuGood modPyr1PL.indNeuGood];
        InitAll.indNeuBad = [modPyr1AL.indNeuBad modPyr1PL.indNeuBad];
        InitAll.idxC2 = [modPyr1AL.idxC2Good modPyr1PL.idxC2Good];
        InitAll.idxC2Bad = [modPyr1AL.idxC2Bad modPyr1PL.idxC2Bad];
        InitAll.avgFRProfile = [modPyr1AL.avgFRProfileGood; modPyr1PL.avgFRProfileGood];
        InitAll.avgFRProfileBad = [modPyr1AL.avgFRProfileBad; modPyr1PL.avgFRProfileBad];
        InitAll.relDepthNeuHDef = [modPyr1AL.relDepthNeuHDefGood modPyr1PL.relDepthNeuHDefGood];
        InitAll.isNeuWithFieldAligned = [modPyr1AL.isNeuWithFieldAlignedGood...
            modPyr1PL.isNeuWithFieldAlignedGood];
        InitAll.taskInt = [modInt1AL.taskGood modInt1PL.taskGood];
        InitAll.taskIntBad = [modInt1AL.taskBad modInt1PL.taskBad];
        InitAll.indRecInt = [modInt1AL.indRecGood modInt1PL.indRecGood];
        InitAll.indRecIntBad = [modInt1AL.indRecBad modInt1PL.indRecBad];
        InitAll.indNeuInt = [modInt1AL.indNeuGood modInt1PL.indNeuGood];
        InitAll.indNeuIntBad = [modInt1AL.indNeuBad modInt1PL.indNeuBad];
        InitAll.idxC2Int = [modInt1AL.idxC2Good modInt1PL.idxC2Good];
        InitAll.idxC2BadInt = [modInt1AL.idxC2Bad modInt1PL.idxC2Bad];
        InitAll.avgFRProfileInt = [modInt1AL.avgFRProfileGood; modInt1PL.avgFRProfileGood];
        InitAll.avgFRProfileBadInt = [modInt1AL.avgFRProfileBad; modInt1PL.avgFRProfileBad];
    elseif(taskSel == 3)
        InitAll.task = modPyr1AL.taskGood;
        InitAll.taskBad = modPyr1AL.taskBad;
        InitAll.indRec = modPyr1AL.indRecGood;
        InitAll.indRecBad = modPyr1AL.indRecBad;
        InitAll.indNeu = modPyr1AL.indNeuGood;
        InitAll.indNeuBad = modPyr1AL.indNeuBad;
        InitAll.idxC2 = modPyr1AL.idxC2Good;
        InitAll.idxC2Bad = modPyr1AL.idxC2Bad;
        InitAll.avgFRProfile = modPyr1AL.avgFRProfileGood;
        InitAll.avgFRProfileBad = modPyr1AL.avgFRProfileBad;
        InitAll.relDepthNeuHDef = modPyr1AL.relDepthNeuHDefGood;
        InitAll.isNeuWithFieldAligned = modPyr1AL.isNeuWithFieldAlignedGood;
        InitAll.taskInt = modInt1AL.taskGood;
        InitAll.taskIntBad = modInt1AL.taskBad;
        InitAll.indRecInt = modInt1AL.indRecGood;
        InitAll.indRecIntBad = modInt1AL.indRecBad;
        InitAll.indNeuInt = modInt1AL.indNeuGood;
        InitAll.indNeuIntBad = modInt1AL.indNeuBad;
        InitAll.idxC2Int = modInt1AL.idxC2Good;
        InitAll.idxC2BadInt = modInt1AL.idxC2Bad;
        InitAll.avgFRProfileInt = modInt1AL.avgFRProfileGood;
        InitAll.avgFRProfileBadInt = modInt1AL.avgFRProfileBad;    
    elseif(taskSel == 4)
        InitAll.task = modPyr1NoCue.taskGood;
        InitAll.taskBad = modPyr1NoCue.taskBad;
        InitAll.indRec = modPyr1NoCue.indRecGood;
        InitAll.indRecBad = modPyr1NoCue.indRecBad;
        InitAll.indNeu = modPyr1NoCue.indNeuGood;
        InitAll.indNeuBad = modPyr1NoCue.indNeuBad;
        InitAll.idxC2 = modPyr1NoCue.idxC2Good;
        InitAll.idxC2Bad = modPyr1NoCue.idxC2Bad;
        InitAll.avgFRProfile = modPyr1NoCue.avgFRProfileGood;
        InitAll.avgFRProfileBad = modPyr1NoCue.avgFRProfileBad;
        InitAll.relDepthNeuHDef = modPyr1NoCue.relDepthNeuHDefGood;
        InitAll.isNeuWithFieldAligned = modPyr1NoCue.isNeuWithFieldAlignedGood;
        InitAll.taskInt = modInt1NoCue.taskGood;
        InitAll.taskIntBad = modInt1NoCue.taskBad;
        InitAll.indRecInt = modInt1NoCue.indRecGood;
        InitAll.indRecIntBad = modInt1NoCue.indRecBad;
        InitAll.indNeuInt = modInt1NoCue.indNeuGood;
        InitAll.indNeuIntBad = modInt1NoCue.indNeuBad;
        InitAll.idxC2Int = modInt1NoCue.idxC2Good;
        InitAll.idxC2BadInt = modInt1NoCue.idxC2Bad;
        InitAll.avgFRProfileInt = modInt1NoCue.avgFRProfileGood;
        InitAll.avgFRProfileBadInt = modInt1NoCue.avgFRProfileBad; 
    else
        ind = modPyr1AL.taskGood == 2 & modPyr1AL.indRecGood > 0 & ...
            modPyr1AL.indRecGood <= 10;
        indBad = modPyr1AL.taskBad == 2 & modPyr1AL.indRecBad > 0 & ...
            modPyr1AL.indRecBad <= 10;
        
        indInt = modInt1AL.taskGood == 2 & modInt1AL.indRecGood > 0 & ...
            modInt1AL.indRecGood <= 10;
        indIntBad = modInt1AL.taskBad == 2 & modInt1AL.indRecBad > 0 & ...
            modInt1AL.indRecBad <= 10;
        
        InitAll.task = modPyr1AL.taskGood(ind);
        InitAll.taskBad = modPyr1AL.taskBad(indBad);
        InitAll.indRec = modPyr1AL.indRecGood(ind);
        InitAll.indRecBad = modPyr1AL.indRecBad(indBad);
        InitAll.indNeu = modPyr1AL.indNeuGood(ind);
        InitAll.indNeuBad = modPyr1AL.indNeuBad(indBad);
        InitAll.idxC2 = modPyr1AL.idxC2Good(ind);
        InitAll.idxC2Bad = modPyr1AL.idxC2Bad(indBad);
        InitAll.avgFRProfile = modPyr1AL.avgFRProfileGood(ind,:);
        InitAll.avgFRProfileBad = modPyr1AL.avgFRProfileBad(indBad,:);
        InitAll.relDepthNeuHDef = modPyr1AL.relDepthNeuHDefGood(ind);
        InitAll.isNeuWithFieldAligned = modPyr1AL.isNeuWithFieldAlignedGood(ind);
        InitAll.taskInt = modInt1AL.taskGood(indInt);
        InitAll.taskIntBad = modInt1AL.taskBad(indIntBad);
        InitAll.indRecInt = modInt1AL.indRecGood(indInt);
        InitAll.indRecIntBad = modInt1AL.indRecBad(indIntBad);
        InitAll.indNeuInt = modInt1AL.indNeuGood(indInt);
        InitAll.indNeuIntBad = modInt1AL.indNeuBad(indIntBad);
        InitAll.idxC2Int = modInt1AL.idxC2Good(indInt);
        InitAll.idxC2BadInt = modInt1AL.idxC2Bad(indIntBad);
        InitAll.avgFRProfileInt = modInt1AL.avgFRProfileGood(indInt,:);
        InitAll.avgFRProfileBadInt = modInt1AL.avgFRProfileBad(indIntBad,:); 
    end
    
    InitAll.idxIsNeuWithFieldAlignedBad = [];
    InitAll.isNeuWithFieldAlignedBad = zeros(1,length(InitAll.taskBad));
    for i = 1:length(InitAll.taskBad)
        idxBad = find(InitAll.task == InitAll.taskBad(i) & InitAll.indRec == InitAll.indRecBad(i)...
            & InitAll.indNeu == InitAll.indNeuBad(i));
        if(length(idxBad) == 1)
            InitAll.isNeuWithFieldAlignedBad(i) = InitAll.isNeuWithFieldAligned(idxBad);
            InitAll.idxIsNeuWithFieldAlignedBad = [InitAll.idxIsNeuWithFieldAlignedBad idxBad];
        end
    end
    
    %% normalization
    indTmp = modPyr1NoCue.timeStepRun >=-1 & modPyr1NoCue.timeStepRun <= 4 ;
    InitAll.avgFRProfileNorm = zeros(size(InitAll.avgFRProfile,1),size(InitAll.avgFRProfile,2));
    for i = 1:size(InitAll.avgFRProfile,1)
        if(max(InitAll.avgFRProfile(i,:)) ~= 0)
            InitAll.avgFRProfileNorm(i,:) = InitAll.avgFRProfile(i,:)/max(InitAll.avgFRProfile(i,:));
        end
    end
    n = 0;
    InitAll.avgFRProfileNorm1 = zeros(size(InitAll.avgFRProfile,1),size(InitAll.avgFRProfile,2));
    for i = 1:size(InitAll.avgFRProfile,1)
        if(max(InitAll.avgFRProfile(i,:)) ~= 0)
            rangeFR = max(InitAll.avgFRProfile(i,indTmp))-min(InitAll.avgFRProfile(i,indTmp));
            if(rangeFR ~= 0)
                InitAll.avgFRProfileNorm1(i,:) = (InitAll.avgFRProfile(i,:)-min(InitAll.avgFRProfile(i,indTmp)))...
                    /rangeFR;
            else
                n = n+1;
            end
        end
    end
    InitAll.avgFRProfileNormBad = zeros(size(InitAll.avgFRProfileBad,1),size(InitAll.avgFRProfileBad,2));
    for i = 1:size(InitAll.avgFRProfileBad,1)
        if(max(InitAll.avgFRProfileBad(i,:)) ~= 0)
            InitAll.avgFRProfileNormBad(i,:) = InitAll.avgFRProfileBad(i,:)/max(InitAll.avgFRProfileBad(i,:));
        end
    end
    
    InitAll.avgFRProfileNormInt = zeros(size(InitAll.avgFRProfileInt,1),size(InitAll.avgFRProfileInt,2));
    for i = 1:size(InitAll.avgFRProfileInt,1)
        if(max(InitAll.avgFRProfileInt(i,:)) ~= 0)
            InitAll.avgFRProfileNormInt(i,:) = InitAll.avgFRProfileInt(i,:)/max(InitAll.avgFRProfileInt(i,:));
        end
    end
    InitAll.avgFRProfileNormInt1 = zeros(size(InitAll.avgFRProfileInt,1),size(InitAll.avgFRProfileInt,2));
    m = 0;
    for i = 1:size(InitAll.avgFRProfileInt,1)
        if(max(InitAll.avgFRProfileInt(i,:)) ~= 0)
            rangeFR = max(InitAll.avgFRProfileInt(i,indTmp))-min(InitAll.avgFRProfileInt(i,indTmp));
            if(rangeFR ~= 0)
                InitAll.avgFRProfileNormInt1(i,:) = (InitAll.avgFRProfileInt(i,:)-min(InitAll.avgFRProfileInt(i,indTmp)))...
                    /rangeFR;
            else
                m = m+1;
            end
        end
    end
    InitAll.avgFRProfileNormBadInt = zeros(size(InitAll.avgFRProfileBadInt,1),size(InitAll.avgFRProfileBadInt,2));
    for i = 1:size(InitAll.avgFRProfileBadInt,1)
        if(max(InitAll.avgFRProfileBadInt(i,:)) ~= 0)
            InitAll.avgFRProfileNormBadInt(i,:) = InitAll.avgFRProfileBadInt(i,:)/max(InitAll.avgFRProfileBadInt(i,:));
        end
    end
end