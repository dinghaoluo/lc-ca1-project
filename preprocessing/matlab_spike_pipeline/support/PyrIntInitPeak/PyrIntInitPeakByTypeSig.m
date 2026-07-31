function InitAll = PyrIntInitPeakByTypeSig(taskSel,modPyr1NoCue,modPyr1AL,modPyr1PL,...
                                modPyr1SigNoCue,modPyr1SigAL,modPyr1SigPL)

    if(taskSel == 1)
        if(sum(modPyr1NoCue.indNeuGood-modPyr1SigNoCue.indNeu) == 0 && ...
                sum(modPyr1AL.indNeuGood-modPyr1SigAL.indNeu) == 0 && ...
                sum(modPyr1PL.indNeuGood-modPyr1SigPL.indNeu) == 0)
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
            InitAll.isPeakNeuArrAll = [modPyr1SigNoCue.isPeakNeuArrAll' modPyr1SigAL.isPeakNeuArrAll' ...
                modPyr1SigPL.isPeakNeuArrAll']; 
        else
            disp('modPyr1 and modPyr1Sig do not have equal number of neurons');
        end
    elseif(taskSel == 2)
        if(sum(modPyr1AL.indNeuGood-modPyr1SigAL.indNeu) == 0 && ...
                sum(modPyr1PL.indNeuGood-modPyr1SigPL.indNeu) == 0)
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
            InitAll.isPeakNeuArrAll = [modPyr1SigAL.isPeakNeuArrAll' ...
                modPyr1SigPL.isPeakNeuArrAll']; 
        else
            disp('modPyr1 and modPyr1Sig do not have equal number of neurons');
        end            
    elseif(taskSel == 3)
        if(sum(modPyr1AL.indNeuGood-modPyr1SigAL.indNeu) == 0)
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
            InitAll.isPeakNeuArrAll = modPyr1SigAL.isPeakNeuArrAll'; 
        else
            disp('modPyr1 and modPyr1Sig do not have equal number of neurons');
        end
    elseif(taskSel == 4)
        if(sum(modPyr1NoCue.indNeuGood-modPyr1SigNoCue.indNeu) == 0)
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
            InitAll.isPeakNeuArrAll = modPyr1SigNoCue.isPeakNeuArrAll'; 
        end
    elseif(taskSel == 5)
        if(sum(modPyr1AL.indNeuGood-modPyr1SigAL.indNeu) == 0)
            ind = modPyr1AL.taskGood == 2 & modPyr1AL.indRecGood > 0 & ...
                modPyr1AL.indRecGood <= 10;
            indBad = modPyr1AL.taskBad == 2 & modPyr1AL.indRecBad > 0 & ...
                modPyr1AL.indRecBad <= 10;
        
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
            InitAll.isPeakNeuArrAll = modPyr1SigAL.isPeakNeuArrAll(ind,:)'; 
        end
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
    InitAll.avgFRProfileNormZ = zscore(InitAll.avgFRProfile,0,2); 
    
    InitAll.avgFRProfileNormBad = zeros(size(InitAll.avgFRProfileBad,1),size(InitAll.avgFRProfileBad,2));
    for i = 1:size(InitAll.avgFRProfileBad,1)
        if(max(InitAll.avgFRProfileBad(i,:)) ~= 0)
            InitAll.avgFRProfileNormBad(i,:) = InitAll.avgFRProfileBad(i,:)/max(InitAll.avgFRProfileBad(i,:));
        end
    end
    InitAll.avgFRProfileNormZBad = zscore(InitAll.avgFRProfileBad,0,2); 
    
end