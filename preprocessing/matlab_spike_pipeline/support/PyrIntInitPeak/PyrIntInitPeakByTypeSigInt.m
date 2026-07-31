function InitAll = PyrIntInitPeakByTypeSigInt(taskSel,modInt1NoCue,modInt1AL,modInt1PL,...
                                modInt1SigNoCue,modInt1SigAL,modInt1SigPL)

    if(taskSel == 1)
        if(sum(modInt1NoCue.indNeuGood-modInt1SigNoCue.indNeu) == 0 && ...
                sum(modInt1AL.indNeuGood-modInt1SigAL.indNeu) == 0 && ...
                sum(modInt1PL.indNeuGood-modInt1SigPL.indNeu) == 0)
            InitAll.task = [modInt1NoCue.taskGood modInt1AL.taskGood modInt1PL.taskGood];
            InitAll.taskBad = [modInt1NoCue.taskBad modInt1AL.taskBad modInt1PL.taskBad];
            InitAll.indRec = [modInt1NoCue.indRecGood modInt1AL.indRecGood modInt1PL.indRecGood];
            InitAll.indRecBad = [modInt1NoCue.indRecBad modInt1AL.indRecBad modInt1PL.indRecBad];
            InitAll.indNeu = [modInt1NoCue.indNeuGood modInt1AL.indNeuGood modInt1PL.indNeuGood];
            InitAll.indNeuBad = [modInt1NoCue.indNeuBad modInt1AL.indNeuBad modInt1PL.indNeuBad];
            InitAll.idxC2 = [modInt1NoCue.idxC2Good modInt1AL.idxC2Good modInt1PL.idxC2Good];
            InitAll.idxC2Bad = [modInt1NoCue.idxC2Bad modInt1AL.idxC2Bad modInt1PL.idxC2Bad];
            InitAll.avgFRProfile = [modInt1NoCue.avgFRProfileGood; modInt1AL.avgFRProfileGood; modInt1PL.avgFRProfileGood];
            InitAll.avgFRProfileBad = [modInt1NoCue.avgFRProfileBad; modInt1AL.avgFRProfileBad; modInt1PL.avgFRProfileBad];
            InitAll.isPeakNeuArrAll = [modInt1SigNoCue.isPeakNeuArrAll' modInt1SigAL.isPeakNeuArrAll' ...
                modInt1SigPL.isPeakNeuArrAll']; 
        else
            disp('modInt1 and modInt1Sig do not have equal number of neurons');
        end
    elseif(taskSel == 2)
        if(sum(modInt1AL.indNeuGood-modInt1SigAL.indNeu) == 0 && ...
                sum(modInt1PL.indNeuGood-modInt1SigPL.indNeu) == 0) 
            InitAll.task = [modInt1AL.taskGood modInt1PL.taskGood];
            InitAll.taskBad = [modInt1AL.taskBad modInt1PL.taskBad];
            InitAll.indRec = [modInt1AL.indRecGood modInt1PL.indRecGood];
            InitAll.indRecBad = [modInt1AL.indRecBad modInt1PL.indRecBad];
            InitAll.indNeu = [modInt1AL.indNeuGood modInt1PL.indNeuGood];
            InitAll.indNeuBad = [modInt1AL.indNeuBad modInt1PL.indNeuBad];
            InitAll.idxC2 = [modInt1AL.idxC2Good modInt1PL.idxC2Good];
            InitAll.idxC2Bad = [modInt1AL.idxC2Bad modInt1PL.idxC2Bad];
            InitAll.avgFRProfile = [modInt1AL.avgFRProfileGood; modInt1PL.avgFRProfileGood];
            InitAll.avgFRProfileBad = [modInt1AL.avgFRProfileBad; modInt1PL.avgFRProfileBad];
            InitAll.isPeakNeuArrAll = [modInt1SigAL.isPeakNeuArrAll' ...
                modInt1SigPL.isPeakNeuArrAll']; 
        else
            disp('modInt1 and modInt1Sig do not have equal number of neurons');
        end
    else
        if(sum(modInt1AL.indNeuGood-modInt1SigAL.indNeu) == 0) 
            InitAll.task = modInt1AL.taskGood;
            InitAll.taskBad = modInt1AL.taskBad;
            InitAll.indRec = modInt1AL.indRecGood;
            InitAll.indRecBad = modInt1AL.indRecBad;
            InitAll.indNeu = modInt1AL.indNeuGood;
            InitAll.indNeuBad = modInt1AL.indNeuBad;
            InitAll.idxC2 = modInt1AL.idxC2Good;
            InitAll.idxC2Bad = modInt1AL.idxC2Bad;
            InitAll.avgFRProfile = modInt1AL.avgFRProfileGood;
            InitAll.avgFRProfileBad = modInt1AL.avgFRProfileBad;      
            InitAll.isPeakNeuArrAll = modInt1SigAL.isPeakNeuArrAll'; 
        else
            disp('modInt1 and modInt1Sig do not have equal number of neurons');
        end
    end
    
    %% normalization
    indTmp = modInt1NoCue.timeStepRun >=-1 & modInt1NoCue.timeStepRun <= 4 ;
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