function PyrIntInitPeakPVStim()
% compare Pyramidal neurons and PV interneurons on their initial peak
    
    pathAnal0 = 'Z:\Yingxue\DataAnalysisRaphi\PyramidalALPVStim\';
    
    if(exist([pathAnal0 'initPeakPyrAllRecStim.mat']))
        load([pathAnal0 'initPeakPyrAllRecStim.mat']);
    end
    if(exist([pathAnal0 'initPeakPyrAllRecStim_km2.mat']))
        load([pathAnal0 'initPeakPyrAllRecStim_km2.mat']);
    end
    
    if(exist([pathAnal0 'initPeakPyrIntAllRecStim.mat']))
        load([pathAnal0 'initPeakPyrIntAllRecStim.mat']);
    end
    
    avgFRProfile = modPyr1AL.avgFRProfile; 
    avgFRProfileStim = modPyr1AL.avgFRProfileStim;
    avgFRProfileStimCtrl = modPyr1AL.avgFRProfileStimCtrl;
    isFieldGoodNonStim = modPyr1AL.isNeuWithFieldAlignedGoodNonStim;
    isFieldGoodNonStimAndStim = (modPyr1AL.isNeuWithFieldAlignedGoodNonStim | ...
        modPyr1AL.isNeuWithFieldAlignedStim);
                
    avgFRProfileNorm = zeros(size(avgFRProfile,1),size(avgFRProfile,2));
    for i = 1:size(avgFRProfile,1)
        if(max(avgFRProfile(i,:)) ~= 0)
            avgFRProfileNorm(i,:) = avgFRProfile(i,:)/max(avgFRProfile(i,:));
        end
    end
    avgFRProfileNormStim = zeros(size(avgFRProfileStim,1),size(avgFRProfileStim,2));
    for i = 1:size(avgFRProfileStim,1)
        if(max(avgFRProfileStim(i,:)) ~= 0)
            avgFRProfileNormStim(i,:) = avgFRProfileStim(i,:)/max(avgFRProfileStim(i,:));
        end
    end
    avgFRProfileNormStimCtrl = zeros(size(avgFRProfileStimCtrl,1),size(avgFRProfileStimCtrl,2));
    for i = 1:size(avgFRProfileStimCtrl,1)
        if(max(avgFRProfileStimCtrl(i,:)) ~= 0)
            avgFRProfileNormStimCtrl(i,:) = avgFRProfileStimCtrl(i,:)/max(avgFRProfileStimCtrl(i,:));
        end
    end
    
    mean0to1 = mean(avgFRProfile(:,FRProfileMeanAll.indFR0to1),2);
    meanBefRun = mean(avgFRProfile(:,FRProfileMeanAll.indFRBefRun),2);
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
    PyrRise.pulseMethod = modPyr1AL.pulseMeth(indOrdTmp);
    PyrRise.actOrInact = modPyr1AL.actOrInact(indOrdTmp);
    PyrRise.indRec = modPyr1AL.indRec(indOrdTmp);
    PyrRise.indNeu = modPyr1AL.indNeu(indOrdTmp);
   
    %% neurons with FR decrease around 0
    indOrdTmp = indOrd(idx1:end);
    idxNanTmp = idxNan(idx1:end);
    indOrdTmp = indOrdTmp(idxNanTmp == 0);
    PyrDown.idxDown = indOrdTmp;
    PyrDown.pulseMethod = modPyr1AL.pulseMeth(indOrdTmp);
    PyrDown.actOrInact = modPyr1AL.actOrInact(indOrdTmp);
    PyrDown.indRec = modPyr1AL.indRec(indOrdTmp);
    PyrDown.indNeu = modPyr1AL.indNeu(indOrdTmp);
    
    actOrInActStr{1} = 'act';
    actOrInActStr{2} = 'inact';
    
    % activation or inactivation
    cond = 0;
    for i = 1: 2
        % different pulse methods
        for j = 1:length(pulseMethod{i})
            cond = cond + 1;
            ind = find(PyrRise.actOrInact == i & PyrRise.pulseMethod == pulseMethod{i}(j));
            FRProfile1{cond}.actOrInact = i;
            FRProfile1{cond}.pulseMethod = pulseMethod{i}(j);
            FRProfile1{cond}.indPyrRise = PyrRise.idxRise(ind);
            FRProfile1{cond}.indRecPyrRise = PyrRise.indRec(ind);
            FRProfile1{cond}.indNeuPyrRise = PyrRise.indNeu(ind);
            ind = find(PyrDown.actOrInact == i & PyrDown.pulseMethod == pulseMethod{i}(j));
            FRProfile1{cond}.indPyrDown = PyrDown.idxDown(ind);
            FRProfile1{cond}.indRecPyrDown = PyrDown.indRec(ind);
            FRProfile1{cond}.indNeuPyrDown = PyrDown.indNeu(ind);
            ind = find(modPyr1AL.actOrInact == i & modPyr1AL.pulseMeth == pulseMethod{i}(j));
            FRProfile1{cond}.ind = ind;
            FRProfile1{cond}.isField = isFieldGoodNonStim(ind);
            FRProfile1{cond}.isFieldComb = isFieldGoodNonStimAndStim(ind);
            FRProfile1{cond}.indField = ind(FRProfile1{cond}.isField == 1);
            FRProfile1{cond}.indFieldComb = ind(FRProfile1{cond}.isFieldComb ==1);
                        
            %% Pyramidal neurons
            FRProfileMean1{cond} = accumMean(avgFRProfile(FRProfile1{cond}.ind,:),modPyr1AL.timeStepRun);
            FRProfileMeanStim1{cond} = accumMean(avgFRProfileStim(FRProfile1{cond}.ind,:),modPyr1AL.timeStepRun);
            FRProfileMeanStimCtrl1{cond} = accumMean(avgFRProfileStimCtrl(FRProfile1{cond}.ind,:),modPyr1AL.timeStepRun);
            
            %% Pyramidal neurons fields
            FRProfileMeanField1{cond} = accumMean(avgFRProfile(FRProfile1{cond}.indField,:),modPyr1AL.timeStepRun);
            FRProfileMeanStimField1{cond} = accumMean(avgFRProfileStim(FRProfile1{cond}.indField,:),modPyr1AL.timeStepRun);
            FRProfileMeanStimCtrlField1{cond} = accumMean(avgFRProfileStimCtrl(FRProfile1{cond}.indField,:),modPyr1AL.timeStepRun);
            
            FRProfileMeanFieldComb{cond} = accumMean(avgFRProfile(FRProfile1{cond}.indFieldComb,:),modPyr1AL.timeStepRun);
            FRProfileMeanStimFieldComb{cond} = accumMean(avgFRProfileStim(FRProfile1{cond}.indFieldComb,:),modPyr1AL.timeStepRun);
            FRProfileMeanStimCtrlFieldComb{cond} = accumMean(avgFRProfileStimCtrl(FRProfile1{cond}.indFieldComb,:),modPyr1AL.timeStepRun);
            
            %% pyramidal neurons rise and down
            FRProfileMeanRise{cond} = accumMean(avgFRProfile(FRProfile1{cond}.indPyrRise,:),modPyr1AL.timeStepRun);
            FRProfileMeanDown{cond} = accumMean(avgFRProfile(FRProfile1{cond}.indPyrDown,:),modPyr1AL.timeStepRun);
            FRProfileMeanRiseStim{cond} = accumMean(avgFRProfileStim(FRProfile1{cond}.indPyrRise,:),modPyr1AL.timeStepRun);
            FRProfileMeanDownStim{cond} = accumMean(avgFRProfileStim(FRProfile1{cond}.indPyrDown,:),modPyr1AL.timeStepRun);
            FRProfileMeanRiseStimCtrl{cond} = accumMean(avgFRProfileStimCtrl(FRProfile1{cond}.indPyrRise,:),modPyr1AL.timeStepRun);
            FRProfileMeanDownStimCtrl{cond} = accumMean(avgFRProfileStimCtrl(FRProfile1{cond}.indPyrDown,:),modPyr1AL.timeStepRun);
            
            %% is field or not
            isNeuWithFieldAlignedRise{cond}.isField = isFieldGoodNonStim(FRProfile1{cond}.indPyrRise);
            isNeuWithFieldAlignedRise{cond}.isFieldComb = isFieldGoodNonStimAndStim(FRProfile1{cond}.indPyrRise);
            isNeuWithFieldAlignedDown{cond}.isField = isFieldGoodNonStim(FRProfile1{cond}.indPyrDown);
            isNeuWithFieldAlignedDown{cond}.isFieldComb = isFieldGoodNonStimAndStim(FRProfile1{cond}.indPyrDown);
            
            isNeuWithFieldAlignedRise{cond}.numField = sum(isNeuWithFieldAlignedRise{cond}.isField);
            isNeuWithFieldAlignedRise{cond}.numFieldComb = sum(isNeuWithFieldAlignedRise{cond}.isFieldComb);
            isNeuWithFieldAlignedDown{cond}.numField = sum(isNeuWithFieldAlignedDown{cond}.isField);
            isNeuWithFieldAlignedDown{cond}.numFieldComb = sum(isNeuWithFieldAlignedDown{cond}.isFieldComb);
            
            isNeuWithFieldAlignedRise{cond}.idxRise = FRProfile1{cond}.indPyrRise(isNeuWithFieldAlignedRise{cond}.isField == 1);
            isNeuWithFieldAlignedRise{cond}.idxRiseComb = FRProfile1{cond}.indPyrRise(isNeuWithFieldAlignedRise{cond}.isFieldComb == 1);
            isNeuWithFieldAlignedDown{cond}.idxDown = FRProfile1{cond}.indPyrDown(isNeuWithFieldAlignedDown{cond}.isField == 1);            
            isNeuWithFieldAlignedDown{cond}.idxDownComb = FRProfile1{cond}.indPyrDown(isNeuWithFieldAlignedDown{cond}.isFieldComb == 1);
            
            %% fields from good non-stim trials
            FRProfileMeanRiseField{cond} = accumMean(avgFRProfile(isNeuWithFieldAlignedRise{cond}.idxRise,:),...
                modPyr1AL.timeStepRun);
            FRProfileMeanDownField{cond} = accumMean(avgFRProfile(isNeuWithFieldAlignedDown{cond}.idxDown,:),...
                modPyr1AL.timeStepRun);
            FRProfileMeanRiseFieldStim{cond} = accumMean(avgFRProfileStim(isNeuWithFieldAlignedRise{cond}.idxRise,:),...
                modPyr1AL.timeStepRun);
            FRProfileMeanDownFieldStim{cond} = accumMean(avgFRProfileStim(isNeuWithFieldAlignedDown{cond}.idxDown,:),...
                modPyr1AL.timeStepRun);
            FRProfileMeanRiseFieldStimCtrl{cond} = accumMean(avgFRProfileStimCtrl(isNeuWithFieldAlignedRise{cond}.idxRise,:),...
                modPyr1AL.timeStepRun);
            FRProfileMeanDownFieldStimCtrl{cond} = accumMean(avgFRProfileStimCtrl(isNeuWithFieldAlignedDown{cond}.idxDown,:),...
                modPyr1AL.timeStepRun);
            
            %% combining the fields from good non-stim and stim trials
            FRProfileMeanRiseFieldComb{cond} = accumMean(avgFRProfile(isNeuWithFieldAlignedRise{cond}.idxRiseComb,:),...
                modPyr1AL.timeStepRun);
            FRProfileMeanDownFieldComb{cond} = accumMean(avgFRProfile(isNeuWithFieldAlignedDown{cond}.idxDownComb,:),...
                modPyr1AL.timeStepRun);
            FRProfileMeanRiseFieldCombStim{cond} = accumMean(avgFRProfileStim(isNeuWithFieldAlignedRise{cond}.idxRiseComb,:),...
                modPyr1AL.timeStepRun);
            FRProfileMeanDownFieldCombStim{cond} = accumMean(avgFRProfileStim(isNeuWithFieldAlignedDown{cond}.idxDownComb,:),...
                modPyr1AL.timeStepRun);
            FRProfileMeanRiseFieldCombStimCtrl{cond} = accumMean(avgFRProfileStimCtrl(isNeuWithFieldAlignedRise{cond}.idxRiseComb,:),...
                modPyr1AL.timeStepRun);
            FRProfileMeanDownFieldCombStimCtrl{cond} = accumMean(avgFRProfileStimCtrl(isNeuWithFieldAlignedDown{cond}.idxDownComb,:),...
                modPyr1AL.timeStepRun);
            
        end
    end
       
    save([pathAnal0 'initPeakPyrIntAllRecStim.mat'],'PyrRise','PyrDown','FRProfile1',...
        'FRProfileMean1','FRProfileMeanStim1','FRProfileMeanStimCtrl1',...
        'FRProfileMeanField1','FRProfileMeanStimField1','FRProfileMeanStimCtrlField1',...
        'FRProfileMeanFieldComb','FRProfileMeanStimFieldComb','FRProfileMeanStimCtrlFieldComb',...    
        'FRProfileMeanRise','FRProfileMeanRiseStim','FRProfileMeanRiseStimCtrl',...
        'FRProfileMeanDown','FRProfileMeanDownStim','FRProfileMeanDownStimCtrl',...
        'isNeuWithFieldAlignedRise','isNeuWithFieldAlignedDown',...
        'FRProfileMeanRiseField','FRProfileMeanRiseFieldStim','FRProfileMeanRiseFieldStimCtrl',...
        'FRProfileMeanRiseFieldComb','FRProfileMeanRiseFieldCombStim','FRProfileMeanRiseFieldCombStimCtrl',...
        'FRProfileMeanDownField','FRProfileMeanDownFieldStim','FRProfileMeanDownFieldStimCtrl',...
        'FRProfileMeanDownFieldComb','FRProfileMeanDownFieldCombStim','FRProfileMeanDownFieldCombStimCtrl'); 
    
    %% pyramidal neurons with fields
    for i = 1:length(FRProfileMean1)
        %% normalized, fields from good non-stim trials
        %% order pyramidal neurons with field based on the peak firing rate after 0
%         plotIndFRProfile(modPyr1AL.timeStepRun,...
%                 avgFRProfileNorm(FRProfile1{i}.indField,:),['Neuron no.'],...
%                 ['Pyr_IndFRProfileNormFRPeakAftRunNeuField_Cond' num2str(i)],...
%                 pathAnal0,[],5,[],[])
%           
        %% order pyramidal neurons with field based on the peak firing rate after 0
        % and compare good non-stim with stim trials
        plotIndFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfileNorm(FRProfile1{i}.indField,:),...
                avgFRProfileNormStim(FRProfile1{i}.indField,:),['Neuron no.'],...
                ['Pyr_IndFRProfileNormFRPeakAftRunNeuFieldCmpGoodNoStimVsStim_Cond' num2str(i)],...
                pathAnal0,[],5,[],[])
            
        
        %% order pyramidal neurons with field based on before and after run FR ratio  
%         plotIndFRProfile(modPyr1AL.timeStepRun,...
%                 avgFRProfileNorm(FRProfile1{i}.indField,:),['Neuron no.'],...
%                 ['Pyr_IndFRProfileNormFR0to1VsBefRunNeuField_Cond' num2str(i)],...
%                 pathAnal0,[],4,FRProfileMeanField1{i}.indFRBefRun,...
%                 FRProfileMeanField1{i}.indFR0to1) % ordered based on -1to0 to 0to1 mean ratio
            
        %% order pyramidal neurons with field based on before and after run FR ratio 
        % and compare good non-stim with stim trials
        plotIndFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfileNorm(FRProfile1{i}.indField,:),...
                avgFRProfileNormStim(FRProfile1{i}.indField,:),['Neuron no.'],...
                ['Pyr_IndFRProfileNormFR0to1VsBefRunNeuFieldCmpGoodNoStimVsStim_Cond' num2str(i)],...
                pathAnal0,[],4,FRProfileMeanField1{i}.indFRBefRun,...
                FRProfileMeanField1{i}.indFR0to1)
    
        %% compare good non stim trials with stim trials, pyr rise
        plotAvgFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfileNorm(isNeuWithFieldAlignedRise{i}.idxRise,:),...
                avgFRProfileNormStim(isNeuWithFieldAlignedRise{i}.idxRise,:),...
                ['FR PyrRise GoodNStim/Stim Field'],...
                ['Pyr_FRProfileNormPyrRiseFieldNoStimGoodVsStim_Cond' num2str(i)],...
                pathAnal0,[0 0.8])
        
        %% compare good non stim trials with stim trials, pyr down
        plotAvgFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfileNorm(isNeuWithFieldAlignedDown{i}.idxDown,:),...
                avgFRProfileNormStim(isNeuWithFieldAlignedDown{i}.idxDown,:),...
                ['FR PyrDown GoodNStim/Stim Field'],...
                ['Pyr_FRProfileNormPyrDownFieldNoStimGoodVsStim_Cond' num2str(i)],...
                pathAnal0,[0 0.8])
            
        %% normalized, fields from both good non-stim trials and stim trials
        %% order pyramidal neurons with field based on the peak firing rate after 0
%         plotIndFRProfile(modPyr1AL.timeStepRun,...
%                 avgFRProfileNorm(FRProfile1{i}.indFieldComb,:),['Neuron no.'],...
%                 ['Pyr_IndFRProfileNormFRPeakAftRunNeuFieldComb_Cond' num2str(i)],...
%                 pathAnal0,[],5,[],[])
          
        %% order pyramidal neurons with field based on the peak firing rate after 0
        % and compare good non-stim with stim trials
        plotIndFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfileNorm(FRProfile1{i}.indFieldComb,:),...
                avgFRProfileNormStim(FRProfile1{i}.indFieldComb,:),['Neuron no.'],...
                ['Pyr_IndFRProfileNormFRPeakAftRunNeuFieldCombCmpGoodNoStimVsStim_Cond' num2str(i)],...
                pathAnal0,[],5,[],[])
            
        
        %% order pyramidal neurons with field based on before and after run FR ratio  
%         plotIndFRProfile(modPyr1AL.timeStepRun,...
%                 avgFRProfileNorm(FRProfile1{i}.indFieldComb,:),['Neuron no.'],...
%                 ['Pyr_IndFRProfileNormFR0to1VsBefRunNeuFieldComb_Cond' num2str(i)],...
%                 pathAnal0,[],4,FRProfileMeanField1{i}.indFRBefRun,...
%                 FRProfileMeanField1{i}.indFR0to1) % ordered based on -1to0 to 0to1 mean ratio
            
        %% order pyramidal neurons with field based on before and after run FR ratio 
        % and compare good non-stim with stim trials
        plotIndFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfileNorm(FRProfile1{i}.indFieldComb,:),...
                avgFRProfileNormStim(FRProfile1{i}.indFieldComb,:),['Neuron no.'],...
                ['Pyr_IndFRProfileNormFR0to1VsBefRunNeuFieldCombCmpGoodNoStimVsStim_Cond' num2str(i)],...
                pathAnal0,[],4,FRProfileMeanField1{i}.indFRBefRun,...
                FRProfileMeanField1{i}.indFR0to1)
    
        %% compare good non stim trials with stim trials, pyr rise
        plotAvgFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfileNorm(isNeuWithFieldAlignedRise{i}.idxRiseComb,:),...
                avgFRProfileNormStim(isNeuWithFieldAlignedRise{i}.idxRiseComb,:),...
                ['FR PyrRise GoodNStim/Stim Field'],...
                ['Pyr_FRProfileNormPyrRiseFieldCombNoStimGoodVsStim_Cond' num2str(i)],...
                pathAnal0,[0 0.8])
        
        %% compare good non stim trials with stim trials, pyr down
        plotAvgFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfileNorm(isNeuWithFieldAlignedDown{i}.idxDownComb,:),...
                avgFRProfileNormStim(isNeuWithFieldAlignedDown{i}.idxDownComb,:),...
                ['FR PyrDown GoodNStim/Stim Field Comb'],...
                ['Pyr_FRProfileNormPyrDownFieldCombNoStimGoodVsStim_Cond' num2str(i)],...
                pathAnal0,[0 0.8])
            
        %% not normalized, fields from good non-stim trials
        %% compare good non stim trials with stim trials, pyr rise
        plotAvgFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfile(isNeuWithFieldAlignedRise{i}.idxRise,:),...
                avgFRProfileStim(isNeuWithFieldAlignedRise{i}.idxRise,:),...
                ['FR PyrRise GoodNStim/Stim Field'],...
                ['Pyr_FRProfilePyrRiseFieldNoStimGoodVsStim_Cond' num2str(i)],...
                pathAnal0,[])
        
        %% compare good non stim trials with stim trials, pyr down
        plotAvgFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfile(isNeuWithFieldAlignedDown{i}.idxDown,:),...
                avgFRProfileStim(isNeuWithFieldAlignedDown{i}.idxDown,:),...
                ['FR PyrDown GoodNStim/Stim Field'],...
                ['Pyr_FRProfilePyrDownFieldNoStimGoodVsStim_Cond' num2str(i)],...
                pathAnal0,[])
            
        %% not  normalized, fields from both good non-stim trials and stim trials
        %% compare good non stim trials with stim trials, pyr rise
        plotAvgFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfile(isNeuWithFieldAlignedRise{i}.idxRiseComb,:),...
                avgFRProfileStim(isNeuWithFieldAlignedRise{i}.idxRiseComb,:),...
                ['FR PyrRise GoodNStim/Stim Field'],...
                ['Pyr_FRProfilePyrRiseFieldCombNoStimGoodVsStim_Cond' num2str(i)],...
                pathAnal0,[])
        
        %% compare good non stim trials with stim trials, pyr down
        plotAvgFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfile(isNeuWithFieldAlignedDown{i}.idxDownComb,:),...
                avgFRProfileStim(isNeuWithFieldAlignedDown{i}.idxDownComb,:),...
                ['FR PyrDown GoodNStim/Stim Field Comb'],...
                ['Pyr_FRProfilePyrDownFieldCombNoStimGoodVsStim_Cond' num2str(i)],...
                pathAnal0,[])
    end
    
    %% all pyramidal neurons
    for i = 1:length(FRProfileMean1)
        %% normalized
        %% order pyramidal neurons based on the peak firing rate after 0
        % and compare good non-stim with stim trials
        plotIndFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfileNorm(FRProfile1{i}.ind,:),...
                avgFRProfileNormStim(FRProfile1{i}.ind,:),['Neuron no.'],...
                ['Pyr_IndFRProfileNormFRPeakAftRunNeuCmpGoodNoStimVsStim_Cond' num2str(i)],...
                pathAnal0,[],5,[],[])
        
        %% order pyramidal neurons based on before and after run FR ratio 
        % and compare good non-stim with stim trials
        plotIndFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfileNorm(FRProfile1{i}.ind,:),...
                avgFRProfileNormStim(FRProfile1{i}.ind,:),['Neuron no.'],...
                ['Pyr_IndFRProfileNormFR0to1VsBefRunNeuCmpGoodNoStimVsStim_Cond' num2str(i)],...
                pathAnal0,[],4,FRProfileMean1{i}.indFRBefRun,...
                FRProfileMean1{i}.indFR0to1)
    
        %% compare good non stim trials with stim trials, pyr rise
        plotAvgFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfileNorm(FRProfile1{i}.indPyrRise,:),...
                avgFRProfileNormStim(FRProfile1{i}.indPyrRise,:),...
                ['FR PyrRise GoodNStim/Stim'],...
                ['Pyr_FRProfileNormPyrRiseNoStimGoodVsStim_Cond' num2str(i)],...
                pathAnal0,[0 0.6])
        
        %% compare good non stim trials with stim trials, pyr down
        plotAvgFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfileNorm(FRProfile1{i}.indPyrDown,:),...
                avgFRProfileNormStim(FRProfile1{i}.indPyrDown,:),...
                ['FR PyrDown GoodNStim/Stim'],...
                ['Pyr_FRProfileNormPyrDownNoStimGoodVsStim_Cond' num2str(i)],...
                pathAnal0,[0 0.6])
 
        %% not normalized, fields from good non-stim trials
        %% compare good non stim trials with stim trials, pyr rise
        plotAvgFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfile(FRProfile1{i}.indPyrRise,:),...
                avgFRProfileStim(FRProfile1{i}.indPyrRise,:),...
                ['FR PyrRise GoodNStim/Stim'],...
                ['Pyr_FRProfilePyrRiseNoStimGoodVsStim_Cond' num2str(i)],...
                pathAnal0,[])
        
        %% compare good non stim trials with stim trials, pyr down
        plotAvgFRProfileCmp(modPyr1AL.timeStepRun,...
                avgFRProfile(FRProfile1{i}.indPyrDown,:),...
                avgFRProfileStim(FRProfile1{i}.indPyrDown,:),...
                ['FR PyrDown GoodNStim/Stim'],...
                ['Pyr_FRProfilePyrDownNoStimGoodVsStim_Cond' num2str(i)],...
                pathAnal0,[])
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

function plotIndFRProfileCmp(timeStepRun,avgFRProfile,avgFRProfile1,yl,fileName,pathAnal,ylimit,ordMethod,indT,indT1)
    if(isempty(avgFRProfile))
        return;
    end
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 450 400]);
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
    subplot(1,2,1)
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
    
    subplot(1,2,2)
    h = imagesc(timeStepRun,1:numNeurons,avgFRProfile1(indOrd,:));
%     set(h,'LineWidth',0.1)
    set(gca,'XLim',[-1 4]);
%     if(~isempty(ylimit))
%         set(gca,'YLim',ylimit);
%     else
%         set(gca,'YLim',[min(avgFRProfile(:)) max(avgFRProfile(:))]);
%     end
    xlabel('Time (s)')
    
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