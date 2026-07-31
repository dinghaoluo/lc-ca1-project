function plotThetaPhaseVsTimeField()

    onlyRun = 1;
    RecordingList;
    pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\Pyramidal\';
    
    if(exist([pathAnal 'autoCorrPyrAllRec.mat']))
        load([pathAnal 'autoCorrPyrAllRec.mat'],...
            'autoCorrPyrNoCue','autoCorrPyrAL','autoCorrPyrPL',...
            'modPyrNoCue','modPyrAL','modPyrPL');
    end
    
    mod.indRec = autoCorrPyrAL.indRec; 
    mod.indNeu = autoCorrPyrAL.indNeu; 
    mod.isNeuWithField = modPyrAL.isNeuWithField; 
    
    disp('Active licking')
    rec1 = unique(mod.indRec);
    for i = 1:length(rec1)
        ind = mod.indRec == rec1(i) & mod.isNeuWithField == 1;
        indN = mod.indNeu(ind);
        if(length(indN) > 0)
            for j = indN
                fullPath = [listRecordingsActiveLickPath(rec1(i),:)...
                    listRecordingsActiveLickFileName(rec1(i),:)...
                    '_PeakFRAligned_msess' num2str(mazeSessionActiveLick(rec1(i)))...
                    '_Run' num2str(onlyRun) '.mat'];    
                if(exist(fullPath) == 0)
                    disp('The _PeakFRAligned file dplotoes not exist');
                    return;
                end
                load(fullPath,'trialNoNonStimGood');
                
                plotSpikeRasterThetaPhaseEg(listRecordingsActiveLickPath(rec1(i),:),...
                    listRecordingsActiveLickFileName(rec1(i),:), onlyRun, ...
                    trialNoNonStimGood, j);
            end
        end
    end
    
    mod.indRec = autoCorrPyrPL.indRec; 
    mod.indNeu = autoCorrPyrPL.indNeu; 
    mod.isNeuWithField = modPyrPL.isNeuWithField; 
    
    disp('Passive licking')
    rec1 = unique(mod.indRec);
    for i = 1:length(rec1)
        ind = mod.indRec == rec1(i) & mod.isNeuWithField == 1;
        indN = mod.indNeu(ind);
        if(length(indN) > 0)
            for j = indN
                fullPath = [listRecordingsPassiveLickPath(rec1(i),:)...
                    listRecordingsPassiveLickFileName(rec1(i),:)...
                    '_PeakFRAligned_msess' num2str(mazeSessionPassiveLick(rec1(i)))...
                    '_Run' num2str(onlyRun) '.mat'];    
                if(exist(fullPath) == 0)
                    disp('The _PeakFRAligned file does not exist');
                    return;
                end
                load(fullPath,'trialNoNonStimGood');
                
                plotSpikeRasterThetaPhaseEg(listRecordingsPassiveLickPath(rec1(i),:),...
                    listRecordingsPassiveLickFileName(rec1(i),:), onlyRun, ...
                    trialNoNonStimGood, j);
            end
        end
    end
end
