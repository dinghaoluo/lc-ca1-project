function PyrIntInitPeakSig(taskSel,pshuffle)
% compare Pyramidal neurons and PV interneurons on their initial peak
% (pshuffle = 1 -> 99.9%, 2 -> 99%, 3 -> 95%)

    pathAnal = 'Z:\Yingxue\Draft\SST\Pyramidal\';
    pathAnalInt0 = 'Z:\Yingxue\Draft\SST\Interneuron\';
    
    if(exist([pathAnal 'initPeakPyrAllRecSig.mat']))
        load([pathAnal 'initPeakPyrAllRec.mat']);
        load([pathAnal 'initPeakPyrAllRecSig.mat']);
    end
    if(exist([pathAnalInt0 'initPeakIntAllRecSig.mat']))
        load([pathAnalInt0 'initPeakIntAllRec.mat']);
        load([pathAnalInt0 'initPeakIntAllRecSig.mat']);
    end
    pathAnalPeak0 = 'Z:\Yingxue\Draft\SST\PyramidalInitPeak\2\';
    if(exist([pathAnalPeak0 'initPeakPyrAllRec_km2.mat']))
        load([pathAnalPeak0 'initPeakPyrAllRec_km2.mat'],'FRProfileMean');
    end
    
    if(taskSel == 1) % including all the neurons
        pathAnal = ['Z:\Yingxue\Draft\SST\PyramidalIntInitPeakSig\' num2str(pshuffle) '\'];
    elseif(taskSel == 2) % including AL and PL neurons
        pathAnal = ['Z:\Yingxue\Draft\SST\PyramidalIntInitPeakSigALPL\' num2str(pshuffle) '\'];
    else
        pathAnal = ['Z:\Yingxue\Draft\SST\PyramidalIntInitPeakSigAL\' num2str(pshuffle) '\'];
    end
    if(exist(pathAnal) == 0)
        mkdir(pathAnal);
    end
    if(exist([pathAnal 'initPeakPyrIntAllRecSig.mat']))
        load([pathAnal 'initPeakPyrIntAllRecSig.mat']);
    end
    
    if(exist('InitAllPyr') == 0)
        InitAllPyr = PyrIntInitPeakByTypeSig(taskSel,modPyr1NoCue,modPyr1AL,modPyr1PL,...
                                    modPyr1SigNoCue,modPyr1SigAL,modPyr1SigPL);

        InitAllInt = PyrIntInitPeakByTypeSigInt(taskSel,modInt1NoCue,modInt1AL,modInt1PL,...
                                    modInt1SigNoCue,modInt1SigAL,modInt1SigPL);

        mean0to1 = mean(InitAllPyr.avgFRProfile(:,FRProfileMean.indFR0to1),2);
        meanBefRun = mean(InitAllPyr.avgFRProfile(:,FRProfileMean.indFRBefRun),2);
        ratio0to1BefRun = mean0to1./meanBefRun;
        idxGood = ~isnan(ratio0to1BefRun) & ~isinf(ratio0to1BefRun);
        idxR = find(InitAllPyr.isPeakNeuArrAll(pshuffle,:) == 1 & idxGood' == 1);
        idxD = find(InitAllPyr.isPeakNeuArrAll(pshuffle,:) == -1 & idxGood' == 1);
        idxO = find(InitAllPyr.isPeakNeuArrAll(pshuffle,:) == 0 & idxGood' == 1);
        tsneEmbed(InitAllPyr.avgFRProfileNorm([idxR,idxO,idxD],FRProfileMean.indFRBefRun(1):FRProfileMean.indFR0to1(end)),...
            length(idxR),length(idxR)+length(idxO)+1,...
            pathAnal,'tsneEmbedPyrRiseDown');
        
        %% pyramidal rise and down
        idx = find(InitAllPyr.isPeakNeuArrAll(pshuffle,:) == 1 & idxGood' == 1); 
        % neurons with FR increase around 0
        PyrRise.idxRise = idx;
        PyrRise.task = InitAllPyr.task(idx);
        PyrRise.indRec = InitAllPyr.indRec(idx);
        PyrRise.indNeu = InitAllPyr.indNeu(idx);
        PyrRise.idxRiseBadBad = [];
        PyrRise.idxRiseBad = [];
        PyrRise.taskBad = [];
        PyrRise.indRecBad = [];
        PyrRise.indNeuBad = [];
        for i = 1:length(idx)
            idxBad = find(InitAllPyr.taskBad == PyrRise.task(i) & InitAllPyr.indRecBad == PyrRise.indRec(i) ...
                & InitAllPyr.indNeuBad == PyrRise.indNeu(i));
            if(length(idxBad) == 1)
                idxBad1 = PyrRise.idxRise(i);
                PyrRise.idxRiseBadBad = [PyrRise.idxRiseBadBad idxBad];
                PyrRise.idxRiseBad = [PyrRise.idxRiseBad idxBad1];
                PyrRise.taskBad = [PyrRise.taskBad InitAllPyr.taskBad(idxBad)];
                PyrRise.indRecBad = [PyrRise.indRecBad InitAllPyr.indRecBad(idxBad)];
                PyrRise.indNeuBad = [PyrRise.indNeuBad InitAllPyr.indNeuBad(idxBad)];
            end
        end

        % neurons with FR decrease around 0
        idx = find(InitAllPyr.isPeakNeuArrAll(pshuffle,:) == -1 & idxGood' == 1); 
        PyrDown.idxDown = idx;
        PyrDown.task = InitAllPyr.task(idx);
        PyrDown.indRec = InitAllPyr.indRec(idx);
        PyrDown.indNeu = InitAllPyr.indNeu(idx);
        PyrDown.idxDownBad = [];
        PyrDown.idxDownBadBad = [];
        PyrDown.taskBad = [];
        PyrDown.indRecBad = [];
        PyrDown.indNeuBad = [];
        for i = 1:length(idx)
            idxBad = find(InitAllPyr.taskBad == PyrDown.task(i) & InitAllPyr.indRecBad == PyrDown.indRec(i) ...
                & InitAllPyr.indNeuBad == PyrDown.indNeu(i));
            if(length(idxBad) == 1)
                idxBad1 = PyrDown.idxDown(i);
                PyrDown.idxDownBadBad = [PyrDown.idxDownBadBad idxBad];
                PyrDown.idxDownBad = [PyrDown.idxDownBad idxBad1];
                PyrDown.taskBad = [PyrDown.taskBad InitAllPyr.taskBad(idxBad)];
                PyrDown.indRecBad = [PyrDown.indRecBad InitAllPyr.indRecBad(idxBad)];
                PyrDown.indNeuBad = [PyrDown.indNeuBad InitAllPyr.indNeuBad(idxBad)];
            end
        end
        
        %% added on 7/8/2022
        % other neurons
        idx = find(InitAllPyr.isPeakNeuArrAll(pshuffle,:) == 0 & idxGood' == 1); 
        PyrOther.idxOther = idx;
        PyrOther.task = InitAllPyr.task(idx);
        PyrOther.indRec = InitAllPyr.indRec(idx);
        PyrOther.indNeu = InitAllPyr.indNeu(idx);
        PyrOther.idxOtherBad = [];
        PyrOther.idxOtherBadBad = [];
        PyrOther.taskBad = [];
        PyrOther.indRecBad = [];
        PyrOther.indNeuBad = [];
        for i = 1:length(idx)
            idxBad = find(InitAllPyr.taskBad == PyrOther.task(i) & InitAllPyr.indRecBad == PyrOther.indRec(i) ...
                & InitAllPyr.indNeuBad == PyrOther.indNeu(i));
            if(length(idxBad) == 1)
                idxBad1 = PyrOther.idxOther(i);
                PyrOther.idxOtherBadBad = [PyrOther.idxOtherBadBad idxBad];
                PyrOther.idxOtherBad = [PyrOther.idxOtherBad idxBad1];
                PyrOther.taskBad = [PyrOther.taskBad InitAllPyr.taskBad(idxBad)];
                PyrOther.indRecBad = [PyrOther.indRecBad InitAllPyr.indRecBad(idxBad)];
                PyrOther.indNeuBad = [PyrOther.indNeuBad InitAllPyr.indNeuBad(idxBad)];
            end
        end

        %% Interneuron rise and down
        for i = 1:max(InitAllInt.idxC2)
            mean0to1Int = mean(InitAllInt.avgFRProfile(InitAllInt.idxC2 == i,FRProfileMean.indFR0to1),2);
            meanBefRunInt = mean(InitAllInt.avgFRProfile(InitAllInt.idxC2 == i,FRProfileMean.indFRBefRun),2);
            ratio0to1BefRunInt = mean0to1Int./meanBefRunInt;
            
            %% neurons with FR increase around 0
            idx = find(InitAllInt.isPeakNeuArrAll(pshuffle,:) == 1 & InitAllInt.idxC2 == i); % significance p = 99%
            IntRise.idxRise{i} = idx;
            IntRise.task{i} = InitAllInt.task(idx);
            IntRise.indRec{i} = InitAllInt.indRec(idx);
            IntRise.indNeu{i} = InitAllInt.indNeu(idx);
            IntRise.idxRiseBad{i} = [];
            IntRise.idxRiseBadBad{i} = [];
            IntRise.taskBad{i} = [];
            IntRise.indRecBad{i} = [];
            IntRise.indNeuBad{i} = [];
            for j = 1:length(idx)
                idxBad = find(InitAllInt.taskBad == IntRise.task{i}(j) & InitAllInt.indRecBad == IntRise.indRec{i}(j) ...
                    & InitAllInt.indNeuBad == IntRise.indNeu{i}(j));
                if(length(idxBad) == 1)
                    idxBad1 = IntRise.idxRise{i}(j);
                    IntRise.idxRiseBadBad{i} = [IntRise.idxRiseBadBad{i} idxBad];
                    IntRise.idxRiseBad{i} = [IntRise.idxRiseBad{i} idxBad1];
                    IntRise.taskBad{i} = [IntRise.taskBad{i} InitAllInt.taskBad(idxBad)];
                    IntRise.indRecBad{i} = [IntRise.indRecBad{i} InitAllInt.indRecBad(idxBad)];
                    IntRise.indNeuBad{i} = [IntRise.indNeuBad{i} InitAllInt.indNeuBad(idxBad)];
                end
            end

            %% neurons with FR decrease around 0
            idx = find(InitAllInt.isPeakNeuArrAll(pshuffle,:) == -1 & InitAllInt.idxC2 == i); % significance p = 99%
            IntDown.idxDown{i} = idx;
            IntDown.task{i} = InitAllInt.task(idx);
            IntDown.indRec{i} = InitAllInt.indRec(idx);
            IntDown.indNeu{i} = InitAllInt.indNeu(idx);
            IntDown.idxDownBad{i} = [];
            IntDown.idxDownBadBad{i} = [];
            IntDown.taskBad{i} = [];
            IntDown.indRecBad{i} = [];
            IntDown.indNeuBad{i} = [];
            for j = 1:length(idx)
                idxBad = find(InitAllInt.taskBad == IntDown.task{i}(j) & InitAllInt.indRecBad == IntDown.indRec{i}(j) ...
                    & InitAllInt.indNeuBad == IntDown.indNeu{i}(j));
                if(length(idxBad) == 1)
                    idxBad1 = IntDown.idxDown{i}(j);
                    IntDown.idxDownBadBad{i} = [IntDown.idxDownBadBad{i} idxBad];
                    IntDown.idxDownBad{i} = [IntDown.idxDownBad{i} idxBad1];
                    IntDown.taskBad{i} = [IntDown.taskBad{i} InitAllInt.taskBad(idxBad)];
                    IntDown.indRecBad{i} = [IntDown.indRecBad{i} InitAllInt.indRecBad(idxBad)];
                    IntDown.indNeuBad{i} = [IntDown.indNeuBad{i} InitAllInt.indNeuBad(idxBad)];
                end
            end
        end

        %% pyramidal neurons
        relDepthNeuHDefPyr.depthRise = InitAllPyr.relDepthNeuHDef(PyrRise.idxRise);
        relDepthNeuHDefPyr.depthDown = InitAllPyr.relDepthNeuHDef(PyrDown.idxDown);
        relDepthNeuHDefPyr.depthOther = InitAllPyr.relDepthNeuHDef(PyrOther.idxOther);% added on 7/8/2022
        relDepthNeuHDefPyr.depthRiseMean = mean(relDepthNeuHDefPyr.depthRise);
        relDepthNeuHDefPyr.depthDownMean = mean(relDepthNeuHDefPyr.depthDown);
        relDepthNeuHDefPyr.depthOtherMean = mean(relDepthNeuHDefPyr.depthOther);% added on 7/8/2022
        relDepthNeuHDefPyr.depthRiseSem = std(relDepthNeuHDefPyr.depthRise)/...
                sqrt(length(relDepthNeuHDefPyr.depthRise));
        relDepthNeuHDefPyr.depthDownSem = std(relDepthNeuHDefPyr.depthDown)/...
                sqrt(length(relDepthNeuHDefPyr.depthDown));  
        relDepthNeuHDefPyr.depthOtherSem = std(relDepthNeuHDefPyr.depthOther)/...
                sqrt(length(relDepthNeuHDefPyr.depthOther));  % added on 7/8/2022
        relDepthNeuHDefPyr.pRSRelDepthNeuHDef = ranksum(relDepthNeuHDefPyr.depthRise,relDepthNeuHDefPyr.depthDown);
        relDepthNeuHDefPyr.pRSRelDepthNeuHDefRO = ranksum(relDepthNeuHDefPyr.depthRise,relDepthNeuHDefPyr.depthOther);% added on 7/8/2022
        relDepthNeuHDefPyr.pRSRelDepthNeuHDefDO = ranksum(relDepthNeuHDefPyr.depthDown,relDepthNeuHDefPyr.depthOther);% added on 7/8/2022

        %% for field
        isNeuWithFieldAligned.isFieldRise = InitAllPyr.isNeuWithFieldAligned(PyrRise.idxRise);
        isNeuWithFieldAligned.isFieldDown = InitAllPyr.isNeuWithFieldAligned(PyrDown.idxDown);
        isNeuWithFieldAligned.isFieldOther = InitAllPyr.isNeuWithFieldAligned(PyrOther.idxOther); % added on 7/8/2022
        isNeuWithFieldAligned.isFieldRiseBad = InitAllPyr.isNeuWithFieldAligned(PyrRise.idxRiseBad);
        isNeuWithFieldAligned.isFieldDownBad = InitAllPyr.isNeuWithFieldAligned(PyrDown.idxDownBad);
        isNeuWithFieldAligned.isFieldOtherBad = InitAllPyr.isNeuWithFieldAligned(PyrOther.idxOtherBad); % added on 7/8/2022
        isNeuWithFieldAligned.numFieldRise = sum(isNeuWithFieldAligned.isFieldRise);
        isNeuWithFieldAligned.numFieldDown = sum(isNeuWithFieldAligned.isFieldDown);
        isNeuWithFieldAligned.numFieldOther = sum(isNeuWithFieldAligned.isFieldOther); % added on 7/8/2022
        isNeuWithFieldAligned.idxRise = PyrRise.idxRise(isNeuWithFieldAligned.isFieldRise == 1);
        isNeuWithFieldAligned.idxDown = PyrDown.idxDown(isNeuWithFieldAligned.isFieldDown == 1);
        isNeuWithFieldAligned.idxOther = PyrOther.idxOther(isNeuWithFieldAligned.isFieldOther == 1); % added on 7/8/2022
        isNeuWithFieldAligned.idxRiseBad = PyrRise.idxRiseBad(isNeuWithFieldAligned.isFieldRiseBad == 1);
        isNeuWithFieldAligned.idxDownBad = PyrDown.idxDownBad(isNeuWithFieldAligned.isFieldDownBad == 1);
        isNeuWithFieldAligned.idxOtherBad = PyrOther.idxOtherBad(isNeuWithFieldAligned.isFieldOtherBad == 1);  % added on 7/8/2022
        isNeuWithFieldAligned.idxRiseBadBad = PyrRise.idxRiseBadBad(isNeuWithFieldAligned.isFieldRiseBad == 1);
        isNeuWithFieldAligned.idxDownBadBad = PyrDown.idxDownBadBad(isNeuWithFieldAligned.isFieldDownBad == 1);
        isNeuWithFieldAligned.idxOtherBadBad = PyrOther.idxOtherBadBad(isNeuWithFieldAligned.isFieldOtherBad == 1); % added on 7/8/2022

        FRProfileMeanPyr.RiseField = accumMeanPVStim(InitAllPyr.avgFRProfile(isNeuWithFieldAligned.idxRise,:),...
            modPyr1NoCue.timeStepRun);
        FRProfileMeanPyr.DownField = accumMeanPVStim(InitAllPyr.avgFRProfile(isNeuWithFieldAligned.idxDown,:),...
            modPyr1NoCue.timeStepRun);
        FRProfileMeanPyr.OtherField = accumMeanPVStim(InitAllPyr.avgFRProfile(isNeuWithFieldAligned.idxOther,:),...
            modPyr1NoCue.timeStepRun); % added on 7/8/2022
        FRProfileMeanPyr.RiseFieldBad = accumMeanPVStim(InitAllPyr.avgFRProfileBad(isNeuWithFieldAligned.idxRiseBadBad,:),...
            modPyr1NoCue.timeStepRun);
        FRProfileMeanPyr.DownFieldBad = accumMeanPVStim(InitAllPyr.avgFRProfileBad(isNeuWithFieldAligned.idxDownBadBad,:),...
            modPyr1NoCue.timeStepRun);
        FRProfileMeanPyr.OtherFieldBad = accumMeanPVStim(InitAllPyr.avgFRProfileBad(isNeuWithFieldAligned.idxOtherBadBad,:),...
            modPyr1NoCue.timeStepRun); % added on 7/8/2022

        FRProfileMeanPyrStat.RiseField = accumMeanPyrIntInitPeakStatC(FRProfileMeanPyr.RiseField);
        FRProfileMeanPyrStat.DownField = accumMeanPyrIntInitPeakStatC(FRProfileMeanPyr.DownField);
        FRProfileMeanPyrStat.OtherField = accumMeanPyrIntInitPeakStatC(FRProfileMeanPyr.OtherField); % added on 7/8/2022

        FRProfileMeanPyrStat.RiseFieldBad = accumMeanPyrIntInitPeakStatC(FRProfileMeanPyr.RiseFieldBad);
        FRProfileMeanPyrStat.DownFieldBad = accumMeanPyrIntInitPeakStatC(FRProfileMeanPyr.DownFieldBad);
        FRProfileMeanPyrStat.OtherFieldBad = accumMeanPyrIntInitPeakStatC(FRProfileMeanPyr.OtherFieldBad); % added on 7/8/2022

        % compare good and bad trials
        FRProfileMeanPyrStat.RiseFieldGoodBad = accumMeanPyrIntInitPeakStatGoodBad(FRProfileMeanPyr.RiseField,FRProfileMeanPyr.RiseFieldBad);
        FRProfileMeanPyrStat.DownFieldGoodBad = accumMeanPyrIntInitPeakStatGoodBad(FRProfileMeanPyr.DownField,FRProfileMeanPyr.DownFieldBad);
        FRProfileMeanPyrStat.OtherFieldGoodBad = accumMeanPyrIntInitPeakStatGoodBad(FRProfileMeanPyr.OtherField,FRProfileMeanPyr.OtherFieldBad); % added on 7/8/2022
        
        %% for all the pyramidal neurons
        FRProfileMeanPyr.Rise = accumMeanPVStim(InitAllPyr.avgFRProfile(PyrRise.idxRise,:),modPyr1NoCue.timeStepRun);
        FRProfileMeanPyr.Down = accumMeanPVStim(InitAllPyr.avgFRProfile(PyrDown.idxDown,:),modPyr1NoCue.timeStepRun);
        FRProfileMeanPyr.Other = accumMeanPVStim(InitAllPyr.avgFRProfile(PyrOther.idxOther,:),modPyr1NoCue.timeStepRun); % added on 7/8/2022

        FRProfileMeanPyr.RiseBad = accumMeanPVStim(InitAllPyr.avgFRProfileBad(PyrRise.idxRiseBadBad,:),modPyr1NoCue.timeStepRun);
        FRProfileMeanPyr.DownBad = accumMeanPVStim(InitAllPyr.avgFRProfileBad(PyrDown.idxDownBadBad,:),modPyr1NoCue.timeStepRun);
        FRProfileMeanPyr.OtherBad = accumMeanPVStim(InitAllPyr.avgFRProfileBad(PyrOther.idxOtherBadBad,:),modPyr1NoCue.timeStepRun); % added on 7/8/2022

        FRProfileMeanPyrStat.Rise = accumMeanPyrIntInitPeakStatC(FRProfileMeanPyr.Rise);
        FRProfileMeanPyrStat.Down = accumMeanPyrIntInitPeakStatC(FRProfileMeanPyr.Down);
        FRProfileMeanPyrStat.Other = accumMeanPyrIntInitPeakStatC(FRProfileMeanPyr.Other); % added on 7/8/2022

        FRProfileMeanPyrStat.RiseBad = accumMeanPyrIntInitPeakStatC(FRProfileMeanPyr.RiseBad);
        FRProfileMeanPyrStat.DownBad = accumMeanPyrIntInitPeakStatC(FRProfileMeanPyr.DownBad);
        FRProfileMeanPyrStat.OtherBad = accumMeanPyrIntInitPeakStatC(FRProfileMeanPyr.OtherBad); % added on 7/8/2022

        % compare good and bad trials
        FRProfileMeanPyrStat.RiseGoodBad = accumMeanPyrIntInitPeakStatGoodBad(FRProfileMeanPyr.Rise,FRProfileMeanPyr.RiseBad);
        FRProfileMeanPyrStat.DownGoodBad = accumMeanPyrIntInitPeakStatGoodBad(FRProfileMeanPyr.Down,FRProfileMeanPyr.DownBad);
        FRProfileMeanPyrStat.OtherGoodBad = accumMeanPyrIntInitPeakStatGoodBad(FRProfileMeanPyr.Other,FRProfileMeanPyr.OtherBad); % added on 7/8/2022

        %% for all the interneurons
        for i = 1:max(InitAllInt.idxC2)
            FRProfileMeanInt.Rise{i} = accumMeanPVStim(InitAllInt.avgFRProfile(IntRise.idxRise{i},:),modInt1NoCue.timeStepRun);
            FRProfileMeanInt.Down{i} = accumMeanPVStim(InitAllInt.avgFRProfile(IntDown.idxDown{i},:),modInt1NoCue.timeStepRun);

            FRProfileMeanInt.RiseBad{i} = accumMeanPVStim(InitAllInt.avgFRProfileBad(IntRise.idxRiseBadBad{i},:),modInt1NoCue.timeStepRun);
            FRProfileMeanInt.DownBad{i} = accumMeanPVStim(InitAllInt.avgFRProfileBad(IntDown.idxDownBadBad{i},:),modInt1NoCue.timeStepRun);

            FRProfileMeanIntStat.Rise{i} = accumMeanPyrIntInitPeakStatC(FRProfileMeanInt.Rise{i});
            FRProfileMeanIntStat.Down{i} = accumMeanPyrIntInitPeakStatC(FRProfileMeanInt.Down{i});

            FRProfileMeanIntStat.RiseBad{i} = accumMeanPyrIntInitPeakStatC(FRProfileMeanInt.RiseBad{i});
            FRProfileMeanIntStat.DownBad{i} = accumMeanPyrIntInitPeakStatC(FRProfileMeanInt.DownBad{i});

            % compare good and bad trials
            FRProfileMeanIntStat.RiseGoodBad{i} = accumMeanPyrIntInitPeakStatGoodBad(FRProfileMeanInt.Rise{i},FRProfileMeanInt.RiseBad{i});
    %         FRProfileMeanIntStat.DownGoodBad{i} = accumMeanPyrIntInitPeakStatGoodBad(FRProfileMeanInt.Down{i},FRProfileMeanInt.DownBad{i});
        end

        save([pathAnal 'initPeakPyrIntAllRecSig.mat'],'InitAllPyr','InitAllInt','PyrRise','PyrDown','PyrOther',...
            'IntRise','IntDown',...
            'FRProfileMeanPyr','FRProfileMeanPyrStat',...
            'FRProfileMeanInt','FRProfileMeanIntStat',...
            'relDepthNeuHDefPyr','isNeuWithFieldAligned','-v7.3'); 
    end
    
%     %% plot depth of rise and down pyramidal neurons
    colorSel = 0;
%     plotBoxPlot(relDepthNeuHDefPyr.depthRise,...
%         relDepthNeuHDefPyr.depthDown,'Depth',...
%         'Pyr_relDepthNeuHDefRiseVsDown',pathAnal,[],relDepthNeuHDefPyr.pRSRelDepthNeuHDef,colorSel,[{'Rise'} {'Down'}]);
%         
%     plotBoxPlot(relDepthNeuHDefPyr.depthRise,...
%         relDepthNeuHDefPyr.depthOther,'Depth',...
%         'Pyr_relDepthNeuHDefRiseVsOther',pathAnal,[],relDepthNeuHDefPyr.pRSRelDepthNeuHDefRO,colorSel,[{'Rise'} {'Other'}]);
%     
%     plotBoxPlot(relDepthNeuHDefPyr.depthDown,...
%         relDepthNeuHDefPyr.depthOther,'Depth',...
%         'Pyr_relDepthNeuHDefDownVsOther',pathAnal,[],relDepthNeuHDefPyr.pRSRelDepthNeuHDefDO,colorSel,[{'Down'} {'Other'}]);
%     
%     plotPyrNeuRiseDownFieldSig(pathAnal,modInt1NoCue.timeStepRun,InitAllPyr,...
%             isNeuWithFieldAligned,FRProfileMean,FRProfileMeanPyr,FRProfileMeanPyrStat,colorSel,[{'Good'} {'Bad'}]);
%     close all;
    
    plotPyrNeuRiseDownSig(pathAnal,modInt1NoCue.timeStepRun,InitAllPyr,PyrRise,PyrDown,PyrOther,...
            FRProfileMean,FRProfileMeanPyr,FRProfileMeanPyrStat,colorSel,[{'Good'} {'Bad'}]);
    close all;
    
    plotIntNeuRiseDownSig(pathAnal,modInt1NoCue.timeStepRun,...
            InitAllInt,IntRise,IntDown,FRProfileMean,FRProfileMeanInt,FRProfileMeanIntStat,colorSel,[{'Good'} {'Bad'}]);
    close all;
    
%     ComparePyrRiseDownIntSig(pathAnal,modInt1NoCue.timeStepRun,InitAllPyr,InitAllInt,PyrRise,PyrDown,IntRise);    
%     close all;
end
