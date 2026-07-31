function [trials, cluList, lapList] = GetTrials_MPFI(filename, varargin)
%%% added delay area wheel running task

% [trials, cluList, lapList] = GetTrials_JF_v4(mazeSection, mazeTrialType, OverwriteTrials, OverwriteAll);
%
% function outputs a cell data structure organized by laps
% input: filename;
% optional input: maze section (see bellow), trial type for maze section 13 (see bellow):
%
% maze section = ID 
%-----------------------------------------
% treadmill -- 1
%
% trial type for maze section between 2-6
% 0 - the whole treadmill run
%
% examples:
%    GetTrials_JF_v4('A132-20111025-01',1,1,0,0);

    [mazeSection, mazeTrialType, OverwriteTrials, OverwriteAll] = ...
        DefaultArgs(varargin,{1, 1, 0, 0});

    if exist([filename '_BehavElectrDataLFP.mat'],'file') == 2 ...
            && ~OverwriteAll
        fprintf('\n Spike, Clu, Track, Laps, xml structures loaded from a file.\n');
        load([filename '_BehavElectrDataLFP.mat'],...
            'Track','Laps','Spike','Clu');
    else
        GenerateBehavElectroDataStructures_MPFIv1(pwd, filename);
        load([filename '_BehavElectrDataLFP.mat'],...
                'Track','Laps','Clu','Spike'); 
    end
        
    load([filename '-param.mat']);    
    lfpFreq = device.RecSyst.lfpSamplRate;
    sampleFreq = device.RecSyst.samplingRate;
    voltNorm = (1/(2^device.RecSyst.ADresolution)) ...
        * device.RecSyst.ADvoltRange ...
        * (1/device.RecSyst.amplification) * 1000;
    
    totNLaps = length(Laps.lapID);
    nTotClu = max(Clu.totClu);
    cluList1 = Clu.totClu;    
   
    eegArtifThresh = 5;%1.5;   % mV
    eegArtifHalfThresh = 7;%2.5;   % mV
    %     eegArtifThresh = 0.9;   % mV
    
    %% check whether the trials are already extracted and stored 
    stimOnlyN = 0;
    if(exist([filename '_DataStructure_mazeSection' num2str(mazeSection) ...
            '_TrialType' num2str(mazeTrialType) '.mat'],'file') == 2)
%         stimOnly=input(['The DataStructure file already exists, do you want to skip' ...
%                   '\nthe trial extraction, Y/N: '],'s');
        stimOnly='Y';
        if(isempty(stimOnly))
            stimOnly = 'N';
        end
        if(~isempty(strfind(stimOnly,'Y')) || ~isempty(strfind(stimOnly,'y')))
           stimOnlyN = 1; 
        end
    end
    
    if(stimOnlyN == 0)
        %% look for the start and stop indices of each section
        lfpIndStartAll = cell(1,totNLaps);
        lfpIndEndAll = cell(1,totNLaps);
        if mazeSection == 1 
            %%%%%%% RUNNING WHEEL %%%%%%%%%%%
            % select start and end point of trials based on whlTrialType variable
            for i = 1:totNLaps                
                ind = find(Laps.lapID == i,1);
                if(isempty(ind) || Laps.behavType(i) ~= 1)
                    continue;
                else
                    lfpIndStartAll{i} = Laps.startLfpInd(i);
                    lfpIndEndAll{i} = Laps.endLfpInd(i);

                    figure(10);
                    plot((Laps.startLfpInd(i):Laps.endLfpInd(i))./lfpFreq,...
                        Track.xMM(Laps.startLfpInd(i):Laps.endLfpInd(i)),'r.')

                    title(['Trial #' num2str(i)]);
                    xlabel('Time (s)')
                    ylabel('Distance (mm)')
%                     pause;
                end
            end
        end
        close(findobj('type','figure','CurrentObject',10));     
            % close fig 10 if it exists
    
    
        %% collect the data according to the start and stop indices
        disp('Assigning values to trials....');
        goodLaps = Laps.behavType;
        totT = 0;
        totCluSp = zeros(nTotClu,1);
        trials = cell(1,totNLaps);
        for nlap=1:totNLaps
            if(i == 194)
                a = 1;
            end
            disp(['Trial no. ' num2str(nlap)])
            trials{nlap} = [];
            if(~isempty(lfpIndStartAll{nlap}))
                totT = totT + lfpIndEndAll{nlap} - lfpIndStartAll{nlap} + 1;
                trials{nlap}.lfpIndStart = lfpIndStartAll{nlap};
                trials{nlap}.lfpIndEnd = lfpIndEndAll{nlap};
                trials{nlap}.Nsamples = lfpIndEndAll{nlap} ...
                    - lfpIndStartAll{nlap} + 1;

                trials{nlap}.corrTrial = Laps.corrChoice(nlap);
                % determine behavior for each lap (good-1)
                trials{nlap}.behavChar = Laps.behavType(nlap);
                trials{nlap}.trackLen = Laps.trackLen(nlap);
                trials{nlap}.lickLfpInd = Laps.lickLfpInd{nlap}...
                    - lfpIndStartAll{nlap}+1;
                trials{nlap}.pumpLfpInd = Laps.pumpLfpInd{nlap}...
                    - lfpIndStartAll{nlap}+1;  
                if(isfield(Laps,'airpuffLfpInd'))
                    trials{nlap}.airpuffLfpInd = Laps.airpuffLfpInd{nlap}...
                        - lfpIndStartAll{nlap}+1;
                else
                    trials{nlap}.airpuffLfpInd = [];
                end
                
                trials{nlap}.movieOnLfpInd = Laps.movieOnLfpInd{nlap}...
                        - lfpIndStartAll{nlap}+1;
                
                %%%%%% added by Yingxue on 03/30/2019
                if(isfield(Laps,'trStartCueLR'))
                    trials{nlap}.movieID = Laps.movieID{nlap};
                    trials{nlap}.movieOnLPulseInd  = Laps.movieOnLPulseLfpInd{nlap}...
                        - lfpIndStartAll{nlap}+1;
                    trials{nlap}.movieOffLPulseInd  = Laps.movieOffLPulseLfpInd{nlap}...
                        - lfpIndStartAll{nlap}+1;
                    trials{nlap}.movieOnRPulseInd  = Laps.movieOnRPulseLfpInd{nlap}...
                        - lfpIndStartAll{nlap}+1;
                    trials{nlap}.movieOffRPulseInd  = Laps.movieOffRPulseLfpInd{nlap}...
                        - lfpIndStartAll{nlap}+1;
                else
                    trials{nlap}.movieOnPulseInd  = Laps.movieOnPulseLfpInd{nlap}...
                        - lfpIndStartAll{nlap}+1;
                    trials{nlap}.movieOffPulseInd  = Laps.movieOffPulseLfpInd{nlap}...
                        - lfpIndStartAll{nlap}+1;
                end
                %%%%%%%%%
                
                %%%%%%% added by Yingxue on 03/11/2020
                if(isfield(Laps,'movieLocation'))
                    trials{nlap}.movieLocation = Laps.movieLocation(nlap);
                end
                %%%%%%
                
                if(isfield(Laps,'lickPeriodLfpInd'))
                    trials{nlap}.lickPeriodLfpInd = Laps.lickPeriodLfpInd{nlap}...
                        - lfpIndStartAll{nlap}+1;
                end
                
                if(isfield(Laps,'stimOnLfpInd'))
                    if(~isempty(Laps.stimOnLfpInd{nlap}))
                        trials{nlap}.stimOnLfpInd = Laps.stimOnLfpInd{nlap}...
                            - lfpIndStartAll{nlap}+1;
                        trials{nlap}.stimPulseLfpInd = Laps.stimPulseLfpInd{nlap}...
                            - lfpIndStartAll{nlap}+1;
                        trials{nlap}.stimPulseWidth = Laps.stimPulseWidth{nlap}/1000*lfpFreq;
                        trials{nlap}.stimDiode = Laps.stimDiode{nlap};
                        trials{nlap}.stimDiodeCurr = Laps.stimDiodeCurr{nlap};
                        trials{nlap}.stimPulseMethod = Laps.stimPulseMethod{nlap};
                        trials{nlap}.stimPulseLoc = Laps.stimPulseLoc{nlap};
                    else
                        trials{nlap}.stimOnLfpInd = [];
                        trials{nlap}.stimPulseLfpInd = [];
                        trials{nlap}.stimPulseWidth = [];
                        trials{nlap}.stimDiode = [];
                        trials{nlap}.stimDiodeCurr = [];
                        trials{nlap}.stimPulseMethod = Laps.stimPulseMethod{nlap};
                        trials{nlap}.stimPulseLoc = [];
                    end
                end
                
                % theta phase
                trials{nlap}.thetaHil = ...
                    Track.thetaPhHilb(lfpIndStartAll{nlap}:lfpIndEndAll{nlap});
                trials{nlap}.cumsumThetaHil = unwrap(trials{nlap}.thetaHil);
                trials{nlap}.thetaLin = ...
                    Track.thetaPhLinInterp(lfpIndStartAll{nlap}:lfpIndEndAll{nlap});
                trials{nlap}.cumsumThetaLin = unwrap(trials{nlap}.thetaLin);
                trials{nlap}.thetaAmpHilb = ...
                    Track.thetaAmpHilb(lfpIndStartAll{nlap}:lfpIndEndAll{nlap});
                trials{nlap}.thetaFreqHilb = ...
                    Track.thetaFreqHilb(lfpIndStartAll{nlap}:lfpIndEndAll{nlap});
                trials{nlap}.eeg = ...
                    Track.eeg(lfpIndStartAll{nlap}:lfpIndEndAll{nlap}) .* voltNorm;
                trials{nlap}.isEegTraceArtifact = ...
                    (sum(trials{nlap}.eeg > eegArtifThresh) > 1 ...
                    && sum(trials{nlap}.eeg < -eegArtifThresh) > 1) ...
                    || (sum(trials{nlap}.eeg > eegArtifHalfThresh) > 1) ...
                    || (sum(trials{nlap}.eeg < -eegArtifHalfThresh) > 1);

                %  detect individual theta cycles and determine the 
                %  amplitude of each (in mV!!!!!!!!!!!!!!!!!!!)
                % theta peak and trough amplitude
                ind = Track.thetaPeak_tAmpl(:,1) >= lfpIndStartAll{nlap} ...
                    & Track.thetaPeak_tAmpl(:,1) <= lfpIndEndAll{nlap};
                trials{nlap}.thetaPeak_tAmpl(:,1) = ...
                    Track.thetaPeak_tAmpl(ind,1) - lfpIndStartAll{nlap} + 1;
                trials{nlap}.thetaPeak_tAmpl(:,2) = ...
                    Track.thetaPeak_tAmpl(ind,2) .* voltNorm;
                
                ind = Track.thetaTrough_tAmpl(:,1) >= lfpIndStartAll{nlap} ...
                    & Track.thetaTrough_tAmpl(:,1) <= lfpIndEndAll{nlap};
                trials{nlap}.thetaTrough_tAmpl(:,1) = ...
                    Track.thetaTrough_tAmpl(ind,1) - lfpIndStartAll{nlap} + 1;
                trials{nlap}.thetaTrough_tAmpl(:,2) = ...
                    Track.thetaTrough_tAmpl(ind,2) .* voltNorm;
                
                ind = Track.thetaPtoTZeros_tAmpl(:,1) >= lfpIndStartAll{nlap} ...
                    & Track.thetaPtoTZeros_tAmpl(:,1) <= lfpIndEndAll{nlap};
                trials{nlap}.thetaPtoTZeros_tAmpl(:,1) = ...
                    Track.thetaPtoTZeros_tAmpl(ind,1) - lfpIndStartAll{nlap} + 1;
                trials{nlap}.thetaPtoTZeros_tAmpl(:,2) = ...
                    Track.thetaPtoTZeros_tAmpl(ind,2) .* voltNorm;
                
                ind = Track.thetaTtoPZeros_tAmpl(:,1) >= lfpIndStartAll{nlap} ...
                    & Track.thetaTtoPZeros_tAmpl(:,1) <= lfpIndEndAll{nlap};
                trials{nlap}.thetaTtoPZeros_tAmpl(:,1) = ...
                    Track.thetaTtoPZeros_tAmpl(ind,1) - lfpIndStartAll{nlap} + 1;
                trials{nlap}.thetaTtoPZeros_tAmpl(:,2) = ...
                    Track.thetaTtoPZeros_tAmpl(ind,2) .* voltNorm;

                % location and distance
                trials{nlap}.xMM = Track.xMM(lfpIndStartAll{nlap}:lfpIndEndAll{nlap});
                trials{nlap}.xMMAll = Track.xMMAll(lfpIndStartAll{nlap}:lfpIndEndAll{nlap});

                % speed and distance
                trials{nlap}.speed = ...
                    Track.speed_MMsec(lfpIndStartAll{nlap}:lfpIndEndAll{nlap});
                if(trials{nlap}.speed(1) < 0)
                    trials{nlap}.speed(1) = trials{nlap}.speed(2);
                end
                trials{nlap}.speedAll = ...
                    Track.speed_MMsecAll(lfpIndStartAll{nlap}:lfpIndEndAll{nlap});
                if(trials{nlap}.speedAll(1) < 0)
                    trials{nlap}.speedAll(1) = trials{nlap}.speedAll(2);
                end

                trials{nlap}.IDShank = zeros(nTotClu,1);
                trials{nlap}.locClu = zeros(nTotClu,1);
                trials{nlap}.spikes = cell(1,nTotClu);
                trials{nlap}.spikes20kHz = cell(1,nTotClu);
                trials{nlap}.spikesCumSumThetaHil = cell(1,nTotClu);
                trials{nlap}.spikesCumSumThetaLin = cell(1,nTotClu);
                trials{nlap}.spikesThetaHil = cell(1,nTotClu);
                trials{nlap}.spikesThetaLin = cell(1,nTotClu);
                trials{nlap}.spikesMM = cell(1,nTotClu);
                trials{nlap}.spikesSpeed = cell(1,nTotClu);
                trials{nlap}.spikesMMAll = cell(1,nTotClu);
                trials{nlap}.spikesSpeedAll = cell(1,nTotClu);

                % times and spikes within each trial
                allSpTrain = [];
                allSpClu = [];
                for nclu = 1 : nTotClu
                    trials{nlap}.IDShank(nclu) = Clu.shank(nclu);
                    trials{nlap}.locClu(nclu) = Clu.localClu(nclu);

                    if ~isempty(find(cluList1 == nclu, 1))
                        ind = Spike.res >= lfpIndStartAll{nlap} ...
                            & Spike.res <= lfpIndEndAll{nlap} ...
                            & Spike.totclu == nclu;
                        totCluSp(nclu) = totCluSp(nclu) + sum(ind);
                        allSp = Spike.res(ind) - lfpIndStartAll{nlap} + 1;
                        trials{nlap}.spikes{nclu} = allSp;
                        % res at 20kHz
                        allSp20kHz = Spike.res20kHz(ind) ...
                            - round(lfpIndStartAll{nlap}/lfpFreq*sampleFreq) + 1;
                        trials{nlap}.spikes20kHz{nclu} = allSp20kHz;

                        % spikes in theta phase 
                        trials{nlap}.spikesCumSumThetaHil{nclu} = ...
                            trials{nlap}.cumsumThetaHil(allSp);
                        trials{nlap}.spikesCumSumThetaLin{nclu} = ...
                            trials{nlap}.cumsumThetaLin(allSp);
                        trials{nlap}.spikesThetaHil{nclu} = ...
                            trials{nlap}.thetaHil(allSp);
                        trials{nlap}.spikesThetaLin{nclu} = ...
                            trials{nlap}.thetaLin(allSp);

                        % spike distance 
                        trials{nlap}.spikesMM{nclu} = ...
                            trials{nlap}.xMM(allSp);
                        trials{nlap}.spikesSpeed{nclu} = ...
                            trials{nlap}.speed(allSp);
                        trials{nlap}.spikesMMAll{nclu} = ...
                            trials{nlap}.xMMAll(allSp);
                        trials{nlap}.spikesSpeedAll{nclu} = ...
                            trials{nlap}.speedAll(allSp);

                        allSpTrain = [allSpTrain; allSp];
                        allSpClu = [allSpClu; nclu*ones(length(allSp),1)];

                    end
                end
                trials{nlap}.allSpTrain = [allSpTrain allSpClu];

%                 if(isempty(Stim))
%                     continue;
%                 end
%                 % whether there is any stimuation within a trial
%                 indStim = find(Stim.startLfpInd >= lfpIndStartAll{nlap} ...
%                     & Stim.startLfpInd < lfpIndEndAll{nlap});
%                 indStim1 = find(Stim.stopLfpInd >= lfpIndStartAll{nlap} ...
%                     & Stim.stopLfpInd < lfpIndEndAll{nlap});
%                 indStim = union(indStim, indStim1,'sorted');
%                 if(isempty(indStim))
%                     trials{nlap}.isStimulated = 0;
%                     trials{nlap}.stimStartInd = [];
%                     trials{nlap}.stimStopInd = [];
%                     trials{nlap}.indStim = [];
%                     trials{nlap}.indDiode = [];
%                     trials{nlap}.indPulseInStim = [];
%                     trials{nlap}.pulsePeriod = [];
%                 else
%                     trials{nlap}.isStimulated = 1;
%                     if(Stim.startLfpInd(indStim(1)) < lfpIndStartAll{nlap})
%                         stimStartInd = [1; ...
%                             Stim.startLfpInd(indStim(2:end)) ...
%                                 - lfpIndStartAll{nlap} + 1];
%                     else
%                         stimStartInd = Stim.startLfpInd(indStim) ...
%                             - lfpIndStartAll{nlap} + 1;
%                     end
%                     trials{nlap}.stimStartInd = stimStartInd;
%                     if(Stim.stopLfpInd(indStim(end)) >= lfpIndEndAll{nlap})
%                         stimStopInd = [Stim.stopLfpInd(indStim(1:end-1)) ...
%                             - lfpIndStartAll{nlap} + 1; trials{nlap}.Nsamples];
%                     else
%                         stimStopInd = Stim.stopLfpInd(indStim) ...
%                             - lfpIndStartAll{nlap} + 1;
%                     end
%                     trials{nlap}.stimStopInd = stimStopInd;
%                     trials{nlap}.indStim = Stim.indStim(indStim);
%                     trials{nlap}.indDiode = Stim.indDiode(indStim);
%                     trials{nlap}.indPulseInStim = Stim.indPulseInStim(indStim);
%                     trials{nlap}.pulsePeriod = Stim.pulsePeriod(indStim);
%                 end  
            end
        end
            
        cluList.all = 1:size(totCluSp,1);
        cluList.firingRate = totCluSp'/(totT/lfpFreq);
        cluList.shank = Clu.shank;
        cluList.localClu = Clu.localClu;
        cluList.spkWidthC = Clu.SpkWidthC;
        cluList.isIntern = Clu.isIntern;
        cluList.refracViolPercent = Clu.RefracViolPercent;
        cluList.mahalDist = Clu.isolDist;
        cluList.spatLocalChan = Clu.SpatLocalChan;
        cluList.spatLocalProbeCh = Clu.SpatLocalProbeCh;
        cluList.spatLocalRelAmpl = Clu.SpatLocalRelAmpl;
        cluList.rightMax = Clu.RightMax;
        cluList.leftMax = Clu.LeftMax;
        cluList.centerMax = Clu.CenterMax;
%         cluList.ccgVal = CCG.ccgVal;
%         cluList.ccgT = CCG.ccgT;

        lapList = goodLaps;
        lap.trackLen = Laps.trackLen;
        lap.mazeSess = Laps.mazeSess;
        lap.mazeType = Laps.mazeType;
        %%%%%% added by Yingxue on 03/30/2019
        if(isfield(Laps,'trStartCueLR'))
            lap.trStartCueLR = Laps.trStartCueLR;
        end
        %%%%%%
    
        if (mazeSection == 1)
            fprintf('\nOutput structures saved: %s\n\n', ...
                [filename '_DataStructure_mazeSection' num2str(mazeSection) ...
                '_TrialType' num2str(mazeTrialType) '.mat']);
            save([filename '_DataStructure_mazeSection' num2str(mazeSection) ...
                '_TrialType' num2str(mazeTrialType) '.mat'], ...
                'trials', 'cluList', 'lapList', 'lap', '-v7.3');
        end
    end
    
    if(isfield(Laps,'stimOnLfpInd'))
        load([filename '_BehavElectrDataLFP.mat'],'Stim');
        
        %% collecting information about the stimulation pulses
        timeBefStim = 3; % s
        timeAftStim = 3; % s
        totStim = length(unique(Stim.indStim));
        indStartStim = min(Stim.indStim);
        stims = cell(1,totStim);

        disp(['Gather information about stimulation'])
        for nStim = 1:totStim
            indCurStim = find(Stim.indStim == nStim+indStartStim-1);
            disp(['Stimulation ' num2str(indCurStim)]);
            
            if(isempty(indCurStim)) %% added 12.28.2015 to account for the stimulation error happened before the pulse definition
                continue;
            end
            nCurStimLen = length(indCurStim);
            disp(['**double check, were there ' num2str(nCurStimLen) ' pulses?**']);
            minStimLen = min(Stim.stopLfpInd(indCurStim) ...
                            - Stim.startLfpInd(indCurStim));
            % stimulation parameters
            stims{nStim}.startLfpInd = Stim.startLfpInd(indCurStim);
            stims{nStim}.stopLfpInd = Stim.stopLfpInd(indCurStim);
            stims{nStim}.indStim = nStim+indStartStim-1;
            stims{nStim}.indDiode = unique(Stim.indDiode{indCurStim(1)});
            stims{nStim}.indPulseInStim = Stim.indPulseInStim(indCurStim);
            stims{nStim}.pulsePeriod = unique(Stim.pulsePeriod(indCurStim));
            numTrain = max(Stim.indPulseInStim(indCurStim));

            stims{nStim}.spikesPulseNo = cell(numTrain,nTotClu);
            stims{nStim}.spikes = cell(numTrain,nTotClu);
            stims{nStim}.spikes20kHz = cell(numTrain,nTotClu);
            stims{nStim}.spikesThetaHil = cell(numTrain,nTotClu);
            stims{nStim}.spikesThetaLin = cell(numTrain,nTotClu);

            stims{nStim}.spikesPulseNoBef = cell(numTrain,nTotClu);
            stims{nStim}.spikesBef = cell(numTrain,nTotClu);
            stims{nStim}.spikesBef20kHz = cell(numTrain,nTotClu);
            stims{nStim}.spikesThetaHilBef = cell(numTrain,nTotClu);
            stims{nStim}.spikesThetaLinBef = cell(numTrain,nTotClu);

            stims{nStim}.spikesPulseNoAft = cell(numTrain,nTotClu);
            stims{nStim}.spikesAft = cell(numTrain,nTotClu);
            stims{nStim}.spikesAft20kHz = cell(numTrain,nTotClu);
            stims{nStim}.spikesThetaHilAft = cell(numTrain,nTotClu);
            stims{nStim}.spikesThetaLinAft = cell(numTrain,nTotClu);

            for nClu = 1:nTotClu
                if ~isempty(find(cluList1 == nClu, 1))
                    ind = Spike.res >= ...
                        Stim.startLfpInd(indCurStim(1)) - timeBefStim * lfpFreq + 1 ...
                        & Spike.res <= ...
                        Stim.stopLfpInd(indCurStim(end)) + timeAftStim * lfpFreq - 1 ...
                        & Spike.totclu == nClu;
                    res = Spike.res(ind);
                    res20k = Spike.res20kHz(ind); % added 20kHz on 4/1/2020
                    for nCurStim = 1:nCurStimLen
                        indStimTemp = indCurStim(nCurStim);
                        indTrain = Stim.indPulseInStim(indStimTemp);
                        % spikes during the stimulation
                        allSpInd = res >= Stim.startLfpInd(indStimTemp) ...
                            & res <= Stim.stopLfpInd(indStimTemp);
                        allSp20kHz = res20k(allSpInd);
                        allSpInd = res(allSpInd);
                        stims{nStim}.spikesPulseNo{indTrain,nClu} = ...
                            [stims{nStim}.spikesPulseNo{indTrain,nClu}; ...
                            Stim.indPulseInStim(indStimTemp)*ones(length(allSpInd),1)];
                        stims{nStim}.spikes{indTrain,nClu} = ...
                            [stims{nStim}.spikes{indTrain,nClu}; ...
                            allSpInd - Stim.startLfpInd(indStimTemp) + 1];
                        stims{nStim}.spikes20kHz{indTrain,nClu} = ...
                            [stims{nStim}.spikes20kHz{indTrain,nClu}; ...
                            allSp20kHz - Stim.startInd(indStimTemp) + 1];
                        stims{nStim}.spikesThetaHil{indTrain,nClu} = ...
                            [stims{nStim}.spikesThetaHil{indTrain,nClu}; ...
                            Track.thetaPhHilb(allSpInd)];
                        stims{nStim}.spikesThetaLin{indTrain,nClu} = ...
                            [stims{nStim}.spikesThetaLin{indTrain,nClu}; ...
                            Track.thetaPhLinInterp(allSpInd)];

                        % spikes before stimulation
                        allSpInd = res >= ...
                            Stim.startLfpInd(indStimTemp) - timeBefStim * lfpFreq + 1 ...
                            & res < Stim.startLfpInd(indStimTemp);
                        allSp20kHz = res20k(allSpInd);
                        allSpInd = res(allSpInd);
                        stims{nStim}.spikesPulseNoBef{indTrain,nClu} = ...
                            [stims{nStim}.spikesPulseNoBef{indTrain,nClu}; ...
                            Stim.indPulseInStim(indStimTemp)*ones(length(allSpInd),1)];
                        stims{nStim}.spikesBef{indTrain,nClu} = ...
                            [stims{nStim}.spikesBef{indTrain,nClu}; ...
                            allSpInd - Stim.startLfpInd(indStimTemp) + 1];
                        stims{nStim}.spikesBef20kHz{indTrain,nClu} = ...
                            [stims{nStim}.spikesBef20kHz{indTrain,nClu}; ...
                            allSp20kHz - Stim.startInd(indStimTemp) + 1];
                        stims{nStim}.spikesThetaHilBef{indTrain,nClu} = ...
                            [stims{nStim}.spikesThetaHilBef{indTrain,nClu}; ...
                            Track.thetaPhHilb(allSpInd)];
                        stims{nStim}.spikesThetaLinBef{indTrain,nClu} = ...
                            [stims{nStim}.spikesThetaLinBef{indTrain,nClu}; ...
                            Track.thetaPhLinInterp(allSpInd)];

                        % spikes after stimulation
                        allSpInd = res > Stim.stopLfpInd(indStimTemp) ...
                            & res <= Stim.stopLfpInd(indStimTemp)  + timeAftStim * lfpFreq - 1;
                        allSp20kHz = res20k(allSpInd);
                        allSpInd = res(allSpInd);
                        stims{nStim}.spikesPulseNoAft{indTrain,nClu} = ...
                            [stims{nStim}.spikesPulseNoAft{indTrain,nClu}; ...
                            Stim.indPulseInStim(indStimTemp)*ones(length(allSpInd),1)];
                        stims{nStim}.spikesAft{indTrain,nClu} = ...
                            [stims{nStim}.spikesAft{indTrain,nClu}; ...
                            allSpInd - Stim.startLfpInd(indStimTemp) + 1];
                        stims{nStim}.spikesAft20kHz{indTrain,nClu} = ...
                            [stims{nStim}.spikesAft20kHz{indTrain,nClu}; ...
                            allSp20kHz - Stim.startInd(indStimTemp) + 1];
                        stims{nStim}.spikesThetaHilAft{indTrain,nClu} = ...
                            [stims{nStim}.spikesThetaHilAft{indTrain,nClu}; ...
                            Track.thetaPhHilb(allSpInd)];
                        stims{nStim}.spikesThetaLinAft{indTrain,nClu} = ...
                            [stims{nStim}.spikesThetaLinAft{indTrain,nClu}; ...
                            Track.thetaPhLinInterp(allSpInd)];
                    end
                end
            end

            for nCurStim = 1:nCurStimLen
                indStimTemp = indCurStim(nCurStim);
                if(nCurStim == 1)
                    % running distance during stimulation
                    stims{nStim}.xMM = ...
                        Track.xMM(Stim.startLfpInd(indStimTemp)...
                            :Stim.startLfpInd(indStimTemp)+minStimLen);
                    stims{nStim}.xMMBef = ...
                        Track.xMM(Stim.startLfpInd(indStimTemp) ...
                            - timeBefStim * lfpFreq + 1 ...
                            :Stim.startLfpInd(indStimTemp)-1);
                    stims{nStim}.xMMAft = ...
                        Track.xMM(Stim.stopLfpInd(indStimTemp)+1 ...
                            :Stim.stopLfpInd(indStimTemp) + timeAftStim * lfpFreq - 1);

                    % running speed during stimulation
                    stims{nStim}.speed = ...
                        Track.speed_MMsec(Stim.startLfpInd(indStimTemp) ...
                            :Stim.startLfpInd(indStimTemp)+minStimLen);
                    stims{nStim}.speedBef = ...
                        Track.speed_MMsec(Stim.startLfpInd(indStimTemp) ...
                            - timeBefStim * lfpFreq + 1:Stim.startLfpInd(indStimTemp)-1);
                    stims{nStim}.speedAft = ...
                        Track.speed_MMsec(Stim.stopLfpInd(indStimTemp)+1 ...
                            :Stim.stopLfpInd(indStimTemp) + timeAftStim * lfpFreq - 1);

                    % hilbert theta phase during stimulation
                    stims{nStim}.thetaHil = ...
                        Track.thetaPhHilb(Stim.startLfpInd(indStimTemp) ...
                            :Stim.startLfpInd(indStimTemp)+minStimLen);
                    stims{nStim}.thetaHilBef = ...
                        Track.thetaPhHilb(Stim.startLfpInd(indStimTemp) ...
                        - timeBefStim * lfpFreq + 1:Stim.startLfpInd(indStimTemp)-1);
                    stims{nStim}.thetaHilAft = ...
                        Track.thetaPhHilb(Stim.stopLfpInd(indStimTemp)+1 ...
                            :Stim.stopLfpInd(indStimTemp) + timeAftStim * lfpFreq - 1);

                    % linear theta phase during stimulation
                    stims{nStim}.thetaLin = ...
                        Track.thetaPhLinInterp(Stim.startLfpInd(indStimTemp) ...
                            :Stim.startLfpInd(indStimTemp)+minStimLen);
                    stims{nStim}.thetaLinBef = ...
                        Track.thetaPhLinInterp(Stim.startLfpInd(indStimTemp) ...
                            - timeBefStim * lfpFreq + 1:Stim.startLfpInd(indStimTemp)-1);
                    stims{nStim}.thetaLinAft = ...
                        Track.thetaPhLinInterp(Stim.stopLfpInd(indStimTemp)+1 ...
                            :Stim.stopLfpInd(indStimTemp) + timeAftStim * lfpFreq - 1);

                    % eeg during stimulation
                    stims{nStim}.eeg = ...
                        Track.eeg(Stim.startLfpInd(indStimTemp) ...
                            :Stim.startLfpInd(indStimTemp)+minStimLen);
                    stims{nStim}.eegBef = ...
                        Track.eeg(Stim.startLfpInd(indStimTemp) ...
                            - timeBefStim * lfpFreq + 1:Stim.startLfpInd(indStimTemp)-1);
                    stims{nStim}.eegAft = ...
                        Track.eeg(Stim.stopLfpInd(indStimTemp)+1 ...
                            :Stim.stopLfpInd(indStimTemp) + timeAftStim * lfpFreq - 1);
                else
                    % running distance during stimulation
                    stims{nStim}.xMM = [stims{nStim}.xMM ...
                        Track.xMM(Stim.startLfpInd(indStimTemp) ...
                            :Stim.startLfpInd(indStimTemp)+minStimLen)];
                    stims{nStim}.xMMBef = [stims{nStim}.xMMBef ...
                        Track.xMM(Stim.startLfpInd(indStimTemp) ...
                            - timeBefStim * lfpFreq + 1:Stim.startLfpInd(indStimTemp)-1)];
                    stims{nStim}.xMMAft = [stims{nStim}.xMMAft ...
                        Track.xMM(Stim.stopLfpInd(indStimTemp)+1 ...
                            :Stim.stopLfpInd(indStimTemp) + timeAftStim * lfpFreq - 1)];

                    % running speed during stimulation
                    stims{nStim}.speed = [stims{nStim}.speed ...
                        Track.speed_MMsec(Stim.startLfpInd(indStimTemp) ...
                            :Stim.startLfpInd(indStimTemp)+minStimLen)];
                    stims{nStim}.speedBef = [stims{nStim}.speedBef ...
                        Track.speed_MMsec(Stim.startLfpInd(indStimTemp) ...
                            - timeBefStim * lfpFreq + 1:Stim.startLfpInd(indStimTemp)-1)];
                    stims{nStim}.speedAft = [stims{nStim}.speedAft ...
                        Track.speed_MMsec(Stim.stopLfpInd(indStimTemp)+1 ...
                            :Stim.stopLfpInd(indStimTemp) + timeAftStim * lfpFreq - 1)];

                    % hilbert theta phase during stimulation
                    stims{nStim}.thetaHil = [stims{nStim}.thetaHil ...
                        Track.thetaPhHilb(Stim.startLfpInd(indStimTemp) ...
                            :Stim.startLfpInd(indStimTemp)+minStimLen)];
                    stims{nStim}.thetaHilBef = [stims{nStim}.thetaHilBef ...
                        Track.thetaPhHilb(Stim.startLfpInd(indStimTemp) ...
                            - timeBefStim * lfpFreq + 1:Stim.startLfpInd(indStimTemp)-1)];
                    stims{nStim}.thetaHilAft = [stims{nStim}.thetaHilAft ...
                        Track.thetaPhHilb(Stim.stopLfpInd(indStimTemp)+1 ...
                            :Stim.stopLfpInd(indStimTemp) + timeAftStim * lfpFreq - 1)];

                    % linear theta phase during stimulation
                    stims{nStim}.thetaLin = [stims{nStim}.thetaLin ...
                        Track.thetaPhLinInterp(Stim.startLfpInd(indStimTemp)...
                            :Stim.startLfpInd(indStimTemp)+minStimLen)];
                    stims{nStim}.thetaLinBef = [stims{nStim}.thetaLinBef ...
                        Track.thetaPhLinInterp(Stim.startLfpInd(indStimTemp) ...
                            - timeBefStim * lfpFreq + 1:Stim.startLfpInd(indStimTemp)-1)];
                    stims{nStim}.thetaLinAft = [stims{nStim}.thetaLinAft ...
                        Track.thetaPhLinInterp(Stim.stopLfpInd(indStimTemp)+1 ...
                            :Stim.stopLfpInd(indStimTemp) + timeAftStim * lfpFreq - 1)];

                    % eeg during stimulation
                    stims{nStim}.eeg = [stims{nStim}.eeg ...
                        Track.eeg(Stim.startLfpInd(indStimTemp) ...
                            :Stim.startLfpInd(indStimTemp)+minStimLen)];
                    stims{nStim}.eegBef = [stims{nStim}.eegBef ...
                        Track.eeg(Stim.startLfpInd(indStimTemp) ...
                            - timeBefStim * lfpFreq + 1:Stim.startLfpInd(indStimTemp)-1)];
                    stims{nStim}.eegAft = [stims{nStim}.eegAft ...
                        Track.eeg(Stim.stopLfpInd(indStimTemp)+1 ...
                            :Stim.stopLfpInd(indStimTemp) + timeAftStim * lfpFreq - 1)];
                end
            end
        end
    
        if (mazeSection == 1)
            fprintf('\nOutput structures saved: %s\n\n', ...
                [filename '_DataStructure_mazeSection' num2str(mazeSection) ...
                '_TrialType' num2str(mazeTrialType) '.mat']);
            save([filename '_DataStructure_mazeSection' num2str(mazeSection) ...
                '_TrialType' num2str(mazeTrialType) '.mat'], 'stims', '-append');
        end
    end
end
