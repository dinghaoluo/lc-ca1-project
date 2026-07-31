function AlignSpikesRunOnset(path, fileName)
% plot spikes according to the run onset
% E.g.: AlignSpikesRunOnset('./','A002-20181005-01_BehavElectrDataLFP')

    GlobalConst;
    paramR.minRunTimePerSeg = 3*sampleFq;
    
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    load([path fileName '.mat'],'Track','Spike','Clu','Laps');
    
    indSpeed = Track.speed_MMsec > minSpeed;
    runSegInd = diff(indSpeed);
    runOn = find(runSegInd == 1);
    runOff = find(runSegInd == -1);
    lenRunOn = length(runOn);
    lenRunOff = length(runOff);
    if(runOn(1) <= runOff(1))
        lenRunSeg = min(lenRunOn,lenRunOff);
        runSegTmp = (runOff(1:lenRunSeg) - runOn(1:lenRunSeg));
        indSel = runSegTmp >= paramR.minRunTimePerSeg;
        runOnSeg = runOn(indSel) + 1;
        runOffSeg = runOff(indSel) + 1;
    else
        lenRunSeg = min(lenRunOn,lenRunOff-1);
        runSegTmp = (runOff(2:lenRunSeg) - runOn(1:lenRunSeg-1));
        indSel = runSegTmp >= paramR.minRunTimePerSeg;
        runOnSeg = runOn(indSel) + 1;
        runOffSeg = runOff(1+indSel) + 1;
    end

    indexFileName = findstr(fileName, '_');
    load([fileName(1:indexFileName(1)-1) '-param.mat']);    
    voltNorm = (1/(2^device.RecSyst.ADresolution)) ...
        * device.RecSyst.ADvoltRange ...
        * (1/device.RecSyst.amplification) * 1000;
    
    nTotClu = max(Clu.totClu);
    cluList1 = Clu.totClu; 
    trackLen = unique(Laps.trackLen);
    
    for i = 1:length(runOnSeg)
        RunSeg{i}.onLfpInd = runOnSeg(i);
        RunSeg{i}.offLfpInd = runOffSeg(i);
        RunSeg{i}.Nsamples = runOffSeg(i) - runOnSeg(i) + 1;
        RunSeg{i}.xMM = Track.xMM(runOnSeg(i):runOffSeg(i));
        diffXMM = diff(RunSeg{i}.xMM);
        zeroCross = find(diffXMM < ...
            -min(trackLen)/2);
        if(~isempty(zeroCross))
            RunSeg{i}.xMMAdj(1:zeroCross(1)) = RunSeg{i}.xMM(1:zeroCross(1));
            for n = 1:length(zeroCross)
                if(n == length(zeroCross))
                    RunSeg{i}.xMMAdj(zeroCross(n)+1:length(RunSeg{i}.xMM))...
                        = RunSeg{i}.xMM(zeroCross(n)+1:end) ...
                        + Laps.trackLen(Track.lapID(runOnSeg(i)+zeroCross(n)-1));
                else
                    RunSeg{i}.xMMAdj(zeroCross(n)+1:zeroCross(n+1))...
                        = RunSeg{i}.xMM(zeroCross(n)+1:zeroCross(n+1)) ...
                          + Laps.trackLen(Track.lapID(runOnSeg(i)+zeroCross(n)-1));
                end
            end
            RunSeg{i}.xMMAdj = RunSeg{i}.xMMAdj - RunSeg{i}.xMMAdj(1);
        else
            RunSeg{i}.xMMAdj = RunSeg{i}.xMM - RunSeg{i}.xMM(1);
        end
        if(length(RunSeg{i}.xMMAdj) ~= length(RunSeg{i}.xMM))
           disp('error'); 
        end
        
        RunSeg{i}.speed = Track.speed_MMsec(runOnSeg(i):runOffSeg(i));
        RunSeg{i}.accel = Track.accel_MMsecSq(runOnSeg(i):runOffSeg(i));
        RunSeg{i}.lapID = Track.lapID(runOnSeg(i):runOffSeg(i));
        RunSeg{i}.behavType = Track.behavType(runOnSeg(i):runOffSeg(i));
        
        RunSeg{i}.thetaHil = Track.thetaPhHilb(runOnSeg(i):runOffSeg(i));
        RunSeg{i}.cumsumThetaHil = unwrap(RunSeg{i}.thetaHil);
        RunSeg{i}.thetaLin = Track.thetaPhLinInterp(runOnSeg(i):runOffSeg(i));
        RunSeg{i}.cumsumThetaLin = unwrap(RunSeg{i}.thetaLin);
        RunSeg{i}.eeg = Track.eeg(runOnSeg(i):runOffSeg(i)) .*voltNorm;
        
        %  detect individual theta cycles and determine the 
        %  amplitude of each (in mV!!!!!!!!!!!!!!!!!!!)
        % theta peak and trough amplitude
        ind = Track.thetaPeak_tAmpl(:,1) >= runOnSeg(i) ...
            & Track.thetaPeak_tAmpl(:,1) <= runOffSeg(i);
        RunSeg{i}.thetaPeak_tAmpl(:,1) = ...
            Track.thetaPeak_tAmpl(ind,1) - runOnSeg(i) + 1;
        RunSeg{i}.thetaPeak_tAmpl(:,2) = ...
            Track.thetaPeak_tAmpl(ind,2) .* voltNorm;

        ind = Track.thetaTrough_tAmpl(:,1) >= runOnSeg(i) ...
            & Track.thetaTrough_tAmpl(:,1) <= runOffSeg(i);
        RunSeg{i}.thetaTrough_tAmpl(:,1) = ...
            Track.thetaTrough_tAmpl(ind,1) - runOnSeg(i) + 1;
        RunSeg{i}.thetaTrough_tAmpl(:,2) = ...
            Track.thetaTrough_tAmpl(ind,2) .* voltNorm;

        ind = Track.thetaPtoTZeros_tAmpl(:,1) >= runOnSeg(i) ...
            & Track.thetaPtoTZeros_tAmpl(:,1) <= runOffSeg(i);
        RunSeg{i}.thetaPtoTZeros_tAmpl(:,1) = ...
            Track.thetaPtoTZeros_tAmpl(ind,1) - runOnSeg(i) + 1;
        RunSeg{i}.thetaPtoTZeros_tAmpl(:,2) = ...
            Track.thetaPtoTZeros_tAmpl(ind,2) .* voltNorm;

        ind = Track.thetaTtoPZeros_tAmpl(:,1) >= runOnSeg(i) ...
            & Track.thetaTtoPZeros_tAmpl(:,1) <= runOffSeg(i);
        RunSeg{i}.thetaTtoPZeros_tAmpl(:,1) = ...
            Track.thetaTtoPZeros_tAmpl(ind,1) - runOnSeg(i) + 1;
        RunSeg{i}.thetaTtoPZeros_tAmpl(:,2) = ...
            Track.thetaTtoPZeros_tAmpl(ind,2) .* voltNorm;
        
        RunSeg{i}.IDShank = zeros(nTotClu,1);
        RunSeg{i}.locClu = zeros(nTotClu,1);
        RunSeg{i}.spikes = cell(1,nTotClu);
        RunSeg{i}.spikes20kHz = cell(1,nTotClu);
        RunSeg{i}.spikesCumSumThetaHil = cell(1,nTotClu);
        RunSeg{i}.spikesCumSumThetaLin = cell(1,nTotClu);
        RunSeg{i}.spikesThetaHil = cell(1,nTotClu);
        RunSeg{i}.spikesThetaLin = cell(1,nTotClu);
        RunSeg{i}.spikesMM = cell(1,nTotClu);
        RunSeg{i}.spikesMMAdj = cell(1,nTotClu);
        RunSeg{i}.spikesSpeed = cell(1,nTotClu);

        % times and spikes within each segment
        allSpTrain = [];
        allSpClu = [];
        totCluSp = zeros(nTotClu,1);
        for nclu = 1 : nTotClu
            RunSeg{i}.IDShank(nclu) = Clu.shank(nclu);
            RunSeg{i}.locClu(nclu) = Clu.localClu(nclu);

            if ~isempty(find(cluList1 == nclu, 1))
                ind = Spike.res >= runOnSeg(i) ...
                    & Spike.res <= runOffSeg(i) ...
                    & Spike.totclu == nclu;
                totCluSp(nclu) = totCluSp(nclu) + sum(ind);
                allSp = Spike.res(ind) - runOnSeg(i) + 1;
                RunSeg{i}.spikes{nclu} = allSp;
                % res at 20kHz
                allSp20kHz = Spike.res20kHz(ind) ...
                    - round(runOnSeg(i)/sampleFqOri*sampleFq) + 1;
                RunSeg{i}.spikes20kHz{nclu} = allSp20kHz;

                % spikes in theta phase 
                RunSeg{i}.spikesCumSumThetaHil{nclu} = ...
                    RunSeg{i}.cumsumThetaHil(allSp);
                RunSeg{i}.spikesCumSumThetaLin{nclu} = ...
                    RunSeg{i}.cumsumThetaLin(allSp);
                RunSeg{i}.spikesThetaHil{nclu} = ...
                    RunSeg{i}.thetaHil(allSp);
                RunSeg{i}.spikesThetaLin{nclu} = ...
                    RunSeg{i}.thetaLin(allSp);

                % spike distance 
                RunSeg{i}.spikesMM{nclu} = ...
                    RunSeg{i}.xMM(allSp);
                RunSeg{i}.spikesMMAdj{nclu} = ...
                    RunSeg{i}.xMMAdj(allSp);
                RunSeg{i}.spikesSpeed{nclu} = ...
                    RunSeg{i}.speed(allSp);

                allSpTrain = [allSpTrain; allSp];
                allSpClu = [allSpClu; nclu*ones(length(allSp),1)];

            end
        end
        RunSeg{i}.allSpTrain = [allSpTrain allSpClu];
    end
    
    trackLen = Laps.trackLen;
    trials = RunSeg;
    
    fprintf('\nOutput structures saved: %s\n\n', ...
        [fileName(1:indexFileName(1)-1) '_DataStructure_mazeSection1'  ...
        '_TrialType_RunOnSet.mat']);
    save([fileName(1:indexFileName(1)-1) '_DataStructure_mazeSection1'  ...
        '_TrialType_RunOnSet.mat'], ...
        'trials', 'paramR', 'trackLen', '-v7.3');
