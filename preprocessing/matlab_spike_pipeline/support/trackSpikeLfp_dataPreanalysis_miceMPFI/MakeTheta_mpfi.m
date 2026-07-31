function MakeTheta_mpfi(filename, SampleRate, lfpSampleRate, ElecGp, nChannels, chSelect)

    % added 07/25/2017, removed _SpikesData.mat file to reduce data 
    % duplication, and also prevented passing large structures as
    % parameters
    if exist([filename '_BehavElectrDataLFP.mat'], 'file') == 2
        load([filename '_BehavElectrDataLFP.mat'], 'Processing');
        if(sum(Processing == 2) == 0)
            fprintf('\nRun getSpikes_v3 first.\n');
            return;
        end
        fprintf(...
            'Load shankList from %s file.\n', ...
            [filename '_BehavElectrDataLFP.mat']);
        load([filename '_BehavElectrDataLFP.mat'], 'shankList');
    end

    if(nargin == 5)
        chSelect = [];
    end

    thetaFiltParam.FreqRange = [4 16]; 
        %originally [4 25], changed to [4 16] on 2016.07.14
    thetaFiltParam.FilterOrd = 4;
    thetaFiltParam.Ripple = 20;
    
    eegOutFileName = ...
        [filename '_eeg_' num2str(lfpSampleRate) 'Hz.mat'];

    nStatus = 0;
    if exist(eegOutFileName,'file') == 2
        vars = whos('-file',eegOutFileName);
        if(ismember('ThetaPhase_linInterp',{vars.name})) 
            return; % has already calculated theta phase
        elseif(ismember('ThetaPhase_hilbert',{vars.name}))
            nStatus = 2;
        elseif(ismember('eeg',{vars.name}))
            nStatus = 1;
        end
        if(nStatus > 0)
            fprintf('\neeg load from a file %s...\n', ...
                    eegOutFileName);
            load(eegOutFileName,'eeg','shankListEEG','chSelect');
        end
    end

    if nStatus == 0 % need to calculate eeg
        % select one channel from each shank 
        if(isempty(chSelect)) % check for ChThetaPerSh.txt file
            fileNameTheta = [filename '-ChForThetaEst.txt'];
            if(exist(fileNameTheta, 'file') == 2)
                fid = fopen(fileNameTheta);
                chSelect = [];
                while ~feof(fid)
                    chSelect = [chSelect str2num(fgetl(fid))];
                end
                fclose(fid);

                % find shank list based on selected channels
                shankListTmp = [];
                chSelectTmp = [];
                for i = 1:length(chSelect)
                    for j = 1:length(shankList)
                        if(sum(ElecGp{shankList(j)} == chSelect(i)) ~= 0)
                            shankListTmp = [shankListTmp, shankList(j)];
                            chSelectTmp = [chSelectTmp chSelect(i)];
                            break;
                        end
                    end
                end
                shankListEEG = shankListTmp;
                chSelect = chSelectTmp;
            else % select the middle channel of each shank
                numShank = length(ElecGp)-1;
                chSelect = zeros(1,numShank);
                for i = 1:numShank
                    % if not - find appropriate channels 
                    chList = ElecGp{shankListEEG(i)};
                    chSelect(i) = chList(round(length(chList)/2));
                end
            end
        end
    
        % load data if not available from mat file (see above)
        datFileName = [filename '.dat'];
        listing = dir(datFileName);
        Nsamples = listing.bytes/2/nChannels;  % sec
        
        locRefFileName = [filename '-locRef.dat'];
        if(exist(locRefFileName) == 2)
            datFileName = locRefFileName;
        end
        if exist(datFileName, 'file') == 2
            fprintf('\n LFP loaded/filtered from: %s  ', datFileName);        
            fid = fopen(datFileName,'r');
            % filter before downsampling
            FreqRange = [0.2 lfpSampleRate];
            FilterOrd = 3;
            Ripple = 20;
            [b,a] = cheby2(FilterOrd, Ripple, FreqRange/(SampleRate/2));
            numCh = length(chSelect);
            eeg = [];
            for i = 1:numCh 
                fprintf('\nLoad channel %d from .dat file.\n', chSelect(i));
                datTmp = LoadDatFile_FL(fid,chSelect(i)+1,...
                            Nsamples,SampleRate,nChannels); 
                % chSelect(i)*Nsamples*2
                fEeg = filtfilt(b,a,datTmp');
                % remove constant term to avoid bias
                fEeg = fEeg - mean(fEeg);
                % downsample
                if(isempty(eeg))
                    eegRaw = zeros(length(fEeg),numCh);
                    eegRaw(:,1) = fEeg;
                    eegTmp = downsample(fEeg, round(SampleRate/lfpSampleRate));                
                    eeg = zeros(length(eegTmp),numCh);
                    eeg(:,i) = eegTmp;
                else
                    eegRaw(:,i) = fEeg;
                    eeg(:,i) = downsample(fEeg, round(SampleRate/lfpSampleRate));
                end                      
            end
            fclose(fid);
        else
            fprintf('\n no %s file was found - theta phase cannot be estimated.\n',...
                    datFileName);
            return;
        end
        save(eegOutFileName, 'eeg', 'eegRaw', 'shankListEEG','chSelect','-v7.3');
    end

    % estimate theta phase (hilbert transform)
    if nStatus < 2 
        ThetaPhase_hilbert = zeros(size(eeg,1),size(eeg,2));
        ThetaAmp_hilbert = zeros(size(eeg,1),size(eeg,2));
        ThetaFreq_hilbert = zeros(size(eeg,1),size(eeg,2));
        Eegf_hilbert = zeros(size(eeg,1),size(eeg,2));
        for i = 1:length(chSelect)
            fprintf('\n theta phase estimation with hilbert transform: ch %d',...
                    chSelect(i));
            [ThetaPhase_hilbert(:,i),ThetaAmp_hilbert(:,i), ~, ...
                ThetaFreqTmp,Eegf_hilbert(:,i)] ...
                = GetThetaPhase_hilbert(eeg(:,i),...
                                        lfpSampleRate, thetaFiltParam);
            ThetaFreq_hilbert(:,i) = [ThetaFreqTmp(1); ThetaFreqTmp];
            
            save(eegOutFileName, '-append', 'ThetaPhase_hilbert', ...
                'ThetaAmp_hilbert','ThetaFreq_hilbert','Eegf_hilbert');
        end
    end

    ThetaPhase_linInterp = zeros(size(eeg,1),size(eeg,2));
    thetaLinear = cell(1,length(chSelect));
    for i = 1:length(chSelect)
        fprintf('\n theta phase estimation with linear interpolation: ch %d',...
                chSelect(i)); 
        thetaLinear{i} = GetThetaPhase_peakValInterp_v2(eeg(1:end,i),...
                        lfpSampleRate, thetaFiltParam, 0);
        ThetaPhase_linInterp(:,i) = thetaLinear{i}.thetaLin;

        save(eegOutFileName, '-append', ...
                'ThetaPhase_linInterp','thetaLinear');   
    end

    return;