function Track = getTheta_mpfi(filename, lfpSampleRate, strExt)
% get theta phase

     if exist([filename '_BehavElectrDataLFP.mat'], 'file') == 2
        load([filename '_BehavElectrDataLFP.mat'], 'Processing');
        if(sum(Processing == 2) == 0)
            fprintf('\nRun getSpikes_v3 first.\n');
            return;
        end
        if(sum(Processing == 5) == 1) 
            return;
        end
        fprintf(...
            '\nLoad shankList and Track from %s file.\n', ...
            [filename '_BehavElectrDataLFP.mat']);
        load([filename '_BehavElectrDataLFP.mat'],'Track','Clu');
    end
    
    if exist([filename '_eeg_' num2str(lfpSampleRate) 'Hz.mat'], 'file') == 2
        load([filename '_eeg_' num2str(lfpSampleRate) 'Hz.mat'], ...
            'ThetaPhase_hilbert','ThetaAmp_hilbert','ThetaFreq_hilbert',...
            'ThetaPhase_linInterp','eeg','eegRaw',...
            'thetaLinear','shankListEEG');
    else
        fprintf('\nRun MakeTheta_v6 first');
    end
    
    numCluPerSh = zeros(1,length(shankListEEG));
    for i = 1:length(shankListEEG)
        numCluPerSh(i) = sum(Clu.shank == shankListEEG(i));
    end
    [~,indSh] = max(numCluPerSh);
    indSh = indSh(1);
%     indSh = find(shankListEEG == maxSh,1);
    
    % theta phase from the shank that has the most spikes     
    Track.(['thetaPhHilb' strExt]) = ThetaPhase_hilbert(:,indSh);
    Track.(['thetaAmpHilb' strExt]) = ThetaAmp_hilbert(:,indSh);
    Track.(['thetaFreqHilb' strExt]) = ThetaFreq_hilbert(:,indSh);
    Track.(['thetaPhLinInterp' strExt]) = ThetaPhase_linInterp(:,indSh);
    Track.(['eeg' strExt]) = eeg(:,indSh);
    Track.(['eegRaw' strExt]) = eegRaw(:,indSh);
    Track.(['thetaPeak_tAmpl' strExt]) = thetaLinear{indSh}.thetaPeak_tAmpl;
    Track.(['thetaTrough_tAmpl' strExt]) = thetaLinear{indSh}.thetaTrough_tAmpl;
    Track.(['thetaPtoTZeros_tAmpl' strExt]) = thetaLinear{indSh}.thetaPtoTZeros_tAmpl;
    Track.(['thetaTtoPZeros_tAmpl' strExt]) = thetaLinear{indSh}.thetaTtoPZeros_tAmpl;
    
    if(sum(Processing == 5) == 0)
        Processing = [Processing 5]; 
        % processing stage five, getting theta phase for Track
    end
    fprintf('Updated Track saved into the structure file: %s....\n',...
            [filename '_BehavElectrDataLFP.mat']);
    save([filename '_BehavElectrDataLFP.mat'], ...
            'Track', 'Processing', 'shankListEEG', 'indSh', '-append');