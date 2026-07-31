function Spike = getSpikeInfo_pix(filename, lfpSampleRate, strExt)
%% get information about spikes from Track and theta phase
%% changed this function to no longer look for xMM and speed variables

    if exist([filename '_BehavElectrDataLFP.mat'], 'file') == 2
        load([filename '_BehavElectrDataLFP.mat'], 'Processing');
        if(sum(Processing == 5) == 0)
            fprintf('\nRun getTheta_VR first.\n');
            return;
        end
        if(sum(Processing == 6) == 1)
            return;
        end
        fprintf(...
            '\nLoad shankList and Track from %s file.\n', ...
            [filename '_BehavElectrDataLFP.mat']);
        load([filename '_BehavElectrDataLFP.mat'], ...
            'shankList','Track','Spike','indSh');
    end
    
    if exist([filename '_eeg_' num2str(lfpSampleRate) 'Hz.mat'], 'file') == 2
        load([filename '_eeg_' num2str(lfpSampleRate) 'Hz.mat'], ...
        'ThetaPhase_hilbert','ThetaPhase_linInterp');
    else
        fprintf('\nRun MakeTheta_v6 first');
    end
    
    % assign behav data to each spike
    
    for sh = 1 : length(shankList)
        IDsh = shankList(sh);
        Spike.(['thPhaseHilb' strExt])(Spike.shank == IDsh, 1) = ...
            ThetaPhase_hilbert(Spike.res(Spike.shank == IDsh),indSh);                   
            % assign theta phase to each spike
        Spike.(['thPhaseInterp' strExt])(Spike.shank == IDsh, 1) = ...
            ThetaPhase_linInterp(Spike.res(Spike.shank == IDsh),indSh);                   
            % assign theta phase to each spike 
    end
    Spike.xMM = Track.xMM(Spike.res);
    Spike.speed_MMsec = Track.speed_MMsec(Spike.res);
    Spike.accel_MMsecSq = Track.accel_MMsecSq(Spike.res);
          
    if(sum(Processing == 6) == 0)
        Processing = [Processing 6]; 
        % processing stage six, getting theta phase and Track information
        % for Spike
    end
    fprintf('\nUpdated Track saved into the structure file: %s....\n',...
            [filename '_BehavElectrDataLFP.mat']);
    save([filename '_BehavElectrDataLFP.mat'], ...
            'Spike', 'Processing', '-append');

end