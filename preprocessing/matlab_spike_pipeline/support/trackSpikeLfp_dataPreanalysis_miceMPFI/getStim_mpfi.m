function  Stim = getStim_mpfi( filename, lfpFreq )
% get information about the optical stimulation

    % added 08/14/2018, removed _SpikesData.mat file to reduce data 
    % duplication, and also prevented passing large structures as
    % parameters
    if exist([filename '_BehavElectrDataLFP.mat'], 'file') == 2
        load([filename '_BehavElectrDataLFP.mat'], 'Processing');
        if(sum(Processing == 7) == 1)
            return;
        end
    end
    
    load([filename 'PTDT.mat']);
    
    for i = 1:size(stimEventsTdt.pulse,1)
        Stim.startLfpInd(i,1) = stimEventsTdt.pulse(i,2);
        Stim.stopLfpInd(i,1) = stimEventsTdt.pulse(i,2)...
                    +round(stimEventsTdt.pulse(i,3)*lfpFreq/1000);
        Stim.indStim(i,1) = stimEventsTdt.pulse(i,7);
        Stim.indTrainInStim(i,1) = stimEventsTdt.pulse(i,8); 
                % train no. within a stimulation
        Stim.indDiode(i,1) = stimEventsTdt.pulse(i,6);
        Stim.indPulseInStim(i,1) = stimEventsTdt.pulse(i,5); 
                % pulse number within a train
        Stim.pulsePeriod(i,1) = stimEventsTdt.pulse(i,3);
    end
    
    if(sum(Processing == 7) == 0)
        Processing = [Processing 7]; 
        % processing stage six, getting theta phase and Track information
        % for Spike
    end
    fprintf('\nStim saved into the structure file: %s....\n',...
            [filename '_BehavElectrDataLFP.mat']);
    save([filename '_BehavElectrDataLFP.mat'], ...
            'Stim', '-append');
end
