function  Stim = getStim_mpfi_smTr( filename, lfpFreq )
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
    
    load([filename 'BTDT.mat']);
    Stim = [];
    
    if(~isfield(behEventsTdt,'stimOn'))
        return;
    end
    
    indStim = find(behEventsTdt.stimOn(:,5) == -1);
    
    indPulse = 1;
    for i = 1:length(indStim)
        startInd = behEventsTdt.stimOn(indStim(i),3);
        stopInd = behEventsTdt.stimOn(indStim(i),3)+...
            behEventsTdt.stimOn(indStim(i),4)-1;
        for j = startInd:stopInd
            Stim.startLfpInd(indPulse) = behEventsTdt.stimPulse(j,2);
            Stim.stopLfpInd(indPulse) = behEventsTdt.stimPulse(j,2)...
                    +behEventsTdt.stimPulse(j,7);
            Stim.startInd(indPulse) = behEventsTdt.stimPulse(j,1);
            Stim.stopInd(indPulse) = behEventsTdt.stimPulse(j,1)...
                    +behEventsTdt.stimPulse(j,6);
            Stim.indStim(indPulse) = behEventsTdt.stimPulse(j,3);
            Stim.indPulseInStim(indPulse) = behEventsTdt.stimPulse(j,5); 
                % pulse no. within a stimulation
            Stim.indDiode{indPulse} = behEventsTdt.stimOnDiode{j}.indDiode;
            Stim.pulsePeriod(indPulse) = behEventsTdt.stimPulse(j,4); %ms
            indPulse = indPulse+ 1;
        end
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
