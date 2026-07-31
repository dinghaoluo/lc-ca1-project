function [p I percentActive] = salt_wrapper(path, filename, neuron_id,baseline_time,time_resolution,window_size)
% baseline_time = amount of time before pulse to collect latency histograms     
% time_resolution = time in between entries of binary matrix to assess if a
% spike has occurred or not
% window_size = windows are non-overlapping, one window corresponds to one
% distribution of latencies used to construct null distribution. should
% have as many as possible before window that occurs after stimulation
% 
% 

    stimEffPath = [path '\' filename '_stimEff.mat'];
    if ~isfile(stimEffPath)
        disp('need to run stim_Effect first');
        return;
    end
    
    load(stimEffPath, 'stimResp');
    load([filename '.mat'], 'cluList');
    shankID = cluList.shank;
    localClu = cluList.localClu;
    GlobalConst;
    
    pulseBefore = stimResp.spPerPulseBef;
    pulseAfter = stimResp.spPerPulseAft;
    
    pBefNeuron = pulseBefore(:,neuron_id);
    pAftNeuron = pulseAfter(:,neuron_id);
    
    % Sometimes we do theta stimulations that gets coded as tagging - need
    % to not count these
    
    stimDur = stimResp.stimDuration;
    short_stims = find(stimDur <= 5);
    first_tag_stim = short_stims(1);
    
    % Create N x M binary matrix, where N is # of pulses(trials)
    % and M is the length of baseline period. 
    
    num_time_steps_baseline = baseline_time/time_resolution;
   
    num_stims = length(short_stims);
    num_pulses_per_stim = length(pBefNeuron{first_tag_stim});
    num_trials = num_stims * num_pulses_per_stim;
    
    spt_baseline = zeros(num_trials, num_time_steps_baseline);   
    spt_test = zeros(num_trials, num_time_steps_baseline);
    
    trial_index = 1;
    total_lookback_time_ms = 3000; 
    
    for s = first_tag_stim:first_tag_stim + num_stims -1
        stim_i_bef = pBefNeuron{s}; 
        stim_i_aft = pAftNeuron{s};
        
        start_seed_time_ms= baseline_time*rand(1); %250; %500*rand(1);
        
        for p = 1:num_pulses_per_stim
            % modified protocol for getting baseline- 
            % we now choose non-overlapping, 0.5s window from time period
            % 3s before the first pulse in the stim train
            
            pulse_i_bef = stim_i_bef{1};
            pulse_i_aft = stim_i_aft{p};
            
            time_bef_window_start_ms = start_seed_time_ms + (baseline_time*(p-1));
            time_bef_window_end_ms = time_bef_window_start_ms + baseline_time;
            
            for i = 1:length(pulse_i_bef)
                time_spike_ms = -1*(pulse_i_bef(i)/1250)*1000;
                abs_time_spike_ms = total_lookback_time_ms - time_spike_ms;
                if abs_time_spike_ms > time_bef_window_start_ms && abs_time_spike_ms < time_bef_window_end_ms
                    index_spike = max(1,round((abs_time_spike_ms - time_bef_window_start_ms)/time_resolution));                   
                    spt_baseline(trial_index, index_spike) = 1;
                end             
            end
            
            for i = 1:length(pulse_i_aft)
                time_spike_ms = (pulse_i_aft(i)/1250)*1000;
                index_spike = round(time_spike_ms/time_resolution);                
                if index_spike > 0 && index_spike < num_time_steps_baseline
                    spt_test(trial_index, index_spike) = 1;
                end
            end
            trial_index = trial_index + 1;
                     
        end
    end

 [p I] = salt(spt_baseline, spt_test,time_resolution/1000,window_size/1000);  
 
 % Determine number of trials where cell was activated by light pulse
 nmbn = round(window_size/time_resolution);
 spt_window = spt_test(:,1:nmbn);
 [ntrials,nbins] = find(spt_window ~= 0);
 ntrialsActive = length(unique(ntrials));
 percentActive = ntrialsActive/size(spt_window,1);
 
 [nrow,ncol] = size(spt_baseline);
 [base_row, base_col] = find(spt_baseline);
 [test_row, test_col] = find(spt_test);
 
%h =  length(findobj('type','figure'));
% h = 1; 
% figure(h);
%  clf;
%  hold on;
%  sz = 10;
%  scatter(base_col, base_row, sz, 'filled');
%  scatter(test_col + ncol,test_row,sz,'filled');
%  xlim([0, ncol + 100]);
%  y1=get(gca,'ylim');
%  for i = 0:window_size/time_resolution : baseline_time
%      plot([i i],y1); 
%  end
end

%                 index_spike = round((baseline_time - time_spike_ms)/time_resolution);
%                 
%                 if index_spike > 0 && index_spike < num_time_steps_baseline
%                     spt_baseline(trial_index, index_spike) = 1;
%                 end
