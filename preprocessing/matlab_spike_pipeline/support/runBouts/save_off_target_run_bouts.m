function run_bout_table = save_off_target_run_bouts(file_path, file_name, msess)
    % 12/14/2022 - Changing parameters to reflect method for identifying 
    % the trial run onset - running period where speed exceeded 10cm/s for
    % at least 0.3 seconds. The onset of the running bout is calculated by
    % tracing back to first time point where speed reached 1 cm/s

    % load behavelfp file 
    behave_lfp_path = [file_path '\' file_name(1:17) '_BehavElectrDataLFP.mat'];
    load(behave_lfp_path, 'Track');
    speed_MMsec = max(Track.speed_MMsec,0);
    
    % load alignRun file
    alignRun_path = [file_path '\' file_name '_alignRun_msess' num2str(msess) '.mat'];
    load(alignRun_path, 'trialsRun');
    startLfpInd = [];
    for i = 1:length(trialsRun.startLfpInd)
        startLfpInd = [startLfpInd trialsRun.startLfpInd(i)];
    end
    
    % running bout params
    run_time = 0.3;
    threshold_speed = 100; % in mm/s
      
    % initialize table entries
    rec_name = [];
    run_start_lfp = [];
    run_length_sec = [];
    mean_speed_run = [];
    acc_run_onset = [];
    precede_pause_length_sec = [];
    xMM_start_run_bout = [];
    
    % identify stim_periods 
    info_path = [file_path '\' file_name '_Info.mat'];
    load(info_path, 'beh');
    stim_lfp_indices = [];
    for i = 1:length(beh.stimPulseLfpInd)
        trial_i_stim_lfp = beh.stimPulseLfpInd{i};
        trial_i_stim_width = beh.stimPulseWidth{i};
        if ~isempty(trial_i_stim_lfp)
            for pulse_i = 1:length(trial_i_stim_lfp)
                stim_lfp_indices = [stim_lfp_indices ...
                                    trial_i_stim_lfp(pulse_i):trial_i_stim_lfp(pulse_i) + round(trial_i_stim_width(pulse_i))];
            end
        end
    end
    
    % identify trial onsets
    % collect lfps around trial running onset(-0.5s->0.5s windows)
    trial_onset_lfp_indices = [];
    for i = 1:length(trialsRun.startLfpInd)
        lfp_ro_i = trialsRun.startLfpInd(i);
        lfp_ro_i_indices = lfp_ro_i - 625 : lfp_ro_i + 625;
        trial_onset_lfp_indices = [trial_onset_lfp_indices lfp_ro_i_indices];
    end
    
    % identify running_bouts
    [run_length, run_start_lfp_indices, precede_pause_length] = find_long_runs(speed_MMsec, run_time, threshold_speed);
    
    for i = 1:length(run_length)
        lfp_indices_run_i = run_start_lfp_indices(i) : run_start_lfp_indices(i) + run_length(i);
        
        good_run = ~any(ismember(lfp_indices_run_i, stim_lfp_indices) | ismember(lfp_indices_run_i, trial_onset_lfp_indices));
        
        % also make sure run bout is not within 50cm of closest trial run
        % onset. 
        good_run_xMM = Track.xMM(run_start_lfp_indices(i));
        [~, lfp_ind_closest_trialRO] = min(abs(run_start_lfp_indices(i) - startLfpInd));
        closest_trialRO_xMM = Track.xMM(max(1,startLfpInd(lfp_ind_closest_trialRO)));
        
        if abs(good_run_xMM - closest_trialRO_xMM) < 500
            good_run = false; 
        end
        
        if good_run_xMM < 500
            good_run = false; % make sure run bout did occur around the beginning of the trial 
        end
        
        if good_run
                      
            rec_name = [rec_name ; string(file_name(1:17))];
            run_start_lfp = [run_start_lfp ; run_start_lfp_indices(i)];
            run_length_sec = [run_length_sec ; run_length(i)/1250];           
            mean_speed_run = [mean_speed_run ; mean(speed_MMsec(lfp_indices_run_i))];            
            acc_run_onset = [acc_run_onset; mean(Track.accel_MMsecSq(lfp_indices_run_i(1:min(625, length(lfp_indices_run_i)))))];
            precede_pause_length_sec = [precede_pause_length_sec; precede_pause_length(i)/1250];
            xMM_start_run_bout = [xMM_start_run_bout ; good_run_xMM];
        end
        
    end
    
    run_bout_table = table(rec_name, run_start_lfp, run_length_sec, mean_speed_run, acc_run_onset, precede_pause_length_sec, xMM_start_run_bout);
    
    save_path = 'Z:\Dinghao\code\lc-ca1-project\data\run_bouts\';
    save_path_py = [save_path file_name(1:17) '_run_bouts_py.csv'];
    save_path = [save_path file_name(1:17) '_run_bouts.mat'];
    save(save_path, 'run_bout_table');
    writetable(run_bout_table, save_path_py)
end

function [run_length, run_start_lfp_indices, precede_pause_length] = find_long_runs(speed_MMsec, run_time, threshold_speed)
      
      % find run bouts where speed was greater than 1cm/s and time spent running
      % was over run_time
      isRunning1 = speed_MMsec > 10;
      [which_run1, length_run1, start_indices1] = RunLength_M(isRunning1);
      which_run1 = which_run1';
      
      run_intervals = which_run1 == 1;
      long_intervals = length_run1/1250 > run_time;
      
      long_runs = find(run_intervals & long_intervals);
      
      run_start_lfp_indices_candidates = start_indices1(long_runs);
      run_length_candidates = length_run1(long_runs);
      
      % count consecutive bins that mouse speed exceeded threshold speed and reset sum 
      % whenever speed dipped below threshold_speed
      
      SpeedOverThresholdCumTime = rcumsum(speed_MMsec > threshold_speed);
      
      run_start_lfp_indices = [];
      run_length = [];
      precede_pause_length = [];
      
      % for each potential run bout, check if speed within the run bout
      % ever exceeded threshold speed for more 
      min_threshold_time = run_time*1250;
      
      for i = 1:length(run_start_lfp_indices_candidates)
          SpeedOverThresholdCumTime_runbout_i = SpeedOverThresholdCumTime(run_start_lfp_indices_candidates(i):...
                                                                          run_start_lfp_indices_candidates(i) + run_length_candidates(i));
                                                                      
          if any(SpeedOverThresholdCumTime_runbout_i > min_threshold_time)
              run_start_lfp_indices = [run_start_lfp_indices run_start_lfp_indices_candidates(i)];
              run_length = [run_length run_length_candidates(i)];
              if long_runs(i) == 1
                    precede_pause_length = [precede_pause_length  0];                   
              else
                    precede_pause_length = [precede_pause_length  length_run1(long_runs(i)-1)];
              end
          end
      end
end
