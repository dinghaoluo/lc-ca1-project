function verify_running_bouts(file_path, file_name)
    run_bout_path = 'Z:\Raphael_tests\Code\matlabAnalysisRaphi\RiseDownOffTargetRunBouts\run_bouts\';
    run_bout_path = [run_bout_path file_name(1:16) '_run_bouts.mat'];
    load(run_bout_path, 'run_bout_table');
    
    behave_lfp_path = [file_path '\' file_name(1:16) '_BehavElectrDataLFP.mat'];
    load(behave_lfp_path, 'Track', 'Laps');
    pumpLfp = downsample(cell2mat(Laps.pumpLfpInd),2);
    speed_MMsec = max(Track.speed_MMsec,0);
    
    save_path_base = 'Z:\Raphael_tests\Code\matlabAnalysisRaphi\RiseDownOffTargetRunBouts\verify_run_bout_plots\';
    save_path_base = [save_path_base file_name(1:16) '\'];
    if ~exist(save_path_base, 'dir')
        mkdir(save_path_base);
    end
    
    run_start_lfp_indices = run_bout_table.run_start_lfp;
    
    for i = 1:size(run_bout_table, 1)      
        run_lfp_indices = run_bout_table.run_start_lfp(i) - 5*1250 : run_bout_table.run_start_lfp(i) + 5*1250;
        
        pump_run_bout_i = pumpLfp(ismember(pumpLfp, run_lfp_indices));
        close_run_bouts = run_start_lfp_indices(ismember(run_start_lfp_indices, run_lfp_indices) & run_start_lfp_indices ~= run_bout_table.run_start_lfp(i));
              
        acc_ro = run_bout_table.acc_run_onset(i);
        prec_pause = run_bout_table.precede_pause_length_sec(i);
        mean_speed = run_bout_table.mean_speed_run(i);
        
        CreateFig();
        hold on;
        plot(run_lfp_indices/1250, speed_MMsec(run_lfp_indices));
        start_run = (run_lfp_indices(1) + 5*1250)/1250;     
        plot([start_run, start_run], [0 max(speed_MMsec(run_lfp_indices))], 'r--');         
        pump_y = 10+max(speed_MMsec(run_lfp_indices));
        scatter(pump_run_bout_i/1250, pump_y*ones(size(pump_run_bout_i)), 'b.');
        for j = 1:size(close_run_bouts)
            plot([close_run_bouts(j), close_run_bouts(j)]/1250, [0 max(speed_MMsec(run_lfp_indices))], 'm--');   
        end
        
        title_string = string(file_name(1:16)) +...
                       "\nmean speed: " + num2str(mean_speed) + ...
                       "\nprec pause: " + num2str(prec_pause) + ...
                       "\nacc_ro: " + num2str(acc_ro);
        title(char(compose(title_string)));
        
        save_path = [save_path_base file_name(1:16) '_run_bout' num2str(i) '.png'];
        saveas(gcf , save_path);
        close all;
    end
end