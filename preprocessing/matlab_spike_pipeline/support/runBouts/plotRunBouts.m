function plotRunBouts(rec_index)
    RecordingList_dinghao;
    
    save_path_base = ['Z:\Dinghao\code\lc-ca1-project\data\run_bouts\fsa_run_bouts_plots\'];
    if ~exist(save_path_base, 'dir')
        mkdir(save_path_base);
    end
    
    run_bout_path = 'Z:\Dinghao\code\lc-ca1-project\data\run_bouts\A029r-20220623-03_run_bouts.mat';
    load(run_bout_path, 'run_bout_table');
    
    behavelfp_path = [listRecordingsActiveLickPath(rec_index,:) '\' listRecordingsActiveLickFileName(rec_index,1:17) '_BehavElectrDataLFP.mat'];
    load(behavelfp_path, 'Track', 'Laps');
    lickLFP = cell2mat(Laps.lickLfpInd);
    speed_MMsec = Track.speed_MMsecAll;
    speed_MMsec(speed_MMsec < 0) = 0;
    
    % load alignRun file
    alignRun_path = [listRecordingsActiveLickPath(rec_index,:) '\'  listRecordingsActiveLickFileName(rec_index,:) '_alignRun_msess' num2str(mazeSessionActiveLick(rec_index)) '.mat'];
    load(alignRun_path, 'trialsRun');
    startLfpInd = [];
    endLfpInd = [];
    for i = 1:length(trialsRun.startLfpInd)
        startLfpInd = [startLfpInd trialsRun.startLfpInd(i)];
        endLfpInd = [endLfpInd trialsRun.endLfpInd(i)];
    end
    
    for t = 2:3:length(endLfpInd) - 2
        CreateFig();
        hold on;
        lfp_indices_t = startLfpInd(t)-1250: min(endLfpInd(t+2) + 1250, length(speed_MMsec));
        lap_start = lfp_indices_t(1);
        plot(speed_MMsec(lfp_indices_t));
        
        startLfpInd_t = startLfpInd(ismember(startLfpInd, lfp_indices_t));
        
        for s = 1:length(startLfpInd_t)
            plot([startLfpInd_t(s) startLfpInd_t(s)] - lap_start, ylim, 'r--');
        end
        
        run_bout_t = run_bout_table.run_start_lfp(ismember(run_bout_table.run_start_lfp, lfp_indices_t));
        
        for r = 1:length(run_bout_t)
            plot([run_bout_t(r) run_bout_t(r)] - lap_start, ylim, 'g--');
        end
        
        licks_t = lickLFP(ismember(lickLFP, lfp_indices_t));
        ylimits  = ylim; 
        for q = 1:length(licks_t)
            plot([licks_t(q) licks_t(q)]  - lap_start, [ylimits(2) - 30, ylimits(2)], 'm-');
        end
        
        set(gcf, 'Position', [-86 527 1444 251]);
        save_path = [save_path_base 'trials' num2str(t) 'to' num2str(t+2) '.png'];
        saveas(gcf, save_path);
        
        if t == 152
            disp('here');
        end
        
        close all;
        
    end
    
end
