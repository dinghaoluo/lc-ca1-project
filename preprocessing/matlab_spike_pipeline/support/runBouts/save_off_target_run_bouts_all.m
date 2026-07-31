function save_off_target_run_bouts_all
    RecordingList_dinghao;
        
    run_bout_table_all = [];
    for i = 1:size(listRecordingsActiveLickPath,1)
        disp(listRecordingsActiveLickFileName(i,1:17));
        % recompute each recording so the combined table comes from one pass
        run_bout_i = save_off_target_run_bouts(listRecordingsActiveLickPath(i,:), listRecordingsActiveLickFileName(i,:), mazeSessionActiveLick(i));
        run_bout_table_all = [run_bout_table_all ; run_bout_i];
        
    end
    
    save_path = 'Z:\Dinghao\code\lc-ca1-project\data\run_bouts\overall_run_bouts.mat';
    save(save_path, 'run_bout_table_all');
    
end
