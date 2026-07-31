function HowManyRunBouts
    % creates a table indicating how many off-target run bouts each
    % recording has
    
    RecordingList_dinghao;
    run_bout_path = 'Z:\Dinghao\code\lc-ca1-project\data\run_bouts\overall_run_bouts.mat';
    load(run_bout_path, 'run_bout_table_all');
    
    rec_index = [];
    num_run_bouts = [];
    
    for i = 1:size(listRecordingsActiveLickFileName,1)
        file_name_i = listRecordingsActiveLickFileName(i,1:17);
        num_run_bouts_i = sum(strcmp(run_bout_table_all.rec_name, string(file_name_i)) & run_bout_table_all.xMM_start_run_bout > 500);
        rec_index = [rec_index ; i];
        num_run_bouts = [num_run_bouts ; num_run_bouts_i];     
    end
    
    num_run_bouts_table = table(rec_index, num_run_bouts);
    
    save_path = 'Z:\Dinghao\code\lc-ca1-project\data\run_bouts\num_run_bouts_per_rec.mat';
    save_path_py = 'Z:\Dinghao\code\lc-ca1-project\data\run_bouts\num_run_bouts_per_rec_py.csv';
    
    save(save_path, 'num_run_bouts_table');
    writetable(num_run_bouts_table, save_path_py); 
end