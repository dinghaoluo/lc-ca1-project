function IdentifyRiseDownRunBouts_all(rerun)
    RecordingList_dinghao;
    RiseDownTable_all = [];
    
    km2_path = 'Z:\Yingxue\Draft\PV\PyramidalInitPeak\2\initPeakPyrAllRec_km2.mat';
    load(km2_path, 'FRProfileMean');
    
    num_run_bout_per_rec_path =  'Z:\Dinghao\code\lc-ca1-project\data\run_bouts\num_run_bouts_per_rec.mat';
    load(num_run_bout_per_rec_path, 'num_run_bouts_table');
    
    for i = 1:size(listRecordingsActiveLickPath,1)
        if num_run_bouts_table.num_run_bouts(num_run_bouts_table.rec_index == i) < 15
            continue;
        end
        
        if str2num(listRecordingsActiveLickFileName(i,2:4)) > 56
            continue;
        end
        
        disp(listRecordingsActiveLickFileName(i,1:17));
        table_path = ['Z:\Raphael_tests\Code\matlabAnalysisRaphi\RiseDownOffTargetRunBouts\RiseDownIDRunBout\'...
                       listRecordingsActiveLickFileName(i,1:17) '_RiseDownID_RunBout.mat'];
        
        if rerun || ~exist(table_path, 'file')
            RiseDownTable = IdentifyRiseDownRunBouts(listRecordingsActiveLickFileName(i,1:17), FRProfileMean.indFR0to1,FRProfileMean.indFRBefRun);
        else
           load(table_path, 'RiseDownTable');
        end
       
        RiseDownTable_all = [RiseDownTable_all;RiseDownTable];
        
    end
    
    save_path = 'Z:\Raphael_tests\Code\matlabAnalysisRaphi\RiseDownOffTargetRunBouts\RiseDownIDRunBout\';
    save_path = [save_path 'RiseDownID_RunBout_all.mat'];
    save(save_path, 'RiseDownTable_all');
end