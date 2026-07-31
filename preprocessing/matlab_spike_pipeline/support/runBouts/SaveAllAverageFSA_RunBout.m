function SaveAllAverageFSA_RunBout
    RecordingList_dinghao;
    
    indRec = [];
    indNeu = [];
    meanFR_RunBout = [];
    
    fsa_base_path = 'Z:\Dinghao\code\lc-ca1-project\data\run_bouts\fsa_run_bouts\';
    
    for rec_i = 1:size(listRecordingsActiveLickPath,1)
        disp(listRecordingsActiveLickFileName(rec_i,1:17));
        fsa_path_i = [fsa_base_path listRecordingsActiveLickFileName(rec_i,1:17) '_BefRunBout0.mat'];
        load(fsa_path_i, 'filteredSpikeArrayRunBoutOnSet','timeStepRun');
        
        for neu_j = 1:length(filteredSpikeArrayRunBoutOnSet)
            fsa_neu_j = filteredSpikeArrayRunBoutOnSet{neu_j};
            
            indRec = [indRec ; rec_i];
            indNeu = [indNeu ; neu_j];
            meanFR_RunBout = [meanFR_RunBout ; mean(fsa_neu_j,1)];
        end
    end
    
    save_path= [fsa_base_path 'BefRunBout0_all.mat'];
    save(save_path, 'indRec', 'indNeu', 'meanFR_RunBout');
    
end