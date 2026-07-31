function RiseDownTable = IdentifyRiseDownRunBouts(filename, indFR0to1,indFRBefRun)
    fsa_path = 'Z:\Raphael_tests\Code\matlabAnalysisRaphi\RiseDownOffTargetRunBouts\fsa_run_bouts\';
    fsa_path = [fsa_path filename '_BefRunBout0.mat'];
    load(fsa_path, 'filteredSpikeArrayRunBoutOnSet','paramC','timeStepRun');
    
    GlobalConstFq;
    p = 95; % overwrite GlobalConstFq to just use p value of 95
    rec_name = [];
    neu_id = [];
    ratio0to1BefRun = [];
    isRise = [];
    isDown = [];
    
    for i = 1:length(filteredSpikeArrayRunBoutOnSet)
        fsa_i = filteredSpikeArrayRunBoutOnSet{i};
        mean_fsa_i = mean(fsa_i,1);
        
        ratio0to1BefRun_i = mean(mean_fsa_i(indFR0to1))/mean(mean_fsa_i(indFRBefRun));        
        ratioAftBefShuf = neuActivityShuffle_runbout(fsa_i,indFR0to1,indFRBefRun,numShuffle);
        sigShuf = prctile(ratioAftBefShuf,[p 100-p]);
        
        rec_name = [rec_name ; string(filename)];
        neu_id = [neu_id ; i];
        ratio0to1BefRun = [ratio0to1BefRun ; ratio0to1BefRun_i];
        isRise = [isRise ; ratio0to1BefRun_i >= sigShuf(1) & ~isinf(ratio0to1BefRun_i)];
        isDown = [isDown ; ratio0to1BefRun_i <= sigShuf(2) & ~isinf(ratio0to1BefRun_i)];
    end
    
    RiseDownTable = table(rec_name, neu_id, ratio0to1BefRun, isRise, isDown);
    
    save_path = 'Z:\Raphael_tests\Code\matlabAnalysisRaphi\RiseDownOffTargetRunBouts\RiseDownIDRunBout\';
    save_path = [save_path filename '_RiseDownID_RunBout.mat'];
    save(save_path, 'RiseDownTable');
end