function verify_ConvSpikeTrain_AlignedRunBout(file_name)
    load_path = 'Z:\Raphael_tests\Code\matlabAnalysisRaphi\RiseDownOffTargetRunBouts\fsa_run_bouts\';
    load_path = [load_path file_name(1:16) '_BefRunBout0.mat'];
    load(load_path, 'filteredSpikeArrayRunBoutOnSet','paramC','timeStepRun'); 
    
    for i = 1:length(filteredSpikeArrayRunBoutOnSet)
        fsa_i = filteredSpikeArrayRunBoutOnSet{i};
        mean_fsa_i = mean(fsa_i,1);
        CreateFig();
        plot(timeStepRun, mean_fsa_i);
        ylabel('mean Fr');
        xlabel('time aligned to run bout');
        title_string = string(file_name(1:16)) + "\nneu " + num2str(i);
        title(char(compose(title_string)));
        save_path = 'Z:\Raphael_tests\Code\matlabAnalysisRaphi\RiseDownOffTargetRunBouts\verify_fsa_run_bouts\';
        save_path = [save_path file_name(1:16) '\'];
        if ~exist(save_path, 'dir')
            mkdir(save_path);
        end
        save_path = [save_path file_name(1:16) '_neu_' num2str(i) '.png'];
        saveas(gcf, save_path);
        close all;
        
    end
end