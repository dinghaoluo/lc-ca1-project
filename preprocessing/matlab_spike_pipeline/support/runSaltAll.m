function runSaltAll_LC_VTA(path, fileName,baseline_time, time_resolution, window_size)
    load([fileName '.mat'], 'cluList');
    num_neurons = length(cluList.all);
    
    salt_p_vals = [];
    salt_I_vals = [];
    tagged = [];
    
    % Salt parameters(ms)
    % baseline_time = 500;
    % time_resolution = 1;
    % window_size = 5; % was 5
    
    for i = 1:num_neurons
        disp("on neuron");
        disp(i);               
        [p I percentActive] = salt_wrapper(path, fileName, i, baseline_time, time_resolution, window_size);
        salt_p_vals = [salt_p_vals p];
        salt_I_vals = [salt_I_vals I];
        if p < 0.01/length(num_neurons) && percentActive > 0.5
            tagged = [tagged 1];
        else
            tagged = [tagged 0];
        end
    end
    
    stimEffPath = [path '\' fileName '_stimEff.mat'];
    load(stimEffPath, 'stimResp');
    load(stimEffPath, 'param');
    
    stimResp.salt_p = salt_p_vals;
    stimResp.salt_I = salt_I_vals;
    stimResp.tagged = tagged;
    
    save(stimEffPath, 'stimResp', 'param');
end