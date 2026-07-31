function plotStimRasterWrapper(path, filename, neuron_id)
    
    stimEffPath = [path '\' filename '_stimEff.mat'];
    if ~isfile(stimEffPath)
        disp('need to run stim_Effect first');
        return;
    end
    
    load(stimEffPath, 'stimResp');
    load([filename '.mat'], 'cluList');
    shankID = cluList.shank;
    localClu = cluList.localClu;
    
    GlobalConst;
    
    if neuron_id == -1 % if neuron_id is -1, make these plots for all neurons
        num_neurons = length(stimResp.indNeurons);
        for i = 1:num_neurons
            plotStimRasterWrapper(path,filename,i);
        end
        return;
    end
    
    pulseBefore = stimResp.spPerPulseBef;
    pulseAfter = stimResp.spPerPulseAft;
    
    pBefNeuron = pulseBefore(:,neuron_id);
    pAftNeuron = pulseAfter(:,neuron_id);
    
    figure(1);
    clf;
    hold on;
    
    numStims = length(pBefNeuron);
    trial_count = 1;
    
    % some sessions have theta stimulation coded as tagging pulses -
    % need to make sure we only take in tagging pulses
    stimDur = stimResp.stimDuration;
    short_pulses = find(stimDur <= 5);
    first_tag_stim = short_pulses(1);
    
    for i = first_tag_stim:numStims
        stim_i_bef = pBefNeuron{i};
        stim_i_aft = pAftNeuron{i};
        pulses_stim_i = length(stim_i_bef);
        for j = 1:pulses_stim_i
            h = plot(stim_i_bef{j}/sampleFq*1000, trial_count * ones(length(stim_i_bef{j})),'k.');
            set(h, 'MarkerSize', 10);
            h = plot(stim_i_aft{j}/sampleFq*1000, trial_count * ones(length(stim_i_aft{j})),'k.');
            set(h, 'MarkerSize', 10);
            trial_count = trial_count+1;
        end     
    end
    PulsePeriod = stimResp.stimDuration(1);
    plot([0 0], [0 trial_count+1], 'c-');
    plot([PulsePeriod PulsePeriod],[0 trial_count+1], 'c-');
    
    % set Xlimit to 20ms before/after pulse for now
    xlim = [-20 20];
    set(gca, 'Xlim', xlim, 'Ylim', [0 trial_count+1]);
    
    xlabel('Time (ms)');
    ylabel('Pulse no.');
    
    p_vals = stimResp.salt_p;
    
    % single shank
    title(['cluster ' num2str(localClu(neuron_id)) ' (p = ' num2str(p_vals(neuron_id)) ')']);

%     % multi-shank
%     title(['neu ' num2str(neuron_id) '(' num2str(shankID(neuron_id)) ' ' num2str(localClu(neuron_id)) ')' 'p = ' num2str(p_vals(neuron_id))]);
    
    RasterDir = [path '\' filename '_Tagging_Raster_Plots\'];
    if ~exist(RasterDir,'dir')
        mkdir(RasterDir);
    end
    
    plot_file_name = [RasterDir 'neuron' num2str(localClu(neuron_id))  '.png'];
    saveas(gcf, plot_file_name);
end
