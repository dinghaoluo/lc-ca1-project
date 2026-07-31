function plotTaggedNeurons(path, filename)
    close all;

    stimEffPath = [path '\' filename '_stimEff.mat'];
        if ~isfile(stimEffPath)
            disp('need to run stim_Effect first');
            return;
        end

    load(stimEffPath, 'stimResp');
    tagged_indices = find(stimResp.tagged == 1);

    load([filename '.mat'], 'cluList');
    shankID = cluList.shank;
    localClu = cluList.localClu;

    load([filename '_Depth.mat']);
    depth = depthNeu.relDepthNeu;
    
    for i = 1:numel(tagged_indices)
        neuron_id = tagged_indices(i);
        shank = shankID(neuron_id);
        clu = localClu(neuron_id);

        plotStimRasterWrapper(path, filename, neuron_id);
%         plotCCG(path,filename(1:16), neuron_id, neuron_id);
%         plotWaveShape(path, filename(1:16), shank, clu, neuron_id);

        figure(4);
        ax = zeros(3,1);
        for j = 1:3
            ax(j) = subplot(3,1,j);
        end

%         for k = 1:3
%             figure(k);
%             h = get(gcf, 'Children');
%             newh = copyobj(h,4);
%             possub = get(ax(k), 'Position');
%             set(newh(1), 'Position', possub);
%             delete(ax(k));
%         end
        figure(4);
        savedir = [path '\' filename '_Tagging_Raster_Plots\Tagged_Clusters\'];
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
        text = (['Depth: ' num2str(depth(neuron_id))]);
        annotation('textbox', [.9 .4 .1 .2], 'String', text, 'EdgeColor', 'none');
        plot_file_name = [savedir 'neuron_tagged_' num2str(neuron_id) '.png'];
        saveas(gcf, plot_file_name);
        % also save to summary plots folder
        savedir2 = ['Z:\Dinghao\LC_recording_summ\tagged_cells\' filename(1:16) '_neuron_' num2str(neuron_id) '.png'];
        saveas(gcf, savedir2);
        close all;
    end

end