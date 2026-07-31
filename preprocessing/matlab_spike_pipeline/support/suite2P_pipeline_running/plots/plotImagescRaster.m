function plotImagescRaster(path,fileName,mazeSess,onlyRun,neuronNo,trialNo)
% plot the raster of convoluted spikes over trials, and plot the averaged
% firing on the top
% e.g.: plotImagescRaster('./','A577-20210817-05_DataStructure_mazeSection1_TrialType1',1,1,[3 4 8],[6:105])

    GlobalConst2P;

    fileName = [fileName '_convSpikesAligned_msess' num2str(mazeSess) '_BefRun' num2str(onlyRun) '.mat']; 
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath, 'filteredSpikeArrayRunOnset');
    
    numTr = length(trialNo);
    
    for n = neuronNo
        spikeArray = filteredSpikeArrayRunOnset{n}(trialNo,:);
        for i = 1:length(trialNo)
            if(sum(spikeArray(i,:)) > 0)
                spikeArray(i,:) = spikeArray(i,:)/max(spikeArray(i,:));
            end
        end
        
        figure
        imagesc((-nSampBef:size(spikeArray,2)-nSampBef-1)/sampleFq,...
            1:numTr,spikeArray);
        colormap pink
        colorbar

        ylabel('Trial no.')
        xlabel('Time (s)')

        yyaxis right
        h = plot((-nSampBef:size(spikeArray,2)-nSampBef-1)/sampleFq,...
            mean(spikeArray));
        set(h,'LineWidth',0.5,'Color',[1 1 1]);
        ylabel('deconvoluted firing rate');        
        set(gca, 'XLim', [-1 5],'FontSize',12);       
        figTitle = ['Neu ' num2str(n)];
        title(figTitle);
        
        ind = strfind(fileName,'_');
        fullpath = [path fileName(1:ind(1)-1) '_msess' num2str(mazeSess)...
            '_Run' num2str(onlyRun) '_DeConvRaster_Neu' num2str(n)];
        print('-painters', '-dpdf', fullpath, '-r600')
        savefig([fullpath '.fig']);
    end
end
