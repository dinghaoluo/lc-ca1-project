function plotImagescRasterDistAlignedRun(path,fileName,onlyRun,mazeSess,neuronNo,trialNo)
% plot the raster of convoluted spikes over trials, and plot the averaged
% firing on the top
% e.g. plotImagescRasterDistAlignedRun('./','A577-20210817-05_DataStructure_mazeSection1_TrialType1',1,1,[3 4 8],[6:105])

    GlobalConst2P;

    fileName = [fileName '_convSpikesDistAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath,'filteredSpikeDistNormTArrayRun');
    
    numTr = length(trialNo);
    
    for n = neuronNo
        spikeArray = filteredSpikeDistNormTArrayRun{n}(trialNo,:);
        for i = 1:length(trialNo)
            if(sum(spikeArray(i,:)) > 0)
                spikeArray(i,:) = spikeArray(i,:)/max(spikeArray(i,:));
            end
        end

        figure
        imagesc(0:0.1:size(spikeArray,2)/10-0.1,1:length(trialNo),spikeArray);
        colormap pink
        colorbar

        ylabel('Trial no.')
        xlabel('Distance (mm)')

        yyaxis right
        h = plot(0:0.1:size(spikeArray,2)/10-0.1,mean(spikeArray),'m-');
        set(h,'LineWidth',0.5,'Color',[1 1 1]);
        ylabel('deconvoluted firing rate');        
        set(gca, 'XLim', [0 180],'FontSize',12);       
        figTitle = ['Neu ' num2str(n)];
        title(figTitle);
        
        ind = strfind(fileName,'_');
        fullpath = [path fileName(1:ind(1)-1) '_msess' num2str(mazeSess)...
            '_Run' num2str(onlyRun) 'AlignedRun_DeConvRaster_Neu' num2str(n)];
        print('-painters', '-dpdf', fullpath, '-r600')
        savefig([fullpath '.fig']);
    end
end
