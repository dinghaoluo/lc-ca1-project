function plotImagescRasterDist(path,fileName,onlyRun,spaceBin,neuronNo,trialNo)
% plot the raster of convoluted spikes over trials, and plot the averaged
% firing on the top
% e.g. plotImagescRaster('./','A577-20210817-05_DataStructure_mazeSection1_TrialType1',1,20,4:21,[6:100])

    GlobalConst2P;

    fileName = [fileName '_convSpikesDist' num2str(spaceBin) 'mm_Run' num2str(onlyRun) '.mat']; 
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath, 'filteredSpikeArrayNormTNormAmp');
    
    numTr = length(trialNo);
    
    for n = neuronNo
        spikeArray = zeros(numTr,size(filteredSpikeArrayNormTNormAmp{end},2));
        for i = trialNo
            spikeArray(i,:) = filteredSpikeArrayNormTNormAmp{i}(n,:);
        end

        figure
        imagesc(spikeArray);
        colormap hot

        ylabel('Trial no.')
        xlabel('Distance (mm)')

        yyaxis right
        h = plot(mean(spikeArray),'m-');
        set(h,'LineWidth',0.5);
        ylabel('FR (Hz)');        
        set(gca, 'XLim', [-1 size(spikeArray,2)]);       
        figTitle = ['Neu ' num2str(n)];
        title(figTitle);
    end
end
