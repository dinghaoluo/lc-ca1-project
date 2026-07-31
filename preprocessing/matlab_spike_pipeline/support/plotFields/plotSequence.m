function plotSequence(pathAnal,spaceBin,onlyRun)
% plot neuronal sequences
%e.g. plotSequence('Z:\Raphael_tests\mice_expdata\Analysis\',20,1)

    RecordingList_dinghao;
    recListPath = listRecordingsActiveLickPathHPCLCOpt;
    recListFileName = listRecordingsActiveLickFileNameHPCLCOpt;
%     recListSess = mazeSessionActiveLickHPCLCOpt;  % unused for now
    
    recListLen = size(recListPath,1);
    
    load([pathAnal 'popActivityCorr_RecList' num2str(recType) '_Run1.mat']);
    
    for i = 1:recListLen
        path = recListPath(i,:);
        fileName = recListFileName(i,:);
        fullPath = [path fileName ...
            '_convSpikesDist' num2str(spaceBin) ...
                    'mm_Run' num2str(onlyRun) '.mat'];
        if(exist(fullPath,'file') == 0)
            disp(['The firing profile file does not exist. Try to run the',...
                    'function again with fileState = 0.']);
            return;
        end
        load(fullPath,'filteredSpikeArrayNormT');
                        
        indNeu = recListPopActivityCorr.indPyr{i};
        indLaps = recListPopActivityCorr.indLaps{i};
        
        filteredSpikeArray = zeros(length(indNeu),...
            size(filteredSpikeArrayNormT{2},2));
        filteredSpikeArrayNorm = zeros(length(indNeu),...
            size(filteredSpikeArrayNormT{2},2));
        indPeakArr = zeros(1,length(indNeu));
        
        if(length(indNeu) < 3)
            continue;
        end
            
        for n = 1:length(indNeu)
            filteredSpikeArrayTmp = zeros(1,size(filteredSpikeArrayNormT{2},2));
            for j = 1:length(indLaps)
                if(isempty(filteredSpikeArrayNormT{indLaps(j)}))
                    continue;
                end
                filteredSpikeArrayTmp = filteredSpikeArrayTmp + ...
                    filteredSpikeArrayNormT{indLaps(j)}(indNeu(n),:);
            end
            filteredSpikeArrayTmp = filteredSpikeArrayTmp/length(indLaps);
            [peakValue,indPeak] = max(filteredSpikeArrayTmp);
            indPeakArr(n) = indPeak(1);
            
            filteredSpikeArray(n,:) = filteredSpikeArrayTmp;
            filteredSpikeArrayNorm(n,:) = filteredSpikeArrayTmp/peakValue;
        end  
        
        [~,indSortPeak] = sort(indPeakArr);

        filteredSpikeArrayOrdered = filteredSpikeArray(indSortPeak,:);
        filteredSpikeArrayNormOrdered = filteredSpikeArrayNorm(indSortPeak,:);

%         figure
%         imagesc(1:length(filteredSpikeArrayTmp),1:length(indNeu),...
%             filteredSpikeArrayOrdered);
%         set(gca,'FontSize',14);
%         xlabel('Dist (mm)');
%         ylabel('Neuron no.')
%         title(fileName(1:13));

        fig = figure;
        set(0,'Units','pixels') 
        set(figure(fig),'OuterPosition',...
            [500 500 280 280])
        imagesc(1:length(filteredSpikeArrayTmp),1:length(indNeu),...
            filteredSpikeArrayNormOrdered);
        colormap('jet');
        colorbar;
        set(gca,'FontSize',14,'YDir','normal');
        xlabel('Dist (mm)');
        ylabel('Neuron no.')
        title(fileName(1:13));
        print('-painters','-dpdf',[pathAnal 'Sequence_' fileName(1:13) '_RecType' ...
            num2str(recType) '_Run' num2str(onlyRun)],'-r600');
    end
end
