function plotdFFFields(path, fileName, mazeSess, onlyRun)
% plot df/f for selected cells
% e.g. plotdFFFields('./','A576-20210914-02_DataStructure_mazeSection1_TrialType1',10,1,1)
    
    GlobalConst2P;

    fileNameFW = [fileName '_FieldSpCorrAligned_Run' num2str(mazeSess) ...
                            '_Run' num2str(onlyRun) '.mat'];
    fullPath = [path fileNameFW];
    if(exist(fullPath) == 0)
        disp('The _FieldSpCorrAligned_Run file does not exist');
        return;
    end
    load(fullPath,'fieldSpCorrSessNonStimGood');
    cellNo = fieldSpCorrSessNonStimGood.indNeuron;
    filenameEnd = '_Field';
    
    indexFileName = findstr(fileName, '.mat');
    if(isempty(indexFileName))
        fileName = [fileName '.mat'];
    end 
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath, 'trials');
    ind = strfind(fileName,'_');
     
    figure
    for trialNo = 1:length(trials)
        dFFTmp = trials{trialNo}.dFF(:,cellNo);
        [dFFMax,maxInd] = max(dFFTmp);
        maxDiffDFF = ceil(max(dFFTmp(:)) - min(dFFTmp(:)));

        strTitle = [];
        for i = 1:length(cellNo)
            if(i == 2)
                hold on
            end
            plot((1:size(dFFTmp,1))/sampleFq,dFFTmp(:,i)'+(i-1)*(maxDiffDFF+1),'k-');
            strTitle = ['Tr' num2str(trialNo)];
%             strNeu = [strNeu '_' num2str(cellNo(i))];
%             strTitle = [strTitle '-' num2str(cellNo(i))];
        end
        set(gca,'FontSize',12);
        xlabel('Time (s)')
        ylabel('dF/F')
        title([fileName(1:ind(1)-1) ' ' strTitle])
        
        hold off;
        pause;
        
        fullpath = [path fileName(1:ind(1)-1) '_dFF' filenameEnd '_Tr' num2str(trialNo)];
        print('-painters', '-dpdf', fullpath, '-r600')
        savefig([fullpath '.fig']);
    end
        
    trialNo = 33;
    cellNoInd = [2 3 4 5 7 11 12 17 21];
    
    dFFTmp = trials{trialNo}.dFF(:,cellNo);
    [dFFMax,maxInd] = max(dFFTmp);
    [a,indSorted] = sort(maxInd(cellNoInd));
    dFFTmp1 = dFFTmp(:,cellNoInd(indSorted));
    maxDiffDFF1 = min(2,ceil(max(dFFTmp1(:)) - min(dFFTmp1(:))));

    figure;
    strTitle = [];
    for i = 1:length(cellNoInd)
        if(i == 2)
            hold on
        end
        plot((1:size(dFFTmp1,1))/sampleFq,dFFTmp1(:,i)'+(i-1)*(maxDiffDFF1+1),'k-');
        strTitle = ['Tr' num2str(trialNo)];
    end
    set(gca,'FontSize',12);
    xlabel('Time (s)')
    ylabel('dF/F')

    title([fileName(1:ind(1)-1) ' ' strTitle])
    hold off;
    pause;

    fullpath = [path fileName(1:ind(1)-1) '_dFF' filenameEnd '_Tr' num2str(trialNo)];
    print('-painters', '-dpdf', fullpath, '-r600')
    savefig([fullpath '.fig']);
    
    trialNo = 70;
    cellNoInd = [2 5 7 8 17 18 19];
    
    dFFTmp = trials{trialNo}.dFF(:,cellNo);
    [dFFMax,maxInd] = max(dFFTmp);
    [a,indSorted] = sort(maxInd(cellNoInd));
    dFFTmp1 = dFFTmp(:,cellNoInd(indSorted));
    maxDiffDFF1 = min(2,ceil(max(dFFTmp1(:)) - min(dFFTmp1(:))));

    figure;
    strTitle = [];
    for i = 1:length(cellNoInd)
        if(i == 2)
            hold on
        end
        plot((1:size(dFFTmp1,1))/sampleFq,dFFTmp1(:,i)'+(i-1)*(maxDiffDFF1+1),'k-');
        strTitle = ['Tr' num2str(trialNo)];
    end
    set(gca,'FontSize',12);
    xlabel('Time (s)')
    ylabel('dF/F')

    title([fileName(1:ind(1)-1) ' ' strTitle])
    hold off;
    pause;

    fullpath = [path fileName(1:ind(1)-1) '_dFF' filenameEnd '_Tr' num2str(trialNo)];
    print('-painters', '-dpdf', fullpath, '-r600')
    savefig([fullpath '.fig']);
    
end
