function plotdFF(path, fileName, cellNo, trialNo)
% plot df/f for selected cells
% e.g. plotdFF('./','A577-20210817-05_DataStructure_mazeSection1_TrialType1',[1,4],58)
    
    GlobalConst2P;

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
    
    dFFTmp = trials{trialNo}.dFF(:,cellNo);
    maxDiffDFF = ceil(max(dFFTmp(:)) - min(dFFTmp(:)));
    
    figure
    hold on
    strNeu = [];
    strTitle = [];
    for i = 1:length(cellNo)
        plot((1:size(dFFTmp,1))/sampleFq,dFFTmp(:,i)'+(i-1)*(maxDiffDFF+1),'k-');
        strNeu = [strNeu '_' num2str(cellNo(i))];
        strTitle = [strTitle '-' num2str(cellNo(i))];
    end
    set(gca,'FontSize',12);
    xlabel('Time (s)')
    ylabel('dF/F')
    
    ind = strfind(fileName,'_');
    title([fileName(1:ind(1)-1) strTitle])
        
    fullpath = [path fileName(1:ind(1)-1) '_dFF_Tr' num2str(trialNo) '_Neu' strNeu];
    print('-painters', '-dpdf', fullpath, '-r600')
    savefig([fullpath '.fig']);
end
