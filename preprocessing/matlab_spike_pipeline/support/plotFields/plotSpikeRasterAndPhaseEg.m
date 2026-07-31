function plotSpikeRasterAndPhaseEg()

    onlyRun = 1;
%     neuNo = 143;
%     path = 'Z:\Raphael_tests\mice_expdata\ANM011\A011-20190219\A011-20190219-01\';
%     fileName = 'A011-20190219-01_DataStructure_mazeSection1_TrialType1';
    neuNo = 26;
    path = 'Z:\Raphael_tests\mice_expdata\ANM012\A012-20190224\A012-20190224-01\';
    fileName = 'A012-20190224-01_DataStructure_mazeSection1_TrialType1';
    mazeSess = 1;
    trInd = [20:40, 60:90];
    
    fullPath = [path fileName...
                '_PeakFRAligned_msess' num2str(mazeSess)...
                '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood');
    if(length(trialNoNonStimGood) >= 50)
        trialNo = trialNoNonStimGood(trInd);
    else
        trialNo = trialNoNonStimGood;
    end
    
    plotSpikeRasterThetaPhaseEg(path, fileName, onlyRun, ...
        trialNo, neuNo);
    
    plotSpikeRaster(path, fileName, 3, trialNo, neuNo);  
    ind = strfind(fileName,'_');
    fileName1 = ['Z:\Yingxue\DataAnalysisRaphi\' ...
        fileName(1:ind(1)) 'Neu' num2str(neuNo) '_SpikeRaster'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end
