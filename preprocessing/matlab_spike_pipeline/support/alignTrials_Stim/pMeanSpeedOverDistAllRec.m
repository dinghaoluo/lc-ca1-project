function pMeanSpeedOverDistAllRec()

    RecordingList;
    onlyRun = 1;

%     path = listRecordingsActiveLickPath(1,:);
%     fileName = listRecordingsActiveLickFileName(1,:);
% 
%     fullPath = [path,fileName, '_lickDist.mat'];
%     load(fullPath,'param');

    fullPath = ['Z:\Yingxue\DataAnalysisRaphi\ALRunSpeedDistMeanAllRec_Run' num2str(onlyRun) '.mat'];
    load(fullPath,'speedOverDistMeanAllArr');
    
    meanSpeedAL = mean(speedOverDistMeanAllArr,2);
    mMeanSpeedAL = mean(meanSpeedAL);
    semMeanSpeedAL = std(meanSpeedAL)/sqrt(length(meanSpeedAL));
    
    fullPath = ['Z:\Yingxue\DataAnalysisRaphi\NCRunSpeedDistMeanAllRec_Run' num2str(onlyRun) '.mat'];
    load(fullPath,'speedOverDistMeanAllArr');
    
    meanSpeedNC = mean(speedOverDistMeanAllArr,2);
    mMeanSpeedNC = mean(meanSpeedNC);
    semMeanSpeedNC = std(meanSpeedNC)/sqrt(length(meanSpeedNC));
    
    pMeanSpeed = ranksum(meanSpeedAL,meanSpeedNC);
end