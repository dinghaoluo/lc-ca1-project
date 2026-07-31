function pMeanLickAlignedOverDistAllRec()

    RecordingList;

%     path = listRecordingsActiveLickPath(1,:);
%     fileName = listRecordingsActiveLickFileName(1,:);
% 
%     fullPath = [path,fileName, '_lickDist.mat'];
%     load(fullPath,'param');

    fullPath = ['Z:\Yingxue\DataAnalysisRaphi\ALLickAlignedDistMeanAllRec.mat'];
    load(fullPath,'lickAlignedOverDistMeanAllArr','param');
    
    indBefRew = param.spaceStepsAligned < 0;
    meanPredLickAL = sum(lickAlignedOverDistMeanAllArr(:,indBefRew),2);
    mMeanPredLickAL = mean(meanPredLickAL);
    semMeanPredLickAL = std(meanPredLickAL)/sqrt(length(meanPredLickAL));
    
    fullPath = ['Z:\Yingxue\DataAnalysisRaphi\NCLickAlignedDistMeanAllRec.mat'];
    load(fullPath,'lickAlignedOverDistMeanAllArr','param');
    
    meanPredLickNC = sum(lickAlignedOverDistMeanAllArr(:,indBefRew),2);
    mMeanPredLickNC = mean(meanPredLickNC);
    semMeanPredLickNC = std(meanPredLickNC)/sqrt(length(meanPredLickNC));
    
    pMeanLickAligned = ranksum(meanPredLickAL,meanPredLickNC);
end