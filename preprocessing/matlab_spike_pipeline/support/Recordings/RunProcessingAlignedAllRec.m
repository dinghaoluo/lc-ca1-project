RecordingList;
% for i = 1:length(listRecordingsNoCuePath)
%     ProcessingMice_smTr(listRecordingsNoCuePath(i,:),listRecordingsNoCueFileName(i,:),1);    
% end
% 
% onlyRun = 1;
% for i = 1:size(listRecordingsActiveLickPath,1)
%     path = listRecordingsActiveLickPath(i,:);   
%     fileName = listRecordingsActiveLickFileName(i,:);
%     mazeSess = mazeSessionActiveLick(i);
%     disp(fileName);
% %     ProcessingAlignedAll(path,fileName,onlyRun,mazeSess); 
% 
% %     disp('Calculate running speed over distance')
% %     RunSpeedOverDist(path, fileName, onlyRun, mazeSess);
% 
%     disp('Calculate lick over distance')
%     LickOverDist(path, fileName, mazeSess);
% end
% 
% for i = 1:size(listRecordingsPassiveLickPath,1)
%     path = listRecordingsPassiveLickPath(i,:);   
%     fileName = listRecordingsPassiveLickFileName(i,:);
%     mazeSess = mazeSessionPassiveLick(i);
%     disp(fileName);
% %     ProcessingAlignedAll(path,fileName,onlyRun,mazeSess); 
% 
% %     disp('Calculate running speed over distance')
% %     RunSpeedOverDist(path, fileName, onlyRun, mazeSess);    
% 
%     disp('Calculate lick over distance')
%     LickOverDist(path, fileName, mazeSess);
% end

for i = 1:size(listRecordingsNoCuePath,1)
    path = listRecordingsNoCuePath(i,:);   
    fileName = listRecordingsNoCueFileName(i,:);
    mazeSess = mazeSessionNoCue(i);
    disp(fileName);
    
%     disp('Calculate running speed over distance')
%     RunSpeedOverDist(path, fileName, onlyRun, mazeSess);    

    disp('Calculate lick over distance')
    LickOverDistNoCue(path, fileName, mazeSess);
end

% spikeTrainSimilarityAllRec(1)
% 
% popSimilarityAllRec(1)

% meanLickOverDistAllRec

meanLickAlignedOverDistAllRec

meanRunSpeedOverDistAllRec(1)
