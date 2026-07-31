function ProcessingActivityCorrRec()

    RecordingList;  
    spacebin = 20; %mm
    timebin = 0.0096; %s
    pathAnal = 'Z:\Raphael_tests\mice_expdata\Analysis\';
    onlyRun = 1;
    
    % passive
    recListPath = listRecordingsPassiveLickPath;
    recListName = listRecordingsPassiveLickFileName;
    recListSess = mazeSessionPassiveLick;
    recListLen = size(recListPath,1);
    recType = 1;
    
%     for i = 1:recListLen
%         activityCorrNeurons(recListPath(i,:),recListName(i,:),spacebin,1,recListSess(i));
%         
%         activityCorrTNeurons(recListPath(i,:),recListName(i,:),timebin,1,recListSess(i));
%         activitySimTNeurons(recListPath(i,:),recListName(i,:),timebin,1,recListSess(i));
%     end
    
    probDistrActivityCorr(pathAnal,recType,onlyRun);
    probDistrPopActivityCorr(pathAnal,recType,spacebin,onlyRun)
    
    % active
    recListPath = listRecordingsActiveLickPath;
    recListName = listRecordingsActiveLickFileName;
    recListSess = mazeSessionActiveLick;
    recListLen = size(recListPath,1);
    recType = 2;    
    
%     for i = 1:recListLen
%         activityCorrNeurons(recListPath(i,:),recListName(i,:),spacebin,1,recListSess(i));
%         activityCorrTNeurons(recListPath(i,:),recListName(i,:),timebin,1,recListSess(i));
%         activitySimTNeurons(recListPath(i,:),recListName(i,:),timebin,1,recListSess(i));
%     end
    
    probDistrActivityCorr(pathAnal,recType,onlyRun);
    probDistrPopActivityCorr(pathAnal,recType,spacebin,onlyRun)

    % passive no cue
    recListPath = listRecordingsNoCuePath;
    recListName = listRecordingsNoCueFileName;
    recListSess = mazeSessionNoCue;
    recListLen = size(recListPath,1);
    recType = 3;   
    
%     for i = 1:recListLen
%         activityCorrNeurons(recListPath(i,:),recListName(i,:),spacebin,1,recListSess(i));
%         activityCorrTNeurons(recListPath(i,:),recListName(i,:),timebin,1,recListSess(i));
%         activitySimTNeurons(recListPath(i,:),recListName(i,:),timebin,1,recListSess(i));
%     end
    
    probDistrActivityCorr(pathAnal,recType,onlyRun);
    probDistrPopActivityCorr(pathAnal,recType,spacebin,onlyRun)
    
    meanActivityCorrDiffBeh(pathAnal,onlyRun,3,2);
    meanPopActivityCorrDiffBeh(pathAnal,onlyRun,3,2);