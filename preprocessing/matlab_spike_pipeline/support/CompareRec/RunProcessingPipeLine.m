RecordingList;

onlyRun = 1;
numRec = size(listRecordingsPath,1);
for i = 1:numRec
    disp(listRecordingsPath(i,:))
    cd(listRecordingsPath(i,:));
    ProcessingMice_smTrSpeedLick(listRecordingsPath(i,:),listRecordingsFileName(i,:),onlyRun,mazeSession(i));
    ProcessingMice_smTr(listRecordingsPath(i,:),listRecordingsFileName(i,:),onlyRun);
    disp(listRecordingsPath(i,:))
    ProcessingMice_smTr_GoodTr(listRecordingsPath(i,:),listRecordingsFileName(i,:),onlyRun,mazeSession(i));
    disp(listRecordingsPath(i,:))
    ProcessingMice_smTrCtrlOnly(listRecordingsPath(i,:),listRecordingsFileName(i,:),onlyRun,mazeSession(i));
    FieldWidthLR_GoodTr(listRecordingsPath(i,:),listRecordingsFileName(i,:),20,1,2,0,...
        onlyRun,mazeSession(i));
    pause;
    close all;
    disp(listRecordingsPath(i,:))
    ProcessingAligned(listRecordingsPath(i,:),listRecordingsFileName(i,:),onlyRun,mazeSession(i),1);
    disp(listRecordingsPath(i,:))
    ProcessingAlignedRun0(listRecordingsPath(i,:),listRecordingsFileName(i,:),mazeSession(i),1);
    close all;
    disp(listRecordingsPath(i,:))
    FieldWidthLRAligned(listRecordingsPath(i,:),listRecordingsFileName(i,:),1,0,onlyRun,mazeSession(i));
    ProcessingAlignedCtrlOnly(listRecordingsPath(i,:),listRecordingsFileName(i,:),onlyRun,mazeSession(i));
end