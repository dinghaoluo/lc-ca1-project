RecordingList;

onlyRun = 1;

numRec = size(listRecordingsActiveLickPath,1);
for i = numRec-1:numRec
    disp(listRecordingsActiveLickPath(i,:))
    cd(listRecordingsActiveLickPath(i,:));
    ProcessingMice_smTr(listRecordingsActiveLickPath(i,:),listRecordingsActiveLickFileName(i,:),1);
    ProcessingAligned(listRecordingsActiveLickPath(i,:),listRecordingsActiveLickFileName(i,:),1,mazeSessionActiveLick(i),4);
    close all;
    FieldWidthLRAligned(listRecordingsActiveLickPath(i,:),listRecordingsActiveLickFileName(i,:),1,2,1,mazeSessionActiveLick(i)); 
    ProcessingAlignedInclNonRun(listRecordingsActiveLickPath(i,:),listRecordingsActiveLickFileName(i,:),mazeSessionActiveLick(i));
%     close all;
end
