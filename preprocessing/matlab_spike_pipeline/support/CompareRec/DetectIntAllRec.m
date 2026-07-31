RecordingList;

numRec = size(listRecordingsNoCuePath,1);
for i = 1:numRec
    disp(listRecordingsNoCuePath(i,:))
    DetectInt(listRecordingsNoCuePath(i,:),listRecordingsNoCueFileName(i,:),1);
end

numRec = size(listRecordingsActiveLickPath,1);
for i = 1:numRec
    disp(listRecordingsActiveLickPath(i,:))
    DetectInt(listRecordingsActiveLickPath(i,:),listRecordingsActiveLickFileName(i,:),1);
end

numRec = size(listRecordingsPassiveLickPath,1);
for i = 1:numRec
    disp(listRecordingsPassiveLickPath(i,:))
    DetectInt(listRecordingsPassiveLickPath(i,:),listRecordingsPassiveLickFileName(i,:),1);
end