RecordingList;
% for i = 2:size(listRecordingsNoCuePath,1)
%     disp(listRecordingsNoCueFileName(i,:));
%     ProcessingMice_smTr(listRecordingsNoCuePath(i,:),listRecordingsNoCueFileName(i,:),1);    
% end

for i = 1:size(listRecordingsPassiveLickPath,1)
    disp(listRecordingsPassiveLickPath(i,:));
    ProcessingMice_smTr(listRecordingsPassiveLickPath(i,:),listRecordingsPassiveLickFileName(i,:),1);    
end

for i = 1:size(listRecordingsActiveLickPath,1)
    disp(listRecordingsActiveLickFileName(i,:));
    ProcessingMice_smTr(listRecordingsActiveLickPath(i,:),listRecordingsActiveLickFileName(i,:),1);    
end