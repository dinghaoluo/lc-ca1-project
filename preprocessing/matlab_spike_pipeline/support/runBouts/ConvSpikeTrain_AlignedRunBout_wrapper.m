function ConvSpikeTrain_AlignedRunBout_wrapper
    RecordingList_dinghao;
    
    for i = 1:size(listRecordingsActiveLickPath,1)
        disp(listRecordingsActiveLickFileName(i,1:17));
        % rerun each recording; the old FSA files are not reused here
        ConvSpikeTrain_AlignedRunBout(listRecordingsActiveLickPath(i,:), listRecordingsActiveLickFileName(i,:));
    end
end
