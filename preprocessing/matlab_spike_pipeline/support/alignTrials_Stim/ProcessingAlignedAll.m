function ProcessingAlignedAll(path,fileName,onlyRun,mazeSess,cond)
% no cue passive --- condition 1
% 60 cm cue passive --- condition 2
% 0.5s start cue passive --- condition 3
% 0.5s start cue active --- condition 4

    ProcessingAligned(path,fileName,onlyRun,mazeSess,cond)
    
%     if(cond == 2)
%         ProcessingAligned_CueOff(path,fileName,onlyRun,mazeSess,cond)
%     end
    
    ProcessingAligned_Corr(path,fileName,onlyRun,mazeSess)
    
end