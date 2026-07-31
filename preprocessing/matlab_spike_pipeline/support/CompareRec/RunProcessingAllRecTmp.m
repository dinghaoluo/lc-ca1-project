RecordingList;

onlyRun = 1;
% numRec = size(listRecordingsNoCuePath,1);
% for i = 7:numRec
%     disp(listRecordingsNoCuePath(i,:))
%     cd(listRecordingsNoCuePath(i,:));
%     ProcessingMice_smTr(listRecordingsNoCuePath(i,:),listRecordingsNoCueFileName(i,:),1);
% %     ProcessingMice_smTr_GoodTr(listRecordingsNoCuePath(i,:),listRecordingsNoCueFileName(i,:),1,mazeSessionNoCue(i));
% % %     FieldWidthLR_GoodTr(listRecordingsNoCuePath(i,:),listRecordingsNoCueFileName(i,:),20,1,2,0,...
% % %         onlyRun,mazeSessionNoCue(i));
% % %     pause;
% % %     close all;
%     ProcessingAligned(listRecordingsNoCuePath(i,:),listRecordingsNoCueFileName(i,:),1,mazeSessionNoCue(i),1);
%     ProcessingAlignedRun0(listRecordingsNoCuePath(i,:),listRecordingsNoCueFileName(i,:),mazeSessionNoCue(i),1);
% %     close all;
% % %     ProcessingAligned_CueOff(listRecordingsNoCuePath(i,:),listRecordingsNoCueFileName(i,:),1,mazeSessionNoCue(i));
% % %     close all;
% % %     FieldWidthLRAligned(listRecordingsNoCuePath(i,:),listRecordingsNoCueFileName(i,:),1,0,1,mazeSessionNoCue(i));
% % %     pause;
% % %     movefile([listRecordingsNoCueFileName(i,:) '_runSpeedDist_Run' num2str(onlyRun) '.mat'], ...
% % %         [listRecordingsNoCueFileName(i,:) '_runSpeedDist_msess' num2str(mazeSessionNoCue(i)) '_Run' num2str(onlyRun) '.mat']);
% % %     movefile([listRecordingsNoCueFileName(i,:) '_lickDist.mat'], ...
% % %         [listRecordingsNoCueFileName(i,:) '_lickDist_msess' num2str(mazeSessionNoCue(i)) '.mat']);
%     close all;
% %     backUpFileAligned(listRecordingsNoCuePath(i,:),listRecordingsNoCueFileName(i,:),mazeSessionNoCue(i),onlyRun);
%     changeNameRun(listRecordingsNoCuePath(i,:),listRecordingsNoCueFileName(i,:),mazeSessionNoCue(i),onlyRun);
% end

numRec = size(listRecordingsPassiveLickPath,1);
for i = 1:numRec
    disp(listRecordingsPassiveLickPath(i,:))
    cd(listRecordingsPassiveLickPath(i,:)); 
% %     ProcessingMice_smTr(listRecordingsPassiveLickPath(i,:),listRecordingsPassiveLickFileName(i,:),1);
% %     ProcessingMice_smTr_GoodTr(listRecordingsPassiveLickPath(i,:),listRecordingsPassiveLickFileName(i,:),1,mazeSessionPassiveLick(i));
% % %     FieldWidthLR_GoodTr(listRecordingsPassiveLickPath(i,:),listRecordingsPassiveLickFileName(i,:),20,1,2,0,...
% % %         onlyRun,mazeSessionPassiveLick(i));
% % %     pause;
% % %     close all;
%     ProcessingAligned(listRecordingsPassiveLickPath(i,:),listRecordingsPassiveLickFileName(i,:),1,mazeSessionPassiveLick(i),3);
%     ProcessingAlignedRun0(listRecordingsPassiveLickPath(i,:),listRecordingsPassiveLickFileName(i,:),mazeSessionPassiveLick(i),3);
% % % %     close all;
% % % %     ProcessingAligned_CueOff(listRecordingsPassiveLickPath(i,:),listRecordingsPassiveLickFileName(i,:),1,mazeSessionPassiveLick(i));
% % % %     close all;
% % %     FieldWidthLRAligned(listRecordingsPassiveLickPath(i,:),listRecordingsPassiveLickFileName(i,:),1,0,1,mazeSessionPassiveLick(i));
% % % %     pause;
% % % %     movefile([listRecordingsPassiveLickFileName(i,:) '_lickDist.mat'], ...
% % % %         [listRecordingsPassiveLickFileName(i,:) '_lickDist_msess' num2str(mazeSessionPassiveLick(i)) '.mat']);
% % % %     movefile([listRecordingsPassiveLickFileName(i,:) '_runSpeedDist_Run' num2str(onlyRun) '.mat'], ...
% % % %         [listRecordingsPassiveLickFileName(i,:) '_runSpeedDist_msess' num2str(mazeSessionPassiveLick(i)) '_Run' num2str(onlyRun) '.mat']);
%     close all;
% %     backUpFileAligned(listRecordingsPassiveLickPath(i,:),listRecordingsPassiveLickFileName(i,:),mazeSessionPassiveLick(i),onlyRun);
%     changeNameRun(listRecordingsPassiveLickPath(i,:),listRecordingsPassiveLickFileName(i,:),mazeSessionPassiveLick(i),onlyRun);
end

numRec = size(listRecordingsActiveLickPath,1);
for i = 1:numRec
%     disp(listRecordingsActiveLickPath(i,:))
    cd(listRecordingsActiveLickPath(i,:));
    load([listRecordingsActiveLickFileName(i,:) '_behPar_msess' num2str(mazeSessionActiveLick(i)) '.mat']);
    if(sum(behPar.rewarded == -1) > 0)
        disp([num2str(i) '  ' listRecordingsActiveLickFileName(i,:) ' ' num2str(sum(behPar.rewarded == -1))]);
    end
%     ProcessingMice_smTr(listRecordingsActiveLickPath(i,:),listRecordingsActiveLickFileName(i,:),1);
%     ProcessingMice_smTr_GoodTr(listRecordingsActiveLickPath(i,:),listRecordingsActiveLickFileName(i,:),1,mazeSessionActiveLick(i));
%     FieldWidthLR_GoodTr(listRecordingsActiveLickPath(i,:),listRecordingsActiveLickFileName(i,:),20,1,2,0,...
%         onlyRun,mazeSessionActiveLick(i));
%     pause;
%     close all;
%     ProcessingAligned(listRecordingsActiveLickPath(i,:),listRecordingsActiveLickFileName(i,:),1,mazeSessionActiveLick(i),4);
%     ProcessingAlignedRun0(listRecordingsActiveLickPath(i,:),listRecordingsActiveLickFileName(i,:),mazeSessionActiveLick(i),4);
%     close all;
% %     ProcessingAligned_CueOff(listRecordingsActiveLickPath(i,:),listRecordingsActiveLickFileName(i,:),1,mazeSessionActiveLick(i));
% %     close all;
%     FieldWidthLRAligned(listRecordingsActiveLickPath(i,:),listRecordingsActiveLickFileName(i,:),1,2,1,mazeSessionActiveLick(i)); 
% %     pause;
% %     movefile([listRecordingsActiveLickFileName(i,:) '_lickDist.mat'], ...
% %         [listRecordingsActiveLickFileName(i,:) '_lickDist_msess' num2str(mazeSessionActiveLick(i)) '.mat']);
% %     movefile([listRecordingsActiveLickFileName(i,:) '_runSpeedDist_Run' num2str(onlyRun) '.mat'], ...
% %         [listRecordingsActiveLickFileName(i,:) '_runSpeedDist_msess' num2str(mazeSessionActiveLick(i)) '_Run' num2str(onlyRun) '.mat']);
%     pause;
%     close all;
%     backUpFileAligned(listRecordingsActiveLickPath(i,:),listRecordingsActiveLickFileName(i,:),mazeSessionActiveLick(i),onlyRun);
%     changeNameRun(listRecordingsActiveLickPath(i,:),listRecordingsActiveLickFileName(i,:),mazeSessionActiveLick(i),onlyRun);
end

% numRec = size(listRecordingsNoCuePath,1);
% for i = 1:numRec %numRec
%     disp(listRecordingsNoCuePath(i,:))
%     cd(listRecordingsNoCuePath(i,:));
% 
%     fileNameFWOld = [listRecordingsNoCueFileName(i,:) '_FieldWidthAligned_'...
%                         'mSess*.mat'];
%     fileNameFW = [listRecordingsNoCueFileName(i,:) '_FieldWidthAligned_'...
%                         'mSess' num2str(mazeSessionNoCue(i)) '_L_Run1_Param.mat'];
%     listing = dir(fileNameFWOld);
%     if(length(listing)~= 0)
%     	movefile(listing.name,fileNameFW);
%     else
%         disp('param file does not exist'); 
%     end
% end

function backUpFileAligned(path, fileName, mazeSess, onlyRun)
    formatOut = 'mm/dd/yy';
    dateS = datestr(now,formatOut);
    dateS = [dateS(1:2) dateS(4:5) dateS(7:8)];
    backupFolderN = [path 'backup_AlignedRun' num2str(onlyRun) '_' dateS];
    if(~exist(backupFolderN,'dir'))
        mkdir(backupFolderN);
    end
    movefile([path fileName '_behPar_msess' num2str(mazeSess) '.mat'], ...
        [backupFolderN '\' fileName '_behPar_msess' num2str(mazeSess) '.mat']);
    
    movefile([path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_convSpikesAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_convSpikesAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_thetaPhaseOverTimeligned_msess' num2str(mazeSess) '.mat'], ...
        [backupFolderN '\' fileName '_thetaPhaseOverTimeligned_msess' num2str(mazeSess) '.mat']);
    
    movefile([path fileName '_lickDist_msess' num2str(mazeSess) '.mat'], ...
        [backupFolderN '\' fileName '_lickDist_msess' num2str(mazeSess) '.mat']);
    
    movefile([path fileName '_runSpeedDist_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_runSpeedDist_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_FRAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_FRAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_ThetaPhaseLAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_ThetaPhaseLAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_ThetaPhaseLAlignedSeg_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_ThetaPhaseLAlignedSeg_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_ThetaPhaseHAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_ThetaPhaseHAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_ThetaPhaseHAlignedSeg_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_ThetaPhaseHAlignedSeg_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_convSpikesAligned_msess' num2str(mazeSess) '_BefRun' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_convSpikesAligned_msess' num2str(mazeSess) '_BefRun' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_PeakFR_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_PeakFR_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_ThetaPhaseH_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_ThetaPhaseH_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_ThetaAlignRun_msess' num2str(mazeSess) '.mat'], ...
        [backupFolderN '\' fileName '_ThetaAlignRun_msess' num2str(mazeSess) '.mat']);
    
    movefile([path fileName '_ConcatspAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_ConcatspAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_CCGAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_CCGAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_ThetaModAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_ThetaModAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_burstAllAlignedRun_THL_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_burstAllAlignedRun_THL_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_burstAllAlignedRun_THH_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_burstAllAlignedRun_THH_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_timePerDistBin_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_timePerDistBin_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_convSpikesDistAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'], ...
        [backupFolderN '\' fileName '_convSpikesDistAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_spikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat'], ...
        [backupFolderN '\' fileName '_spikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat']);
    
    movefile([path fileName '_meanSpikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat'], ...
        [backupFolderN '\' fileName '_meanSpikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat']);
    
%     movefile([path fileName '_spikesCorrDistAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intD1800.mat'], ...
%         [backupFolderN '\' fileName '_spikesCorrDistAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intD1800.mat']);
%     
%     movefile([path fileName '_meanSpikesCorrDistAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intD1800.mat'], ...
%         [backupFolderN '\' fileName '_meanSpikesCorrDistAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intD1800.mat']);
    
    movefile([path fileName '_spikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat'], ...
        [backupFolderN '\' fileName '_spikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat']);
    
    movefile([path fileName '_meanSpikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat'], ...
        [backupFolderN '\' fileName '_meanSpikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat']);
    
    movefile([path fileName '_spikeTrainSimVP_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_q0p0032_intT20.mat'], ...
        [backupFolderN '\' fileName '_spikeTrainSimVP_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_q0p0032_intT20.mat']);
    
    movefile([path fileName '_meanSpikeTrainSimVP_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_q0p0032_intT20.mat'], ...
        [backupFolderN '\' fileName '_meanSpikeTrainSimVP_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_q0p0032_intT20.mat']);
    
    movefile([path fileName '_spikeTrainSimiVanR_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_tc500_intT20.mat'], ...
        [backupFolderN '\' fileName '_spikeTrainSimiVanR_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_tc500_intT20.mat']);
    
    movefile([path fileName '_mean_spikeTrainSimiVanR_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_tc500_intT20.mat'], ...
        [backupFolderN '\' fileName '_mean_spikeTrainSimiVanR_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_tc500_intT20.mat']);
    
    movefile([path fileName '_popCorrT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat'], ...
        [backupFolderN '\' fileName '_popCorrT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat']);
    
    movefile([path fileName '_meanPopCorrT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat'], ...
        [backupFolderN '\' fileName '_meanPopCorrT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat']);
    
    movefile([path fileName '_popSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat'], ...
        [backupFolderN '\' fileName '_popSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat']);
    
    movefile([path fileName '_meanPopSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat'], ...
        [backupFolderN '\' fileName '_meanPopSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT20.mat']);
end

function changeNameRun(path, fileName, mazeSess, onlyRun)
    movefile([path fileName '_ThetaAlignRun_msess' num2str(mazeSess) '.mat'], ...
        [path fileName '_ThetaAlignRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_lickDist_msess' num2str(mazeSess) '.mat'], ...
        [path fileName '_lickDist_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    movefile([path fileName '_thetaPhaseOverTimeligned_msess' num2str(mazeSess) '.mat'], ...
        [path '\' fileName '_thetaPhaseOverTimeligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']);
    
    files = dir([path fileName '_Raster*']);
    for ii = 1:length(files)
        fn = files(ii).name;
        idx = strfind(fn,'_');
        fnnew = [fn(1:idx(end)-1) '_Run' num2str(onlyRun) fn(idx(end):end)];
        movefile([path fn], [path fnnew]);
    end
end
