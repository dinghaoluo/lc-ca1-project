function ProcessAllRecBeh()
% process all the recordings

    RecordingListNT;
    
%     % muscimol active licking
%     nRec = length(recSessionsALMusc);
%     for i = 1:nRec
%         ind = findstr(activeLickMuscPath(i,:),'\');
%         recName = activeLickMuscPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickMuscPath(i,:),recName,recSessionsALMusc{i},...
%             ALMuscMazeSession{i},ALMuscSessions{i}(1)-1);
%     end
%     
%     % muscimol passive licking
%     nRec = length(recSessionsPLMusc);
%     for i = 1:nRec
%         ind = findstr(passiveLickMuscPath(i,:),'\');
%         recName = passiveLickMuscPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(passiveLickMuscPath(i,:),recName,recSessionsPLMusc{i},...
%             PLMuscMazeSession{i},PLMuscSessions{i}(1)-1);
%     end
%     
%     % scopolamine active licking
%     nRec = length(recSessionsALScop);
%     for i = 1:nRec
%         ind = findstr(activeLickScopPath(i,:),'\');
%         recName = activeLickScopPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickScopPath(i,:),recName,recSessionsALScop{i},...
%             ALScopMazeSession{i},ALScopSessions{i}(1)-1);
%     end
%     
%     % pirenzepine active licking
%     nRec = length(recSessionsALPire);
%     for i = 1:nRec
%         ind = findstr(activeLickPirePath(i,:),'\');
%         recName = activeLickPirePath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickPirePath(i,:),recName,recSessionsALPire{i},...
%             ALPireMazeSession{i},ALPireSessions{i}(1)-1);
%     end
%     
%     % propranolol active licking
%     nRec = length(recSessionsALProp);
%     for i = 1:nRec
%         ind = findstr(activeLickPropPath(i,:),'\');
%         recName = activeLickPropPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickPropPath(i,:),recName,recSessionsALProp{i},...
%             ALPropMazeSession{i},ALPropSessions{i}(1)-1);
%     end
%     
%     % SCH23390 active licking
%     nRec = length(recSessionsALSCH);
%     for i = 1:nRec
%         ind = findstr(activeLickSCHPath(i,:),'\');
%         recName = activeLickSCHPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickSCHPath(i,:),recName,recSessionsALSCH{i},...
%             ALSCHMazeSession{i},ALSCHSessions{i}(1)-1);
%     end
%     
%     %% saline active licking
%     nRec = length(recSessionsALSal);
%     for i = 1:nRec
%         ind = findstr(activeLickSalPath(i,:),'\');
%         recName = activeLickSalPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickSalPath(i,:),recName,recSessionsALSal{i},...
%             ALSalMazeSession{i},ALSalSessions{i}(1)-1);
%     end
%     
%     %% CNO Som-cre active licking
%     nRec = length(recSessionsALSOMCNO);
%     for i = 1:nRec
%         ind = findstr(activeLickSOMCNOPath(i,:),'\');
%         recName = activeLickSOMCNOPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickSOMCNOPath(i,:),recName,recSessionsALSOMCNO{i},...
%             ALSOMCNOMazeSession{i},ALSOMCNOSessions{i}(1)-1);
%     end
%     
%     %% Saline Som-cre active licking
%     nRec = length(recSessionsALSOMSal);
%     for i = 1:nRec
%         ind = findstr(activeLickSOMSalPath(i,:),'\');
%         recName = activeLickSOMSalPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickSOMSalPath(i,:),recName,recSessionsALSOMSal{i},...
%             ALSOMSalMazeSession{i},ALSOMSalSessions{i}(1)-1);
%     end
%     
%     %% CNO PV-cre active licking
%     nRec = length(recSessionsALPVCNO);
%     for i = 1:nRec
%         ind = findstr(activeLickPVCNOPath(i,:),'\');
%         recName = activeLickPVCNOPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickPVCNOPath(i,:),recName,recSessionsALPVCNO{i},...
%             ALPVCNOMazeSession{i},ALPVCNOSessions{i}(1)-1);
%     end
%     
%     %% Saline PV-cre active licking
%     nRec = length(recSessionsALPVSal);
%     for i = 1:nRec
%         ind = findstr(activeLickPVSalPath(i,:),'\');
%         recName = activeLickPVSalPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickPVSalPath(i,:),recName,recSessionsALPVSal{i},...
%             ALPVSalMazeSession{i},ALPVSalSessions{i}(1)-1);
%     end
%     
%     %% CNO control animals active licking
%     nRec = length(recSessionsALCNO);
%     for i = 1:nRec
%         ind = findstr(activeLickCNOPath(i,:),'\');
%         recName = activeLickCNOPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickCNOPath(i,:),recName,recSessionsALCNO{i},...
%             ALCNOMazeSession{i},ALCNOSessions{i}(1)-1);
%     end
    
    accumulateRecData();
    accumRecDataRunStat();
    accumRecDataCueStat();
    
    plotAccumulateRecData();
    plotAccumulateRecDataCue();
    plotAccumulateRecDataCorr();
end
