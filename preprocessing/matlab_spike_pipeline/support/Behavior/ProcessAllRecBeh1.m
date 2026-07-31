function ProcessAllRecBeh1()
% process all the recordings
    RecordingListNT1;
    
%    
%     % muscimol active licking
%     nRec = length(recSessionsALMusc);
%     for i = 1:nRec
%         ind = findstr(activeLickMuscPath(i,:),'\');
%         recName = activeLickMuscPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickMuscPath(i,:),recName,recSessionsALMusc{i},...
%             ALMuscMazeSession{i},ALMuscSessions{i}(1)-1);
%     end
% %     
% %     % muscimol passive licking
% %     nRec = length(recSessionsPLMusc);
% %     for i = 1:nRec
% %         ind = findstr(passiveLickMuscPath(i,:),'\');
% %         recName = passiveLickMuscPath(i,ind(end-1)+1:ind(end)-1);
% %         disp(recName);
% %         ProcessRecordingBeh(passiveLickMuscPath(i,:),recName,recSessionsPLMusc{i},...
% %             PLMuscMazeSession{i},PLMuscSessions{i}(1)-1);
% %     end
% %     
%     % scopolamine active licking
%     nRec = length(recSessionsALScop);
%     for i = 1:nRec
%         ind = findstr(activeLickScopPath(i,:),'\');
%         recName = activeLickScopPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickScopPath(i,:),recName,recSessionsALScop{i},...
%             ALScopMazeSession{i},ALScopSessions{i}(1)-1);
%     end
    
% %     % pirenzepine active licking
% %     nRec = length(recSessionsALPire);
% %     for i = 1:nRec
% %         ind = findstr(activeLickPirePath(i,:),'\');
% %         recName = activeLickPirePath(i,ind(end-1)+1:ind(end)-1);
% %         disp(recName);
% %         ProcessRecordingBeh(activeLickPirePath(i,:),recName,recSessionsALPire{i},...
% %             ALPireMazeSession{i},ALPireSessions{i}(1)-1);
% %     end
% %     
% %     % propranolol active licking
% %     nRec = length(recSessionsALProp);
% %     for i = 1:nRec
% %         ind = findstr(activeLickPropPath(i,:),'\');
% %         recName = activeLickPropPath(i,ind(end-1)+1:ind(end)-1);
% %         disp(recName);
% %         ProcessRecordingBeh(activeLickPropPath(i,:),recName,recSessionsALProp{i},...
% %             ALPropMazeSession{i},ALPropSessions{i}(1)-1);
% %     end
%     
% %     % SCH23390 active licking
% %     nRec = length(recSessionsALSCH);
% %     for i = 1:3
% %         ind = findstr(activeLickSCHPath(i,:),'\');
% %         recName = activeLickSCHPath(i,ind(end-1)+1:ind(end)-1);
% %         disp(recName);
% %         ProcessRecordingBeh(activeLickSCHPath(i,:),recName,recSessionsALSCH{i},...
% %             ALSCHMazeSession{i},ALSCHSessions{i}(1)-1);
% %     end
%     
% %     %% saline active licking
% %     nRec = length(recSessionsALSal);
% %     for i = 1:nRec
% %         ind = findstr(activeLickSalPath(i,:),'\');
% %         recName = activeLickSalPath(i,ind(end-1)+1:ind(end)-1);
% %         disp(recName);
% %         ProcessRecordingBeh(activeLickSalPath(i,:),recName,recSessionsALSal{i},...
% %             ALSalMazeSession{i},ALSalSessions{i}(1)-1);
% %     end
% %     
%   %% CNO Som-cre active licking (1mM)
%     nRec = length(recSessionsALSOMCNO1mM);
%     for i = nRec
%         ind = findstr(activeLickSOMCNO1mMPath(i,:),'\');
%         recName = activeLickSOMCNO1mMPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickSOMCNO1mMPath(i,:),recName,recSessionsALSOMCNO1mM{i},...
%             ALSOMCNO1mMMazeSession{i},ALSOMCNO1mMSessions{i}(1)-1);
%     end
% %     
%     %% CNO Som-cre active licking (3mM)
%     nRec = length(recSessionsALSOMCNO3mM);
%     for i = 1:nRec
%         ind = findstr(activeLickSOMCNO3mMPath(i,:),'\');
%         recName = activeLickSOMCNO3mMPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickSOMCNO3mMPath(i,:),recName,recSessionsALSOMCNO3mM{i},...
%             ALSOMCNO3mMMazeSession{i},ALSOMCNO3mMSessions{i}(1)-1);
%     end
% %     
% %   Saline Som-cre active licking
%     nRec = length(recSessionsALSOMSal);
%     for i = nRec
%         ind = findstr(activeLickSOMSalPath(i,:),'\');
%         recName = activeLickSOMSalPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickSOMSalPath(i,:),recName,recSessionsALSOMSal{i},...
%             ALSOMSalMazeSession{i},ALSOMSalSessions{i}(1)-1);
%     end
% % %     
% %    CNO PV-cre active licking (1mM)
%     nRec = length(recSessionsALPVCNO1mM);
%     for i = nRec
%         ind = findstr(activeLickPVCNO1mMPath(i,:),'\');
%         recName = activeLickPVCNO1mMPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickPVCNO1mMPath(i,:),recName,recSessionsALPVCNO1mM{i},...
%             ALPVCNO1mMMazeSession{i},ALPVCNO1mMSessions{i}(1)-1);
%     end
% % %     
%     % CNO PV-cre active licking (3mM)
%     nRec = length(recSessionsALPVCNO3mM);
%     for i = nRec
%         ind = findstr(activeLickPVCNO3mMPath(i,:),'\');
%         recName = activeLickPVCNO3mMPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickPVCNO3mMPath(i,:),recName,recSessionsALPVCNO3mM{i},...
%             ALPVCNO3mMMazeSession{i},ALPVCNO3mMSessions{i}(1)-1);
%     end
% %     
%     %% Saline PV-cre active licking
%     nRec = length(recSessionsALPVSal);
%     for i = nRec
%         ind = findstr(activeLickPVSalPath(i,:),'\');
%         recName = activeLickPVSalPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickPVSalPath(i,:),recName,recSessionsALPVSal{i},...
%             ALPVSalMazeSession{i},ALPVSalSessions{i}(1)-1);
%     end
    
%     %% CNO control animals active licking
%     nRec = length(recSessionsALCNO);
%     for i = nRec
%         ind = findstr(activeLickCNOPath(i,:),'\');
%         recName = activeLickCNOPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickCNOPath(i,:),recName,recSessionsALCNO{i},...
%             ALCNOMazeSession{i},ALCNOSessions{i}(1)-1);
%     end
    
%     %% atropine-DART active licking
%     nRec = length(recSessionsALAtroDART);
%     for i =  nRec
%         ind = findstr(activeLickAtroDARTPath(i,:),'\');
%         recName = activeLickAtroDARTPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickAtroDARTPath(i,:),recName,recSessionsALAtroDART{i},...
%             ALAtroDARTMazeSession{i},ALAtroDARTSessions{i}(1)-1);
%     end
% %      
%     % atropine-DART control active licking
%     nRec = length(recSessionsALAtroDARTCtrl);
%     for i = nRec
%         ind = findstr(activeLickAtroDARTCtrlPath(i,:),'\');
%         recName = activeLickAtroDARTCtrlPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickAtroDARTCtrlPath(i,:),recName,recSessionsALAtroDARTCtrl{i},...
%             ALAtroDARTCtrlMazeSession{i},ALAtroDARTCtrlSessions{i}(1)-1);
%     end
% % 
%     % CNO Som-cre active 3mM Kori
%     nRec = length(recSessionsALSOMCNO3mMKori);
%     for i = nRec
%         ind = findstr(activeLickSOMCNO3mMKoriPath(i,:),'\');
%         recName = activeLickSOMCNO3mMKoriPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickSOMCNO3mMKoriPath(i,:),recName,recSessionsALSOMCNO3mMKori{i},...
%             ALSOMCNO3mMKoriMazeSession{i},ALSOMCNO3mMKoriSessions{i}(1)-1);
%    end
%     % CNO Som-cre active 1mM Kori
%     nRec = length(recSessionsALSOMCNO1mMKori);
%     for i = 1:nRec
%         ind = findstr(activeLickSOMCNO1mMKoriPath(i,:),'\');
%         recName = activeLickSOMCNO1mMKoriPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickSOMCNO1mMKoriPath(i,:),recName,recSessionsALSOMCNO1mMKori{i},...
%             ALSOMCNO1mMKoriMazeSession{i},ALSOMCNO1mMKoriSessions{i}(1)-1);
%       end
    % CNO Som-cre Saline active  Kori
%     nRec = length(recSessionsALSOMSalineKori);
%     for i = 1:nRec
%         ind = findstr(activeLickSOMSalineKoriPath(i,:),'\');
%         recName = activeLickSOMSalineKoriPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickSOMSalineKoriPath(i,:),recName,recSessionsALSOMSalineKori{i},...
%             ALSOMSalineKoriMazeSession{i},ALSOMSalineKoriSessions{i}(1)-1);
%     end
% % 
%     % CNO PV-cre active 3mM Kori
%     nRec = length(recSessionsALPVCNO3mMKori);
%     for i = nRec
%         ind = findstr(activeLickPVCNO3mMKoriPath(i,:),'\');
%         recName = activeLickPVCNO3mMKoriPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickPVCNO3mMKoriPath(i,:),recName,recSessionsALPVCNO3mMKori{i},...
%             ALPVCNO3mMKoriMazeSession{i},ALPVCNO3mMKoriSessions{i}(1)-1);
%     end

%     % CNO PV-cre active 1mM Kori
%     nRec = length(recSessionsALPVCNO1mMKori);
%     for i = 1:nRec
%         ind = findstr(activeLickPVCNO1mMKoriPath(i,:),'\');
%         recName = activeLickPVCNO1mMKoriPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickPVCNO1mMKoriPath(i,:),recName,recSessionsALPVCNO1mMKori{i},...
%             ALPVCNO1mMKoriMazeSession{i},ALPVCNO1mMKoriSessions{i}(1)-1);
%      end

%     % CNO PV-cre Saline active  Kori
%     nRec = length(recSessionsALPVSalineKori);
%     for i = nRec
%         ind = findstr(activeLickPVSalineKoriPath(i,:),'\');
%         recName = activeLickPVSalineKoriPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickPVSalineKoriPath(i,:),recName,recSessionsALPVSalineKori{i},...
%             ALPVSalineKoriMazeSession{i},ALPVSalineKoriSessions{i}(1)-1);
%     end

% %     % Compound 21 Som-cre active 3mM 
%     nRec = length(recSessionsALSOMC213mM);
%     for i = nRec
%         ind = findstr(activeLickSOMC213mMPath(i,:),'\');
%         recName = activeLickSOMC213mMPath(i,ind(end-1)+1:ind(end)-1);
%         disp(recName);
%         ProcessRecordingBeh(activeLickSOMC213mMPath(i,:),recName,recSessionsALSOMC213mM{i},...
%             ALSOMC213mMMazeSession{i},ALSOMC213mMSessions{i}(1)-1);
%    end
% 
     accumulateRecData();
     accumRecDataRunStat();
%      accumRecDataCueStat();
% % % % % % % % %     
       plotAccumulateRecData();
%         plotAccumulateRecDataCue();
%     plotAccumulateRecDataCorr();
end
