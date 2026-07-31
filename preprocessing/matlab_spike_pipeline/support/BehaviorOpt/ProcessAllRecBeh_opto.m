function ProcessAllRecBeh_opto(mode, exp)
% process all the recordings
% modified to fit experiments; 2nd input parameter added to improve user exp. (Dinghao, 27 Aug. 2021)

%     if(mode == 1)
%         RecordingListNT_opto1; % including the continous stimulation
%     elseif(mode == 2)
%         RecordingListNT_opto2; % does not include the continous stimulation
%     else
%         RecordingListNT_opto;
%     end
%     
%     %% Dbh active licking
%     if(exp == 2)
%         disp('Dbh active licking');
%         nRec = length(ALRecSessionsOptoDbh);
%         for i = nRec-1:nRec
%             ind = findstr(activeLickOptoDbhPath(i,:),'\');
%             recName = activeLickOptoDbhPath(i,ind(end-1)+1:ind(end)-1);
%             disp(recName);
%             ProcessRecordingBeh_opto(activeLickOptoDbhPath(i,:),recName,ALRecSessionsOptoDbh{i},...
%                 ALMazeSessionOptoDbh{i},ALMazeTypeOptoDbh{i});
%         end
%     end
% 
%     %% Dbh active licking activation
%     if(exp == 6)
%         disp('Dbh active licking activation');
%         nRec = length(ALRecSessionsOptoDbhAct);
%         for i = nRec:nRec
%             ind = findstr(activeLickOptoDbhActPath(i,:),'\');
%             recName = activeLickOptoDbhActPath(i,ind(end-1)+1:ind(end)-1);
%             disp(recName);
%             ProcessRecordingBeh_opto(activeLickOptoDbhActPath(i,:),recName,ALRecSessionsOptoDbhAct{i},...
%                 ALMazeSessionOptoDbhAct{i},ALMazeTypeOptoDbhAct{i});
%         end
%     end
% 
%     %% Dbh active licking activation (EPHYS - NEURONEXUS) THIS DOESNT WORK YET
%     if(exp == 7)
%         disp('Dbh active licking activation (EPHYS - NEURONEXUS)');
%         nRec = length(ALRecSessionsOptoDbhActEphys);
%         for i = nRec:nRec
%             ind = findstr(activeLickOptoDbhActEphysPath(i,:),'\');
%             recName = activeLickOptoDbhActEphysPath(i,ind(end-1)+1:ind(end)-1);
%             disp(recName);
%             ProcessRecordingBeh_opto(activeLickOptoDbhActEphysPath(i,:),recName,ALRecSessionsOptoDbhActEphys{i},...
%                 ALMazeSessionOptoDbhActEphys{i},ALMazeTypeOptoDbhActEphys{i});
%         end
%     end
%     
%     %% Sulpiride
%     if(exp == 1)
%         disp('SUL active licking');
%         nRec = length(ALRecSessionsSul);
%         for i = nRec:nRec
%             ind = findstr(activeLickSulPath(i,:),'\');
%             recName = activeLickSulPath(i,ind(end-1)+1:ind(end)-1);
%             disp(recName);
%             ProcessRecordingBeh_opto(activeLickSulPath(i,:),recName,ALRecSessionsSul{i},...
%                 ALMazeSessionSul{i},ALMazeTypeSul{i});
%         end
%     end
% 
%     %% SCH
%     if(exp == 3)
%         disp('SCH active licking');
%         nRec = length(ALRecSessionsSCH);
%         for i = nRec:nRec
%             ind = findstr(activeLickSCHPath(i,:),'\');
%             recName = activeLickSCHPath(i,ind(end-1)+1:ind(end)-1);
%             disp(recName);
%             ProcessRecordingBeh_opto(activeLickSCHPath(i,:),recName,ALRecSessionsSCH{i},...
%                 ALMazeSessionSCH{i},ALMazeTypeSCH{i});
%         end
%     end
% 
%     %% Saline
%     if(exp == 4)
%         disp('Saline active licking');
%         nRec = length(ALRecSessionsSal);
%         for i = nRec:nRec
%             ind = findstr(activeLickSalPath(i,:),'\');
%             recName = activeLickSalPath(i,ind(end-1)+1:ind(end)-1);
%             disp(recName);
%             ProcessRecordingBeh_opto(activeLickSalPath(i,:),recName,ALRecSessionsSal{i},...
%                 ALMazeSessionSal{i},ALMazeTypeSal{i});
%         end
%     end
%     
%     %% Dat
%     if(exp == 5)
%         disp('Dat active licking');
%         nRec = length(ALRecSessionsOptoDat);
%         for i = nRec:nRec
%             ind = findstr(activeLickOptoDatPath(i,:),'\');
%             recName = activeLickOptoDatPath(i,ind(end-1)+1:ind(end)-1);
%             disp(recName);
%             ProcessRecordingBeh_opto(activeLickOptoDatPath(i,:),recName,ALRecSessionsOptoDat{i},...
%                 ALMazeSessionOptoDat{i},ALMazeTypeOptoDat{i});
%         end
%     end
%     
%     %% Training practice
%     if(exp == 0)
%         nRec = length(ALRecSessionsOptoPractice);
%         for i = nRec:nRec
%             ind = findstr(PracticePath(i,:),'\');
%             recName = PracticePath(i,ind(end-1)+1:ind(end)-1);
%             disp(recName);
%             ProcessRecordingBeh_opto(PracticePath(i,:),recName,ALRecSessionsOptoPractice{i},...
%                ALOptoMazeSessionPractice{i},ALOptoMazeTypePractice{i});
%         end
%     end

    summExp = input('Which experiment to summarise? (enter 0 to break out)\n'); 
                         % prompt for the experiment to summarise
                         % 1. SCH-23390
                         % 2. Sulpiride
                         % 3. Saline
                         % 4. Dbh-inhibition
                         % 5. Dat-inhibition
                         % 6. Dbh-activation
    if(summExp==0)
        error('No summary.');
    end
    optoStim = input('Which stimulation condition to summarise? (enter 0 to break out)\n'); 
                         % prompt for the stim condition to summarise
                         % 2. run-onset stim
                         % 3. mid-trial stim
                         % 4. reward-delivery stim
    if(optoStim==0)
        error('No summary.');
    end
    accumulateRecData_opto(mode, summExp, optoStim);
    accumRecDataRunStat_opto();
    accumRecDataCueStat_opto();
    accumRecDataRewStat_opto();
    
    % all the trials during opt
    plotAccumulateRecData_opto(mode, 0);
%     plotAccumulateRecData_opto(mode,1);
    plotAccumulateRecData_optovsCtrl(mode);
    plotAccumulateRecDataCue_opto(mode, 0);
%     plotAccumulateRecDataCue_opto(mode,1);
    plotAccumulateRecDataCue_optovsCtrl(mode);
    plotAccumulateRecDataRew_opto(mode, 0);
%     plotAccumulateRecDataRew_opto(mode,1);
    plotAccumulateRecDataRew_optovsCtrl(mode);
end
