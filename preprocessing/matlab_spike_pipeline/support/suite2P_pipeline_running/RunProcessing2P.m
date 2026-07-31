function RunProcessing2P(BehaviorPath, FileNameBase,onlyRun,mazeSess)
%% second stage processing
%% e.g.: RunProcessing2P('Z:\Kori\2P_free_moving\Kori_data_2P\20220203\ANMC105\A105-20220203-01T\', 'A105-20220203-01', 1, 1)

disp('***Final processing....')
ProcessingMice_smTr2P(BehaviorPath,[FileNameBase '_DataStructure_mazeSection1_TrialType1'],onlyRun);

disp('***Calculate licking and speed')
ProcessingMice_smTrSpeedLick2P(BehaviorPath,[FileNameBase '_DataStructure_mazeSection1_TrialType1'],onlyRun,mazeSess);

disp('***Align to run onset, during run')
ProcessingAligned2P(BehaviorPath,[FileNameBase '_DataStructure_mazeSection1_TrialType1'],1,mazeSess);

disp('***Align to run onset, including immobile period')
ProcessingAligned2P(BehaviorPath,[FileNameBase '_DataStructure_mazeSection1_TrialType1'],0,mazeSess);
