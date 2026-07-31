function RunSuite2pPipeline()
    % built on Jingyu and Yingxue's code
    % runs all suite2p pipeline functions
    % Dinghao 8 Jan, 2024
    
    % DIRECTORY ONLY WORKS WITH DINGHAO'S RECORDINGS
    % CHANGE DIRECTORY IF ERROR

    % in order to read .npy files...
    addpath('Z:\Dinghao\code\npy-matlab-master\npy-matlab')

    % read session list
    RecordingList2P_dinghao;
    paths = listRecordingActiveBehPathdLight;

    % define processing mode
    mode = 'axons_v2.0'; % the .npy file (mode) with which suite2p was run

    % get session count
    totSess = length(paths);

    % extract info from each string
    % merge coresponding running data into images' folder
    for i = 1:totSess

        % string manipulations
        anmlDate = paths(i,1:13);  % e.g. A054-2023031
        anmlNo = paths(i,1:4);  % e.g. A054
        session = paths(i,15:16);  % e.g. 03
        filebase = [anmlDate '-' session];  % e.g. A054-20230310-03
        
        % announce 
        print(filebase)

        % create folders for 2chan
        chan1Folder = 'Z:\Dinghao\2p_recording\'+anmlNo+'\'+anmlDate+'\'+session+'\'+mode+'\' + 'Channel1';
        chan2Folder = 'Z:\Dinghao\2p_recording\'+anmlNo+'\'+anmlDate+'\'+session+'\'+mode+'\' + 'Channel2';
        mkdir(chan1Folder); 
        mkdir(chan2Folder);

        % get running log files
        if anmlNo=='A054' || anmlNo=='A058' || anmlNo=='A074'
            copyRelatedFiles(('Z:\Jingyu\mice-expdata\' + anmlNo), extractBetween(anmlDate, 6, 13), session, chan1Folder);
            copyRelatedFiles(('Z:\Jingyu\mice-expdata\' + anmlNo), extractBetween(anmlDate, 6, 13), session, chan2Folder);
        else
            print('no action was allocated to other animals at the moment');
%             copyRelatedFiles(('Z:\Jingyu\mice-expdata\' + anmlNo), extractBetween(anmlDate, 6, 13), session, chan1Folder);
%             copyRelatedFiles(('Z:\Jingyu\mice-expdata\' + anmlNo), extractBetween(anmlDate, 6, 13), session, chan2Folder);
        end

        % change folder and read channel 2 data
        dataFolder = 'Z:\Dinghao\2p_recording\'+anmlNo+'\'+anmlDate+'\'+session+'\'+mode+ '\suite2p\plane0\';
        cd(dataFolder);
        load('Fall.mat');
        F = readNPY('F_chan2.npy');
        Fneu = readNPY('Fneu_chan2.npy');
        save('Fall_chan2.mat');

        %imageFolder  = targetFolder;
        chan1matFile =  dataFolder + 'Fall.mat';
        chan2matFile = dataFolder + 'Fall_chan2.mat';
        %beginPipeline = anmlDate + '-' + session;

        % Run pipeline for each session
        cd Z:\Dinghao\code\suite2P_pipeline_running;

        disp(paths(i,:)+'GenerateBehav2PDataStructures--------Channel1');
        tic
        GenerateBehav2PDataStructures_smTrMPFI2Pv1(char(chan1matFile),char(chan1Folder),char(filebase),0);
        toc
        disp(paths(i,:)+'Finished_GenerateBehav2PDataStructures--------Channel1');

        disp(paths(i,:)+'GenerateBehav2PDataStructures--------Channel2');
        tic
        GenerateBehav2PDataStructures_smTrMPFI2Pv1(char(chan2matFile),char(chan2Folder),char(filebase),0);
        toc
        disp(paths(i,:)+'Finished_GenerateBehav2PDataStructures--------Channel2');

    end

    %% run pipeline step2

    for i = 1:totSess

        % Extract info
        anmlDate = extractBetween(paths(i,:), 1, 13);
        anmlNo = extractBetween(paths(i,:), 1, 4);
        session = extractBetween(paths(i,:), 15, 16);
        filebase = anmlDate +  '-' + session;

        chan1Folder = 'Z:\Dinghao\2p_recording\'+anmlNo+'\'+anmlDate+'\'+session+'\'+mode+'\' + 'Channel1';
        chan2Folder = 'Z:\Dinghao\2p_recording\'+anmlNo+'\'+anmlDate+'\'+session+'\'+mode+'\' + 'Channel2';

        disp(paths(i,:)+'RunProcessing--------Channel1');
        tic
        RunProcessing2P(char(chan1Folder+'\'),char(filebase),1, 1);
        toc
        disp(paths(i,:)+'Finished_RunProcessing--------Channel1');

        disp(paths(i,:)+'RunProcessing--------Channel2');
        tic
        RunProcessing2P(char(chan2Folder+'\'),char(filebase),1,1);
        toc
        disp(paths(i,:)+'Finished_RunProcessing--------Channel2');

    end
end

%%
function copyRelatedFiles(sourceFolder, searchString, searchSession, destinationFolder)
    % sourceFolder: The folder where you want to search for related files
    % searchString: The string you want to search for in file names
    % destinationFolder: The folder where you want to copy the related files

    % Find all files (including in subfolders) in the source folder
    cd (sourceFolder);
    files = dir('*.txt');

    % Initialize an empty list to store related file names
    relatedFiles = {};

    % Loop through all the files and check if they contain the searchString
    for i = 1:length(files)
        if contains(files(i).name, (searchString + "-" + searchSession), 'IgnoreCase', true)
            relatedFiles{end+1} = files(i).name;
        end
    end

    % If any related files are found, copy them to the destination folder
    if ~isempty(relatedFiles)
        for i = 1:length(relatedFiles)
            sourceFile = fullfile(files(i).folder, relatedFiles{i});

            destinationFile = fullfile(destinationFolder, relatedFiles{i});
            copyfile(sourceFile, destinationFile);
        end
        disp(searchString + "-" + searchSession + ', related files copied successfully!');
    else
        disp('No related files found.');
    end
end
