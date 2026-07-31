function GenerateBehavDataStructures_smallTMPFI(path, baseFileName)
% This function reads behavior file and produce data-structures that
% align the behavioral data 
% ouput structures:
%
%
% e.g.: GenerateBehavDataStructures_smallTMPFI
%  ('Z:\Raphael_tests\mice_expdata\ANM001\A001-20180928', 'A001-20180928-01')

    cd(path);
    GlobalConst;
      
    if exist([baseFileName 'BLFP.mat'], 'file') == 2
        disp('BTDT file already exists.')
    else
        fullNameB = [baseFileName 'B.mat'];
        if(exist(fullNameB,'file') ~= 0)
            load(fullNameB);
        else
            disp('parsing the behavioral file.');
            behEvents = LoadBehMazeFile_smTr([baseFileName 'T.txt']);
                      
            save([baseFileName 'B.mat'], 'behEvents');
        end
        Arduino2LFPtime_smTr(baseFileName,sampleFq);
    end
    
    % interpolate tracking data
    fprintf('\nInterpolate tracking data...\n')
    InterpolateWheel_smTr(baseFileName,sampleFq);
    
    % get tracking related data
    fprintf('\nGet Track and Laps...\n')
    SortTrials_Beh_smTr(baseFileName, sampleFq);

end
