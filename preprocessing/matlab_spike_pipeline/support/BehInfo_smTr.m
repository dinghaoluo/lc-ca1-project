function BehInfo_smTr(path,fileName)
% This is to extract the behavior information from the session 
% Input arguments:
% path:         the path of the recording file
% fileName:     name of the recording file
%
% e.g.: BehInfo_smTr('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1')

    %%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin > 3
        disp('Too many arguments');        
        return;
    end
    
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(isempty(indexFileName))
        fileNameMinLen = [fileName '_BehInfo.mat'];
        fileName = [fileName '.mat'];
    else
        fileNameMinLen = [fileName(1:indexFileName(end)-1) '_BehInfo.mat'];
    end 
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath);
    
    GlobalConst;
    
    for i = 1:size(lapList,1)
        if(lapList(i,1) ~= -1 && ~isempty(trials{i}))
            %% figure out stop times
            stopInd = trials{i}.speed < minSpeed;
            s = sprintf('%d',stopInd);
            [t,p]=textscan(s,'%s','delimiter','0','multipleDelimsAsOne',1);
            d = t{:};
            dp = p{:};
            data = zeros(length(d));
            for j = 1:length(d)
                data(k) = length(d{j});
            end
            [behInfo.stop_number_times{i} behInfo.stop_length{i}] = ...
                hist(data, [1:max(data)]);
            [behInfo.max_stop_len(i),ind] = max(behInfo.stop_length{i});
            behInfo.max_stop_xMMAll = trials{i}.xMMAll(dp{ind});
            behInfo.max_stop_xMM = trials{i}.xMM(dp{ind});
            
            %% figure out lick times
            behInfo.lick_ind{i} = trials{i}.lickLfpInd;
            behInfo.lick_xMMAll{i} = trials{i}.xMMAll(trials{i}.lickLfpInd); 
            behInfo.lick_speedAll{i} = trials{i}.speedAll(trials{i}.lickLfpInd);
            
        end
    end
    
end
