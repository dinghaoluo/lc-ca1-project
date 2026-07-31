function [fieldArr] = getRecField(trials,fieldName,lapList)
% This function is to get the EEG from the trial data structure
% Inputs:
% trials:       the trial data structure
% fieldName:    field name string
% lapList:      the array recording the valid trials
%
% Outputs:
% fieldArr:     the requested field

fieldArr = [];
fieldExist = 0;

for i = 1:length(lapList)
    if(~isempty(trials{lapList(i)}))
        if(fieldExist == 0)
            if(~isfield(trials{lapList(i)},fieldName))
                disp(['Field ' fieldName ' does not exist.']);
                fieldArr{i} = {};
                return;
            else
                fieldExist = 1;
            end
        end
        fieldArr{i} = trials{lapList(i)}.(genvarname(fieldName));
    else
        fieldArr{i} = [];
    end
end
