function RunAllFunctions
    % Runs all of the code necessary for identifying run bouts, aligning
    % neural activity to run bouts, identifying Rise/Down neurons based off
    % of run bouts, and making summary plots
    
    % Global 
    
    addpath('Z:\Yingxue\code\codePVPaper\'); % path to GlobalConstFq
    
    % Identify run bouts for each recording
    save_off_target_run_bouts_all;
    
    % Save table with # of Run Bouts per recording
    HowManyRunBouts;
    
    % Align FSA to run bouts for each recording     
    ConvSpikeTrain_AlignedRunBout_wrapper;
    
    % Save table of averaged fsa aligned to runbout for each recording
    SaveAllAverageFSA_RunBout;
    
%     % Identify Rise and Down Neurons based on run bouts
%     IdentifyRiseDownRunBouts_all(1);
    
%     % Plot average firing rate profile for Rise and Down neurons
%     FR_Profile_Overall('Rise',1);
%     FR_Profile_Overall('Down',1);
%     FR_Profile_Overall('Rise',0);
%     FR_Profile_Overall('Down',0);
%     % Pie chart of proportion of old rise/down neurons that also have
%     % strong responses to run bouts
%     PieChartPropRecruitedOldRiseDown;
%     
%     % Pie chart of proportion of neurons that have responses to the trial
%     % ro and the run bout onset
%     PieChartPropRiseDown;
%     
%     % Calculates the ratio0to1BefRun for the old Rise and Down neurons to
%     % the trial run onset and the run bout onset
%     FR_Ratio_RiseDownRunBouts('Rise');
%     FR_Ratio_RiseDownRunBouts('Down');
%     
%     % Create a bargraph comparing proportion of rise and down neurons per
%     % recording and condition
%     
%     BarGraphPropRiseDowntrialROvsRunBout;
end
