function [trials, cluList, lapList] = GetTrials_smTrMPFI(TwoPhotonPath, BehaviorPath, filename)
%%% added delay area wheel running task

% [trials, cluList, lapList] = GetTrials_JF_v4(mazeSection, mazeTrialType, OverwriteTrials, OverwriteAll);
%
% function outputs a cell data structure organized by laps
% input: filename;
% optional input: maze section (see bellow), trial type for maze section 13 (see bellow):
%
% maze section = ID 
%-----------------------------------------
% treadmill -- 1
%
% trial type for maze section between 2-6
% 0 - the whole treadmill run
%
% examples:
%    GetTrials_JF_v4('A132-20111025-01',1,1,0,0);
    cd(BehaviorPath); 

    if exist([filename '_Behav2PDataLFP.mat'],'file') == 2 
        fprintf('\n Track, Laps, Clu structures loaded from a file.\n');
        load([filename '_Behav2PDataLFP.mat'],...
            'Track','Laps','Clu','lfpFreq');
    else
        GenerateBehav2PDataStructures_smTrMPFI2Pv1(TwoPhotonPath, BehaviorPath, filename);
        load([filename '_Behav2PDataLFP.mat'],...
                'Track','Laps','Clu','lfpFreq'); 
    end
               
    totNLaps = length(Laps.lapID);
    nTotClu = length(Clu.localClu);   
   
    %% check whether the trials are already extracted and stored 
    stimOnlyN = 0;
    mazeSection = 1;
    mazeTrialType = 1;
    if(exist([filename '_DataStructure_mazeSection' num2str(mazeSection) ...
            '_TrialType' num2str(mazeTrialType) '.mat'],'file') == 2)
        stimOnly=input(['The DataStructure file already exists, do you want to skip' ...
                  '\nthe trial extraction, Y/N: '],'s');
        if(isempty(stimOnly))
            stimOnly = 'N';
        end
        if(~isempty(strfind(stimOnly,'Y')) || ~isempty(strfind(stimOnly,'y')))
           stimOnlyN = 1; 
        end
    end
    
    if(stimOnlyN == 0)
        %% look for the start and stop indices of each section
        lfpIndStartAll = cell(1,totNLaps);
        lfpIndEndAll = cell(1,totNLaps);
        if mazeSection == 1 
            %%%%%%% RUNNING WHEEL %%%%%%%%%%%
            % select start and end point of trials based on whlTrialType variable
            for i = 1:totNLaps                
                ind = find(Laps.lapID == i,1);
                if(isempty(ind) || Laps.behavType(i) ~= 1)
                    continue;
                else
                    lfpIndStartAll{i} = Laps.startLfpInd(i);
                    lfpIndEndAll{i} = Laps.endLfpInd(i);

                    figure(10);
                    plot((Laps.startLfpInd(i):Laps.endLfpInd(i))./lfpFreq,...
                        Track.xMM(Laps.startLfpInd(i):Laps.endLfpInd(i)),'r.')

                    title(['Trial #' num2str(i)]);
                    xlabel('Time (s)')
                    ylabel('Distance (mm)')
%                     pause;
                end
            end
        end
        close(findobj('type','figure','CurrentObject',10));     
            % close fig 10 if it exists
    
    
        %% collect the data according to the start and stop indices
        disp('Assigning values to trials....');
        goodLaps = Laps.behavType;
        totT = 0;
        trials = cell(1,totNLaps);
        for nlap=1:totNLaps
%             if(i == 194)
%                 a = 1;
%             end
            disp(['Trial no. ' num2str(nlap)])
            trials{nlap} = [];
            if(~isempty(lfpIndStartAll{nlap}))
                totT = totT + lfpIndEndAll{nlap} - lfpIndStartAll{nlap} + 1;
                trials{nlap}.lfpIndStart = lfpIndStartAll{nlap};
                trials{nlap}.lfpIndEnd = lfpIndEndAll{nlap};
                trials{nlap}.Nsamples = lfpIndEndAll{nlap} ...
                    - lfpIndStartAll{nlap} + 1;

                trials{nlap}.corrTrial = Laps.corrChoice(nlap);
                % determine behavior for each lap (good-1)
                trials{nlap}.behavChar = Laps.behavType(nlap);
                trials{nlap}.trackLen = Laps.trackLen(nlap);
                trials{nlap}.lickLfpInd = Laps.lickLfpInd{nlap}...
                    - lfpIndStartAll{nlap}+1;
                trials{nlap}.pumpLfpInd = Laps.pumpLfpInd{nlap}...
                    - lfpIndStartAll{nlap}+1;  
                if(isfield(Laps,'airpuffLfpInd'))
                    trials{nlap}.airpuffLfpInd = Laps.airpuffLfpInd{nlap}...
                        - lfpIndStartAll{nlap}+1;
                else
                    trials{nlap}.airpuffLfpInd = [];
                end
                
                trials{nlap}.movieOnLfpInd = Laps.movieOnLfpInd{nlap}...
                        - lfpIndStartAll{nlap}+1;
                                
                %%%%%%% added by Yingxue on 03/11/2020
                if(isfield(Laps,'movieLocation'))
                    trials{nlap}.movieLocation = Laps.movieLocation(nlap);
                end
                %%%%%%
                
                if(isfield(Laps,'lickPeriodLfpInd'))
                    trials{nlap}.lickPeriodLfpInd = Laps.lickPeriodLfpInd{nlap}...
                        - lfpIndStartAll{nlap}+1;
                end
                
                if(isfield(Laps,'stimOnLfpInd'))
                    if(~isempty(Laps.stimOnLfpInd{nlap}))
                        trials{nlap}.stimOnLfpInd = Laps.stimOnLfpInd{nlap}...
                            - lfpIndStartAll{nlap}+1;
                        trials{nlap}.stimPulseLfpInd = Laps.stimPulseLfpInd{nlap}...
                            - lfpIndStartAll{nlap}+1;
                        trials{nlap}.stimPulseWidth = Laps.stimPulseWidth{nlap}/1000*lfpFreq;
                        trials{nlap}.stimDiode = Laps.stimDiode{nlap};
                        trials{nlap}.stimDiodeCurr = Laps.stimDiodeCurr{nlap};
                        trials{nlap}.stimPulseMethod = Laps.stimPulseMethod{nlap};
                        trials{nlap}.stimPulseLoc = Laps.stimPulseLoc{nlap};
                    else
                        trials{nlap}.stimOnLfpInd = [];
                        trials{nlap}.stimPulseLfpInd = [];
                        trials{nlap}.stimPulseWidth = [];
                        trials{nlap}.stimDiode = [];
                        trials{nlap}.stimDiodeCurr = [];
                        trials{nlap}.stimPulseMethod = Laps.stimPulseMethod{nlap};
                        trials{nlap}.stimPulseLoc = [];
                    end
                end
                
                % location and distance
                trials{nlap}.xMM = Track.xMM(lfpIndStartAll{nlap}:lfpIndEndAll{nlap});
                trials{nlap}.xMMAll = Track.xMMAll(lfpIndStartAll{nlap}:lfpIndEndAll{nlap});
                
                %  neuro data
                trials{nlap}.dFF = Track.dFF(lfpIndStartAll{nlap}:lfpIndEndAll{nlap}, :);
                trials{nlap}.F = Track.F(lfpIndStartAll{nlap}:lfpIndEndAll{nlap}, :);
                trials{nlap}.Fneu = Track.Fneu(lfpIndStartAll{nlap}:lfpIndEndAll{nlap}, :);
                trials{nlap}.spikes = Track.spks(lfpIndStartAll{nlap}:lfpIndEndAll{nlap}, :);
                trials{nlap}.dFFSM = Track.dFFSM(lfpIndStartAll{nlap}:lfpIndEndAll{nlap}, :);
                trials{nlap}.dFFGF = Track.dFFGF(lfpIndStartAll{nlap}:lfpIndEndAll{nlap}, :);
                trials{nlap}.spikesSM = Track.spikesSM(lfpIndStartAll{nlap}:lfpIndEndAll{nlap}, :);
                trials{nlap}.indNoise = Track.indNoise(lfpIndStartAll{nlap}:lfpIndEndAll{nlap}, :);

                % speed and distance
                trials{nlap}.speed = ...
                    Track.speed_MMsec(lfpIndStartAll{nlap}:lfpIndEndAll{nlap});
                if(trials{nlap}.speed(1) < 0)
                    trials{nlap}.speed(1) = trials{nlap}.speed(2);
                end
                trials{nlap}.speedAll = ...
                    Track.speed_MMsecAll(lfpIndStartAll{nlap}:lfpIndEndAll{nlap});
                if(trials{nlap}.speedAll(1) < 0)
                    trials{nlap}.speedAll(1) = trials{nlap}.speedAll(2);
                end

            end
        end
            
        cluList = Clu;
        
        lapList = goodLaps;
        lap.trackLen = Laps.trackLen;
        lap.mazeSess = Laps.mazeSess;
        lap.mazeType = Laps.mazeType;
            
        if (mazeSection == 1)
            fprintf('\nOutput structures saved: %s\n\n', ...
                [filename '_DataStructure_mazeSection' num2str(mazeSection) ...
                '_TrialType' num2str(mazeTrialType) '.mat']);
            save([filename '_DataStructure_mazeSection' num2str(mazeSection) ...
                '_TrialType' num2str(mazeTrialType) '.mat'], ...
                'trials', 'cluList', 'lapList', 'lap', '-v7.3');
        end
    end
end
