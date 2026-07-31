function BasicInfo_smTr2P(path,fileName)
% This is to extract the basic information from the session (heating and cooling)
% Input arguments:
% path:         the path of the recording file
% fileName:     name of the recording file
%
% e.g.: BasicInfo_smTr('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1')

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
        fileNameMinLen = [fileName '_Info.mat'];
        fileName = [fileName '.mat'];
    else
        fileNameMinLen = [fileName(1:indexFileName(end)-1) '_Info.mat'];
    end 
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath);
    
    GlobalConst2P;

    rec.numNeurons = size(trials{1}.F,2);
    
    %%%%%%%%% collect the basic information from a session    
    beh.trackLen = lap.trackLen';
    beh.mazeSess = lap.mazeSess;
    beh.mazeType = lap.mazeType;
    beh.mazeSessAll = unique(lap.mazeSess);
    if(sum(beh.mazeSessAll~=0) ~= 0)
        beh.mazeSessAll = beh.mazeSessAll(beh.mazeSessAll~=0);
    end
    beh.numSamples = floor(lap.trackLen')+1; 
    beh.numTrials = size(lapList,1);
    beh.lenTrials = zeros(1,beh.numTrials); 
    beh.indGoodLap = [];
    beh.indErrLap = [];
    beh.indCorrLap = [];
    beh.indStimLap = zeros(1,beh.numTrials); 
    beh.pulseMethod = zeros(1,beh.numTrials);
    beh.stimPulseLfpInd = cell(1,beh.numTrials);
    beh.stimPulseWidth = cell(1,beh.numTrials);
    beh.stimDiode = cell(1,beh.numTrials);
    beh.stimDiodeCurr = cell(1,beh.numTrials); 
    beh.stimPulseLoc = cell(1,beh.numTrials);
    %%% added by Yingxue on 12/14/2020
    beh.numLickBef100cm = zeros(1,beh.numTrials);
    beh.numLick100to180cm = zeros(1,beh.numTrials);
    beh.numLickAft180cm = zeros(1,beh.numTrials);
    beh.numLickAftRew = zeros(1,beh.numTrials);
    beh.indTrCtrl = [];
    beh.indGoodTrCtrl = [];
    beh.indBadTrCtrl = [];
    beh.numRun = zeros(1,beh.numTrials);
    beh.medDistFirst5Licks = zeros(1,beh.numTrials);
    beh.meanDistFirst5Licks = zeros(1,beh.numTrials);
    beh.maxRunSegStartDist = zeros(1,beh.numTrials);
    beh.maxRunSegStopDist = zeros(1,beh.numTrials);
    beh.totStopTime = zeros(1,beh.numTrials); 
    %%%
    
    noStimTr = 0;
    smoothSpan = 100;  
    for i = 1:beh.numTrials
        if(lapList(i,1) ~= -1 && ~isempty(trials{i}))
            beh.indRunInLap{i} = trials{i}.speed > minSpeed;
            beh.lenTrials(i) = trials{i}.Nsamples; % the length of each trial in time
            %% licks; added by Yingxue on 12/14/2020
            beh.numLickBef100cm(i) = sum(trials{i}.xMM(trials{i}.lickLfpInd) < 1000 & ...
                trials{i}.xMM(trials{i}.lickLfpInd) >= 300);
            beh.numLick100to180cm(i) = sum(trials{i}.xMM(trials{i}.lickLfpInd) >= 1000 & ...
                trials{i}.xMM(trials{i}.lickLfpInd) < 1800);
            beh.numLickAft180cm(i) = sum(trials{i}.xMM(trials{i}.lickLfpInd) >= 1800);
            if(~isempty(trials{i}.pumpLfpInd))
                beh.numLickAftRew(i) = sum(trials{i}.lickLfpInd >= trials{i}.pumpLfpInd(1));
            else
                disp(['Trial no. ' num2str(i) ' is not rewarded. ']);
            end
            indAft30cm = find(trials{i}.xMM(trials{i}.lickLfpInd) > 300);
            numLicks = length(indAft30cm);
            if(numLicks >= 5)
                licksTmp = trials{i}.lickLfpInd(indAft30cm(1:5));
            else
                licksTmp = trials{i}.lickLfpInd(indAft30cm);
            end
            beh.medDistFirst5Licks(i) = median(trials{i}.xMM(licksTmp));
            beh.meanDistFirst5Licks(i) = mean(trials{i}.xMM(licksTmp));
            %% 
            
            %% number of runs per trial
            speed = trials{i}.speed;
            speedSM = smooth(speed,smoothSpan);
            minSpeed1 = 100; % mm/s
            indSpeed = speedSM >= minSpeed1;
            [continuousRun,stopRun] = numOfConsecutiveOnes(indSpeed);
            beh.totStopTime(i) = sum(stopRun)/sampleFq;
            beh.numRun(i) = sum(continuousRun > 250);
            [maxRunSeg,maxRunInd] = max(continuousRun);
            if(indSpeed(1) == 0)
                if(maxRunInd(1) < 2)
                    indMaxRunStart = stopRun(1)+1;
                else
                    indMaxRunStart = sum(stopRun(1:maxRunInd(1))) + ...
                        sum(continuousRun(1:maxRunInd(1)-1))+1;
                end
            else
                if(maxRunInd(1) < 2)
                    indMaxRunStart = 1;
                else
                    indMaxRunStart = sum(stopRun(1:maxRunInd(1)-1)) + ...
                        sum(continuousRun(1:maxRunInd(1)-1))+1;
                end
            end
            indMaxRunEnd = indMaxRunStart + maxRunSeg - 1;
%             if(indMaxRunEnd>length(trials{i}.xMM))
%                 a = 1;
%             end
            beh.maxRunSegStartDist(i) = trials{i}.xMM(indMaxRunStart);
            beh.maxRunSegEndDist(i) = trials{i}.xMM(indMaxRunEnd);
            %%
            
            if(beh.lenTrials(i) < sampleFq*100)
                beh.indGoodLap = [beh.indGoodLap i]; 
            else
                disp(['Trial no. ' num2str(i) ' is longer than 100 s']);
            end
            
            if(isfield(trials{i},'stimOnLfpInd'))      
                if(~isempty(trials{i}.stimOnLfpInd))
                    noStimTr = noStimTr + 1; %% added by Yingxue on 12/14/2020
                    beh.pulseMethod(i) = trials{i}.stimPulseMethod;
                    beh.indStimLap(i) = 1;
                    beh.stimPulseLfpInd{i} = trials{i}.stimPulseLfpInd;
                    beh.stimPulseWidth{i} = trials{i}.stimPulseWidth;
                    beh.stimDiode{i} = trials{i}.stimDiode;
                    beh.stimDiodeCurr{i} = trials{i}.stimDiodeCurr;
                    beh.stimPulseLoc{i} = trials{i}.stimPulseLoc;
                end
            end
            %%% added by Yingxue on 12/14/2020
            if(noStimTr == 0) % control trials before the first stimulation trial
                if(beh.lenTrials(i) < sampleFq*100)
                    beh.indTrCtrl = [beh.indTrCtrl i]; 
                else
                    disp(['Trial no. ' num2str(i) ' is longer than 100 s']);
                end
            end
            %%% added by Yingxue on 12/14/2020
        end
    end
    
    %%% added by Yingxue on 12/14/2020
    medNumSample = median(beh.lenTrials(beh.indTrCtrl));
    figure
    for i = beh.indTrCtrl
%         if(i == 70 || i == 71)
%             a = 1;
%         end
        hold off;
        plot(trials{i}.speed);
        hold on;
        h = plot(trials{i}.lickLfpInd,600*ones(1,length(trials{i}.lickLfpInd)),'m.');
        set(h,'MarkerSize',8);
        if(~isempty(trials{i}.pumpLfpInd))
            plot([trials{i}.pumpLfpInd(1), trials{i}.pumpLfpInd(1)],[0 600],'g-');
        end
        if(~isempty(trials{i}.pumpLfpInd) && beh.numLickBef100cm(i) == 0 && ...
                beh.lenTrials(i) < medNumSample*3 && ...
                (beh.numRun(i) <= 2 || (beh.maxRunSegEndDist(i)-beh.maxRunSegStartDist(i))>1500) && ... %  
                mean(trials{i}.speed(1:125)) < minSpeed1)
            beh.indGoodTrCtrl = [beh.indGoodTrCtrl i];
            title(['trial no. ' num2str(i) ' Good Tr, MaxDist = ' ...
                num2str(beh.maxRunSegEndDist(i)-beh.maxRunSegStartDist(i)) ', numLickBef100cm = ' ...
                num2str(beh.numLickBef100cm(i)) ', NumRun = ' ...
                num2str(beh.numRun(i))]);
        elseif(isempty(trials{i}.pumpLfpInd) || beh.numLickBef100cm(i) > 2 || ...
                (beh.lenTrials(i) > medNumSample*3 && ...
                beh.numRun(i) > 2) || (beh.maxRunSegEndDist(i)-beh.maxRunSegStartDist(i))<1200 || ...
                mean(trials{i}.speed(1:250)) > minSpeed1)
            beh.indBadTrCtrl = [beh.indBadTrCtrl i];
            title(['trial no. ' num2str(i) ' Bad Tr, MaxDist = ' ...
                num2str(beh.maxRunSegEndDist(i)-beh.maxRunSegStartDist(i)) ', numLickBef100cm = ' ...
                num2str(beh.numLickBef100cm(i)) ', NumRun = ' ...
                num2str(beh.numRun(i))]);
        else
            title(['trial no. ' num2str(i) ', MaxDist = ' ...
                num2str(beh.maxRunSegEndDist(i)-beh.maxRunSegStartDist(i)) ', numLickBef100cm = ' ...
                num2str(beh.numLickBef100cm(i)) ', NumRun = ' ...
                num2str(beh.numRun(i))]);
        end
%         pause;
    end
    disp(['Number of good trials in the control: ' num2str(length(beh.indGoodTrCtrl))]);
    close;
   
    fullPath = [path fileNameMinLen];
    save(fullPath, 'rec', 'beh');
    return;
end

function [data,data1] = numOfConsecutiveOnes(arr)
    data = [];
    data1 = [];
    s = sprintf('%d', arr);
    %Reading the consequences of 1's from the string by using 0's as delimiters
    t1=textscan(s,'%s','delimiter','0','multipleDelimsAsOne',1);
    % Converting cell array of cell into a single cell array
    d = t1{:};
    % Computing the length of each run by going through the array and assigning
    % it into 
    for k = 1:length(d)
          data(k) = length(d{k});
    end
    
    t2=textscan(s,'%s','delimiter','1','multipleDelimsAsOne',1);
    % Converting cell array of cell into a single cell array
    f = t2{:};
    % Computing the length of each run by going through the array and assigning
    % it into 
    for k = 1:length(f)
          data1(k) = length(f{k});
    end
end
