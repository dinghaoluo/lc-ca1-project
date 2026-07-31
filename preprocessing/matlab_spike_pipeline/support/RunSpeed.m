function RunSpeed(path,fileName,figState)
% running speed (speed > certain threshold)
%
% e.g.: RunSpeedVR('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1')
%
% by Yingxue, 08/24/2017

    %%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin == 2
        figState = 0;
    elseif nargin > 3
        disp('Too many arguments');        
        return;
    end
    
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(isempty(indexFileName))
        fileNameInfo = [fileName '_Info.mat'];
        fileNameSpeed = [fileName '_runSpeed.mat'];
        fileName = [fileName '.mat'];
    else
        fileNameInfo = [fileName(1:indexFileName(end)-1) '_Info.mat'];
        fileNameSpeed = [fileName(1:indexFileName(end)-1) '_runSpeed.mat'];
    end 
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath,'trials','lapList');
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
    
    GlobalConst;
    
    speed = getRecField(trials,'speed',1:length(lapList));
    lick = getRecField(trials,'lickLfpInd',1:length(lapList));
    xmm = getRecField(trials,'xMM',1:length(lapList));
    
    param.minSpeed = minSpeed;
    param.sampleFq = sampleFq;
    param.minSpeedSeg = 100;
    param.spectroWin = floor(1.5*sampleFq);
    param.spectroOverlap = floor(1.2*sampleFq);
    param.sampleFq = sampleFq;
    
    % calculate the running speed and lick events for the whole session
    disp(['Calculate the running speed and count lick events for the whole'...
          ' recording session']);
    [runSpeed,lickEvents] = ...
        speedandLick(speed,lick,beh.indGoodLap,beh.numTrials,param);
    
    disp('Find the longest continuous run in each trial');
    runSegments = runSegment(speed,xmm,beh.indGoodLap,beh.numTrials,param);
    
    disp('Calculate the power spectrogram of the running speed');
    speedSpectro = speedSpectrogram(speed,beh.indGoodLap,beh.numTrials,param);
                
    disp(['Calculate the running speed and count lick events for each'...
          ' subsession']);
    runSpeedSess = cell(length(mazeSess),1);
    lickEventsSess = cell(length(mazeSess),1);
    if(length(mazeSess)>1)
        for i = 1:length(mazeSess) 
            fprintf('\nSession %d\n',i);
            indLaps = find(beh.mazeSess == mazeSess(i)); 
            %%% calculate mean firing rate for each neuron over specified trials
            [runSpeedSess{i},lickEventsSess{i}] = ...
                speedandLick(speed,lick,indLaps,beh.numTrials,param);
        end
    end 
            
    fullPath = [path fileNameSpeed];
    save(fullPath, 'runSpeed','runSegments','speedSpectro','lickEvents',...
                'runSpeedSess','lickEventsSess','param','-v7.3');
    
    if(figState == 1)
        barPlot(1,meanSpeed,stdSpeed,...
                'Session','Mean speed (cm/sec)','Mean speed per session');
        barPlot(1,meanRest,stdRest,...
                'Session','Mean low speed time (sec)',...
                'Mean low speed time per session');
        barPlot(1,meanLick,stdLick,...
                'Session','Mean number of licks',...
                'Mean number of licks per session');
    end
end

function [runSpeed, lickEvents] = ...
                speedandLick(speed,lick,indLaps,totNumTrials,param)
    runSpeed = struct('lapList', indLaps',... trials involved 
                      'meanSpeedPerTr',zeros(1,length(indLaps)),... % mean running speed per trial
                      'timeLowSpeedPerTr',zeros(1,length(indLaps)),... % the amont of time per trial where the speed is lower than threshold 
                      'meanSpeed',0,... % mean running speed
                      'stdSpeed', 0,...   % std running speed
                      'meanTimeLowSpeed',0,... % mean of the amount of time with low speed
                      'stdTimeLowSpeed',0); % std of the amount of time with low speed

    lickEvents = struct('lapList',indLaps',... trial involved
                        'numLickPerTr',zeros(1,length(indLaps)),... % number of lick events per trial
                        'meanNumLick',0,... % mean number of licks
                        'stdNumLick',0);    % std number of licks
                        
    tr = 0;              
    for i = 1:totNumTrials
        if(sum(indLaps == i)>0)
            ind = speed{i} > param.minSpeed;
            tr = tr+1;
            runSpeed.meanSpeedPerTr(tr) = mean(speed{i}(ind));
            runSpeed.timeLowSpeedPerTr(tr) = sum(ind == 0)/param.sampleFq;
            
            lickEvents.numLickPerTr(tr) = length(lick{i});         
        end
    end
    
    runSpeed.meanSpeed = mean(runSpeed.meanSpeedPerTr);
    runSpeed.stdSpeed = std(runSpeed.meanSpeedPerTr);
    
    runSpeed.meanTimeLowSpeed = mean(runSpeed.timeLowSpeedPerTr);
    runSpeed.stdTimeLowSpeed = std(runSpeed.timeLowSpeedPerTr);
    
    lickEvents.meanNumLick = mean(lickEvents.numLickPerTr);
    lickEvents.stdNumLick = std(lickEvents.numLickPerTr);
    
end

function runSegments = runSegment(speed,xmm,indLaps,totNumTrials,param)
    runSegments = struct('lapList', indLaps',... trials involved 
                      'indStartRunSegments',[],... % indices of the start point of each running segment
                      'indStopRunSegments',[],... % indices of the stop point of each running segment
                      'numRunSegments',zeros(1,length(indLaps)),... % number of running segments per trial
                      'distRunSegments',[],... % distance of each running segment
                      'speedRunSegments',[],... % speed of each running segment
                      'meanSpRunSegments',[],... % average speed of each running segment
                      'dist1stRunSegment',zeros(1,length(indLaps)),... % distance of the first run segment
                      'speed1stRunSegment',zeros(1,length(indLaps)),... % mean speed of the first run segment
                      'distMaxRunSegment',zeros(1,length(indLaps)),... % distance of the longest run segment
                      'speedMaxRunSegment',zeros(1,length(indLaps))); % mean speed of the longest run segment
      
    tr = 0;
    for i = 1:totNumTrials
        if(sum(indLaps == i)>0)
%             disp(num2str(tr));
            tr = tr+1;
            speedAvg = movmean(speed{i},10);
            ind = speedAvg>param.minSpeedSeg;
            ii1 = strfind([0 ind' 0],[0 1]);
            ii2 = strfind([0 ind' 0],[1 0]);
            ii = (ii2-ii1+1) >= 4;
            ii1 = ii1(ii);
            ii2 = ii2(ii);
            if(ii2(end) > length(speedAvg))
                ii2(end) = length(speedAvg);
            end
    
            runSegments.indStartRunSegments{tr} = ii1;
            runSegments.indStopRunSegments{tr} = ii2;
            runSegments.numRunSegments(tr) = length(ii1);
            runSegments.distRunSegments{tr} = ...
                xmm{i}(ii2)-xmm{i}(ii1);
            runSegments.dist1stRunSegment(tr) = ...
                runSegments.distRunSegments{tr}(1);
            runSegments.distMaxRunSegment(tr) = ...
                max(runSegments.distRunSegments{tr});
            
            runSegments.speedRunSegments{tr} = ...
                arrayfun(@(x,y) speedAvg(x:y),ii1,ii2,'un',0);           
            for j = 1:length(runSegments.speedRunSegments{tr})
                runSegments.meanSpRunSegments{tr}(j) = ...
                    mean(runSegments.speedRunSegments{tr}{j});
                if(j == 1)
                    runSegments.speed1stRunSegment(tr) = ...
                        runSegments.meanSpRunSegments{tr}(j);
                end
            end
            runSegments.speedMaxRunSegment(tr) = ...
                        max(runSegments.meanSpRunSegments{tr});
%             out = arrayfun(@(x,y) speedAvg(x:y),ii1(ii),ii2(ii),'un',0);
%             celldisp(out)
        else
            continue;
        end
    end
end
