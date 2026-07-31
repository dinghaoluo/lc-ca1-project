function thetaCycleStruct = ThetaCycle1(theta, indLapList, minSamplePerCycle, thetaPhaseJump)
% this function is to extract each theta cycle in each valid trial
%
% Inputs:
% thetaPeak:            theta phase structure including index of peak theta 
%                       phase within each cycle of each valid trial
% indLapList:           the array containing the index numbers of the valid
%                       trials (Caution: since theta strcuture records the 
%                       theta traces from all the valid trials in 
%                       the task, indLapList here includes the indices of 
%                       all to be selected cells in the theta structure, 
%                       rather than referring back to the 
%                       original trial no.)
% minSamplePerCycle:    the minimum number of samples a valid theta cycle
%                       should contain
% thetaPhaseJump:       the min sudden change in theta phase which indicates
%                       the start of a new cycle
% 
% Outputs: thetaCycleStruct which includes following fields:
% indLapList:           
% indCycleStart:        indices of the starting sample number of each theta
%                       cycle (indCycleStart{i} contains all the indices 
%                       from trial i)
% numCyclesPerTrial:    number of cycles per trial
% meanCycleLenPerTrial: mean cycle length per trial
% stdCycleLenPerTrial:  std cycle length per trial
% cycleLenPerCycle:     theta cycle length per cycle
% meanCycleLenPerCycle: mean theta cycle length per cycle
% stdCycleLenPerCycle:  std theta cycle length per cycle
% numCyclesMin:         the minimum number of cycles over all valid trials
% meanCycleLen:         mean cycle length over all valid trials
% stdCycleLen:          the std of cycle length over all valid trials

numTrials = length(indLapList);
if(numTrials == 0)
    thetaCycleStruct = [];
    return;
end
    
thetaCycleStruct = struct('indLapList', indLapList,...
                          'indCycleStart',[],...
                          'numCyclesPerTrial',zeros(1,numTrials),...
                          'cycleLenPerTrial',[],...
                          'meanCycleLenPerTrial',zeros(1,numTrials),...
                          'stdCycleLenPerTrial',zeros(1,numTrials),...
                          'cycleLenPerCycle',[],...
                          'meanCycleLenPerCycle',[],...
                          'stdCycleLenPerCycle',[],...
                          'numCyclesMin',0,...
                          'meanCycleLen',0,...
                          'stdCycleLen',0);
                      
% extract the start point of each theta cycle 
lenPerTrial = zeros(1,numTrials);
allCycleLen = [];
for i = 1:numTrials
    if(i ==3)
        a = 1;
    end
    lenPerTrial(i) = length(theta{indLapList(i)});
    
    % find the points where the phase jumps from -180 to 180
    diffPhase = diff(theta{indLapList(i)});
    indTmp = find(diffPhase < thetaPhaseJump);
    
    if(isempty(indTmp))
       continue; 
    end
    
    % find points where the phase crosses 0
    indZeros = [];
    if(indTmp(1)) > 0 % the first cycle
        [ind, ind0, s0] = detectCrossing_yx(theta{indLapList(i)}...
                                        (1:indTmp(1))',1:indTmp(1),0);
        if(length(s0) > 1)
            indUniq0 = find(abs(s0) == min(abs(s0)),1);
            indZeros = [indZeros, ind0(indUniq0)];
        else
            indZeros = [indZeros, ind0];
        end
    end
    for j = 1:length(indTmp)-1 % cycles in between
        if(j == 236)
            a = 1;
        end
        [ind, ind0, s0] = detectCrossing_yx(theta{indLapList(i)}...
                    (indTmp(j)+1:indTmp(j+1)),indTmp(j)+1:indTmp(j+1),0);
        if(length(s0) > 1)
            indUniq0 = find(abs(s0) == min(abs(s0)),1);
            indZeros = [indZeros, ind0(indUniq0)];
        else
            indZeros = [indZeros, ind0];
        end
    end
    lenThetaPhase = length(theta{indLapList(i)});
    if(indTmp(end) < lenThetaPhase) % the last cycle
        [ind, ind0, s0] = detectCrossing_yx(theta{indLapList(i)}...
                    (indTmp(end)+1:end),indTmp(end)+1:lenThetaPhase,0);
        if(length(s0) > 1)
            indUniq0 = find(abs(s0) == min(abs(s0)),1);
            indZeros = [indZeros, ind0(indUniq0)];
        else
            indZeros = [indZeros, ind0];
        end
    end
    
    % find the length of each theta cycle
    diffInd = diff(indZeros);   
    indChk = find(diffInd > minSamplePerCycle); 
                    % guarantee that every theta cycle is >
                    % minSampleThetaCycle, otherwise combine the cycle with
                    % the following cycle
    % searching for long theta cycles and segment them into several cycles
    diffIndChk = diff(indChk);                
    indNonCont = find(diffIndChk > 2);
    
    for nInd = 1:length(indNonCont)
        indInsert = [];
        startInd = indChk(indNonCont(nInd))+1;
        while(startInd < indChk(indNonCont(nInd)+1))
            cumTmp = cumsum(diffInd(startInd:indChk(indNonCont(nInd)+1)));
            indTmp = find(cumTmp > minSamplePerCycle,1);
            indInsert = [indInsert startInd+indTmp-1];
            startInd = startInd + indTmp;
        end
        indChk = [indChk indInsert];
    end
    indChk = sort(indChk);
    thetaCycleStruct.indCycleStart{i} = indZeros(indChk); 
        % the starting index of each theta cycle within each trial
    thetaCycleStruct.numCyclesPerTrial(i) = ...
        length(thetaCycleStruct.indCycleStart{i}); 
        % number of theta cycles in a single trial
    cycleLenPerTrial = diffInd(indChk);
    thetaCycleStruct.cycleLenPerTrial{i} = cycleLenPerTrial;
    allCycleLen = [allCycleLen, cycleLenPerTrial];
    thetaCycleStruct.meanCycleLenPerTrial(i) = mean(cycleLenPerTrial); 
        % mean length of a theta cycle over a single trial
    thetaCycleStruct.stdCycleLenPerTrial(i) = std(cycleLenPerTrial); 
        % std length of a theta cycle over a single trial
end
thetaCycleStruct.numCyclesMin = min(thetaCycleStruct.numCyclesPerTrial...
    (thetaCycleStruct.numCyclesPerTrial > 0)); 
    % the minimum number of cycles over all the valid trials
thetaCycleStruct.meanCycleLen = mean(allCycleLen); 
    % the mean cycle length over all the valid trials
thetaCycleStruct.stdCycleLen = std(allCycleLen); 
    % the std of cycle length over all the valid trials

% calculate the averaged theta cycle length per cycle
thetaCycleStruct.cycleLenPerCycle = ...
    zeros(numTrials,thetaCycleStruct.numCyclesMin);
for i = 1:numTrials
    cycleLenPerCycleTmp = thetaCycleStruct.cycleLenPerTrial{i};
    if(isempty(cycleLenPerCycleTmp))
        continue;
    end
    thetaCycleStruct.cycleLenPerCycle(i,:) = ...
        cycleLenPerCycleTmp(1:thetaCycleStruct.numCyclesMin);
end
thetaCycleStruct.meanCycleLenPerCycle = ...
    mean(thetaCycleStruct.cycleLenPerCycle);
thetaCycleStruct.stdCycleLenPerCycle = ...
    std(thetaCycleStruct.cycleLenPerCycle);
