function plotField(path,fileName,spaceBin,methodTheta,onlyRun,sortMethod)
% plot field according to licking time within a trial
% sortMethod: 1: based on the first lick position
%             2: based on the running distance of the first segment
%             3: based on the running speed of the first segment

    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    fileNameConv = [fileName '_convSpikesDist' num2str(spaceBin)...
                        'mm_Run' num2str(onlyRun) '.mat'];     
    fileNameSpInfo = [fileName '_SpInfo_Run' num2str(onlyRun) '.mat'];
    fileNameInfo = [fileName '_Info.mat'];
    fileNameRun = [fileName '_ext.mat'];
    fileNameFR = [fileName '_FR_Run' num2str(onlyRun) '.mat'];
    fileNameDepth = [fileName '_Depth.mat'];
    fileNameSpeed = [fileName '_runSpeed.mat'];
    fileNameRec = [fileName '.mat'];
    
    fullPath = [path fileNameConv];
    if(exist(fullPath) == 0)
        disp(['The firing profile file does not exist. Please call ',...
                'function "ConvSpikeTrainDistParVR" first.']);
        return;
    end
    load(fullPath,'filteredSpikeArrayNormT','paramC');
    
    fullPath = [path fileNameSpInfo];
    if(exist(fullPath) == 0)
        disp(['The spatial information file does not exist. Please call ',...
                'function "GetFRMapInfo" first.']);
        return;
    end
    load(fullPath,'spatialInfoSess');  
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath,'file') == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
    numSamples = zeros(1,length(unique(beh.trackLen)));
    for i = 1:length(numSamples)
        numSamples(i) = length(paramC.spaceSteps{i});
    end
    maxNumSamples = max(numSamples);
    
    fullPath = [path fileNameRun];
    if(exist(fullPath) == 0)
        disp(['The spikes during run file does not exist. Please call ',...
                'function "SpikeDuringRun" first.']);
        return;
    end
    load(fullPath);
    
    fullPath = [path fileNameFR];
    if(exist(fullPath) == 0)
        disp(['The neuron depth file does not exist. Please call ',...
                'function "GetNeuRelativeDepth" first.']);
        return;
    end
    load(fullPath,'mFRStruct','mFRStructSess');
    
    fullPath = [path fileNameDepth];
    if(exist(fullPath) == 0)
        disp(['The mean firing rate file does not exist. Please call ',...
                'function "MeanFiringRateVR" first.']);
        return;
    end
    load(fullPath);
    
    fullPath = [path fileNameSpeed];
    if(exist(fullPath) == 0)
        disp(['The run speed file does not exist. Please call ',...
                'function "RunSpeed" first.']);
        return;
    end
    load(fullPath,'runSegments','speedSpectro');
    
    fullPath = [path fileNameRec];
    if(exist(fullPath) == 0)
        disp('The recording file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath,'trials','cluList');
    
    indNeurons = 1:rec.numNeurons;
    indGoodLap = beh.indGoodLap;
    %% sort trials 
    mazeSessTr = beh.mazeSess(indGoodLap);
    if(sortMethod == 1) % based on the first lick position
        % find first lick position
        firstLickPos = zeros(1,length(mazeSessTr));
        tr = 0;
        for i = indGoodLap
            tr = tr+1;
            lickDist = trials{i}.xMM(trials{i}.lickLfpInd);
            ind = find(lickDist > 300, 1, 'first'); 
            firstLickPos(tr) = mean(lickDist(ind));
        end
        
%         ind = [1:119]';
%         [orderedParam{1}, indGoodLapOrdered{1}] = ...
%                 sort(firstLickPos(ind));
%         indGoodLapOrdered{1} = ind(indGoodLapOrdered{1});

        % order trials by where the lick occurs
        for i = mazeSess'
            ind = find(mazeSessTr == i);
            indLaps = indGoodLap(ind);
            [orderedParam{i}, indGoodLapOrdered{i}] = ...
                sort(firstLickPos(ind));
            indGoodLapOrdered{i} = indLaps(indGoodLapOrdered{i});
        end
    elseif(sortMethod == 2) % based on the running distance of the longest segment
        for i = mazeSess'
            ind = find(mazeSessTr == i);
            [orderedParam{i}, indGoodLapOrdered{i}] = ...
                sort(runSegments.distMaxRunSegment(ind));
            indGoodLapOrdered{i} = ind(indGoodLapOrdered{i});
        end
    elseif(sortMethod == 3) % based on the running speed of the longest segment
        for i = mazeSess'
            ind = find(mazeSessTr == i);
            [orderedParam{i}, indGoodLapOrdered{i}] = ...
                sort(runSegments.speedMaxRunSegment(ind));
            indGoodLapOrdered{i} = ind(indGoodLapOrdered{i});
        end
    elseif(sortMethod == 4)
        % distance ran during the start cue
        distDuringStartCue = [0.0 0.08 0.6 0.44 0.28 0.16 0.88 0.24 6.72 ...
            0.0 0.72 0.44 9.16 0.28 0.68 0.28 0.12 1.56 0.44 0.68 1.0 0.0 ...
            4.44 0.48 0.4 0.24 0.28 0.16 0.92 1.2 0.0 0.2 0.44 0.28 0.08 ...
            0.4 0.72 0.08 0.52 0.0 0.2 0.0 0.32 0.0 5.32 0.84 0.48 0.68 ...
            1.04 0.88 0.2 0.84 0.44 0.88 0.16 0.44 1.12 0.48 0.16 0.68 ...
            6.44 0.6 0.48 0.32 1.08 0.0 0.28 0.28 0.24 0.48 0.84 1.04 ...
            0.6 0.44 8.68 4.24 0.08 0.68 0.56 0.44 0.96 0.2 0.4 1.56 ...
            0.12 0.84 0.68 1.36 0.48 1.2 0.2 3.8 9.04 0.8 0.88 3.72 0.56 ...
            1.68 1.0 0.96 1.4 1.88 2.6 0.72 0.2 0.32 1.04 1.0 0.52 2.08 ...
            0.84 1.0 0.4 2.08 0.6 0.36 0.72 1.2 0.0 11.2 1.28 5.32 0.68 ...
            0.12 1.8 1.08 0.64 0.24 0.92 0.92 0.16 2.24 2.0 0.64 1.28 0.16 ...
            0.68 0.28 0.16 10.52 2.88 0.8 0.44 0.32 0.2 4.2 0.52 0.08 0.6 ...
            0.68 2.08 0.0 0.2 0.32 0.92 0.28 0.68 1.0 0.64 1.0 0.96 0.68 ...
            1.92 0.36 3.96 0.96 0.2 0.2 0.48 0.36 0.04 0.48 0.44 0.64 1.2 ...
            3.92 0.68 1.12 0.56 0.36 0.16 0.64 1.04 0.28 1.24 0.92 0.6 ...
            0.44 0.32 0.52 5.8 0.88 0.16 0.24 0.28 1.52 0.68 0.56 0.6 0.12 ...
            0.28 0.56 1.04 0.04 1.24 1.36 2.28 0.68 1.2 1.4 0.68 1.0 1.0 ...
            0.72 1.0 0.32 0.24 2.6 0.28 0.52 1.2 0.56 0.44 0.4 0.6 0.28 ...
            1.68 0.56 0.36 1.12 8.48 1.48 1.04 0.08 0.52 0.68 0.72 1.12 ...
            0.84 1.2 0.32 0.96 0.56 0.24 0.48 0.6 0.28 0.44 0.16 0.72 0.32 ...
            0.08 0.2 0.56 0.24 0.28];
        for i = mazeSess'
            ind = find(mazeSessTr == i);
            [orderedParam{i}, indGoodLapOrdered{i}] = ...
                sort(distDuringStartCue(ind));
            indGoodLapOrdered{i} = ind(indGoodLapOrdered{i});
        end
    elseif(sortMethod == 5)
        % Time between start cue offset and speed exceeding 10cm/s (in seconds)
        timeStartCueOffToSpeedUp = [1.519 1.191 5.125 3.632 1.039 13.503 6.059 4.673 0.017 11.545 14.885 7.43 0.003 0.325 1.849 0.642 0.392 4.203 2.564 1.845 2.34 2.983 0.93 3.345 0.377 1.081 1.043 5.004 0.323 4.969 0.443 6.943 1.729 2.924 6.643 2.962 0.395 2.767 0.563 0.763 2.018 1.136 1.343 2.096 0.092 0.824 1.457 2.119 1.941 1.163 5.294 2.418 2.383 7.317 2.183 1.387 1.277 1.863 0.8 1.555 0.003 0.083 1.123 7.083 4.697 1.133 1.26 4.323 4.888 2.662 1.394 0.635 3.347 3.403 0.003 0.007 4.687 0.523 1.294 3.867 1.41 3.028 0.956 0.423 0.228 1.021 1.954 1.168 0.943 1.268 2.144 3.049 0.001 0.048 0.622 0.083 0.877 0.019 0.326 1.223 0.692 0.759 0.383 0.617 2.763 5.828 1.788 0.103 1.363 3.163 1.129 0.863 0.367 2.203 0.661 1.718 1.123 0.143 0.347 0.003 0.126 0.006 5.863 2.607 0.863 0.783 1.274 0.611 3.704 0.945 1.119 1.252 1.329 1.421 0.736 0.643 2.176 0.621 0.323 0.003 0.016 2.304 0.582 1.143 0.636 0.423 0.923 1.252 2.442 1.343 0.028 1.452 0.707 0.622 0.223 2.799 1.103 0.816 0.99 1.359 1.123 0.863 0.832 1.175 0.015 0.015 2.844 0.28 2.123 0.603 0.723 1.683 1.396 1.398 0.521 0.003 1.095 1.123 2.769 1.129 0.331 0.814 0.416 1.144 2.283 0.683 0.182 0.756 1.443 1.151 0.003 0.627 0.743 0.243 0.141 0.024 1.001 0.411 2.243 1.093 0.063 0.343 0.356 0.968 0.503 0.736 0.203 3.031 0.117 0.883 0.583 0.076 0.001 0.61 1.483 0.503 0.708 0.019 2.853 0.361 0.083 0.984 0.141 0.126 1.403 3.027 2.165 2.323 0.637 0.668 0.003 0.243 0.419 0.762 0.686 0.042 0.352 0.569 0.457 1.009 0.715 0.882 0.661 0.961 1.668 0.008 1.883 0.962 1.863 0.723 1.78 1.751 1.777 0.669 1.205 2.96];
        for i = mazeSess'
            ind = find(mazeSessTr == i);
            [orderedParam{i}, indGoodLapOrdered{i}] = ...
                sort(timeStartCueOffToSpeedUp(ind));
            indGoodLapOrdered{i} = ind(indGoodLapOrdered{i});
        end
    elseif(sortMethod == 6)
        % # Counts < 10cm/s between period of running onset and running offset after the first lick
        countLT10cm = [111 0 68 15 0 0 0 0 4 0 0 0 22 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 0 0 0 0 0 0 0 0 0 0 0 0 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 8 0 0 0 0 0 0 0 0 0 0 0 126 35 49 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 87 2 0 4 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 4 0 3 0 0 7 0 0 0 0 2 0 0 0 0 2 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 1 0 0 0 0 0 8 0 0 0 0 0 43 0 2 0 0 11 0 0 0 0 0 0 0 2 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 29 0 2 0 9];
        for i = mazeSess'
            ind = find(mazeSessTr == i);
            [orderedParam{i}, indGoodLapOrdered{i}] = ...
                sort(countLT10cm(ind));
            indGoodLapOrdered{i} = ind(indGoodLapOrdered{i});
        end
    elseif(sortMethod == 7) % speed power spectrogram percentage of samples with high power, between 3-5Hz 
        for i = mazeSess'
            ind = find(mazeSessTr == i);
            [orderedParam{i}, indGoodLapOrdered{i}] = ...
                sort(speedSpectro.percHighPower3_5(ind));
            indGoodLapOrdered{i} = ind(indGoodLapOrdered{i});
        end
    elseif(sortMethod == 8) % speed power spectrogram mean power between 6-10Hz
        for i = mazeSess'
            ind = find(mazeSessTr == i);
            [orderedParam{i}, indGoodLapOrdered{i}] = ...
                sort(speedSpectro.percHighPower6_10(ind));
            indGoodLapOrdered{i} = ind(indGoodLapOrdered{i});
        end
    elseif(sortMethod == 9) %* speed power spectrogram mean power between 3-5Hz
        for i = mazeSess'
            ind = find(mazeSessTr == i);
            [orderedParam{i}, indGoodLapOrdered{i}] = ...
                sort(speedSpectro.mmeanAmpSpectro3_5(ind));
            indGoodLapOrdered{i} = ind(indGoodLapOrdered{i});
        end
    elseif(sortMethod == 10)
        % Time between 1st lick after 20cm and reward delivery:
        time1stLickToRew = [15.24 2.682 7.383 5.989 3.105 1.175 2.195 1.072 0.203 4.146 3.698 1.7 2.283 0.47 0.922 0.792 1.403 2.129 1.526 1.562 17.956 4.527 2.45 1.758 0.016 1.317 3.762 2.088 2.231 2.847 2.386 2.056 2.031 1.81 0.739 0.504 0.657 2.452 1.351 3.853 0.24 2.839 0.726 2.457 2.217 0.567 1.917 1.37 0.721 1.155 1.506 1.396 1.238 2.142 1.656 0.016 1.469 0.905 2.075 0.208 1.419 2.168 0.032 4.116 1.332 1.682 0.764 2.35 1.494 1.454 1.298 1.129 1.525 1.212 1.482 0.016 2.139 1.63 1.787 2.205 1.836 1.528 1.128 2.779 0.73 2.172 1.659 0.976 1.161 2.052 1.321 2.511 1.685 1.979 1.81 0.578 1.828 0.89 2.09 1.486 1.155 1.066 0.73 0.689 2.335 2.416 2.538 1.061 0.842 2.047 1.194 1.765 0.847 1.253 1.518 0.962 0.978 1.006 1.298 13.066 9.053 3.029 2.12 1.49 1.303 1.184 1.234 2.085 0.542 2.064 0.553 0.488 1.455 1.018 0.636 1.324 0.663 3.137 1.773 1.403 0.898 1.148 0.763 0.767 0.639 2.937 0.628 2.693 0.703 0.651 0.378 2.183 0.016 0.455 0.849 1.733 0.915 0.686 0.732 0.751 1.299 0.52 17.818 4.134 1.54 0.016 3.039 0.951 1.555 1.149 0.849 0.652 1.37 11.302 3.8 2.239 2.223 1.952 3.94 0.966 1.928 0.016 1.531 1.453 2.491 1.132 1.454 0.461 1.672 1.716 1.21 0.433 2.225 1.01 0.79 1.351 1.719 0.886 0.809 0.886 0.629 1.942 0.505 1.608 0.894 0.674 0.812 5.051 0.297 0.918 0.454 2.906 0.752 1.296 0.781 2.544 0.639 1.703 4.404 0.058 1.099 1.638 0.016 1.086 1.185 1.04 0.868 0.949 0.49 1.141 0.843 0.452 1.293 0.683 1.699 0.654 0.65 1.859 1.455 1.094 0.763 1.045 0.597 1.103 1.648 0.956 1.195 0.415 0.723 1.934 1.387 2.556 1.751 0.803 1.131 2.152 ];
        for i = mazeSess'
            ind = find(mazeSessTr == i);
            [orderedParam{i}, indGoodLapOrdered{i}] = ...
                sort(time1stLickToRew(ind));
            indGoodLapOrdered{i} = ind(indGoodLapOrdered{i});
        end
    elseif(sortMethod == 11) %*
        % mean time between licks that occur after 1st 20cm of trial
        meanLickInterval = [0.7 0.207 0.491 0.27 0.203 0.163 0.19 0.173 0.126 0.225 0.238 0.193 0.224 0.138 0.162 0.163 0.183 0.181 0.174 0.179 0.419 0.219 0.232 0.189 0.134 0.163 0.195 0.192 0.183 0.181 0.278 0.23 0.284 0.191 0.201 0.157 0.141 0.178 0.254 0 0.187 0.182 0.156 0.178 0.207 0.133 0.173 0.172 0.141 0.157 0.17 0.17 0.163 0.19 0.189 0.163 0.158 0.164 0.196 0.127 0.233 0.192 0.123 0.179 0.173 0.216 0.154 0.24 0.16 0.185 0.167 0.186 0.158 0.167 0.207 0.124 0.154 0.215 0.159 0.231 0.161 0.158 0.176 0.195 0.192 0.194 0.183 0.163 0.211 0.165 0.179 0.206 0.19 0.2 0.216 0.158 0.165 0.217 0.189 0.157 0.19 0.152 0.168 0.148 0.173 0.2 0.182 0.151 0.15 0.164 0.158 0.163 0.157 0.15 0.185 0.145 0.158 0.144 0.172 1.048 0.593 0.2 0.199 0.204 0.177 0.205 0.179 0.168 0.142 0.155 0.154 0.181 0.147 0.178 0.141 0.148 0.155 0.195 0.275 0.165 0.15 0.156 0.151 0.152 0.145 0.171 0.144 0.203 0.149 0.159 0.138 0.165 0.124 0.145 0.145 0.162 0.14 0.155 0.146 0.144 0.147 0.144 0.863 0.203 0.19 0.118 0.184 0.154 0.192 0.154 0.143 0.156 0.153 0.239 0.221 0.207 0.177 0.191 0.187 0.159 0.149 0.134 0.158 0.156 0.182 0.162 0.159 0.134 0.144 0.156 0.15 0.14 0.155 0.166 0.16 0.16 0.195 0.146 0.146 0.155 0.151 0.157 0.141 0.169 0.159 0.137 0.147 0.25 0.126 0.141 0.147 0.18 0.131 0.15 0.135 0.17 0.155 0.157 0.199 0.118 0.156 0.184 0.12 0.156 0.15 0.165 0.148 0.143 0.147 0.155 0.145 0.141 0.144 0.158 0.159 0.144 0.144 0.179 0.178 0.16 0.153 0.152 0.133 0.151 0.175 0.142 0.155 0.128 0.132 0.156 0.159 0.164 0.157 0.139 0.147 0.161];
        for i = mazeSess'
            ind = find(mazeSessTr == i);
            [orderedParam{i}, indGoodLapOrdered{i}] = ...
                sort(meanLickInterval(ind));
            indGoodLapOrdered{i} = ind(indGoodLapOrdered{i});
        end
    elseif(sortMethod == 12)
        % standard deviation of time between licks that occur after 1st 20cm of trial
        stdLickInterval = [ 1.563 0.181 0.542 0.267 0.11 0.107 0.1 0.073 0.047 0.096 0.112 0.059 0.086 0.064 0.053 0.088 0.092 0.063 0.083 0.076 0.472 0.09 0.18 0.095 0.014 0.084 0.078 0.066 0.047 0.047 0.117 0.088 0.231 0.077 0.112 0.055 0.055 0.03 0.13 0 0.098 0.05 0.06 0.034 0.08 0.072 0.045 0.06 0.028 0.033 0.041 0.064 0.045 0.056 0.07 0.081 0.05 0.053 0.084 0.024 0.113 0.061 0.01 0.034 0.062 0.1 0.051 0.103 0.056 0.067 0.082 0.072 0.037 0.065 0.058 0.017 0.054 0.102 0.037 0.115 0.062 0.064 0.058 0.063 0.136 0.073 0.076 0.065 0.092 0.055 0.069 0.051 0.083 0.074 0.089 0.071 0.062 0.113 0.081 0.055 0.066 0.054 0.086 0.046 0.044 0.061 0.055 0.042 0.039 0.054 0.054 0.04 0.044 0.048 0.077 0.031 0.052 0.037 0.048 1.869 0.599 0.077 0.108 0.074 0.059 0.099 0.066 0.034 0.036 0.04 0.047 0.084 0.036 0.065 0.04 0.036 0.052 0.063 0.228 0.081 0.036 0.035 0.052 0.039 0.034 0.044 0.034 0.069 0.055 0.072 0.03 0.051 0.014 0.04 0.049 0.044 0.024 0.048 0.029 0.028 0.029 0.039 1.989 0.064 0.063 0.022 0.059 0.044 0.074 0.052 0.037 0.064 0.047 0.174 0.133 0.112 0.051 0.073 0.044 0.058 0.057 0.044 0.043 0.048 0.056 0.073 0.059 0.031 0.028 0.028 0.05 0.035 0.031 0.053 0.045 0.044 0.091 0.033 0.036 0.059 0.042 0.05 0.039 0.056 0.052 0.038 0.041 0.289 0.02 0.037 0.051 0.076 0.039 0.032 0.038 0.049 0.047 0.038 0.055 0.01 0.056 0.074 0.011 0.053 0.044 0.067 0.049 0.027 0.038 0.044 0.037 0.047 0.023 0.053 0.029 0.052 0.055 0.057 0.05 0.056 0.07 0.045 0.038 0.035 0.058 0.033 0.044 0.028 0.034 0.046 0.041 0.048 0.046 0.042 0.031 0.048];
        for i = mazeSess'
            ind = find(mazeSessTr == i);
            indGoodLapOrdered{i} = ind;
            orderedParam{i} = ones(1,length(ind));    
        end
    elseif(sortMethod == 0)
        % original field
%         for i = mazeSess'
        for i = 1
%             ind = find(mazeSessTr == i);
%             ind = [45:53 72:88 92 93 106:108 117:123];
            % A013-20190504-01
            ind = [2.0, 3.0, 4.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 17.0, 18.0, 19.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 34.0, 35.0, 36.0, 37.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 64.0, 66.0, 72.0, 74.0, 79.0, 82.0, 83.0, 84.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 93.0, 94.0, 95.0, 97.0, 98.0, 101.0, 102.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 119.0, 120.0, 121.0, 122.0, 123.0, 124.0, 125.0, 126.0, 127.0, 131.0, 132.0, 135.0, 138.0, 139.0, 140.0, 141.0, 143.0, 146.0, 147.0, 148.0, 158.0, 159.0, 160.0, 162.0, 170.0, 171.0, 172.0, 173.0, 174.0, 175.0];
            ind = intersect(ind,beh.indGoodLap);
            orderedParam{i} = zeros(1,length(ind));
            indGoodLapOrdered{i} = ind;
        end
    elseif(sortMethod == 13)
        % kmean clustering result based on behavior
%         % A012-20190224
%         % 3 features
%         kmeanClu = [1 1 1 1 1 2 1 2 2 1 1 2 2 0 2 0 2 2 2 2 1 1 1 2 0 2 1 2 1 2 2 2 2 2 2 0 0 2 2 2 2 1 0 2 1 0 2 2 0 0 0 2 0 2 2 2 0 0 2 0 0 2 0 1 0 2 0 1 2 2 2 2 0 0 2 0 2 2 0 2 0 0 0 1 0 2 0 0 0 2 0 1 2 1 0 0 0 0 0 2 0 0 0 0 2 1 2 0 0 2 0 2 0 0 0 0 0 0 0 1 1 1 2 0 0 0 2 2 0 2 0 0 0 0 0 0 0 1 2 0 0 0 0 0 0 1 0 2 0 0 0 2 0 0 0 2 0 0 0 0 2 0 1 1 0 0 1 0 2 0 0 0 0 0 0 1];
%         % 5 features
%         kmeanClu = [2 2 2 2 2 2 2 2 0 2 2 2 0 1 1 1 2 2 1 1 2 2 0 2 1 1 2 2 2 2 1 2 2 2 1 1 1 2 1 1 1 2 1 1 0 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 0 1 1 2 1 1 1 2 2 1 1 1 1 1 0 1 2 1 1 2 1 1 1 2 1 1 1 1 1 1 1 2 0 2 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 2 1 1 1 1 1 1 1 1 1 0 2 0 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 0 1 1 1 1 1 0 1 2 1 1 1 2 1 1 1 1 1 1 1 1 1 1 2 2 1 1 2 1 1 1 1 1 1 2 0 0];
%         % 8 features
%         kmeanClu = [2 0 0 0 0 3 3 3 4 3 3 3 4 5 5 1 0 3 1 1 0 0 0 3 5 1 1 3 0 3 1 3 0 1 3 5 5 1 1 1 1 0 5 1 0 5 1 1 5 5 3 1 1 3 1 5 1 1 1 5 4 1 5 3 5 1 5 3 3 1 1 1 5 5 4 5 1 1 1 3 5 5 1 0 1 1 1 5 1 1 5 0 4 1 1 1 1 5 1 1 5 5 5 5 1 3 1 5 1 1 1 1 5 5 1 5 5 5 5 2 0 0 3 5 1 1 1 1 5 1 5 1 1 5 5 1 5 0 1 4 5 5 5 5 1 0 5 0 5 5 5 0 5 5 5 5 5 5 5 5 1 1 2 0 1 5 0 5 1 5 5 5 5 0 ];

% %       A012-20190221
%         kmeanClu = [0, 4, 4, 1, 4, 4, 4, 1, 0, 4, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0,...
%            4, 0, 4, 0, 0, 0, 0, 4, 1, 4, 1, 1, 1, 1, 1, 1, 4, 1, 4, 0, 0, 1,...
%            0, 1, 1, 0, 0, 0, 1, 4, 1, 1, 1, 4, 4, 4, 1, 0, 4, 0, 0, 4, 4, 4,...
%            0, 0, 4, 4, 1, 4, 4, 1, 1, 1, 4, 0, 1, 3, 0, 1, 4, 3, 4, 1, 4, 3,...
%            1, 1, 1, 1, 4, 4, 3, 4, 0, 1, 4, 1, 1, 1, 4, 4, 1, 1, 4, 0, 1, 2,...
%            4, 4, 1, 1, 4, 4, 4, 3, 3, 1, 4, 4, 4, 4, 3, 3, 4, 1, 3, 1, 1, 4,...
%            3, 4, 4, 1, 4, 4, 3, 4, 1, 4, 4, 0, 1, 3, 3, 4, 0, 0, 4, 4, 1, 4,...
%            4, 4, 4, 4, 1, 1, 5, 4, 1, 3, 4, 1, 1, 3, 1, 3, 4, 1, 0, 4, 4];
%        kmeanClu = [1, 2, 2, 0, 2, 2, 2, 0, 1, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,...
%        2, 1, 2, 1, 1, 1, 1, 2, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 2, 1, 1, 0,...
%        1, 0, 0, 1, 1, 1, 0, 2, 0, 0, 0, 2, 2, 2, 0, 1, 2, 1, 1, 2, 2, 2,...
%        1, 1, 2, 2, 0, 2, 2, 2, 2, 0, 2, 1, 0, 4, 1, 0, 2, 4, 2, 0, 2, 4,...
%        0, 0, 0, 0, 2, 2, 4, 2, 1, 2, 2, 0, 0, 0, 2, 2, 0, 0, 2, 1, 0, 3,...
%        2, 2, 0, 0, 2, 2, 2, 4, 4, 0, 2, 2, 2, 2, 4, 4, 2, 0, 4, 0, 0, 2,...
%        4, 2, 2, 0, 2, 2, 4, 2, 0, 2, 2, 1, 0, 4, 4, 2, 1, 1, 2, 2, 0, 2,...
%        2, 2, 2, 2, 0, 0, 5, 2, 0, 4, 2, 0, 0, 4, 2, 4, 2, 2, 1, 2, 2];
   
%        % A011-20190218
%        kmeanClu = [0, 4, 0, 1, 4, 0, 0, 1, 4, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 2, 0, 0,...
%        2, 0, 0, 2, 2, 2, 2, 5, 5, 2, 2, 2, 2, 5, 2, 2, 2, 2, 5, 1, 0, 4,...
%        0, 2, 2, 2, 2, 2, 5, 0, 3, 3, 0, 2, 2, 5, 2, 5, 5, 2, 2, 5, 5, 2,...
%        5, 5, 5, 5, 5, 5, 2, 2, 2, 5, 5, 5, 4, 3, 2, 2, 5, 5, 2, 2, 0, 2,...
%        2, 3, 3];    
%         ind =  [72:114 116:161 163:164];
%         [cluTmp, indGoodLapOrderedTmp] = ...
%                 sort(kmeanClu);
%         clus = unique(kmeanClu);
%         for i = 1:length(clus)
%             ind1 = find(cluTmp == clus(i));
%             orderedParam{i} = zeros(1,length(ind1));    
%             indGoodLapOrdered{i} = ind(indGoodLapOrderedTmp(ind1));
%         end    

%         % A011-20190219
        kmeanClu = [0, 1, 1, 1, 1, 1, 4, 4, 4, 1, 4, 4, 1, 4, 5, 4, 5, 5, 4, 4, 4, 0,...
       0, 1, 4, 4, 5, 1, 4, 1, 4, 4, 4, 1, 1, 1, 1, 4, 4, 4, 4, 5, 4, 4,...
       4, 4, 1, 1, 1, 4, 4, 1, 4, 1, 4, 1, 5, 4, 4, 4, 4, 4, 1, 4, 1, 4,...
       4, 4, 4, 4, 1, 4, 5, 4, 4, 4, 4, 5, 4, 4, 0, 4, 4, 4, 4, 4, 5, 4,...
       4, 4, 4, 4, 1, 4, 1, 4, 4, 1, 4, 0, 4, 4, 1, 4, 5, 5, 4, 1, 1, 3,...
       4, 0, 4, 4, 4, 1, 1, 4, 4, 5, 5, 0, 0, 1, 1, 4, 1, 0, 2, 0, 3, 4,...
       1, 1, 0, 0, 5, 0, 4, 5, 0, 0, 0, 2, 5, 2, 0, 0, 0, 0, 4, 0, 4, 1,...
       1, 4, 1, 4, 4, 1, 4, 4, 0, 1, 5, 0, 4, 4, 0, 0, 0, 4, 4, 4, 4, 4,...
       4, 4, 4, 4, 4, 4, 4, 0, 0, 5, 5, 4, 1, 0];

        ind = find(mazeSessTr == 1);
        [cluTmp, indGoodLapOrderedTmp] = ...
                sort(kmeanClu(ind));
        clus = unique(kmeanClu);
        for i = 1:length(clus)
            ind = find(cluTmp == clus(i));
            orderedParam{i} = zeros(1,length(ind));    
            indGoodLapOrdered{i} = indGoodLapOrderedTmp(ind);
        end
    elseif(sortMethod == 14) % ordered by activity correlation
            activityCorr = [...
           255   243   222   240   245   242   250   211   236   220   246   192   254   195   216   217 ...
           248   210   232   253   241   218   244   225   228   191   219   212   224   194   230   229 ...
           238   182   226   235   249   215   198   251   234   209   221   203   231   214   227   256 ...
           199   252   204   247   201   213   205   193   206   233   189   239   196   207   180   223 ...
           185   181   237   187   197   183   200   208   188   202   190   174   186   179   178   184 ...
           177   175   176   152   121   170   171   155   153   165   120   166   154     4   169   168 ...
           122    93   119   124   123   147     1   141   173   136     3   128   162   149   101   159 ...
           167   163   148   160     8   140   104   164   138    99     9   172   114   161   102   132 ...
           145    76   150   151   158    95     5   139   135   125    61   143   156   131   146   133 ...
           118    63   100   157   103    96   137    98     2   115   113    94    60   142    81    83 ...
           107   134   116   106   126    75    58    82   105    53   144     6   127     7    59    46 ...
            97    14    62    54    92   108   129    89    85    44    47    15    73    45    16    11 ...
            52    86    19    50    87    91   117   109    71   130    24    21    78   112    13   110 ...
            84    90    20    43    17    33    37    66    40    10    26    57    88    42    72    29 ...
            18    25    12    56    79    55    51   111    67    70    74    49    65    77    27    41 ...
            68    36    31    69    34    28    80    23    32    64    22    48    39    30    35    38]; 
        for i = mazeSess'
            ind = find(mazeSessTr == i);
            orderedParam{i} = zeros(1,length(ind)); 
            [~,indGoodLapOrdered{i}] = intersect(activityCorr,ind,'stable');
            indGoodLapOrdered{i} = activityCorr(indGoodLapOrdered{i});
        end
    elseif(sortMethod == 15) % ordered by left/right trials
        for i = mazeSess'
            ind = find(mazeSessTr == i);
            ind = intersect(ind, beh.indGoodLap);
            orderedParam{i} = zeros(1,length(ind)); 
%             [~,indGoodLapOrdered{i}] = intersect(beh.trStartCueLR,ind,'stable');
            [orderedParam{i}, indGoodLapOrdered{i}] = ...
                sort(beh.trStartCueLR(ind)');
            indGoodLapOrdered{i} = ind(indGoodLapOrdered{i});
        end
    end
    
    count = 0;
    for i = indNeurons
        count = count + 1;

        if(mod(count-1,16) == 0)
            [figNew,pos] = CreateFig();
            set(0,'Units','pixels') 
            figTitle = 'All the trials';
            set(figure(figNew),'OuterPosition',...
                [pos(1) pos(2) pos(3)*2 pos(4)*2.2],'Name',figTitle)
        end
        
        subplot(4,4,mod(count-1,16)+1)
        
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ') D' ...
                    num2str(depthNeu.relDepthNeuHDef(i)) ' Sp' ...
                    num2str(spatialInfoSess{1}.adaptSpatialInfo(i))];
                
%         plotFRProfIndNeuronIndTrial(gca,filteredSpikeArrayNormT,maxNumSamples,i,...
%                    orderedParam,indGoodLapOrdered,figTitle);
        plotSpikesIndNeuronIndTrial(gca,trialsExt,i,orderedParam,...
            indGoodLapOrdered,figTitle);
               
    end
end

function plotSpikesIndNeuronIndTrial(handle,spikes,neuronNo,paramSess,indLapSess,figTitle)
    numTrTot = 0;
    numTrSess = [];
    paramArr = [];
    for sess = 1:length(paramSess)
        numTrTot = numTrTot + length(paramSess{sess});
        numTrSess = [numTrSess numTrTot];
        paramArr = [paramArr paramSess{sess}];
    end
    
    numTr = 0;
    hold on;
    for sess = 1:length(paramSess)
        for i = 1:length(indLapSess{sess})
            numTr = numTr+1;
            h = plot(spikes{indLapSess{sess}(i)}.spikesMM{neuronNo},...
                numTr*ones(1,length(spikes{indLapSess{sess}(i)}.spikesMM{neuronNo})),...
                'k.');
            set(h,'MarkerSize',5);
        end
    end
    for i = 1:length(numTrSess)-1
        h = plot([0 1800],numTrSess(i)*ones(1,2),'r');
        set(h,'LineWidth',1);
    end
    set(gca,'XLim',[0 1800],'YLim',[0 numTrTot],'Ydir','reverse');
    xlabel('Dist (mm)');
    ylabel('Trial No.');
    title(figTitle);
end

function plotFRProfIndNeuronIndTrial...
            (handle,filteredSpikeArrayNormT,numSamples,neuronNo,paramSess,indLapSess,figTitle)
    numTrTot = 0;
    numTrSess = [];
    paramArr = [];
    for sess = 1:length(paramSess)
        numTrTot = numTrTot + length(paramSess{sess});
        numTrSess = [numTrSess numTrTot];
        paramArr = [paramArr paramSess{sess}];
    end
    FRProfilePerTrial = zeros(numTrTot,numSamples);
    
    numTr = 0;
    
    for sess = 1:length(paramSess)
        for i = 1:length(indLapSess{sess})
            numTr = numTr+1;
            szArr = size(filteredSpikeArrayNormT{indLapSess{sess}(i)},2);
            FRProfilePerTrial(numTr,1:szArr) = filteredSpikeArrayNormT{indLapSess{sess}(i)}(neuronNo,:); 
            %./max(filteredSpikeArrayNormT{indLaps(i)}(neuronNo,:));
        end
    end  
    FRProfilePerTrial = FRProfilePerTrial/max(FRProfilePerTrial(:));
    h = imagesc(0:numSamples-1,1:numTr,FRProfilePerTrial);
    hold on
    for i = 1:length(numTrSess)-1
        h = plot([0 numSamples],numTrSess(i)*ones(1,2),'r');
        set(h,'LineWidth',1);
    end
    for i = 1:length(paramArr)
        if(~isnan(paramArr(i)))
            h = plot([paramArr(i) paramArr(i)],[i-1, i],'r');
            set(h,'LineWidth',1);
        end
    end
            
    set(gca,'FontSize',8.0,'Box','on','XLim',[0 numSamples-1],'YLim',[0 numTr]);
    xlabel('Dist (mm)');
    ylabel('Trial No.');
    title(figTitle);
end