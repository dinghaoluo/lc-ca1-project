function InterpolateWheel_smTr(file,lfpFreq)
% interpolate the wheel events to sampling frequency

if exist([file '-whl.mat'], 'file') == 2
    disp('whl file already exists.')
    return;
else
    fullNameBLFP = [file 'BLFP.mat'];
    if(exist(fullNameBLFP,'file')~=0)
        load(fullNameBLFP);
    else return;
    end
    
    DatLength=(behEventsTdt.taskDescr(end,1)-behEventsTdt.taskDescr(1,1)+1); %recording length
    
    %% initiate variables
    timeStepLfp = 1/lfpFreq;
    whlDataLfp(:,1) = (0:timeStepLfp:(DatLength-1)/lfpFreq)*1000; % in ms
    
    behType = zeros(max(behEventsTdt.trialDescr(:,2)),1);
    trStartLfpIndArr = zeros(max(behEventsTdt.trialDescr(:,2)),1);
    trEndLfpIndArr = zeros(max(behEventsTdt.trialDescr(:,2)),1);
    trackLenArr = zeros(max(behEventsTdt.trialDescr(:,2)),1);
    
    % remove any treadmill event before the beginning of the first trial     
    indTmp = find(behEventsTdt.wheel(:,1) >= behEventsTdt.taskDescr(1,1) & ...
        behEventsTdt.wheel(:,1) <= behEventsTdt.taskDescr(end,1));
    wheelDist = behEventsTdt.wheel(indTmp,3);
    wheelTimeLfp = behEventsTdt.wheel(indTmp,1);
    indErr = diff(wheelTimeLfp) == 0;
    if(sum(indErr)~=0)
        indGood = [true; ~indErr];
        wheelDist = wheelDist(indGood);
        wheelTimeLfp = wheelTimeLfp(indGood);
    end
        
    % interpolate the wheel distance between each pair of wheel event
    interpTspDistLfp1 = zeros(length(whlDataLfp),1);
    
    EncoderToDist = 0.4; % mm- conversion factor is 0.04cm/click
    for i = 1:length(wheelTimeLfp)-1
        if(wheelDist(i+1) >= wheelDist(i))
            xWheelLfp = interp1([wheelTimeLfp(i) wheelTimeLfp(i+1)],...
                [wheelDist(i) wheelDist(i+1)]*EncoderToDist,...
                wheelTimeLfp(i):wheelTimeLfp(i+1)-1,'linear');
            interpTspDistLfp1(wheelTimeLfp(i)-behEventsTdt.taskDescr(1,1)+1:...
                wheelTimeLfp(i+1)-behEventsTdt.taskDescr(1,1)) = ...
                xWheelLfp;
        else % if there is a distance reset, then start interpolation from 0
            xWheelLfp = interp1([wheelTimeLfp(i)+1 wheelTimeLfp(i+1)],...
                    [0 wheelDist(i+1)]*EncoderToDist,...
                    wheelTimeLfp(i)+1:wheelTimeLfp(i+1)-1,'linear');
            interpTspDistLfp1(wheelTimeLfp(i)-behEventsTdt.taskDescr(1,1)+1:...
                wheelTimeLfp(i+1)-behEventsTdt.taskDescr(1,1)) = ...
                    [wheelDist(i)*EncoderToDist xWheelLfp];
        end
        if(i == length(wheelTimeLfp)-1)    
            interpTspDistLfp1(wheelTimeLfp(i+1)-behEventsTdt.taskDescr(1,1)+1) = ...
                wheelDist(i+1)*EncoderToDist;
        end
        
    end
     
    %% interpolate the running distance for each trial
    for tr = 1:max(behEventsTdt.trialDescr(:,2))
        
        indTaskDescr = find(behEventsTdt.taskDescr(:,1) < ...
                behEventsTdt.trialDescr(tr,1),1,'last');
        trackLenArr(tr) = behEventsTdt.taskDescr(indTaskDescr,2)*10; % mm
        behType(tr) = 1;
        
        if(str2num(file(2:4)) < 11) % A001-A010
            %%% changed on 1/21/2019 the trial start is defined by TR, while
            %%% trial end is defined by the first pump event in the trial
            indPump = find(behEventsTdt.pump(:,1) > ...
                behEventsTdt.taskDescr(indTaskDescr,1),1,'first');
            trStartLfpIndArr(tr) = behEventsTdt.taskDescr(tr,1) - ...
                behEventsTdt.taskDescr(1,1) + 1;
            trEndLfpIndArr(tr) = behEventsTdt.pump(indPump,1) - ...
                behEventsTdt.taskDescr(1,1) + 1;
        else % changed on 2/16/2019 for the active licking task
            trStartLfpIndArr(tr) = behEventsTdt.taskDescr(indTaskDescr,1) - ...
                behEventsTdt.taskDescr(1,1) + 1;           
            trEndLfpIndArr(tr) = behEventsTdt.trialDescr(tr,1) - ...
                behEventsTdt.taskDescr(1,1) + 1;
        end
        
    end
    
    whlDataLfp(:,2) = interpTspDistLfp1;
        
    disp(['Total number of samples: ' int2str(DatLength)]);
    disp(['Lines in whlDataLfp:' int2str(size(whlDataLfp,1))]);

    %save it
    save([file '-whl.mat'], 'whlDataLfp','behType',...
         'trStartLfpIndArr','trEndLfpIndArr', 'trackLenArr','-v7.3');
end
