function aligntsp2dat_mpfi_pix(file,sampleFreq,lfpFreq,ChanNum)
% also, intan does not create meta file, had to change how DatLength is
% calculated
if exist([file '-whl.mat'], 'file') == 2
    disp('whl file already exists.')
    return;
else
   
    listing = dir([file '.dat']);
    DatLength = listing.bytes/(ChanNum*2);
       
    %% initiate variables
    timeStep = 1/sampleFreq;
    timeStepLfp = 1/lfpFreq;
    interpTsp(:,1) = (0:timeStep:(DatLength-1)/sampleFreq)*1000; % in ms
    whlDataLfp(:,1) = (0:timeStepLfp:(DatLength-1)/sampleFreq)*1000; % in ms
      
    %% get timestamps and positions from BTDT file
    load([file 'BTDT.mat']);
    
    behType = zeros(max(behEventsTdt.trialDescr(:,3)),1);
    trStartLfpIndArr = zeros(max(behEventsTdt.trialDescr(:,3)),1);
    trEndLfpIndArr = zeros(max(behEventsTdt.trialDescr(:,3)),1);
    trackLenArr = zeros(max(behEventsTdt.trialDescr(:,3)),1);
      
    % remove any treadmill event before the beginning of the first trial     
    indTmp = find(behEventsTdt.wheel(:,1) >= behEventsTdt.taskDescr(1,1));
    wheelDist = behEventsTdt.wheel(indTmp,4);
    wheelTime = behEventsTdt.wheel(indTmp,1);
    wheelTimeLfp = behEventsTdt.wheel(indTmp,2);
    
    % added to avoid 0 starting index 
    if(wheelTime(1)==0)
        wheelTime = wheelTime+1;
    end
    if(wheelTimeLfp(1)==0)
        wheelTimeLfp = wheelTimeLfp+1;
    end
    
    indErr = diff(wheelTimeLfp) == 0;
    if(sum(indErr)~=0)
        indGood = [true; ~indErr];
        wheelDist = wheelDist(indGood);
        wheelTime = wheelTime(indGood);
        wheelTimeLfp = wheelTimeLfp(indGood);
    end
    
    % interpolate the wheel distance from start of each trial (TR) to the
    % end of each trial (NT)
    interpTspDist = zeros(length(interpTsp),1);
    interpTspDistLfp = zeros(length(whlDataLfp),1);
    
    % interpolate the wheel distance between each pair of wheel event
    interpTspDist1 = zeros(length(interpTsp),1);
    interpTspDistLfp1 = zeros(length(whlDataLfp),1);
    
    EncoderToDist = 0.4; % mm- conversion factor is 0.04cm/click
    for i = 1:length(wheelTime)-1
        if(wheelDist(i+1) >= wheelDist(i))
            
            xWheel = interp1([wheelTime(i) wheelTime(i+1)],...
                    [wheelDist(i) wheelDist(i+1)]*EncoderToDist,...
                    wheelTime(i):wheelTime(i+1)-1,'linear');
            interpTspDist1(wheelTime(i):wheelTime(i+1)-1) = xWheel;
            xWheelLfp = interp1([wheelTimeLfp(i) wheelTimeLfp(i+1)],...
                [wheelDist(i) wheelDist(i+1)]*EncoderToDist,...
                wheelTimeLfp(i):wheelTimeLfp(i+1)-1,'linear');
            interpTspDistLfp1(wheelTimeLfp(i):wheelTimeLfp(i+1)-1) = ...
                xWheelLfp;
        else % if there is a distance reset, then start interpolation from 0
            xWheel = interp1([wheelTime(i)+1 wheelTime(i+1)],...
                    [0 wheelDist(i+1)]*EncoderToDist,...
                    wheelTime(i)+1:wheelTime(i+1)-1,'linear');           
            interpTspDist1(wheelTime(i):wheelTime(i+1)-1) = ...
                    [wheelDist(i)*EncoderToDist xWheel];    
            xWheelLfp = interp1([wheelTimeLfp(i)+1 wheelTimeLfp(i+1)],...
                    [0 wheelDist(i+1)]*EncoderToDist,...
                    wheelTimeLfp(i)+1:wheelTimeLfp(i+1)-1,'linear');
            interpTspDistLfp1(wheelTimeLfp(i):wheelTimeLfp(i+1)-1) = ...
                    [wheelDist(i)*EncoderToDist xWheelLfp];
        end
        if(i == length(wheelTime)-1)    
            interpTspDist1(wheelTime(i+1)) = wheelDist(i+1)*EncoderToDist;
            interpTspDistLfp1(wheelTimeLfp(i+1)) = wheelDist(i+1)*EncoderToDist;
        end
        
    end
    
    %% interpolate the running distance for each trial
    for tr = 1:max(behEventsTdt.trialDescr(:,3))
        
        indTaskDescr = find(behEventsTdt.taskDescr(:,1) < ...
                behEventsTdt.trialDescr(tr,1),1,'last');
        trackLenArr(tr) = behEventsTdt.taskDescr(indTaskDescr,3)*10; % mm
        behType(tr) = 1;
        
%         trStartInd = behEventsTdt.trialDescr(tr,1)+1;
%         trEndInd = behEventsTdt.trialDescr(tr+1,1);
%         
%         trStartLfpIndArr(tr) = behEventsTdt.trialDescr(tr,2)+1;
%         trEndLfpIndArr(tr) = behEventsTdt.trialDescr(tr+1,2);
        
        if(str2num(file(2:4)) < 11) % A001-A010
            %%% changed on 1/21/2019 the trial start is defined by TR, while
            %%% trial end is defined by the first pump event in the trial
            indPump = find(behEventsTdt.pump(:,1) > ...
                behEventsTdt.taskDescr(indTaskDescr,1),1,'first');
            trStartInd = behEventsTdt.taskDescr(tr,1);
%             trEndInd = behEventsTdt.trialDescr(tr,1); % changed on 1/22/2020
            trEndInd = behEventsTdt.pump(indPump);
            
            trStartLfpIndArr(tr) = behEventsTdt.taskDescr(tr,2);
            trEndLfpIndArr(tr) = behEventsTdt.pump(indPump,2);
        else % changed on 2/16/2019 for the active licking task
            trStartInd = behEventsTdt.taskDescr(indTaskDescr,1);
            trEndInd = behEventsTdt.trialDescr(tr,1);

            trStartLfpIndArr(tr) = behEventsTdt.taskDescr(indTaskDescr,2);           
            trEndLfpIndArr(tr) = behEventsTdt.trialDescr(tr,2);
        end
        
        % added to avoid 0 staring index
        if(trStartInd(1)==0)
            trStartInd = 1;
        end
        if(trStartLfpIndArr(1)==0)
            trStartLfpIndArr = 1;
        end
        
        % added to avoid negative staring index, 15 Feb 2024 Dinghao 
        if(trStartInd(1)<0)
            trStartInd = 1;
        end
        if(trStartLfpIndArr(1)<0)
            trStartLfpIndArr = 1;
        end
        
        disp(tr);
        
        %% interpolate the distance of the treadmill events based 
        % on the track length
        indTmp = find(wheelTime >= trStartInd ...
                        & wheelTime <= trEndInd);
            % find the treadmill events between two trials
        wheelDistTr = wheelDist(indTmp);
        diffDist = diff(wheelDistTr);
        ind = find(diffDist < -10,1);
            % check whether the distance is monotonically increasing 
        if(~isempty(ind))
            indTmp = indTmp(ind+1:end);
            wheelDistTr = wheelDist(indTmp);
        end % adjust the start point of the trial
        
        %%        
        wheelTimeTr = wheelTime(indTmp);
        wheelTimeTrLfp = wheelTimeLfp(indTmp);
        distStart = 0;
        distEnd = wheelDistTr(end) + ...
            (wheelDistTr(end) - wheelDistTr(end-1))...
            *(trEndInd - wheelTimeTr(end)) ...
            /(wheelTimeTr(end) - wheelTimeTr(end-1));
        %% changed by yingxue on 6/1/2019
        distEnd = distEnd*EncoderToDist;
        %%

        %%%%% Changed on 4/8/2019 by Yingxue 
        %%% the trial end distance can by larger than trackLen
        xBeh = [distStart; wheelDistTr*EncoderToDist; distEnd];
%         xBeh = interp1([distStart distEnd],[0 trackLenArr(tr)],...
%                 [distStart; wheelDistTr; distEnd],'linear');
            % adjust the running distance according to the actual track 
            % length
            
%       disp(['Trial ' num2str(tr) ' beam ' num2str(i)])
        % interpolate running distance
        if(wheelTimeTr(1) == trStartInd)
            wheelTimeTr(1) = wheelTimeTr(1) + 0.0001;
        end
        if(wheelTimeTr(end) == trEndInd)
            wheelTimeTr(end) = wheelTimeTr(end) - 0.0001;
        end
        xBeh2Time = interp1(...
           [trStartInd;wheelTimeTr;trEndInd],...
           xBeh,trStartInd:trEndInd,'linear');
                % interpolate the running distance to sampling frequency
        if(trStartLfpIndArr(tr) == wheelTimeTrLfp(1))
            wheelTimeTrLfp(1) = wheelTimeTrLfp(1) ...
                                + lfpFreq/sampleFreq*0.0001;
        end
        if(trEndLfpIndArr(tr) == wheelTimeTrLfp(end))
            wheelTimeTrLfp(end) = wheelTimeTrLfp(end) ...
                                - lfpFreq/sampleFreq*0.0001;
        end
        xBeh2TimeLfp = interp1([trStartLfpIndArr(tr); wheelTimeTrLfp;...
            trEndLfpIndArr(tr)],xBeh,...
            trStartLfpIndArr(tr):trEndLfpIndArr(tr),'linear');
        
        interpTspDist(trStartInd:trEndInd) = xBeh2Time;
        interpTspDistLfp(trStartLfpIndArr(tr):trEndLfpIndArr(tr)) = ...
            xBeh2TimeLfp;
    end
    
    interpTsp(:,2) = interpTspDist;
    whlDataLfp(:,2) = interpTspDistLfp;
    
    interpTsp(:,3) = interpTspDist1;
    whlDataLfp(:,3) = interpTspDistLfp1;
        
    disp(['Samples in .dat file per channel: ' int2str(DatLength)]);
    disp(['Lines in interpTsp:' int2str(size(interpTsp,1))]);
    disp(['Lines in whlDataLfp:' int2str(size(whlDataLfp,1))]);

    %save it
    save([file '-whl.mat'], 'whlDataLfp','interpTsp','behType',...
         'trStartLfpIndArr','trEndLfpIndArr', 'trackLenArr','-v7.3');
   
end
