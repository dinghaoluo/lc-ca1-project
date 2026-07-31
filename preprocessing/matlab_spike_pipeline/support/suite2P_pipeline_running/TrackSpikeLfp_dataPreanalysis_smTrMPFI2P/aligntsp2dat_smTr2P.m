function aligntsp2dat_mpfi(file,sampleFreq,lfpFreq,ChanNum, data_2p)

if exist([file '-whl.mat'], 'file') == 2
    disp('whl file already exists.')
    return;
else    
    %% get timestamps and positions from BTDT file
    load([file 'BTDT.mat']);
    
    timeStep = 1000/lfpFreq;
    timeStepSf = 1000/sampleFreq;
    
%     interpTsp(:,1) = 0:timeStepSf:behEventsTdt.TDTsyncMsecLfp(end); % in ms
%     whlDataLfp(:,1) = 0:timeStep:behEventsTdt.TDTsyncMsecLfp(end); % in ms
    
    interpTsp(:,1) = behEventsTdt.ArdSyncMsec; % in ms
    whlDataLfp(:,1) = 0:timeStep:behEventsTdt.ArdSyncMsecLfp(end)/lfpFreq*1000; % in ms
    
    behType = zeros(max(behEventsTdt.trialDescr(:,5)),1);
    trStartLfpIndArr = zeros(max(behEventsTdt.trialDescr(:,5)),1);
    trEndLfpIndArr = zeros(max(behEventsTdt.trialDescr(:,5)),1);
    trackLenArr = zeros(max(behEventsTdt.trialDescr(:,5)),1);
    
    % remove any treadmill event before the beginning of the first trial     
    indTmp = find(behEventsTdt.wheel(:,4) >= behEventsTdt.taskDescr(1,4));
    wheelDist = behEventsTdt.wheel(indTmp,6);
    wheelTime = behEventsTdt.wheel(indTmp,3);
    wheelTimeLfp = behEventsTdt.wheel(indTmp,4);
    
    twoPTime = behEventsTdt.ArdSyncMsec;
    twoPTimeLfp = behEventsTdt.ArdSyncMsecLfp;
    
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
    
    % interpolate the 2p data
    load([file '_corrFluo.mat'],'Clu','dFF');
    cluNo = Clu.localClu;
    interpTspF = zeros(length(cluNo), length(whlDataLfp));
    interpTspFneu = zeros(length(cluNo),length(whlDataLfp));
    interpTspSpks = zeros(length(cluNo), length(whlDataLfp));
    interpTspdFF = zeros(length(cluNo), length(whlDataLfp));
    
    %% interpolate wheel events (sample freq)
    disp('interpolate wheel events (sample freq)')
    EncoderToDist = 0.4; % mm- conversion factor is 0.04cm/click
    [wheelTimeTmp,indTime] = unique(wheelTime);
    wheelDistTmp = wheelDist(indTime);
    for i = 1:length(wheelTimeTmp)-1
        if(wheelDistTmp(i+1) >= wheelDistTmp(i))
            if(wheelTimeTmp(i+1) > wheelTimeTmp(i)+1)
                xWheel = interp1([wheelTimeTmp(i) wheelTimeTmp(i+1)],...
                    [wheelDistTmp(i) wheelDistTmp(i+1)]*EncoderToDist,...
                    wheelTimeTmp(i):wheelTimeTmp(i+1)-1,'linear');
                interpTspDist1(wheelTimeTmp(i):wheelTimeTmp(i+1)-1) = xWheel;
            elseif(wheelTimeTmp(i+1) == wheelTimeTmp(i)+1)
                interpTspDist1(wheelTimeTmp(i)) = wheelDistTmp(i)*EncoderToDist;
            end
                        
        else % if there is a distance reset, then start interpolation from 0
            if(wheelTimeTmp(i+1) > wheelTimeTmp(i)+2)
                xWheel = interp1([wheelTimeTmp(i)+1 wheelTimeTmp(i+1)],...
                        [0 wheelDistTmp(i+1)]*EncoderToDist,...
                        wheelTimeTmp(i)+1:wheelTimeTmp(i+1)-1,'linear');   
                interpTspDist1(wheelTimeTmp(i):wheelTimeTmp(i+1)-1) = ...
                    [wheelDistTmp(i)*EncoderToDist xWheel];
            elseif(wheelTimeTmp(i+1) == wheelTimeTmp(i)+2)
                interpTspDist1(wheelTimeTmp(i):wheelTimeTmp(i+1)-1) = ...
                    [wheelDistTmp(i)*EncoderToDist 0];
            else
                interpTspDist1(wheelTimeTmp(i)) = 0;
            end
        end
        if(i == length(wheelTimeTmp)-1)    
            interpTspDist1(wheelTimeTmp(i+1)) = wheelDistTmp(i+1)*EncoderToDist;
        end        
    end
    
    %% interpolate wheel events (lfp freq)
    disp('interpolate wheel events (lfp freq)')
    for i = 1:length(wheelTimeLfp)-1
        if(wheelDist(i+1) >= wheelDist(i))       
            xWheelLfp = interp1([wheelTimeLfp(i) wheelTimeLfp(i+1)],...
                [wheelDist(i) wheelDist(i+1)]*EncoderToDist,...
                wheelTimeLfp(i):wheelTimeLfp(i+1)-1,'linear');
            interpTspDistLfp1(wheelTimeLfp(i):wheelTimeLfp(i+1)-1) = ...
                xWheelLfp;
        else % if there is a distance reset, then start interpolation from 0               
            xWheelLfp = interp1([wheelTimeLfp(i)+1 wheelTimeLfp(i+1)],...
                    [0 wheelDist(i+1)]*EncoderToDist,...
                    wheelTimeLfp(i)+1:wheelTimeLfp(i+1)-1,'linear');
            interpTspDistLfp1(wheelTimeLfp(i):wheelTimeLfp(i+1)-1) = ...
                    [wheelDist(i)*EncoderToDist xWheelLfp];
        end
        if(i == length(wheelTimeLfp)-1)    
            interpTspDistLfp1(wheelTimeLfp(i+1)) = wheelDist(i+1)*EncoderToDist;
        end
    end
    
    %% interpolate the running distance for each trial
    disp('interpolate the running distance for each trial')
    for tr = 1:max(behEventsTdt.trialDescr(:,5))
        indTaskDescr = find(behEventsTdt.taskDescr(:,4) < ...
                behEventsTdt.trialDescr(tr,4),1,'last');
        trackLenArr(tr) = behEventsTdt.taskDescr(indTaskDescr,5)*10; % mm
        behType(tr) = 1;
        
        %% for the active licking task
        trStartInd = behEventsTdt.taskDescr(indTaskDescr,3);
        trEndInd = behEventsTdt.trialDescr(tr,3);

        trStartLfpIndArr(tr) = behEventsTdt.taskDescr(indTaskDescr,4);           
        trEndLfpIndArr(tr) = behEventsTdt.trialDescr(tr,4);
        
        disp(tr);
        
        %% interpolate the distance of the treadmill events based 
        % on the track length
        indTmp = find(wheelTimeLfp >= trStartLfpIndArr(tr) ...
                        & wheelTimeLfp <= trEndLfpIndArr(tr));
            % find the treadmill events between two trials
        wheelDistTr = wheelDist(indTmp);
        diffDist = diff(wheelDistTr);
        ind = find(diffDist < -10,1);
            % check whether the distance is monotonically increasing 
        if(~isempty(ind))
            indTmp = indTmp(ind+1:end);
            wheelDistTr = wheelDist(indTmp);
        end % adjust the start point of the trial
        
        %% for sample frequency wheel distance       
        [wheelTimeTr,indWhT] = unique(wheelTime(indTmp));  
        distStart = 0;
        distEnd = wheelDistTr(end) + ...
            (wheelDistTr(end) - wheelDistTr(end-1))...
            *(trEndInd - wheelTimeTr(end)) ...
            /(wheelTimeTr(end) - wheelTimeTr(end-1));
        distEnd = distEnd*EncoderToDist;
        %%% the trial end distance can by larger than trackLen
        xBeh = [distStart; wheelDistTr(indWhT)*EncoderToDist; distEnd];
        %%
        
        %% for LFP frequency wheel distance   
        wheelTimeTrLfp = wheelTimeLfp(indTmp);
        distEnd = wheelDistTr(end) + ...
            (wheelDistTr(end) - wheelDistTr(end-1))...
            *(trEndLfpIndArr(tr) - wheelTimeTrLfp(end)) ...
            /(wheelTimeTrLfp(end) - wheelTimeTrLfp(end-1));
        distEnd = distEnd*EncoderToDist;
        %%% the trial end distance can by larger than trackLen
        xBehLfp = [distStart; wheelDistTr*EncoderToDist; distEnd];
            % adjust the running distance according to the actual track 
            % length
        %%
            
%       disp(['Trial ' num2str(tr) ' beam ' num2str(i)])
        % interpolate running distance to sampling frequency
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
                
        % interpolate running distance to LFP frequency
        if(trStartLfpIndArr(tr) == wheelTimeTrLfp(1))
            wheelTimeTrLfp(1) = wheelTimeTrLfp(1) ...
                                + 0.0001;
        end
        if(trEndLfpIndArr(tr) == wheelTimeTrLfp(end))
            wheelTimeTrLfp(end) = wheelTimeTrLfp(end) ...
                                - 0.0001;
        end
        xBeh2TimeLfp = interp1([trStartLfpIndArr(tr); wheelTimeTrLfp;...
            trEndLfpIndArr(tr)],xBehLfp,...
            trStartLfpIndArr(tr):trEndLfpIndArr(tr),'linear');
                % interpolate the running distance to lfp frequency
        interpTspDist(trStartInd:trEndInd) = xBeh2Time;
        interpTspDistLfp(trStartLfpIndArr(tr):trEndLfpIndArr(tr)) = ...
            xBeh2TimeLfp;
    end
           
    disp('interpolate the neural data')
    for n = 1:length(cluNo) % range all the cells
        i = cluNo(n);
        for j = 1:(length(twoPTimeLfp)-1)
%           range of arduinoarcsync

%           interpolate in between
            interpF = interp1([twoPTimeLfp(j) (twoPTimeLfp(j+1))],...
                    [data_2p.F(i, j) data_2p.F(i, j+1)],...
                    twoPTimeLfp(j)+1:(twoPTimeLfp(j+1)),'linear');
            interpTspF(n, twoPTimeLfp(j)+1:twoPTimeLfp(j+1)) = interpF;

            interpFneu = interp1([twoPTimeLfp(j) (twoPTimeLfp(j+1))],...
                    [data_2p.Fneu(i, j) data_2p.Fneu(i, j+1)],...
                    twoPTimeLfp(j)+1:(twoPTimeLfp(j+1)),'linear');
            interpTspFneu(n, twoPTimeLfp(j)+1:twoPTimeLfp(j+1)) = interpFneu;

            interpSpks = interp1([twoPTimeLfp(j) ( twoPTimeLfp(j+1))],...
                    [data_2p.spks(i, j) data_2p.spks(i, j+1)],...
                    twoPTimeLfp(j)+1:(twoPTimeLfp(j+1)),'linear');
            interpTspSpks(n, twoPTimeLfp(j)+1:twoPTimeLfp(j+1)) = interpSpks;
            
            interpdFF = interp1([twoPTimeLfp(j) (twoPTimeLfp(j+1))],...
                    [dFF(n, j) dFF(n, j+1)],...
                    twoPTimeLfp(j)+1:(twoPTimeLfp(j+1)),'linear');
            interpTspdFF(n, twoPTimeLfp(j)+1:twoPTimeLfp(j+1)) = interpdFF;

%             now repeat for Fneu and spikes (and also get the total num of
%             frames between the first and last recorded/saved FMs for
%             step_2p
        end  
    end
    
    interpTsp(:,2) = interpTspDist;
    whlDataLfp(:,2) = interpTspDistLfp;
    
    interpTsp(:,3) = interpTspDist1;
    whlDataLfp(:,3) = interpTspDistLfp1;
        
%     disp(['Samples in .dat file per channel: ' int2str(DatLength)]);
    disp(['Lines in interpTsp:' int2str(size(interpTsp,1))]);
    disp(['Lines in whlDataLfp:' int2str(size(whlDataLfp,1))]);

    %save it
    disp('save data')
    save([file '-whl.mat'], 'whlDataLfp','interpTsp','behType',...
         'trStartLfpIndArr','trEndLfpIndArr', 'trackLenArr', ...
         'interpTspF','interpTspFneu','interpTspSpks','interpTspdFF','-v7.3');
end
