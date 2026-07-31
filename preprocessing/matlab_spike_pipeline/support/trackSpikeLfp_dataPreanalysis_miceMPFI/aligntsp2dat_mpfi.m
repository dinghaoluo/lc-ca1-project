function aligntsp2dat_mpfi(file,sampleFreq,lfpFreq,ChanNum)

if exist([file '-whl.mat'], 'file') == 2
    disp('whl file already exists.')
    return;
else
    %% load the length of each segment of the treadmill belt
    if exist([file '-segLen.mat'], 'file') == 2
        load([file '-segLen.mat']);
    else
        segLenDef=input('Do you want to use the default inter-beambreak distance values? Y/N [Y]: ','s');
        if(isempty(segLenDef))
            segLenDef = 'Y';
        end
        if(~isempty(strfind(segLenDef,'Y')) || ~isempty(strfind(segLenDef,'y')))
            % belt 1 (with dense marks)
%             seg1Len = 612.775; %mm
%             seg2Len = 608.0125; %mm
%             seg3Len = 587.375; %mm
            % belt 2 (with sparse marks) 
            seg1Len = 609.60; %mm
            seg2Len = 609.60; %mm
            seg3Len = 609.60; %mm
        else
            inputs = input('Please enter three inter-beambreak distance values [a b c]: ');
            seg1Len = inputs(1);
            seg2Len = inputs(2);
            seg3Len = inputs(3);
        end
        save([file '-segLen.mat'],'seg1Len','seg2Len','seg3Len');
    end
    
    %% get start and end timestamp of dat file
    fid=fopen([file '.meta']);
    tline= fgetl(fid);
    while ischar(tline)
        if startsWith(tline,'TimeStamp of the end')
            tline=tline(59:end);
            EndTimestamp=sscanf(tline,'%d',1);
        end
        if startsWith(tline,'TimeStamp of the start')
            tline=tline(61:end);
            StartTimestamp=sscanf(tline,'%d',1);
        end
        if startsWith(tline,'File size')
            tline=tline(21:end);
            DatSize=sscanf(tline,'%lu',1);
            DatSize=double(DatSize);
        end

        tline= fgetl(fid);
    end
    fclose(fid);
    
    DatLength=DatSize/(ChanNum*2); %Dat file size in ms
    
    %% initiate variables
    if EndTimestamp<0 & StartTimestamp<0
        disp('META file time stamps are corrupted.');
        timeStep = 1/sampleFreq;
        timeStepLfp = 1/lfpFreq;
        interpTsp(:,1) = (0:timeStep:(DatLength-1)/sampleFreq)*1000; % in ms
        whlDataLfp(:,1) = (0:timeStepLfp:(DatLength-1)/sampleFreq)*1000; % in ms
        StartTimestamp = 0;
    else
        timeStep = (EndTimestamp-StartTimestamp)/(DatLength-1); % in ms
        timeStepLfp = timeStep*sampleFreq/lfpFreq;
        interpTsp(:,1) = 0:timeStep:EndTimestamp-StartTimestamp;
        whlDataLfp(:,1) = 0:timeStepLfp:EndTimestamp-StartTimestamp;
    end
    interpTspDist = zeros(length(interpTsp),1);
    interpTspDistLfp = zeros(length(whlDataLfp),1);
    
    %% get timestamps and positions from BTDT file
    load([file 'BTDT.mat']);
    
    behType = zeros(max(behEventsTdt.trialDescr(:,3)),1);
    trStartLfpIndArr = zeros(max(behEventsTdt.trialDescr(:,3)),1);
    trEndLfpIndArr = zeros(max(behEventsTdt.trialDescr(:,3)),1);
    
    % remove any treadmill event before the beginning of the first trial     
    indTmp = find(behEventsTdt.wheel(:,1) >= behEventsTdt.trialDescr(1,1));
    wheelDist = behEventsTdt.wheel(indTmp,4);
    wheelTime = behEventsTdt.wheel(indTmp,1);
    wheelTimeLfp = behEventsTdt.wheel(indTmp,2);
    
    % find all the beam break events starting immediately before the
    % beginning of the first trial
    offsetBeam = 0.01*sampleFreq; % no. samples
    indBeam = behEventsTdt.beam(:,1) > behEventsTdt.trialDescr(1,1) ...
                - offsetBeam & behEventsTdt.beam(:,4) == 0; 
    beam = behEventsTdt.beam(indBeam,:);
    
    %% interpolate the running distance for each trial
    for tr = 1:max(behEventsTdt.trialDescr(:,3))
        indBeamCurTrTmp = find(beam(:,1) > behEventsTdt.trialDescr(tr,1) ...
            - offsetBeam & beam(:,1) < behEventsTdt.trialDescr(tr+1,1) + offsetBeam);
        % get the last breaking time of the first beam break
        diffBeam = [0; diff(beam(indBeamCurTrTmp,3))];
        ind = find(cumsum(diffBeam) ~= 0,1);
        indBeamCurTrTmp = indBeamCurTrTmp(ind-1:end);

        if(isempty(indBeamCurTrTmp)) 
            % run backwards and lead to a trial without correct beam break
            disp(['Trial no. ' num2str(tr) ' has ' num2str(length(indBeamCurTrTmp))...
                ' beam breaks. Therefore, it is considered as a bad trial']);
            behType(tr) = -1; %% wrong number of beam breaks
            trStartLfpIndArr(tr) = behEventsTdt.trialDescr(tr,2);
            trEndLfpIndArr(tr) = behEventsTdt.trialDescr(tr+1,2)-1;
            continue;
        else
            % check whether there are multiple triggering of the same beam break
            indBeamCurTr = indBeamCurTrTmp(1);
            trStartLfpIndArr(tr) = beam(indBeamCurTr,2);
            trEndLfpIndArr(tr) = behEventsTdt.trialDescr(tr+1,2)-1;
            for i = 1:length(indBeamCurTrTmp)-1
                if(beam(indBeamCurTrTmp(i+1),3) ~= beam(indBeamCurTrTmp(i),3))
                    indBeamCurTr = [indBeamCurTr indBeamCurTrTmp(i+1)];
                end
            end
        end
        
        disp(tr);
        if(length(indBeamCurTr) == 4 ...
                & sum(beam(indBeamCurTr,3) == [1; 2; 3; 1]) == 4)
            % changed by Yingxue on 08/21/2018, since the order of the beam
            % break should have already been corrected in the end of
            % Arduino2TDTtime_mpfi, here we should only consider the
            % correct beam break order (1 2 3)
            behType(tr) = 1;
            for i = 1:length(indBeamCurTr)-1
                
                %% find the treadmill events between two beam breaks
                indTmp = find(wheelTime > beam(indBeamCurTr(i),1) ...
                    & wheelTime < beam(indBeamCurTr(i+1),1));
                beamEndT = beam(indBeamCurTr(i+1),1);
                beamStartT = beam(indBeamCurTr(i),1);
                beamEndTLfp = beam(indBeamCurTr(i+1),2);
                beamStartTLfp = beam(indBeamCurTr(i),2);

                %% interpolate the distance of the treadmill events based 
                % on the distance between two beam breaks
                wheelSegDist = wheelDist(indTmp);
                % check whether the distance is monotonically increasing
                diffDist = diff(wheelSegDist);
                ind = find(diffDist < -500,1);
                if(~isempty(ind))
                    indTmp = indTmp(ind+1:end);
                    wheelSegDist = wheelDist(indTmp);
                end
                wheelSegTime = wheelTime(indTmp);
                wheelSegTimeLfp = wheelTimeLfp(indTmp);
                if(beam(indBeamCurTr(i),3) == 1) % the first segment of the treadmill
                    distStart = 0;
                    distEnd = wheelSegDist(end) + (wheelSegDist(end) ...
                        - wheelSegDist(end-1))*(beamEndT - wheelSegTime(end))...
                        /(wheelSegTime(end) - wheelSegTime(end-1));
                    xBeh = interp1([distStart distEnd],[0 seg1Len],...
                        [distStart; wheelSegDist; distEnd],'linear');
                elseif(beam(indBeamCurTr(i),3) == 2) % the second segment of the treadmill    
                    distStart = wheelSegDist(1) - (wheelSegDist(2) ...
                        - wheelSegDist(1))*(wheelSegTime(1) - beamStartT)...
                        /(wheelSegTime(2) - wheelSegTime(1));
                    distEnd = wheelSegDist(end) + (wheelSegDist(end) ...
                        - wheelSegDist(end-1))*(beamEndT - wheelSegTime(end))...
                        /(wheelSegTime(end) - wheelSegTime(end-1));
                    xBeh = interp1([distStart distEnd],[seg1Len seg1Len+seg2Len],...
                        [distStart; wheelSegDist; distEnd],'linear');
                elseif(beam(indBeamCurTr(i),3) == 3) % the third segment of the treadmill
                    distStart = wheelSegDist(1) - (wheelSegDist(2) ...
                        - wheelSegDist(1))*(wheelSegTime(1) - beamStartT)...
                        /(wheelSegTime(2) - wheelSegTime(1));
                    distEnd = wheelSegDist(end) + (wheelSegDist(end) ...
                        - wheelSegDist(end-1))*(beamEndT - wheelSegTime(end))...
                        /(wheelSegTime(end) - wheelSegTime(end-1));
                    xBeh = interp1([distStart distEnd],...
                        [seg1Len+seg2Len seg1Len+seg2Len+seg3Len],...
                        [distStart; wheelSegDist; distEnd],'linear');
                end
%                 disp(['Trial ' num2str(tr) ' beam ' num2str(i)])
                xBeh2Time = interp1([beamStartT; wheelSegTime; beamEndT],...
                    xBeh,beamStartT:beamEndT,'linear');
                if(beamStartTLfp == wheelSegTimeLfp(1))
                    wheelSegTimeLfp(1) = wheelSegTimeLfp(1) + 0.0001;
                end
                if(beamEndTLfp == wheelSegTimeLfp(end))
                    wheelSegTimeLfp(end) = wheelSegTimeLfp(end) - 0.0001;
                end
                xBeh2TimeLfp = interp1([beamStartTLfp; wheelSegTimeLfp; beamEndTLfp],...
                    xBeh,beamStartTLfp:beamEndTLfp,'linear');

                interpTspDist(beamStartT:beamEndT-1) = xBeh2Time(1:end-1);
                interpTspDistLfp(beamStartTLfp:beamEndTLfp-1) = xBeh2TimeLfp(1:end-1);
            end
        elseif(length(indBeamCurTr) == 3 ...
                & (sum(beam(indBeamCurTr,3) == [1; 2; 1]) == 3 ...
                | sum(beam(indBeamCurTr,3) == [1; 3; 1]) == 3)) 
            % added by Yingxue on 08/21/2018, to account for the situation
            % where one beam break is consistently missing
            behType(tr) = 1;
            for i = 1:length(indBeamCurTr)-1
                %% find the treadmill events between two beam breaks
                indTmp = find(wheelTime > beam(indBeamCurTr(i),1) ...
                    & wheelTime < beam(indBeamCurTr(i+1),1));
                beamEndT = beam(indBeamCurTr(i+1),1);
                beamStartT = beam(indBeamCurTr(i),1);
                beamEndTLfp = beam(indBeamCurTr(i+1),2);
                beamStartTLfp = beam(indBeamCurTr(i),2);

                %% interpolate the distance of the treadmill events based 
                % on the distance between two beam breaks
                wheelSegDist = wheelDist(indTmp);
                % check whether the distance is monotonically increasing
                diffDist = diff(wheelSegDist);
                ind = find(diffDist < -500,1);
                if(~isempty(ind))
                    indTmp = indTmp(ind+1:end);
                    wheelSegDist = wheelDist(indTmp);
                end
                wheelSegTime = wheelTime(indTmp);
                wheelSegTimeLfp = wheelTimeLfp(indTmp);
                if(beam(indBeamCurTr(i),3) == 1) % the first segment of the treadmill
                    distStart = 0;
                    distEnd = wheelSegDist(end) + (wheelSegDist(end) ...
                        - wheelSegDist(end-1))*(beamEndT - wheelSegTime(end))...
                        /(wheelSegTime(end) - wheelSegTime(end-1));
                    if(~isempty(beam(indBeamCurTr,3) == 2))
                        distEndRef = seg1Len;
                    elseif(~isempty(beam(indBeamCurTr,3) == 3))
                        distEndRef = seg1Len+seg2Len;
                    end
                    xBeh = interp1([distStart distEnd],[0 distEndRef],...
                        [distStart; wheelSegDist; distEnd],'linear');
                else % the second segment of the treadmill    
                    distStart = wheelSegDist(1) - (wheelSegDist(2) ...
                        - wheelSegDist(1))*(wheelSegTime(1) - beamStartT)...
                        /(wheelSegTime(2) - wheelSegTime(1));
                    distEnd = wheelSegDist(end) + (wheelSegDist(end) ...
                        - wheelSegDist(end-1))*(beamEndT - wheelSegTime(end))...
                        /(wheelSegTime(end) - wheelSegTime(end-1));
                    if(~isempty(beam(indBeamCurTr,3) == 2))
                        distStartRef = seg1Len;
                    elseif(~isempty(beam(indBeamCurTr,3) == 3))
                        distStartRef = seg1Len+seg2Len;
                    end
                    xBeh = interp1([distStart distEnd],...
                        [distStartRef seg1Len+seg2Len+seg3Len],...
                        [distStart; wheelSegDist; distEnd],'linear');

                end
%                 disp(['Trial ' num2str(tr) ' beam ' num2str(i)])
                xBeh2Time = interp1([beamStartT; wheelSegTime; beamEndT],...
                    xBeh,beamStartT:beamEndT,'linear');
                if(beamStartTLfp == wheelSegTimeLfp(1))
                    wheelSegTimeLfp(1) = wheelSegTimeLfp(1) + 0.0001;
                end
                if(beamEndTLfp == wheelSegTimeLfp(end))
                    wheelSegTimeLfp(end) = wheelSegTimeLfp(end) - 0.0001;
                end
                xBeh2TimeLfp = interp1([beamStartTLfp; wheelSegTimeLfp; beamEndTLfp],...
                    xBeh,beamStartTLfp:beamEndTLfp,'linear');

                interpTspDist(beamStartT:beamEndT-1) = xBeh2Time(1:end-1);
                interpTspDistLfp(beamStartTLfp:beamEndTLfp-1) = xBeh2TimeLfp(1:end-1);
            end
        else
            disp(['Trial no. ' num2str(tr) ' has ' num2str(length(indBeamCurTr))...
                ' beam breaks. Therefore, it is considered as a bad trial']);
            behType(tr) = -1; %% wrong number of beam breaks
        end
    end
    
    interpTsp(:,2) = interpTspDist;
    whlDataLfp(:,2) = interpTspDistLfp;
        
    disp(['Samples in .dat file per channel: ' int2str(DatLength)]);
    disp(['Lines in interpTsp:' int2str(size(interpTsp,1))]);
    disp(['Lines in whlDataLfp:' int2str(size(whlDataLfp,1))]);

    %save it
    save([file '-whl.mat'], 'whlDataLfp','interpTsp','behType',...
         'trStartLfpIndArr','trEndLfpIndArr');
end
