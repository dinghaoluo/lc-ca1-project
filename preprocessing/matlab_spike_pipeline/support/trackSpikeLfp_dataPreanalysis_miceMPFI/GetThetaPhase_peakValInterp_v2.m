function thetaLinear = GetThetaPhase_peakValInterp_v2(eeg, lfpSampleRate, thetaFiltParam, numStartSmp)
% determine theta phase from the eeg trace
% eeg:  the eeg trace
% lfpSampleRate: lfp sample frequency
% thetaFiltParam:
% thetaFiltParam.FreqRange = [4 15];
% thetaFiltParam.FilterOrd = 4;
% thetaFiltParam.Ripple = 20;
% numStartSmp: number of samples which should be added to the beginning of
% eeg

    % filter theta band and find local minimum
    [b a] = cheby2(thetaFiltParam.FilterOrd, thetaFiltParam.Ripple, ...
                thetaFiltParam.FreqRange/(lfpSampleRate/2),'bandpass');
    eegf = filtfilt(b,a,eeg);
    eegf = eegf - mean(eegf);
    
    [minX, minY] = LocalMaxima(eegf .* (-1));       
        % theta peaks and valleys
    minY = minY.* (-1);
    [maxX, maxY] = LocalMaxima(eegf);
    
    % remove peaks/vallyes with too small amplitude => 
    % these do appear as very short cycles
    %  use both, peaks and vallyes  
    minCycleLen = 100;      % number of samples per half cycle
    cycleLen = 150;
    pp = sortrows([[minX; maxX] [minY; maxY]], 1);            
    minPP = abs(min(pp(:,2)));
    pp(:,2) = pp(:,2) + minPP;
    minCycAmpl = std(pp(:,2))*0.3;      % bytes
    smallPP = find(abs(diff(pp(:,2))) <  minCycAmpl);
    if(~isempty(smallPP)) 
        % at this moment, smallPP corresponding to the indices of intervals,
        % not the indices of points, here makes the conversion
        smallPP = [smallPP(1); smallPP+1];
     
        if(smallPP(end) < size(pp,1))
            indLongCycles = pp(smallPP+1,1) - pp(smallPP,1) > minCycleLen; 
                % find the half cycles with small amplitude change but large time span
        else
            indLongCycles = pp(smallPP(1:end-1)+1,1) - pp(smallPP(end-1),1)...
                            > minCycleLen;
        end
        smallPP(indLongCycles) = []; 
        % if the half cycle is long, consider it as a half cycle regardless
        % of the amplitude change

        % collect all the peak points to be removed 
        goodSmall = [];
        n = 2;
        while(n <= size(smallPP,1))
            if(isempty(goodSmall))
                 goodSmall = smallPP(n-1);
            end
            if smallPP(n) ~= goodSmall(end) + 1  
                % if the next short interval is not directly connected with
                % the previous one, then add the corresponding point into
                % to be removed list
                goodSmall = [goodSmall; smallPP(n)];
                n = n + 1;
            else
                % if the next short interval continues to be a short
                % interval
                indNonContPP = find(diff(smallPP(n-1:end)) > 1, 1); 
                    % find the next long interval
                if(~isempty(indNonContPP))  
                    if(n == 2) % the first cloud of points with small intervals
                        smallContPP = smallPP(n-1:n+indNonContPP-2);
                    else
                        % if not the first cloud of points, then including 
                        % the end point of the last long interval into the current cloud
                        smallContPP = [smallPP(n-1)-1; ...
                                       smallPP(n-1:n+indNonContPP-2)];
                    end
                    indLongCycle = 1:length(smallContPP); 
                        % indices of all long cycles
                    while(~isempty(indLongCycle)) 
                        [maxContS,maxContSInd] = max(abs(pp(smallContPP,2)-minPP));
                        % find the peak point among the short interval cloud
                        curPP = smallContPP(maxContSInd(1));    
                        % set the current point as the peak point, 
                        % this is the point to be kept
                        if(isempty(goodSmall))
                            goodSmall = smallContPP(smallContPP < curPP); 
                        else
                            goodSmall = [goodSmall(1:end-1); ...
                                         smallContPP(smallContPP < curPP)];     
                        end
                            % include all the points in the current cloud 
                            % which are before the curPP into to be removed list
                        smallContPP = smallContPP(smallContPP > curPP);     
                            % keep all the point after the curPP 
                        if(~isempty(smallContPP))  
                            % check for long cycle points
                            indLongCycle = find(pp(smallContPP,1) - pp(curPP,1)...
                                                > minCycleLen); 
                                % points further away from the curPP
                            if(isempty(indLongCycle))
                                goodSmall = [goodSmall; smallContPP];  
                                    % if no long cycle points, then include
                                    % all points after the curPP into to be removed list
                            else
                                % else only delete points with small distance from curPP
                                goodSmall = [goodSmall; ...
                                    smallContPP(pp(smallContPP,1) - pp(curPP,1) <= minCycleLen)];
                                if(isempty(goodSmall))
                                    smallContPP = smallContPP(indLongCycle); 
                                    % keep the long distance points
                                else
                                    smallContPP = [goodSmall(end); ...
                                        smallContPP(indLongCycle)]; 
                                    % keep the long distance points
                                end
                            end
                        else
                            indLongCycle = [];
                        end
                    end
                    n = n + indNonContPP - 1;
                else
                    % reach the end of the pp array, add all the small 
                    % interval points to the list
                    indContPP = find(diff(smallPP(n-1:end)));
                    if(~isempty(indContPP))
                        goodSmall = [goodSmall; smallPP(n-1+indContPP)];
                        n = n + length(indContPP);
                    else
                        n = n + 1;
                    end
                end
            end
        end
        
        % delete small interval points
        if ~isempty(goodSmall)
            pp(goodSmall,:) = [];
        end
        
        % check whether two adjacent peak points have the same polarity, if
        % so, either delete the low amplitude one or add in an opposite
        % polarity peak in between
        ppComp = []; 
        for i = size(pp,1)-1:-1:1
            eqSignPP1Amp = pp(i,2)-minPP;
            eqSignPP2Amp = pp(i+1,2)-minPP;
            if(eqSignPP1Amp*eqSignPP2Amp > 0)
                if(pp(i+1,1) - pp(i,1) > cycleLen)
                    if(eqSignPP1Amp > 0)
                        [troughTmp,indTrough] = min(eegf(pp(i,1):pp(i+1,1)));
                        ppComp = [ppComp; pp(i,1)+indTrough-1, troughTmp+minPP];
                    else
                        [peakTmp,indPeak] = max(eegf(pp(i,1):pp(i+1,1)));
                        ppComp = [ppComp; pp(i,1)+indPeak-1, peakTmp+minPP];
                    end
                else
                    if(abs(pp(i,2)-minPP) > abs(pp(i+1,2)-minPP))
                        pp(i+1,:) = [];
                    else
                        pp(i,:) = [];
                    end 
                end
            end
        end
        pp = [pp; ppComp];
        [pp(:,1),indPP] = sort(pp(:,1));
        pp(:,2) = pp(indPP,2);
    end
    
    % find all crossing point of the middle of each half cycle segment
    zeroCross = [];
    t0Close = [];
    for i = 1:size(pp,1)-1
        if length(eegf(pp(i,1):pp(i+1,1))) > 3
        [indCross,t0Close,s0Close] = ...
            detectCrossing(eegf(pp(i,1):pp(i+1,1)),pp(i,1):pp(i+1,1),...
                            (eegf(pp(i,1))+eegf(pp(i+1,1)))/2);
        end
        if(length(t0Close) > 1)
            compMidTime = abs(t0Close - (pp(i,1)+pp(i+1,1))/2);
            indCross = find(compMidTime == min(compMidTime));
            t0Close = t0Close(indCross(1));
            s0Close = s0Close(indCross(1));
        end
        zeroCross(i,:) = [t0Close s0Close];
    end
    
    if(~isempty(pp))
        if(pp(1,2) - minPP <= zeroCross(1,2)) % start with trough
            thetaLinear.thetaPeak_tAmpl(:,1) = pp(2:2:end,1) + numStartSmp;
            thetaLinear.thetaPeak_tAmpl(:,2) = pp(2:2:end,2)- minPP;
            thetaLinear.thetaTrough_tAmpl(:,1) = pp(1:2:end,1) + numStartSmp;
            thetaLinear.thetaTrough_tAmpl(:,2) = pp(1:2:end,2) - minPP;
            thetaLinear.thetaPtoTZeros_tAmpl(:,1) = zeroCross(2:2:end,1) + numStartSmp;
            thetaLinear.thetaPtoTZeros_tAmpl(:,2) = zeroCross(2:2:end,2);
            thetaLinear.thetaTtoPZeros_tAmpl(:,1) = zeroCross(1:2:end,1) + numStartSmp;
            thetaLinear.thetaTtoPZeros_tAmpl(:,2) = zeroCross(1:2:end,2);
        else % start with peak
            thetaLinear.thetaPeak_tAmpl(:,1) = pp(1:2:end,1) + numStartSmp;
            thetaLinear.thetaPeak_tAmpl(:,2) = pp(1:2:end,2) - minPP;
            thetaLinear.thetaTrough_tAmpl(:,1) = pp(2:2:end,1) + numStartSmp;
            thetaLinear.thetaTrough_tAmpl(:,2) = pp(2:2:end,2) -minPP;
            thetaLinear.thetaPtoTZeros_tAmpl(:,1) = zeroCross(1:2:end,1) + numStartSmp;
            thetaLinear.thetaPtoTZeros_tAmpl(:,2) = zeroCross(1:2:end,2);
            thetaLinear.thetaTtoPZeros_tAmpl(:,1) = zeroCross(2:2:end,1) + numStartSmp;
            thetaLinear.thetaTtoPZeros_tAmpl(:,2) = zeroCross(2:2:end,2);
        end
    else
        thetaLinear.thetaPeak_tAmpl = [];
        thetaLinear.thetaTrough_tAmpl = [];
        thetaLinear.thetaPtoTZeros_tAmpl = [];
        thetaLinear.thetaTtoPZeros_tAmpl = [];
    end

    % step-wise interpolation between valeys to match Hilbert (cosyne: -pi->zero->pi)
    if(~isempty(thetaLinear.thetaTrough_tAmpl))
        startInterp = [thetaLinear.thetaTrough_tAmpl(:,1) ...
                    -1*pi*ones(size(thetaLinear.thetaTrough_tAmpl,1),1)]; 
                % theta through corresponding to -pi and pi(-180 and 180deg)
        endInterp = [thetaLinear.thetaTrough_tAmpl(:,1)-1 ...
                pi*ones(size(thetaLinear.thetaTrough_tAmpl,1),1)];
    else
        startInterp = [];
        endInterp = [];
    end
    if(~isempty(thetaLinear.thetaTtoPZeros_tAmpl))
        interp90 = [thetaLinear.thetaTtoPZeros_tAmpl(:,1) ...
            -0.5*pi*ones(size(thetaLinear.thetaTtoPZeros_tAmpl,1),1)]; 
            % trough to peak zero crossing corresponding to -0.5*pi (-90 degree)
    else
        interp90 = [];
    end
    if(~isempty(thetaLinear.thetaPeak_tAmpl))
        interp180 =  [thetaLinear.thetaPeak_tAmpl(:,1) ...
            zeros(size(thetaLinear.thetaPeak_tAmpl,1),1)];     
            % theta peak corresponding to 0 degree         
    else
        interp180 = [];
    end
    if(~isempty(thetaLinear.thetaPtoTZeros_tAmpl))  
        interp270 = [thetaLinear.thetaPtoTZeros_tAmpl(:,1) ...
            0.5*pi*ones(size(thetaLinear.thetaPtoTZeros_tAmpl,1),1)]; 
            % peak to trough zero crossing corresponding to 0.5*pi (90 degree)
    else
        interp270 = [];
    end

    allInterp = [startInterp; interp90; interp180; interp270]; 
    thetaLin = [];

    if(size(allInterp,1) > 1)
        allInterp = sortrows(allInterp, 1);
        if(~isempty(find(diff(allInterp(:,1)) <= 0,1)))
            disp('The intervals to be interpolated are not monotonically increasing');
%             pause;
        end
        [tmp,indTmp,indStInterp] = intersect(startInterp(:,1),allInterp(:,1));
        if(~isempty(indStInterp))
            if(allInterp(1,1) ~= 1) 
                % interpolate a constant value before the first peak/trough
                time = 1:allInterp(1,1);
                phase = allInterp(1,2)*ones(1,length(time));
                thetaLin = phase;
            end
            
            if(indStInterp(1) ~= 1) 
                % interpolate the segment before the first peak
                time = allInterp(1:indStInterp(1),1);
                phase = allInterp(1:indStInterp(1),2);
                phase(end) = phase(end) + 2*pi;
                repT = find(diff(time)==0);
                if ~isempty(repT)
                    time(repT+1) = [];
                    phase(repT+1) = [];
                end
                interpTmp = interp1q(time,phase,(allInterp(1,1):startInterp(1,1))');
                interpTmp(isnan(interpTmp)) = 0;
                thetaLin = [thetaLin(1:end-1) interpTmp'];
            end
            
            for i = 1:length(startInterp(:,1))-1    
                % interpolate between startInterp one by one
                time = allInterp(indStInterp(i):indStInterp(i+1),1);
                phase = allInterp(indStInterp(i):indStInterp(i+1),2);
                phase(end) = phase(end) + 2*pi;
                repT = find(diff(time)==0);
                if ~isempty(repT)
                    time(repT+1) = [];
                    phase(repT+1) = [];
                end
                interpTmp = interp1q(time,phase,(startInterp(i,1):startInterp(i+1,1))');
                interpTmp(isnan(interpTmp)) = 0;
                thetaLin = [thetaLin(1:end-1) interpTmp'];
            end
            
            if(startInterp(end,1) ~= allInterp(end,1))  
                % the last startInterp to the end: interpolate with a constant value
                time = allInterp(indStInterp(end):end,1);
                phase = allInterp(indStInterp(end):end,2);
                repT = find(diff(time)==0);
                if ~isempty(repT)
                    time(repT+1) = [];
                    phase(repT+1) = [];
                end
                interpTmp = ...
                    interp1q(time,phase,(startInterp(end,1):length(eeg)+numStartSmp)');
                interpTmp(isnan(interpTmp)) = 0;
                thetaLin = [thetaLin(1:end-1) interpTmp'];
            end
            
            if(length(thetaLin) ~= length(eeg)+numStartSmp)
                thetaLin = [thetaLin(1:end-1) ...
                    zeros(1,length(eeg) + numStartSmp - length(thetaLin) + 1)];
            end
        else
            thetaLin = interp1(allInterp(:,1),allInterp(:,2),...
                1:length(eeg)+numStartSmp,'linear',mean(allInterp(:,2)));
        end
    else
        thetaLin = zeros(1,length(eeg)+numStartSmp);
    end
    
    thetaLinear.thetaLin = thetaLin';  
    
    return;