function [sp ac] = MazeSpeedAccel_e(whldata, varargin)

% function [speed, accel] = MazeSpeedAccel(whldata,smoothlen)
% Calculates running speed and acceleration from XY position data. 

[SamplRate, maxSpeed, minAccel, maxAccel, viewSpeedMap, viewAccelMap] ...
    = DefaultArgs(varargin,{1250, 100, -200, 200, 0, 0});

smoothlen = round(SamplRate / 3);

% filter with a hanning window
findex = find(whldata(:,1)>0);
hanfilter = hanning(smoothlen);
hanfilter = hanfilter./sum(hanfilter);

whldata(findex,1) = Filter0(hanfilter,whldata(findex,1));
whldata(findex,2) = Filter0(hanfilter,whldata(findex,2));

%calculate speed for values that aren't -1 or distorded by filtering
%---------------------------------
intval = findex(ceil(smoothlen/2):end-ceil(smoothlen/2));
speeddata(:,1) = diff(whldata(intval,1)) + 1i*diff(whldata(intval,2));
Amp = abs(speeddata)*SamplRate;
Phase = mod(angle(speeddata),2*pi);
speed = ones(size(whldata));
speed(intval(1:end-1),:) = [Amp Phase];
sp = speed(:,1);
ac = speed(:,2);

%---------------------------------

if viewSpeedMap == 1        
    ViewParamMap(whldata(findex,1:2), sp(findex), 'gray');        
end
if viewAccelMap == 1    
    ViewParamMap(whldata(findex,1:2), ac(findex), 'gray');
end

return

% findex = find(whldata(:,1)>0);
% [findex_m n] = size(findex);
% beginnotdistorted = ceil(smoothlen/2);
% endnotdistorted = findex_m-ceil(smoothlen/2);
% % speeddata starts one sample later due to loss by diff
% speeddata(findex(beginnotdistorted+1:endnotdistorted), 1) = diff(whldata(findex(beginnotdistorted:endnotdistorted),1));
% speeddata(findex(beginnotdistorted+1:endnotdistorted), 2) = diff(whldata(findex(beginnotdistorted:endnotdistorted),2));
% [whl_m n] = size(whldata);
% speeddata2d = -1*ones(whl_m, 1);
% speeddata2d(findex(beginnotdistorted+1:endnotdistorted)) = ...
% sqrt(speeddata(findex(beginnotdistorted+1:endnotdistorted),1).^2+speeddata(findex(beginnotdistorted+1:endnotdistorted),2).^2);
% findex = find(speeddata2d(:)~=-1);
% speed = -1*ones(whl_m, 1);
% speed(findex) = speeddata2d(findex).*SamplRate;
% speed(speed > maxSpeed) = maxSpeed;
% [findex_m n] = size(findex);
% accel = -1*ones(whl_m, 1);
% accel(findex(1:findex_m-1)) = diff(speed(findex(1:findex_m))).*SamplRate; % account for lost datapoint by shifting findex position back
% accel(accel < minAccel) = minAccel;
% accel(accel > maxAccel) = maxAccel;
