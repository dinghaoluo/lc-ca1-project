% Crossings = SchmittTriggerUpDownMarked(Signal,UpThresh,DownThresh)
%
% finds all those points where the signal crosses UpThresh in the
% upwards direction, but only for the first time after each time
% it was below DownThresh

function [UpCrossings DownCrossings] = TriggerUpDownMarked(Signal,UpThresh,DownThresh)

% this is the sort of thing that would be so much easier not vectorized

if prod(size(Signal))~=max(size(Signal))
	error('Can only take a vector input');
end
Signal = Signal(:);

n = length(Signal);
% PrevVal = [(UpThresh+DownThresh)/2; Signal(1:n-1)];

UpCrossings = find([Signal(1); Signal(1:n-1)]<UpThresh & Signal>=UpThresh);
DownCrossings = find([Signal(1); Signal(1:n-1)]>DownThresh & Signal<=DownThresh);
UpCrossings = find([(UpThresh+DownThresh)/2; Signal(1:n-1)]<UpThresh & Signal>=UpThresh);
DownCrossings = find([(UpThresh+DownThresh)/2; Signal(1:n-1)]>DownThresh & Signal<=DownThresh);

if(UpCrossings(1)~=1)
    UpCrossings = UpCrossings((Signal(UpCrossings)-Signal(UpCrossings-1)) > abs(UpThresh-DownThresh)/2);
else
    UpCrossings = UpCrossings((Signal(UpCrossings(2:end))-Signal(UpCrossings(2:end)-1)) > abs(UpThresh-DownThresh)/2);
    if((Signal(1)-(UpThresh+DownThresh)/2) > abs(UpThresh-DownThresh))
        UpCrossings = [1; UpCrossings];
    end
end

if(DownCrossings(1)~=1)
    DownCrossings = DownCrossings((Signal(DownCrossings-1)-Signal(DownCrossings)) > abs(UpThresh-DownThresh)/2);
else
    DownCrossings = DownCrossings((Signal(DownCrossings(2:end)-1)-Signal(DownCrossings(2:end))) > abs(UpThresh-DownThresh)/2);
    if(((UpThresh+DownThresh)/2-Signal(1)) > abs(UpThresh-DownThresh))
        DownCrossings = [1; DownCrossings];
    end

end