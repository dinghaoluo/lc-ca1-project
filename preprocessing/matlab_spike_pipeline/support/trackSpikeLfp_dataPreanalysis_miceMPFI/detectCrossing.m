function [ind,t0close,s0close] = detectCrossing(S,t,level)
% CROSSING find the crossings of a given level of a signal
%   ind = CROSSING(S) returns an index vector ind, the signal
%   S crosses zero at ind or at between ind and ind+1
%   [ind,t0close] = CROSSING(S,t) additionally returns a time
%   vector t0 of the zero crossings of the signal S. The crossing
%   times are linearly interpolated between the given times t
%   [ind,t0close] = CROSSING(S,t,level) returns the crossings of the
%   given level instead of the zero crossings
%   ind = CROSSING(S,[],level) as above but without time interpolation
%
%	[ind,t0close,s0close] = ... also returns the data vector corresponding to 
%	the t0 values.

% check the number of input arguments
error(nargchk(1,3,nargin));

% check the time vector input for consistency
if nargin < 2 || isempty(t)
	% if no time vector is given, use the index vector as time
    t = 1:length(S);
elseif length(t) ~= length(S)
	% if S and t are not of the same length, throw an error
    error('t and S must be of identical length!');    
end

% check the level input
if nargin < 3
	% set standard value 0, if level is not given
    level = 0;
end

% make row vectors
t = t(:)';
S = S(:)';

% always search for zeros. So if we want the crossing of 
% any other threshold value "level", we subtract it from
% the values and search for zeros.
S   = S - level;

% first look for exact zeros
ind0 = find( S == 0 ); 

% then look for zero crossings between data points
S1 = S(1:end-1) .* S(2:end);
ind1 = find( S1 < 0 );

% bring exact zeros and "in-between" zeros together 
ind = sort([ind0 ind1]);

% and pick the associated time values
t0 = t(ind); 
s0 = S(ind);

% Addition:
% Some people like to get the data points closest to the zero crossing,
% so we return these as well
if ind(1) == 1; 
    [CC,II] = min(abs([S(ind) ; S(ind+1)]),[],1); 
    ind2 = ind + (II-1); %update indices 
else
    [CC,II] = min(abs([S(ind-1) ; S(ind) ; S(ind+1)]),[],1); 
    ind2 = ind + (II-2); %update indices 
end

t0close = t(unique(ind2));
s0close = S(unique(ind2)) + level;