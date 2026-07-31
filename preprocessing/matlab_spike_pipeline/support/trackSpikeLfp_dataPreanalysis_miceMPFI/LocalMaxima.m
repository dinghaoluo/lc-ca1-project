%[ maxsind, maxsvalues] = LocalMaxima(x)
%
% finds positions of all strict local maxima in input array

function [maxsind, maxsvalues] = LocalMaxima(x)

nPoints = length(x);

Middle = x(2:(nPoints-1));
Left = x(1:(nPoints-2));
Right = x(3:nPoints);

maxsind = 1+find(Middle > Left & Middle > Right);
maxsvalues = x(maxsind);