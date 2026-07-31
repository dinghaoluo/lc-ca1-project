% Measure of cluster Mahalanobis distance
% [IsolationDist] = Cluster_MahalDist(Fet, allGoodSpikes, MySpikes, Res)
%
% originally from ClusterQuality_e (use the original fce if you want Mahal for bursting spikes)
%
% Inputs: Fet - N by D array of feature vectors (N spikes, D dimensional feature space)
% MySpikes: list of spikes corresponding to cell whose quality is to be evaluated.
% Res - Spike times
%

function [IsolationDist] = Cluster_MahalDist(Fet, MySpikes, Res)

% check there are enough spikes (but not more than half)
if sum(MySpikes) < size(Fet,2) | sum(MySpikes)>length(Res)/2
	IsolationDist = -1;	
	return
end

m = mahal(Fet, Fet(MySpikes,:));
mNoise = m(~MySpikes);
mClu = m(MySpikes);

% based on Harris et al., 2001 and Schmitzer-Torbert et al., 2005: 
% Isolation Distance = a point where mahalDist of other spikes = n of this cell
nCluSpikes = sum(MySpikes);
nSpikes = length(Res);
if (nCluSpikes < nSpikes/2)
	[sortedNoise order] = sort(mNoise);
    [sortedClu order] = sort(mClu);
    aboveCross = find(sortedClu>sortedNoise(1:length(sortedClu)));
    if length(aboveCross) > 0        
        IsolationDist = sortedNoise(aboveCross(1));    
    else
        IsolationDist = max(mClu);
    end
    
else
	IsolationDist = 0; % If there are more of this cell than every thing else, forget it.
end

return;
