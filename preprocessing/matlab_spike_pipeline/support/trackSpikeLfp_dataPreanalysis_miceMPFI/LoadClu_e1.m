% [nClusters, Clu] = LoadClu_e1(FileName, shank);
%
% A simple matlab function to load a .clu file

function [nClusters, Clu] = LoadClu_e1(FileName, varargin);

[shank] = DefaultArgs(varargin,{0});

if shank == 0
    cluFileName = FileName;
else
    cluFileName = [FileName '.clu.' num2str(shank)];
end

Fp = fopen(cluFileName, 'r');

if Fp==-1
    error(['Could not open file ' FileName]);
end

% first #: total N of clusters
nClusters = fscanf(Fp, '%d', 1);
% one # per one spike in time order as detected (time stamps are in res file)
Clu = fscanf(Fp, '%d');

fclose(Fp);

return;