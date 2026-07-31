% SaveSpk(FileName, Spk)
%
% Saves a .spk file
%
% takes a 3d array Spk(Channel, Sample, Spike Number)

function SaveSpk(FileName, Spk)

nChannels = size(Spk, 1);
SpkSampls = size(Spk, 2);
nSpikes = size(Spk, 3);

Spk = reshape(Spk, [nChannels* SpkSampls* nSpikes, 1]);

bsave(FileName, Spk, 'short');