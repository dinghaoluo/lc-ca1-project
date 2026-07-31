% Spk = LoadSpkPartial(FileName, nChannels, SpkSampls, Spikes2Load)
%
% Loads part of a .spk file
% use this when loading the whole thing will cause theta to crash.
%
% returns a 3d array Spk(Channel, Sample, Spike Number)

function Spk = LoadSpkPartial(FileName, nChannels, SpkSampls, Spikes2Load)

nSpikes2Load = length(Spikes2Load);

Spk = zeros(nChannels, SpkSampls, nSpikes2Load);

fp = fopen(FileName);

for i=1:nSpikes2Load

	% go to correct part of file	
	status = fseek(fp, (Spikes2Load(i)-1)*nChannels*SpkSampls*2, 'bof');
	if (status == -1) ferror(fp); end;

	% load in wave form
	Spk(:,:,i) = fread(fp, [nChannels, SpkSampls], 'short');
end

fclose(fp);
