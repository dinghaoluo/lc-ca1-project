function [eeg] = LoadDatFile(filename, chID, startRead, chunkLength, nChannelsTot)

% set up array for Dat
eeg = zeros(length(chID), chunkLength);

% load dat data	
fp = fopen([filename '.dat'], 'r');
status = fseek(fp, startRead*nChannelsTot*2, 'bof');

buffersize = 2^10; %2^18;

eeg=zeros(length(chID),chunkLength);
N_EL=0;
numelm=0;

while N_EL < chunkLength
    [data,count] = fread(fp,[nChannelsTot,buffersize],'int16');
    numelm = count/nChannelsTot;
    if numelm>0 % Kenji modified 061009.Otherwise if numelm == 0 an error occur.
        eeg(:,N_EL+1:N_EL+numelm) = data(chID,:);
        N_EL = N_EL+numelm;
    else break;
    end 
end
fclose(fp);

eeg = eeg(1:chunkLength);

return;
