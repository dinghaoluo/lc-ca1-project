function [eeg] = LoadDatFile_FL_uint16(fid, chID, Nsamples, samplRate, nChannelsTot)
%%% load the selected channels from the .dat file
%%% filename: file name without extension
%%% chID: channel number, the index starts from 1 instead of 0

%%% had to modify this script to use uint16 for intan file types


%     %% load meta file to get the sample rate and file size
%     fid = fopen([filename '.meta']);
%     metaInfo=textscan(fid, '%s');
%     metaInfo=metaInfo{1};
%     fclose(fid);
%     
%     samplRate = str2double(metaInfo(55));
%     fileSize = str2double(metaInfo(26));    % bytes
%     nChannelsTot = str2double(metaInfo(38));    % total N of ch
    
%     %% calculate the total number of samples recorded in the file
%     if(nargin == 2)
%         listing = dir([filename '.dat']);
%         Nsamples = listing.bytes/2/nChannelsTot;  
%     end
    
    %% load 1 min of data at a time
    nChunks = 1;
    lastChunkSize = 0;
    chunkLength = round(samplRate * 60);       
    if Nsamples > chunkLength    
        nChunks = floor(Nsamples / chunkLength);
        lastChunkSize = floor(Nsamples - (nChunks*chunkLength));
    end
    
     %% load the channel from the .dat file
     startRead = 0;
     eeg = zeros(length(chID),Nsamples);
     for n = 1:nChunks
        disp(['segment #' num2str(n) ' (out of ' num2str(nChunks+1) ')']);    
        eeg(:,(n-1)*chunkLength+1:n*chunkLength) = LoadDatFile(fid, chID, startRead, chunkLength, nChannelsTot); 
        startRead = startRead + chunkLength;
     end
     if lastChunkSize > 0
        disp(['segment #' num2str(nChunks+1) ' (out of ' num2str(nChunks+1) ')']);    
        eeg(:,nChunks*chunkLength+1:end) = LoadDatFile(fid, chID, startRead, lastChunkSize, nChannelsTot); 
     end
end

function [eeg] = LoadDatFile(fid, chID, startRead, chunkLength, nChannelsTot)

    % set up array for Dat
    eeg = zeros(length(chID), chunkLength);

    % load dat data	
    status = fseek(fid, startRead*nChannelsTot*2, 'bof');

    buffersize = 2^10; %2^18;

    eeg=zeros(length(chID),chunkLength);
    N_EL=0;
    numelm=0;

    while N_EL < chunkLength
        [data,count] = fread(fid,[nChannelsTot,buffersize],'uint16');
        numelm = count/nChannelsTot;
        if numelm>0 % Kenji modified 061009.Otherwise if numelm == 0 an error occur.
            eeg(:,N_EL+1:N_EL+numelm) = data(chID,:);
            N_EL = N_EL+numelm;
        else break;
        end 
    end

    eeg = eeg(:,1:chunkLength);

    return;
end