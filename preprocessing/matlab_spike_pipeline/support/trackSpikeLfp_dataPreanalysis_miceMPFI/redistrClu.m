function [nTotSampl] = redistrClu(path, concatFolder, recordFolder)

fprintf(['\nRedistributing concatenated clu, fet, spk and res files',...
        ' into individual recordings...']);

% load the xml file from the CONCATENATED folder
cd([path '/' concatFolder]);
if exist([concatFolder '.xml'], 'file') == 2
    [cXml] = LoadXml_e([concatFolder '.xml']);
else fprintf('\n Cannot find %s - file -> analysis aborted.\n', ...
             [concatFolder '.xml']);
    return;
end
cd ..;

% get the number of total samples and the number of channels 
% for each recording
nTotSampl = zeros(1,length(recordFolder));
for nf = 1 : length(recordFolder)
    FileName = recordFolder{nf};
    if(exist(FileName,'dir') == 0)
        mkdir(FileName);
    end
    if (exist([path '/' FileName '/' FileName '.xml'],'file') ==  0)
        copyfile([path '/' concatFolder '/' concatFolder '.xml'], ...
            [path '/' FileName '/' FileName '.xml']);
    end

    fileinfo = dir([path '/' concatFolder '/' FileName '-locRef.dat']);
    nTotSampl(nf) = ceil(fileinfo(1).bytes / 2 / cXml.nChannels);
%     dat = readmulti([path '/' FileName '/' FileName '.dat'], ...
%                      cXml.nChannels,cXml.SpkGrps(1).Channels(1)+1); 
%     nTotSampl(nf) = length(dat); 
end

% Redistributing concatenated clu, fet, spk and res files into individual 
% recordings
totNspikes = [];

for nsh = 1 : cXml.nElecGps
    cd([path '/' concatFolder]);
    if ~isempty(dir([concatFolder '.clu.' num2str(nsh)]))
        [ncClu, cClu] = LoadClu_e1(concatFolder, nsh); 
                    % load n-th clu file.
        cRes = load([concatFolder '.res.' num2str(nsh)]); 
                    % load n-th res file.
        cSpk = LoadSpk([concatFolder '.spk.' num2str(nsh)], ...
                    length(cXml.ElecGp{nsh}), cXml.SpkGrps(1).nSamples); 
                    % load n-th spk file.
        cFet = LoadFet_e1([concatFolder '.fet.' num2str(nsh)]);
                    % load n-th fet file.
        
        maxClu = max(cClu);
        
        cd ..;
        
        recTStart = 1;     
        
        if ncClu > 0
            % copy sections of clu, res, spk and fet files into subfolders 
            % with original data
            for nf = 1 : length(recordFolder)                
                FileName = recordFolder{nf};
                disp(['shank ID: ' num2str(nsh) ': recording file: ' ...
                      FileName]);
                cd([path '/' FileName]);
                if(exist([FileName '.clu.' num2str(nsh)],'file') == 2)
                    fprintf('File %s already exists.\n', ...
                            [FileName '.clu.' num2str(nsh)]);
                    continue;
                end
                                
                % recTStart; 
                recTEnd = nTotSampl(nf) + recTStart - 1;        
                    % N samples in each file
                
                myRes = cRes>recTStart & cRes<=recTEnd;
                rRes = cRes(myRes,:)-recTStart+1;
                rSpk = cSpk(:,:,myRes);
                rClu = cClu(myRes,:);
                rFet = cFet(myRes,:);
                totNspikes(nsh, nf) = length(rRes);
                
                % add one spike for the last cluster if it is not present
                % in the current recording to make sure the total N of clu
                % is the same in all recordings
                if max(rClu) < maxClu
                    missCluList = setdiff(unique(cClu), unique(rClu));
                    for nMissClu = 1 : length(missCluList)
                        myClu = missCluList(nMissClu);
                        rRes = [1; rRes];
                        rSpk = cat(3, rSpk(:,:,1), rSpk);
                        rClu = [myClu; rClu];                        
                        % add zeros into fet file and time to the end
                        lastFet = [-ones(1,size(rFet,2)-1) rFet(end,end)];
                        rFet = [lastFet;rFet];
                        disp(['One artificial spike inserted at the begining',...
                            ' of fet, clu, res and spk files so that',...
                            ' the maxClu ID is the same in all the recordings',...
                            ' (NOTE: res value is 1).']);
                    end
                end
                nrClu = length(unique(rClu));
                
                %save new clu, spk and res files.
                rCluOut = [nrClu; rClu];   
                    % put nClusters and Clu back together.
                rSpkOut = reshape(rSpk, [numel(rSpk) 1]); 
                    % reshape Spk 3-D array into the spkList vector.
                
                fid1 = fopen([FileName '.clu.' num2str(nsh)], 'w');
                fprintf(fid1, '%d\n', rCluOut); 
                    % save cluFile into new clu.n file.
                fid2 = fopen([FileName '.spk.' num2str(nsh)], 'w');
                fwrite(fid2, rSpkOut, 'short'); 
                    % save spkList into new spk.n file.
                fid3 = fopen([FileName '.res.' num2str(nsh)], 'w');
                fprintf(fid3, '%d\n', rRes);  
                    % save tStamps into new res.n file.
                SaveFetFile([FileName '.fet.' num2str(nsh)], rFet);
                    % save fet file
                fclose(fid1); fclose(fid2); fclose(fid3);
                
                cd ..;
                recTStart = recTEnd + 1;    
                % start the next rec 1 sample after the end of 
                % the previous recording
                
            end
        end
    end
end
