function GetNeuRelativeDepth(path,fileName,onlyRun)
% Get the depth of each neuron relative to the pyramidal layer center
% Input arguments:
% path:         the path of the recording file
% fileName:     name of the recording file
% onlyRun:       1: only consider the time period when the animal is running 
%
% e.g.: BasicInfo('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1',1)

    %%%%%%%% check arguments
    if nargin<3
        disp('At least three arguments are needed for this function.');
        return;
    elseif nargin > 3
        disp('Too many arguments');        
        return;
    end
    
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(isempty(indexFileName))
        fileNameDepth = [fileName '_Depth.mat'];
        fileNameFR = [fileName '_FR_Run' num2str(onlyRun) '.mat'];
        fileName = [fileName '.mat'];  
    else
        fileNameDepth = [fileName(1:indexFileName(end)-1) '_Depth.mat'];
        fileNameFR = [fileName(1:indexFileName(end)-1) '_FR_Run' ...
                        num2str(onlyRun) '.mat'];
    end 
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath,'cluList');
    
    fullPath = [path fileNameFR];
    if(exist(fullPath) == 0)
        disp('Mean firing rate file does not exist.');
        return;
    end
    load(fullPath,'mFRStruct');
    
    if(isunix == 0)
        templatePath = 'Z:\Yingxue\mice_ephys_old\xmlTemplate\64Probe6Shanks.xml';
    else
        templatePath = ...
        '~/Desktop/Network-Drives/Wang_Lab/Yingxue/mice_ephys_old/xmlTemplate/64Probe6Shanks.xml';
    end
    xml_t = LoadXml_e(templatePath(1:end-4));
    ind = strfind(fileName,'_');
    xml = LoadXml_e([path fileName(1:ind(1)-1)]);
    
    % RH 12/08/2020 - Need to use different template with single shank
    % recordings for MS - just use original xml as template
    xml_t = xml;
    
    %%%%%%%%% define pyramidal layer center according to the channel that has the max
    %%%%%%%%% number of clusters
    shanks = unique(cluList.shank);
        
    depthNeu = struct('trueSpatLocalChan', [],...
                      'layerCenter', zeros(1,max(shanks)),...
                      'relDepthNeu', [],...
                      'layerCenterHDef',[],...
                      'relDepthNeuHDef', []);
                  
    %% get the true local spatial position for each clu 
    indNeuPerSh = cell(1,max(shanks));
    rowInd = cell(1,max(shanks));
    numClu = zeros(1,length(shanks));
    for n = shanks
        indNeuPerSh{n} = cluList.shank == n;   
        numClu(n) = sum(indNeuPerSh{n});
        
        % check whether bad channels have been deleted from the shank
        if(length(xml_t.ElecGp{n}) ~= length(xml.ElecGp{n}))
            templateChRep = repmat(xml_t.ElecGp{n}',1,numClu(n));
            cluChRep = repmat(cluList.spatLocalProbeCh(indNeuPerSh{n}),...
                length(xml_t.ElecGp{n}),1);
            mismatchChRep = templateChRep - cluChRep;
            [rowInd{n},colInd] = find(mismatchChRep == 0);    
        else 
            rowInd{n} = cluList.spatLocalChan(indNeuPerSh{n});
        end
        depthNeu.trueSpatLocalChan = [depthNeu.trueSpatLocalChan, rowInd{n}]; 
    end
    
    %% calculate the relative depth using the channel with max number of clus
    %  as the center
    for n = shanks
        numNeuPerCh = zeros(1,sum(indNeuPerSh{n}));
        for m = 1:length(xml_t.ElecGp{n})
            numNeuPerCh(m) = sum(rowInd{n}==m & ...
                mFRStruct.mFR(indNeuPerSh{n}) > 0.15);
        end
        [~,depthNeu.layerCenter(n)] = max(numNeuPerCh);
        depthNeu.relDepthNeu = [depthNeu.relDepthNeu, ...
            rowInd{n}-depthNeu.layerCenter(n)];
    end
    
    %% human defined pyramidal layer center
    indTmp = strfind(fileName,'_');
    layerCenterName = [fileName(1:indTmp(1)-1) '-ChForLayerCenter.txt'];
    if(exist(layerCenterName,'file')==2)
        fid=fopen(layerCenterName);
        layerCenterHDef = [];
        while(~feof(fid))
            layerCenterHDef = [layerCenterHDef, str2num(fgetl(fid))];
        end
        fclose(fid);
        depthNeu.layerCenterHDef = layerCenterHDef;
        for n = shanks
            if(layerCenterHDef(n) >= 0)
                indLocalCenter = find(xml_t.ElecGp{n} ...
                    == layerCenterHDef(n));
            elseif(layerCenterHDef(n) < 0 && ...
                    layerCenterHDef(n) > -20 ) % in SR
                indLocalCenter = layerCenterHDef(n)+1;
            else % in Orien
                indLocalCenter = length(xml_t.ElecGp{n}) + ...
                    abs(layerCenterHDef(n))-20;
            end
            disp(n);
%             depthNeu.relDepthNeuHDef = [depthNeu.relDepthNeuHDef,...
%                 depthNeu.trueSpatLocalChan(indNeuPerSh{n})-indLocalCenter];
            % positive -- SR; negative -- oriens
        end
    end
    
    %%%%%%%%% save data to file     
    fullPath = [path fileNameDepth];
    save(fullPath,'depthNeu');
        
end