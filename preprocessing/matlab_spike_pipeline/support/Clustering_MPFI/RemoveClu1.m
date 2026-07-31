function RemoveClu1(filename, shank)
    % Remove clu 1 from the res, spk, and fet files
    % This is to remove the cluster with artifacts that contaminate and skew 
    % the PCA space
    
    cluFilename = [filename '.clu.' num2str(shank)];
    resFilename = [filename '.res.' num2str(shank)];
    spkFilename = [filename '.spk.' num2str(shank)];
    fetFilename = [filename '.fet.' num2str(shank)];
    xmlFilename = [filename '.xml'];

    if exist(cluFilename, 'file') == 2
        xml = LoadXml_e(xmlFilename(1:end-4));
        SpkSamples = xml.SpkGrps(1).nSamples;

        fprintf('El. group #: %d\n', shank);
        [~,cluList] = LoadClu_e1(cluFilename);   
        indSpikes = (cluList > 1);
        indSpikes0 = (cluList == 0);
        
        if(exist('./bkSortingFile') == 0)
            mkdir bkSortingFile
        end
        copyfile(cluFilename, './bkSortingFile');
        copyfile(resFilename, './bkSortingFile');
        copyfile(spkFilename, './bkSortingFile');
        copyfile(fetFilename, './bkSortingFile');            
                
        if (sum(indSpikes) > 0 && sum(indSpikes0) == 0)
            % if there is some spikes
            disp('Save new clu file.');     
            cluList1 = cluList(indSpikes,:);
            SaveClu(['./' cluFilename],cluList1);
            
            % if there is some spikes
            disp('Save new res file.');
            resList = load(resFilename,'r');      
            resList = resList(indSpikes,:);
            SaveRes(['./' resFilename],resList);
            
            disp('Save new fet file.');
            fetList = LoadFet_e1(fetFilename);
            fetList = fetList(indSpikes,:);
            SaveFet(['./' fetFilename],fetList);

            disp('Save new spk file.');
            spkList = LoadSpk(spkFilename,length(xml.ElecGp{shank}),...
                SpkSamples,length(cluList));
            spkList = spkList(:,:,indSpikes);
            SaveSpk(['./' spkFilename],spkList);
        end
    end
end