function GenCluDuringImmobile(filename, shank)
    % generate new clu, fet, res, and spk files containning the information
    % of all the spikes generated during non-running period
    
    cluFilename = [filename '.clu.' num2str(shank)];
    resFilename = [filename '.res.' num2str(shank)];
    spkFilename = [filename '.spk.' num2str(shank)];
    fetFilename = [filename '.fet.' num2str(shank)];
    xmlFilename = [filename '.xml'];
    
    GlobalConst;
    
    if exist(cluFilename, 'file') == 2
        xml = LoadXml_e(xmlFilename(1:end-4));
        SpkSamples = xml.SpkGrps(1).nSamples;

        fprintf('El. group #: %d\n', shank);
        [~,cluList] = LoadClu_e1(cluFilename);   
        indSpikes = (cluList>1);
        
        if sum(cluList > 1) > 0              % if no spikes
            resList = load(resFilename,'r');      
            resList_eegSamplRate = round(resList / sampleFqOri...
                    * sampleFq);
            resList_eegSamplRate(resList_eegSamplRate == 0) = 1;
            
            if exist([filename '_BehavElectrDataLFP.mat'], 'file') == 2
            	load([filename '_BehavElectrDataLFP.mat'],'Track');
                indSpeed = Track.speed_MMsec(resList_eegSamplRate) ...
                    < minSpeed; 
                
                indClu = indSpikes & indSpeed;
                
                dir1 = 'immobile';
                mkdir(dir1);
                
                disp('Copy .xml file.');
                copyfile(xmlFilename, ['./' dir1 '/' xmlFilename]);
                
                disp('Save new clu file.');
                cluList1 = cluList(indClu);
                SaveClu(['./' dir1 '/' cluFilename],cluList1);
                
                disp('Save new res file. \n');
                resList = resList(indClu);
                SaveRes(['./' dir1 '/' resFilename],resList);
                
                disp('Save new fet file. \n');
                fetList = LoadFet_e1(fetFilename);
                fetList = fetList(indClu,:);
                SaveFet(['./' dir1 '/' fetFilename],fetList);
                
                disp('Save new spk file. \n');
                spkList = LoadSpk(spkFilename,length(xml.ElecGp{shank}),...
                    SpkSamples,length(cluList));
                spkList = spkList(:,:,indClu);
                SaveSpk(['./' dir1 '/' spkFilename],spkList);
            else
                disp('Please run GenerateBehavElectroDataStructures func. first');
            end
        end
    end
end
