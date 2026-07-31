function [Spike,Clu] = getSpikes_mpfi(filename, SampleRate, lfpSampleRate)
    
    % added 08/14/2018, removed _SpikesData.mat file to reduce data 
    % duplication, and also prevented passing large structures as
    % parameters
    if exist([filename '_BehavElectrDataLFP.mat'], 'file') == 2
        load([filename '_BehavElectrDataLFP.mat'], 'Processing');
        if(sum(Processing == 2) == 1)
            return;
        end
        fprintf(...
            '\nLoad Track from %s file.\n', ...
            [filename '_BehavElectrDataLFP.mat']);
        load([filename '_BehavElectrDataLFP.mat'], 'Track');
    end
    
    % there is really no need to make this an input, Dinghao 23 Jan 2023
%     shankList=input('Specify a list of good shanks    [1:6]: ');
%     if(isempty(shankList))
%         shankList=1:6;
%     end    
    shankList=1:6;
    
    numSpStartRaw = 0; %round(tStart/1000 * xml.SampleRate);    
    numSpStart = 0; %round(numSpStartRaw / xml.SampleRate * xml.lfpSampleRate);
    
    Spike.res=[];
    Spike.res20kHz=[];
    Spike.clu=[];
    Spike.shank=[];
    Spike.totclu=[];
    
    Clu.shank = [];
    Clu.localClu = [];
    Clu.totClu = [];
    
    Clu.isIntern = [];
    Clu.SpkWidthC = [];
    Clu.RefracViolPercent = [];
    Clu.FirRateHz = [];
    Clu.isolDist = [];
    Clu.SpatLocalChan = [];
    Clu.SpatLocalRelAmpl = [];
    Clu.RightMax = [];
    Clu.LeftMax = [];
    Clu.CenterMax = [];         
    Clu.SpkWidthR = [];
    Clu.SpkWidthL = [];
    Clu.AmpSym = [];
    Clu.TimeSym = [];
    Clu.IsPositive = [];
        
    nClusters = [];
    
    fprintf('\n    Clu and Res files loaded; \n');
    
    for nsh = 1 : length(shankList);
        IDsh = shankList(nsh);
        cluFilename = [filename '.clu.' num2str(IDsh)];
        resFilename = [filename '.res.' num2str(IDsh)];
        spkFilename = [filename '.spk.' num2str(IDsh)];
        fetFilename = [filename '.fet.' num2str(IDsh)];
        xmlFilename = [filename '.xml'];
        
        if exist(cluFilename, 'file') == 2
            fprintf('El. group #: %d\n', IDsh);
            [~,cluList] = LoadClu_e1(cluFilename);   
            uniqueCluList = unique(cluList(cluList>1));
            Nclu = length(2:max(uniqueCluList));            
            
            if sum(cluList > 1) > 0              % if no spikes
                resList = load(resFilename,'r');       
                    % TimeStamps at the .dat sampling frequency 20.000 Hz
                    % -> change for xml.lfpSampleRate (round!)
                resList_eegSamplRate = round(resList / SampleRate...
                                        * lfpSampleRate);
                resList_eegSamplRate(resList_eegSamplRate == 0) = 1;

                if length(cluList)~=length(resList_eegSamplRate);
                    fprintf(...
                        '\nEl. group %d: The lengths of res and clu do not match',...
                        IDsh);
                end
                
                units = find(cluList > 1);
                cluUnits = cluList(units)-1;         
                    % the first good cluster ==2 -> change to 1
                resUnits_eegSamplRate = resList_eegSamplRate(units);
                resUnits = resList(units);
                shankUnits = ones(length(cluUnits),1) * IDsh;

                % save data into sturctures
                nClusters(IDsh)=Nclu;

                Spike.res = [Spike.res; resUnits_eegSamplRate];
                Spike.res20kHz = [Spike.res20kHz; resUnits];
                Spike.clu = [Spike.clu; cluUnits];
                Spike.shank = [Spike.shank; shankUnits];

                % do not count missing artifact and noise clusters 
                % from each shank
                if ~isempty(Spike.totclu)         % IDsh > 1
                    Spike.totclu = [Spike.totclu; ...
                            cluUnits + sum(nClusters(1:IDsh-1))];
                else Spike.totclu = cluUnits;
                end

                % which cells came from which shank and which channel
                IDTotClu = sum(nClusters(1:IDsh-1)) + (1:Nclu);
                Clu.shank(IDTotClu) = IDsh;
                Clu.localClu(IDTotClu) = 2:max(uniqueCluList); 
                Clu.totClu(IDTotClu) = IDTotClu';

                % get spike-shape characteristics for all clusters (Anton's fce:)
                if exist(spkFilename,'file') == 2                   
                    nq = NeuronQuality_e4(filename, IDsh, cluFilename, ...
                        resFilename, spkFilename, fetFilename, xmlFilename);  
                                % function NeuronQuality_e1(FileBase, ...
                                % Electrodes,Display=0,Batch=0,Overwrite=0)
                else
                    fprintf(['\nno *.spk file found - neuron quality function' ...
                            'is not called.\n']);
                end

                if exist(spkFilename,'file') == 2 && ~isempty(nq)
                    as = nq.AmpSym;    
                        % positive: hyperpolar peak > depolarization peak
                    swr = nq.SpkWidthR;
                    swl = nq.SpkWidthL;
                    tmas = nq.TimeSym;
                    fr = nq.FirRate;
                    indClu = nq.indClu;
                    ispos = nq.IsPositive;
                    rmax = nq.RightMax;
                    
                    Clu.isIntern(IDTotClu) = zeros(1,Nclu);
                    Clu.SpkWidthC(IDTotClu) = zeros(1,Nclu);
                    Clu.RefracViolPercent(IDTotClu) = zeros(1,Nclu);
                    Clu.FirRateHz(IDTotClu) = zeros(1,Nclu);
                    Clu.isolDist(IDTotClu) = zeros(1,Nclu);
                    Clu.SpatLocalChan(IDTotClu) = zeros(1,Nclu);
                    Clu.SpatLocalProbeCh(IDTotClu) = zeros(1,Nclu);
                    Clu.SpatLocalRelAmpl(IDTotClu) = zeros(1,Nclu);
                    Clu.RightMax(IDTotClu) = zeros(1,Nclu);
                    Clu.LeftMax(IDTotClu) = zeros(1,Nclu);
                    Clu.CenterMax(IDTotClu) = zeros(1,Nclu);
                    Clu.SpkWidthR(IDTotClu) = zeros(1,Nclu);
                    Clu.SpkWidthL(IDTotClu) = zeros(1,Nclu);
                    Clu.AmpSym(IDTotClu) = zeros(1,Nclu);
                    Clu.TimeSym(IDTotClu) = zeros(1,Nclu);
                    Clu.IsPositive(IDTotClu(indClu)) = zeros(1,Nclu);
                                       
                    Clu.confidIsIntern(IDTotClu(indClu)) = ...
                        sum([((as(indClu') > -0.3) & (as(indClu') < 0.7))...
                        swr(indClu')<0.5 fr(indClu')>3.5...
                        swl(indClu') > 0.3 tmas(indClu') <3 rmax > 25]');
                    Clu.isIntern(IDTotClu(indClu)) =...
                        Clu.confidIsIntern(IDTotClu(indClu)) > 2 & ...
                        ((as(indClu) > -0.3) & (as(indClu) < 0.7) == 1)' &...
                        (ispos(indClu) == 0)';   
                        % if at least 2 criteria are satisfied
                    Clu.SpkWidthC(IDTotClu(indClu)) = nq.SpkWidthC(indClu');
                    Clu.RefracViolPercent(IDTotClu(indClu)) = ...
                            nq.RefracViolPercent(indClu');
                    Clu.FirRateHz(IDTotClu(indClu)) = nq.FirRate(indClu');
                    Clu.isolDist(IDTotClu(indClu)) = nq.isolDist(indClu');
                    Clu.SpatLocalChan(IDTotClu(indClu)) = nq.SpatLocal(indClu');
                    Clu.SpatLocalProbeCh(IDTotClu(indClu)) = ...
                            nq.SpatLocalProbeCh(indClu)';
                    Clu.SpatLocalRelAmpl(IDTotClu(indClu)) = ...
                            nq.SpatLocalRelAmpl(indClu');
                    Clu.RightMax(IDTotClu(indClu)) = nq.RightMax(indClu');
                    Clu.LeftMax(IDTotClu(indClu)) = nq.LeftMax(indClu');
                    Clu.CenterMax(IDTotClu(indClu)) = nq.CenterMax(indClu');
                    Clu.SpkWidthR(IDTotClu(indClu)) = nq.SpkWidthR(indClu');
                    Clu.SpkWidthL(IDTotClu(indClu)) = nq.SpkWidthL(indClu');
                    Clu.AmpSym(IDTotClu(indClu)) = nq.AmpSym(indClu');
                    Clu.TimeSym(IDTotClu(indClu)) = nq.TimeSym(indClu');
                    Clu.IsPositive(IDTotClu(indClu)) = nq.IsPositive(indClu');
                    
                end
            end
        end
    end
    
    Spike.res20kHz = Spike.res20kHz + numSpStartRaw;
    Spike.res = Spike.res + numSpStart;
    
    if(sum(Processing == 2) == 0)
        Processing = [Processing 2]; 
        % processing stage two, getting Spike and Clu
    end
    fprintf('\nSpike-data saved into the structure file: %s....\n',...
            [filename '_BehavElectrDataLFP.mat']);
    save([filename '_BehavElectrDataLFP.mat'], ...
            'Spike', 'Clu', 'Processing', 'shankList', '-append');

end
