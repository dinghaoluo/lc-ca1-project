function nq = NeuronQuality_e4(FileBase, IDsh, cluFilename, resFilename, spkFilename, fetFilename, xmlFilename)

% give cluster quality measures for all clusters recorded from one
% electrode (shank)

Display = 0;
Overwrite = 0;

xml = LoadXml_e(xmlFilename(1:end-4));

%------------------------------------------------------------------------
% load data from a file if available
nqFileName = [FileBase '.NeuronQuality.mat'];
if exist(nqFileName, 'file') == 2 & ~Overwrite
    inFileShank = who('-file', nqFileName);
    isNqShank = strmatch(['nqShank' num2str(IDsh)], inFileShank, 'exact');
    if length(isNqShank) == 1
        fprintf('\nNeuron quality info for shank # %d is loaded from a file %s...\n', IDsh, nqFileName);
                
        load(nqFileName, ['nqShank' num2str(IDsh)]);
        eval(['nq = nqShank' num2str(IDsh) ';']); 
        return;
    end
end

%------------------------------------------------------------------------
% calcualte values if no stored data
%------------------------------------------------------------------------
SampleRate = 1e6./xml.SampleTime;
SpkSamples = xml.SpkGrps(1).nSamples;

% load all data
[nClusters, Clu] = LoadClu_e1(cluFilename);
Fet = LoadFet_e1(fetFilename);
Res = load(resFilename);
Res = Res ./ SampleRate .* 1000;      % Res in ms (not sampling rate)

nspk = FileLength(spkFilename)/2/length(xml.ElecGp{IDsh})/SpkSamples;
if nspk ~= length(Clu)
    fprintf('%s - Electrode %d - length of spk file does not correspond to clu length\n\n',FileBase,IDsh);
end
uClu = setdiff(unique(Clu),[0 1]);     % find all non-artifact/no-noise spikes
nClu = length(uClu);

%noise - comes from first 10000 spikes
SpkNoise = LoadSpk(spkFilename,length(xml.ElecGp{IDsh}),SpkSamples,10000);
stdSpkNoise =squeeze(std(SpkNoise,0,3));

%------------------------------------------------------------------------
%load only enough spike to sample all clusters
% create represantative spikes sample for good cells
avSpk =[]; stdSpk = [];SpatLocal=[];SpkWidthC=[];SpkWidthL=[];SpkWidthR=[];posSpk=[];FirRate = [];AvSpkAll=[];
leftmax=[]; rightmax=[];troughamp=[];troughSD=[]; indClu=[]; IsolationDist=[];RefracViolPercent=[];SpatLocalRelAmpl=[];
nq = struct;

for cnum=uClu'
    % get spike wavesdhapes and compute SNR
    SampleSize = 1000;
    myClu=find(Clu==cnum);
    
    if length(myClu)>2
        indClu = [indClu cnum];
        SampleSize = min(length(myClu),SampleSize);
        RndSample = sort(myClu(randsample(length(myClu),SampleSize)));
        mySpk = LoadSpkPartial_e1(spkFilename,length(xml.ElecGp{IDsh}),SpkSamples,RndSample);
        
        avSpk(:,:,cnum) = squeeze(mean(mySpk, 3));
        stdSpk(:,:,cnum)  = squeeze(std(mySpk,0,3));% may not need it
        
        %find the channel of largest amp (positive or negative)
        [amps ampch] = max(abs(avSpk(:,:,cnum)),[],2);
        [dummy maxampch] = max(squeeze(amps));
        nch = length(xml.ElecGp{IDsh});
        nonmax = setdiff([1:nch],maxampch);
        %compute spatial localization as ratio of max ch amplitude to the mean over all others.
        if nch>1
            SpatLocal(cnum) = maxampch;
            SpatLocalProbeCh(cnum) = xml.ElecGp{IDsh}(maxampch);
            SpatLocalRelAmpl(cnum) = amps(maxampch)/(mean(amps(nonmax))+eps);
            
        else
            SpatLocal(cnum) = 0;
            SpatLocalProbeCh(cnum) = -1;
            SpatLocalRelAmpl(cnum) = 0;
        end
        myAvSpk = squeeze(avSpk(maxampch,:,cnum)); % largest channel spike wave for that cluster
        
        %now we need to take care of the positive spikes (we reverse them)
        minamp  = abs(min(myAvSpk));
        maxamp  = max(myAvSpk);
        %        if (minamp-maxamp)/minamp < 0.05 %(spike is more positive then negative)
        if maxamp>1.2*minamp %(spike is more positive then        negative)
            myAvSpk = -myAvSpk; %reverse the spike
            posSpk(cnum) = 1;
        else
            posSpk(cnum) = 0;
        end
        
        %now let's upsample the spike waveform
        ResCoef = 10; %                           
        Sample2Msec = 1000./SampleRate/ResCoef; %to get fromnew samplerate to the msec
        myAvSpk = resample(myAvSpk,ResCoef,1);
        
        [troughamp(cnum) troughTime] = min(myAvSpk);
        if troughTime <= 5
            troughTime = 6;
        end
        pts= myAvSpk(troughTime+[-5 0 5]);
        pts=pts(:);
        troughSD(cnum) = pts'*[1 -2 1]';
        amphalf = 0.5*min(myAvSpk);
        both=0;cnt=0;halfAmpTimes=[0 0];
        while both<2
            if halfAmpTimes(1)==0 && myAvSpk(troughTime-cnt)>amphalf
                halfAmpTimes(1)=troughTime-cnt;
                both=both+1;
            end
            if halfAmpTimes(2)==0 && myAvSpk(troughTime+cnt)>amphalf
                halfAmpTimes(2)=troughTime+cnt;
                both=both+1;
            end
            cnt=cnt+1;
            if cnt == troughTime
                break;
            end
        end
        %width
        SpkWidthC(cnum) = diff(halfAmpTimes)*Sample2Msec;        
        
        %dmyAvSpk = diff(myAvSpk);
        SpkPieceR = myAvSpk(troughTime:end);
        [rightmax(cnum) SpkWidthR(cnum)] = max(SpkPieceR); % this is the distance from the trough to the rise peak
        SpkWidthR(cnum)= SpkWidthR(cnum)*Sample2Msec;
        
        SpkPieceL = myAvSpk(1:troughTime);
        SpkPieceL = flipud(SpkPieceL(:)); % to look at the right time lag
        [leftmax(cnum) SpkWidthL(cnum)] = max(SpkPieceL); % this is the distance from the peak to the trough
        SpkWidthL(cnum)= SpkWidthL(cnum)*Sample2Msec;
        
        troughTime = troughTime*Sample2Msec;
        if posSpk(cnum);	myAvSpk = -myAvSpk; end
        AvSpkAll(cnum,:) = myAvSpk;
                
        if Display
            figure(765)
            if cnum==2
                clf;
            end
            subplotfit(cnum-1,length(uClu));
            shift = troughTime;%SpkSamples/2*Sample2Msec;
            plot([1:SpkSamples*ResCoef]*Sample2Msec-shift,myAvSpk, 'k');
            axis tight
            hold on
            Lines(0,[],'g');%trough line
            Lines(halfAmpTimes*Sample2Msec-shift,amphalf,'r');
            Lines(SpkWidthR(cnum),[],'b');
            Lines(-SpkWidthL(cnum),[],'b');
            %Lines([-SpkWidthL SpkWidthR], troughamp,'r');
            mystr = sprintf('El=%d Clu=%d',IDsh,cnum);
            title(mystr);
        end
        
        % firing rate
        myRes = Res(find(Clu==cnum));
        ISI = diff(myRes);
        MeanISI = mean(bootstrp(100,'mean',ISI));
        FirRate(cnum) = (1/MeanISI) * 1000;     % FR per sec (Res is per msec)
        %              check whether the refractory period is clean
        Refrac = 2; %msec
        RefracViolPercent(cnum) = sum(ISI<Refrac)/length(ISI) * 100;
        
        % get Isolation Distance between the cluster and all other good clusters
        IsolationDist(cnum) = Cluster_MahalDist(Fet, Clu==cnum, Res);    % fet for all clu, no-noise&no-artif clu, target clu, res
        
    end
    
end

snr = squeeze(mean(mean(abs(avSpk),1),2))./mean(stdSpkNoise(:));

nq.SNR = snr(2:end);
nq.indClu = indClu-1;
nq.SpkWidthC = SpkWidthC(2:end)';
nq.SpkWidthR = SpkWidthR(2:end)';
nq.SpkWidthL = SpkWidthL(2:end)';
%fix for empty clusters
SpkWidthL(SpkWidthL==0)=1000000;
nq.TimeSym = SpkWidthR(2:end)'./SpkWidthL(2:end)';
nq.isolDist = IsolationDist(2:end)';   % Isolation distance based on Harris et al, 2001 - mahalanobis distance of the n-th noise spike from the center of myClu, n=length(myClu);
nq.RefracViolPercent = RefracViolPercent(2:end)';
nq.Clus = uClu; % list of clusters
nq.ElNum = repmat(IDsh,length(uClu),1);
nq.IsPositive = posSpk(2:end)';
nq.SpatLocal=SpatLocal(2:end)';
nq.SpatLocalProbeCh = SpatLocalProbeCh(2:end);
nq.SpatLocalRelAmpl=SpatLocalRelAmpl(2:end)';
nq.FirRate = FirRate(2:end)';
nq.AvSpk = AvSpkAll(2:end,:);
nq.RightMax= rightmax(2:end)';
nq.LeftMax= leftmax(2:end)';
nq.CenterMax= troughamp(2:end)';
%fix for empty clusters
rightmax(rightmax==0)=1e6;
leftmax(leftmax==0)=1e6;
nq.AmpSym = (abs(rightmax(2:end))'-abs(leftmax(2:end))')./(abs(rightmax(2:end))'+abs(leftmax(2:end))');
nq.troughSD = troughSD(2:end)';

fprintf('\nNeuron Quality file for shank # %d was saved in the file %s.\n\n', IDsh, nqFileName);
variabName = (genvarname(['nqShank' num2str(IDsh)]));
eval([variabName ' = nq;']);

if exist(nqFileName, 'file') == 2
    save(nqFileName, '-append', variabName);
else
    save(nqFileName, variabName);
end

return;
