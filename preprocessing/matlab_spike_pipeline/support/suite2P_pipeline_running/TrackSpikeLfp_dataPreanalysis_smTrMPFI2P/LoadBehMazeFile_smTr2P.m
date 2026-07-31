function [behEvents] = LoadBehMazeFile_smTr2P(FileName,sampleFreq)

fid=fopen(FileName);
tr = 0;
nt = 0;
te = 0;
nw = 0;
nl = 0;
np = 0;
nq = 0;
sy = 0;
sy_2p = 0;
mv = 0;
lp = 0;
mt = 0;
pc = 0;
po = 0;
pp = 0;
ds = 0;
of = 0; % microsecond timer overflow
ofConst = (2^32-1)/1000;
line = 0;
should_start = 0;

behEvents.dropped_fm_count = 0;
behEvents.dropped_idx = [];

skipped_first_tr = 0;

while ~feof(fid)
  line = line + 1;
  descr = fscanf(fid, '%c', 1);
  while(descr ~= '$')
      descr = fscanf(fid, '%c', 1);
  end
  descr = [descr fscanf(fid, '%c', 2)];
  if(isempty(descr))
      fclose(fid);
      return;
  end
  dat = str2num(fgetl(fid));

  if(line == 2)
      timePre = dat(1);
  end
%   disp(['line = ' num2str(line) ' ' descr ' ' dat]);
%   if(line == 417)
%         a = 1;
%   end
  if(line > 2 && timePre - dat(1) > 4e+6)
    of = of+1;
  end

  if(descr(2:3) == "FM" && should_start < 2)
    should_start = should_start + 1;
  end
  
  if(should_start < 2)
      continue;
  end
  
  switch descr(2:3)
    % task description
    case 'TR'  
        if(length(dat) < 3)%
            continue;
        end
%         if(skipped_first_tr == 0)
%            skipped_first_tr = 1;
%            continue;
%         end
        tr = tr + 1;
		behEvents.taskDescr(tr,1:3) = dat(1:3); 
        behEvents.taskDescr(tr,1) = behEvents.taskDescr(tr,1)+of*ofConst;
        if(length(dat) > 3)
            behEvents.movieTDescr{tr} = dat(4:end);
        end
    % trial description   
    case 'NT'  
        if(length(dat) ~= 2 && length(dat) ~= 3)%
            continue;
        end
        nt = nt + 1;
		behEvents.trialDescr(nt,:) = dat; 
        behEvents.trialDescr(nt,3) = nt; 

        behEvents.trialDescr(nt,1) = behEvents.trialDescr(nt,1)+of*ofConst;
    % trial start/end   
    case 'TE'  
        if(length(dat) ~= 2)%
            continue;
        end
        te = te + 1;
		behEvents.trialT(te,:) = dat; 
        behEvents.trialT(te,1) = behEvents.trialT(te,1)+of*ofConst;
    % treadmill  
    case 'WE'  
        if(length(dat) ~= 3)%
            continue;
        end
        nw = nw + 1;
		behEvents.wheel(nw,:) = dat; 
        behEvents.wheel(nw,1) = behEvents.wheel(nw,1)+of*ofConst;
    % lick port    
    case 'LE'  
        if(length(dat) ~= 3)
            continue;
        end
        nl = nl + 1;
		behEvents.lick(nl,:) = dat; 
        behEvents.lick(nl,1) = behEvents.lick(nl,1)+of*ofConst;
    % pump    
    case 'PE'  
        if(length(dat) ~= 3)
            continue;
        end
        np = np + 1;
		behEvents.pump(np,:) = dat; 
        behEvents.pump(np,1) = behEvents.pump(np,1)+of*ofConst;
    % tone   
    case 'TN'  
        if(length(dat) ~= 2)%
            continue;
        end
        nq = nq + 1;
		behEvents.tone(nq,:) = dat;
        behEvents.tone(nq,1) = behEvents.tone(nq,1)+of*ofConst;
    % airpuff   
    case 'AE'  
        if(length(dat) ~= 2)%
            continue;
        end
        nq = nq + 1;
		behEvents.airpuff(nq,:) = dat;
        behEvents.airpuff(nq,1) = behEvents.airpuff(nq,1)+of*ofConst;
    % SYNC pulse
    case 'SY'  
        if(length(dat) ~= 1)%
            continue;
        end
        sy = sy + 1;
		behEvents.ArdSyncMsec(sy,:) = dat; 
        behEvents.ArdSyncMsec(sy,1) = behEvents.ArdSyncMsec(sy,1)+of*ofConst;
    % 2p SYNC pulse
    case 'FM'
        if(length(dat) ~= 2)%
            continue;
        end
        if(dat(:, 2) == 0)
            
            sy_2p = sy_2p + 1;
            behEvents.TwoPSyncMsec(sy_2p,:) = dat(:, 1); 
            behEvents.TwoPSyncMsec(sy_2p,1) = behEvents.TwoPSyncMsec(sy_2p,1)+of*ofConst;
            
%             if the interval between two FMs is greater than 50, then add
%             a FM in between. this is to alleviate the issue of the 2p
%             system not always sending an FM signal per frame it takes
            if(sy_2p > 1 && behEvents.TwoPSyncMsec(sy_2p,1) - behEvents.TwoPSyncMsec(sy_2p-1,1) > 1000/sampleFreq*1.7)
                behEvents.dropped_idx(end+1) = sy_2p - behEvents.dropped_fm_count;
                numDroppedFm = round((behEvents.TwoPSyncMsec(sy_2p,1) - behEvents.TwoPSyncMsec(sy_2p-1,1))/(1000/sampleFreq));
                avgTsp = (behEvents.TwoPSyncMsec(sy_2p,1) - behEvents.TwoPSyncMsec(sy_2p-1,1))/numDroppedFm;
                behEvents.TwoPSyncMsec(sy_2p+numDroppedFm-1,1) = behEvents.TwoPSyncMsec(sy_2p,1);
                for i = 1:numDroppedFm-1
                    behEvents.TwoPSyncMsec(sy_2p,1) = behEvents.TwoPSyncMsec(sy_2p-1,1) + avgTsp;
                    sy_2p = sy_2p + 1;
                end    
                
                behEvents.dropped_fm_count = behEvents.dropped_fm_count + numDroppedFm-1;
            end
            
            if(sy_2p == 1)
               behEvents.firstFM = behEvents.TwoPSyncMsec(sy_2p,1);
            end
            behEvents.lastFM = behEvents.TwoPSyncMsec(sy_2p,1);
        end
    % movie
    case 'MV'
        if(length(dat) ~= 2 && length(dat) ~= 1)%
            continue;
        end
        mv = mv + 1;
		behEvents.movieOn(mv,:) = dat; 
        behEvents.movieOn(mv,1) = behEvents.movieOn(mv,1)+of*ofConst;
    % lick period
    case 'LP'
        if(length(dat) ~= 2)%
            continue;
        end
        lp = lp + 1;
		behEvents.lickPeriod(lp,:) = dat; 
        behEvents.lickPeriod(lp,1) = behEvents.lickPeriod(lp,1)+of*ofConst;
    % brake
    case 'MT'
        if(length(dat) ~= 2)
            continue;
        end
        mt = mt + 1;
		behEvents.brakeOn(mt,:) = dat; 
        behEvents.brakeOn(mt,1) = behEvents.brakeOn(mt,1)+of*ofConst;
    % optogenetics stimulation pulse
    case 'PC'
        if(length(dat) ~= 1)
            continue;
        end
        pc = pc + 1;
		behEvents.stimOn(pc,:) = dat; 
        behEvents.stimOn(pc,1) = behEvents.stimOn(pc,1)+of*ofConst;
    % optogenetics stimulation pulse, stim off
    case 'PO'
        if(length(dat) ~= 1)
            continue;
        end
        po = po + 1;
		behEvents.stimOff(po,:) = dat; 
        behEvents.stimOff(po,1) = behEvents.stimOff(po,1)+of*ofConst;
    % pulse parameter
    case 'PP'
        if(length(dat) < 8)
            continue;
        end
        pp = pp + 1;
		behEvents.pulsePar(pp,:) = dat; 
        behEvents.pulsePar(pp,1) = behEvents.pulsePar(pp,1)+of*ofConst;
    % diode parameter
    case 'DS'
        if(length(dat) ~= 7)
            continue;
        end
        ds = ds + 1;
		behEvents.diodePar(ds,:) = dat; 
        behEvents.diodePar(ds,1) = behEvents.diodePar(ds,1)+of*ofConst;
  end
  
  if(line > 2)
    timePre = dat(1);
  end
end

fclose(fid);

return