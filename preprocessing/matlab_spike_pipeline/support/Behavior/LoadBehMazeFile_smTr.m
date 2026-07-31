function [behEvents] = LoadBehMazeFile_smTr(FileName)

fid=fopen(FileName);
tr = 0;
nt = 0;
te = 0;
nw = 0;
nl = 0;
np = 0;
nq = 0;
sy = 0;
mv = 0;
lp = 0;
of = 0; % microsecond timer overflow
ofConst = (2^32-1)/1000;
line = 0;

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
  if(line > 2 && timePre - dat(1) > 4e+6)
    of = of+1;
  end

  switch descr(2:3)
    % task description
    case 'TR';  
        if(length(dat) < 4)%
            continue;
        end
        tr = tr + 1;
		behEvents.taskDescr(tr,1:3) = dat(1:3); 
        behEvents.taskDescr(tr,1) = behEvents.taskDescr(tr,1)+of*ofConst;
        behEvents.movieTDescr{tr} = dat(4:end);
     % trial description   
    case 'NT';  
        if(length(dat) ~= 2)%
            continue;
        end
        nt = nt + 1;
		behEvents.trialDescr(nt,:) = dat; 
        behEvents.trialDescr(nt,1) = behEvents.trialDescr(nt,1)+of*ofConst;
    % trial start/end   
    case 'TE';  
        if(length(dat) ~= 2)%
            continue;
        end
        te = te + 1;
		behEvents.trialT(te,:) = dat; 
        behEvents.trialT(te,1) = behEvents.trialT(te,1)+of*ofConst;
    % treadmill  
    case 'WE';  
        if(length(dat) ~= 3)%
            continue;
        end
        nw = nw + 1;
		behEvents.wheel(nw,:) = dat; 
        behEvents.wheel(nw,1) = behEvents.wheel(nw,1)+of*ofConst;
    % lick port    
    case 'LE';  
        if(length(dat) ~= 3)
            continue;
        end
        nl = nl + 1;
		behEvents.lick(nl,:) = dat; 
        behEvents.lick(nl,1) = behEvents.lick(nl,1)+of*ofConst;
    % pump    
    case 'PE';  
        if(length(dat) ~= 3)
            continue;
        end
        np = np + 1;
		behEvents.pump(np,:) = dat; 
        behEvents.pump(np,1) = behEvents.pump(np,1)+of*ofConst;
    % tone   
    case 'TN';  
        if(length(dat) ~= 2)%
            continue;
        end
        nq = nq + 1;
		behEvents.tone(nq,:) = dat;
        behEvents.tone(nq,1) = behEvents.tone(nq,1)+of*ofConst;
    % airpuff   
    case 'AE';  
        if(length(dat) ~= 2)%
            continue;
        end
        nq = nq + 1;
		behEvents.airpuff(nq,:) = dat;
        behEvents.airpuff(nq,1) = behEvents.airpuff(nq,1)+of*ofConst;
    % SYNC pulse
    case 'SY';  
        if(length(dat) ~= 1)%
            continue;
        end
        sy = sy + 1;
		behEvents.ArdSyncMsec(sy,:) = dat; 
        behEvents.ArdSyncMsec(sy,1) = behEvents.ArdSyncMsec(sy,1)+of*ofConst;
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
  end
  
  if(line > 2)
    timePre = dat(1);
  end
end

fclose(fid);

return