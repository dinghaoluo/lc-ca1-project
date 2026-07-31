function [behEvents] = LoadBehMazeFile(FileName)

fid=fopen(FileName);
nt = 0;
te = 0;
nw = 0;
nd = 0;
nb = 0;
nl = 0;
np = 0;
nq = 0;
sy = 0;
br = 0;
ht = 0;
of = 0; % microsecond timer overflow
ofConst = (2^32-1)/1000;
line = 0;

while ~feof(fid)
  line = line + 1;
  descr = fscanf(fid, '%c', 3);
  dat = str2num(fgetl(fid));
  if(line == 2)
      timePre = dat(1);
  end
  if(line > 2 && timePre - dat(1) > 4e+6)
    of = of+1;
  end

  switch descr(2:3)
     % trial description   
    case 'NW';  
        if(length(dat) ~= 5)%
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
    % treadmill  
    case 'TM';  
        if(length(dat) ~= 4)%
            continue;
        end
        nw = nw + 1;
		behEvents.treadmill(nw,:) = dat; 
        behEvents.treadmill(nw,1) = behEvents.treadmill(nw,1)+of*ofConst;
    % light beam    
    case 'BE';  %
        if(length(dat) ~= 3)
            continue;
        end
        nb = nb + 1;
		behEvents.beam(nb,:) = dat; 
        behEvents.beam(nb,1) = behEvents.beam(nb,1)+of*ofConst;
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
  end
  
  if(line > 2)
    timePre = dat(1);
  end
end

fclose(fid);

return