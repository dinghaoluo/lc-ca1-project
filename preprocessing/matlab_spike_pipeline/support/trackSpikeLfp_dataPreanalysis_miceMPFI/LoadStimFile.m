function [stimEvents] = LoadStimFile(FileName)

fid=fopen(FileName);
pu = 0;
st = 0;
ru = 0;
sto = 0;

stimEvents = [];

while ~feof(fid)
  line = fgetl(fid);
  ind = strfind(line,' ');
  if(length(ind) == 1)
      ind(2) = length(line)+1;
  end
  descr = line(ind(1)+1:ind(2)-1);
  dat = str2num(line(ind(2):end));
  time = str2double(line(1:ind(1)-1));

  switch descr(1:3)
     % pulse definition   
    case 'PUL';  
        if(length(dat) ~= 5)%
            continue;
        end
        pu = pu + 1;
		stimEvents.pulse(pu,:) = [time dat]; 
    % set diodes   
    case 'SET';  
        if(length(dat) ~= 6)%
            continue;
        end
        indNonZero = find(dat ~= 0,1);
        if(~isempty(indNonZero))
            st = st + 1;
            stimEvents.diode(st,:) = [time dat]; 
        end
    % run stimulation
    case 'RUN';  
        if(~isempty(dat))%
            continue;
        end
        ru = ru + 1;
		stimEvents.stim(ru,:) = time; 
    % stop stimulation
    case 'STOP';
        if(~isempty(dat))%
            continue;
        end
        sto = sto + 1;
		stimEvents.stop(sto,:) = time; 
  end
end

fclose(fid);

return