function [UpCrossings, DownCrossings] = getMovieOnTimeSmTr(baseFileName,...
                nChannelsTot,sampleRate,numMovies)

 % movie:  are we trying to detect movie pulses? 1: yes
 % 21 - file length (sec)
 % 26 - file size (bytes)
 
    %movieCh = str2double(metaInfo(38));    % the last ch is SYNC
    % ask for ID of the SYNC pulse
    resp2 = 'n';
    while strcmp(resp2, 'y') ~= 1
        resp1 = input(['\nWhich ch is the Movie signal ch? Remember to',...
                       ' report the neuronscope ch # plus 1. [67]']);
        if isempty(resp1)
            movieCh = 67;
        else movieCh = resp1;
        end
        resp2 = input(['\nIs Movie signal ch # ' num2str(movieCh) '? [y/n]'], 's');
        if isempty(resp2)
            resp2 = 'y';
        end
    end
    UpCrossings = [];
    DownCrossings = [];
    
    % changed 3/2/2017, line 22-44, reading the channel directly from .dat,
    % removed the dependence on .sev file
    datFileName = [baseFileName '.dat'];
    if exist(datFileName, 'file') == 2
        listing = dir(datFileName);
        Nsamples = listing.bytes/2/nChannelsTot;  % sec
        
        fprintf('\n Loading Movie signal from: %s  ', datFileName);        
        fid = fopen(datFileName,'r');
         % filter before downsampling
        movie = LoadDatFile_FL_uint16(fid, movieCh, Nsamples, sampleRate, nChannelsTot);
        movie = movie';
    end
    
    % the width of the visual stimulus
    if(~isempty(strfind(baseFileName,'A002-20181016-01')))
        numFrame = [1 6];
    elseif(~isempty(strfind(baseFileName,'A001')) || ...
            ~isempty(strfind(baseFileName,'A002')))
        numFrame = 10; 
    elseif(~isempty(strfind(baseFileName,'A004-20181030-01')))
        numFrame = 31;  
    elseif(~isempty(strfind(baseFileName,'A004')))
        numFrames = 13;
    elseif(~isempty(strfind(baseFileName,'A010')) ...
            || ~isempty(strfind(baseFileName,'A007')) ...
            || ~isempty(strfind(baseFileName,'A009')))
        numFrame = 6; 
    elseif(~isempty(strfind(baseFileName,'A011')) || ...
            ~isempty(strfind(baseFileName,'A012')))
        numFrame = 3;
    else
        numFrame = 3;
    end
    moviePulseLen = 1/60*[numFrame]; %s (10 frames 60 Hz refresh rate by default)
    
    % use a wide band pass filter to identify the final peak or trough
    [b,a]=butter(2,[15 900]/(sampleRate/2));
    lfMovie = filtfilt(b,a,movie);
    lfMovie = lfMovie - lfMovie(1);
    
    % use a narrow band pass filter to identify the range where the peak
    % or trough is located
    [b,a]=butter(2,[10 50]/(sampleRate/2));
    lfMovie1 = filtfilt(b,a,movie);
    lfMovie1 = lfMovie1 - lfMovie1(1);
    
    % movie pulse goes first down and then up!!!  
    % changed by Yingxue 20181113
    % identify peak and trough segment through thresholding
    ampN = min(lfMovie1); 
    ampP = max(lfMovie1);
    stdlfMovie = std(lfMovie1);
    thrPeakInit = 15; % initial threshold for peak detection
    thrTroughInit = 8; % initial threshold for trough detection
    thrPeakCur = thrPeakInit; % current threshould for peak detection
    thrTroughCur = thrTroughInit; % current threshold for trough detection
    maxIterPeak = 70; % max number of iterations for peak detection
    maxIterTrough = 55; % max number of iterations for trough detection
    iterStep = 0.2; % threshold change step size in each iteration
    peakStatus = 0; % whether correct number of peaks has been detected
    troughStatus = 0; % whether correct number of troughs has been detected
    
    for i = 1:maxIterPeak 
        if(abs(ampP) > thrPeakCur*stdlfMovie) % peak segments
            indP = lfMovie1 > thrPeakCur*stdlfMovie;
        else
            thrPeakCur = thrPeakCur - iterStep;
            if(thrPeakCur <= 0)
                disp('Movie detection peak threshold <= 0');
                break;
            end            
            continue;
        end
    
        % identify the start of the stimulus by calculating the peak location
        % of the wide band passed data. This is because the detector has a shorter
        % falling edge than rising edge, the absolute peak of a lfMovie1 pulse 
        % corresponds to the real onset of the stimulus.
        indPStart = find(diff([0; indP]) == 1);
        indPEnd = find(diff([0; indP]) == -1);
        lenIndPStart = length(indPStart);
        if(lenIndPStart == length(indPEnd))
            [peaks,indPeaks] = arrayfun(@(x,y) max(lfMovie(x:y)),...
                        indPStart,indPEnd);
            indPeaks = indPeaks + indPStart - 1;
            indTmp = [1; (diff(indPeaks) > max(moviePulseLen)*sampleRate*1.5)];
            indPeaks = indPeaks(indTmp == 1);
            disp(['MovieDetection ' num2str(i) ' numPeaks ' num2str(length(indPeaks))...
                [' PeakThre ' num2str(thrPeakCur)]]);
        end
        
        % modify detection threshold if less peaks are detected   
        if(length(indPeaks) < numMovies)
            thrPeakCur = thrPeakCur - iterStep;
            if(thrPeakCur <= 0)
                disp('Movie detection peak threshold <= 0');
                break;
            end                
        elseif(length(indPeaks == numMovies))
            peakStatus = 1; % found all the peaks
            break;
        else
            thrPeakCur = thrPeakCur + iterStep;
        end
    end
    
    if(peakStatus == 0)
        disp(['Can not find the same number of movie onset time points' ...
            ' as reported in the behavior file']);
        return;
    end

    for i = 1:maxIterTrough
        if(abs(ampN) > thrTroughCur*stdlfMovie) % trough segments
            indN = lfMovie1 < -1*thrTroughCur*stdlfMovie;
        else
            thrTroughCur = thrTroughCur - iterStep;
            continue;
        end
        
        % idenfify when the stimulus stops
        % first find the next trough location that is moviePulseLen * 0.9 away
        % from the identified peak
        % then, use the same method as above to calculate the offset of the
        % stimulus
        indNStart = find(diff([0; indN]) == 1);
        indNEnd = find(diff([0; indN]) == -1);
        lenIndNStart = length(indNStart);
        meanLen = mean(indNEnd-indNStart)/2;
        if(lenIndNStart == length(indNEnd))
            % deal with the "UniformOutput" error, which means no 
            % trough detected for the last peak
            try
                if(length(moviePulseLen) > 1)
                    selTroughStart = [];
                    selTroughEnd = [];
                    for n = 1:length(indPeaks)
                        selTroughStartTmp = find(indNStart >...
                            indPeaks(n), 2);
                        diffToPeak = indNStart(selTroughStartTmp)-indPeaks(n);
                        if(diffToPeak(1) < min(moviePulseLen)*0.7*sampleRate...
                            && diffToPeak(2) < max(moviePulseLen)*1.5*sampleRate...
                            && thrTroughCur > 1)
                            selTroughStart = [selTroughStart selTroughStartTmp(2)];
                            diffToP = diffToPeak(2);
                        else
                            selTroughStart = [selTroughStart selTroughStartTmp(1)];
                            diffToP = diffToPeak(1);
                        end
                        selTroughendTmp =  ...
                            find(indNEnd > indPeaks(n) + diffToP, 1);
                        selTroughEnd = [selTroughEnd selTroughendTmp];
                    end
                else
                    minPulseWidth = (moviePulseLen * 0.7)*sampleRate;
                    if(indNStart(end) - indPeaks(end) < minPulseWidth)
                        disp('Check the last movie pulse');
                        indNStart(end) = indNStart(end) + minPulseWidth;
                        indNEnd(end) = indNEnd(end) + minPulseWidth;
                    end
                    selTroughStart = arrayfun(@(y) ...
                        find(indNStart > y + minPulseWidth, 1),...
                        indPeaks);
                    
                    selTroughEnd = arrayfun(@(y) ...
                        find(indNEnd > ...
                        y + minPulseWidth + meanLen, 1),...
                        indPeaks);
                end
                
                indNSelStart = indNStart(selTroughStart);
                indNSelEnd = indNEnd(selTroughEnd);
                [trough,indTroughs] = arrayfun(@(x,y) min(lfMovie(x:y)),...
                                            indNSelStart,indNSelEnd);
                indTroughs = indTroughs + indNSelStart - 1;
                indTmp = [1; (diff(indTroughs) > max(moviePulseLen)*sampleRate)];
                indTroughs = indTroughs(indTmp == 1);
                
                % decrease trough detection threshold until detect the correct number 
                % of movie pulses
                if(length(indTroughs) == length(indPeaks))
                    pulseWidth = indTroughs - indPeaks + 1;
                    for nPulLen = 1:length(moviePulseLen)
                        if(nPulLen == 1)
                            indIncorr = find(pulseWidth > moviePulseLen(nPulLen) ...
                                * sampleRate * 1.2 | ...
                                pulseWidth < moviePulseLen(nPulLen) * sampleRate * 0.8);  
                        else
                            indIncorrTmp = find(pulseWidth(indIncorr) > ...
                                moviePulseLen(nPulLen) * sampleRate * 1.2 | ...
                                pulseWidth(indIncorr) < moviePulseLen(nPulLen) ...
                                * sampleRate * 0.8); 
                            indIncorr = indIncorr(indIncorrTmp);
                        end
                    end

                    if(isempty(indIncorr))
                        troughStatus = 1;
                        break;
                    else
                        thrTroughCur = thrTroughCur - iterStep;
                    end
                else
                    indIncorr = [];
                    thrTroughCur = thrTroughCur - iterStep;    
                end
                disp(['Movie detection ' num2str(i) ...
                        ' numIncorr ' num2str(length(indIncorr))...
                        ' numTroughs ' num2str(length(indTroughs))...
                        [' TroughThre ' num2str(thrTroughCur)]]);
                if(length(indTroughs) == 1133)
                    a = 1;
                end
                if(thrTroughCur <= 0)
                    disp('Movie detection trough threshold <= 0');
                    break;
                end
                
            catch ME
                thrTroughCur = thrTroughCur - iterStep; 
                disp(['TroughThre ' num2str(thrTroughCur)]);
                if(thrTroughCur <= 0)
                    disp('Movie detection trough threshold <= 0');
                    return;
                end   
                continue;
            end
        end
    end
    
    if(troughStatus == 0)
        disp(['Can not find the same number of movie offset time points' ...
            ' to match with the movie onsets']);
        return;
    end
            
    DownCrossings = indPeaks;
    UpCrossings = indTroughs;

end
