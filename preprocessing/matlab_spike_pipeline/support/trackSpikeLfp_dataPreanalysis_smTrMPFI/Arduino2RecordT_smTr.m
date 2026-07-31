function check = Arduino2RecordT_smTr(baseFileName,sampleFreq,lfpFreq,nChannelsTot)
% Convert arduino time to recording time
    
    check = 0; 
    
    if exist([baseFileName 'BTDT.mat'], 'file') == 2
        disp('BTDT file already exists.')
        return;
    else
        fullNameB = [baseFileName 'B.mat'];
        if(exist(fullNameB,'file') ~= 0)
            load(fullNameB);
        else
            disp('parsing the behavioral file.');
            behEvents = LoadBehMazeFile_smTr([baseFileName 'T.txt']);
                      
            save([baseFileName 'B.mat'], 'behEvents');
        end
        
        if(~isfield(behEvents,'TDTsyncInd'))
            % get the recording time of the sync pulses            
            [UpCrossings] = getSyncPulseTimeSmTr(baseFileName,nChannelsTot,...
                sampleFreq); 
          
            % store the sync pulse recording time into Arduino behavior file
            behEvents.TDTsyncInd = UpCrossings;
            behEvents.TDTsyncMsec = UpCrossings / sampleFreq * 1000;
            
            if(length(UpCrossings) ~= length(behEvents.ArdSyncMsec(:,1)))
                disp(['Unmatching numbers of sync pulses between the '...
                    'behavior file and the recording. Please check.']);
                check = 1;
                return;
            end
            save(fullNameB,'behEvents');
        end
        
        %% stimulation pulses 
        %%%% added by yingxue on 11/1/2019
        if(~isfield(behEvents,'TDTstimInd') && isfield(behEvents,'pulsePar'))
            % get the recording time of the sync pulses            
            [UpCrossings, DownCrossings] = getStimPulseTime(baseFileName,nChannelsTot,...
                sampleFreq); 
          
            % store the sync pulse recording time into Arduino behavior file
            behEvents.TDTstimInd = UpCrossings;
            behEvents.TDTstimMsec = UpCrossings / sampleFreq * 1000;
            behEvents.TDTpulseWidthMsec = (DownCrossings - UpCrossings) / ...
                sampleFreq * 1000;
            
            [behEvents,check] = CheckStimPulses(behEvents,UpCrossings,DownCrossings,...
                                                sampleFreq);
            
            if(check == 1)
                return;
            end
            save(fullNameB,'behEvents');
        end
        
        
        %%
        %%%%%% added by Yingxue on 03/30/3019
        if(~isfield(behEvents,'TDTMovieOnInd') && isfield(behEvents,'movieOn'))        
            % get the recording time of actual movie on time
              anmNum = str2num(baseFileName(2:4));
              if(anmNum >= 13) 
                  %% MV = 1 grey; MV = 2 black out; MV = 3 left start cue; 
                  % MV = 4 right start cue
                  % count the number of movies displayed on the left and
                  % right screens
                  behEvents.movieNumL = 0;
                  behEvents.movieNumR = 0;
                  behEvents.startCueEndDist = zeros(size(behEvents.taskDescr,1),1);
                  behEvents.movieID = zeros(size(behEvents.movieOn,1),1);
                  behEvents.trStartCueLR = zeros(size(behEvents.taskDescr,1),1);
                  numTr = 1;
                  for i = 1:size(behEvents.movieOn,1)
                        if(behEvents.movieOn(i,2) == 1 || ...
                                behEvents.movieOn(i,2) == 2)
                            behEvents.movieNumL = behEvents.movieNumL + 1;
                            behEvents.movieNumR = behEvents.movieNumR + 1;
                            if(behEvents.movieOn(i,2) == 1)
                                behEvents.movieID(i) = 2; % grey
                                %%%% the distance at which the start cue ends
                                %%%% added by Yingxue on 4/11/2019
                                if(behEvents.movieOn(i,1) >= behEvents.taskDescr(1,1))
                                    indDist = find(behEvents.wheel(:,1) >= ...
                                        behEvents.movieOn(i,1),1,'first');
%                                     indDistLast = find(behEvents.wheel(:,1) < ...
%                                         behEvents.movieOn(i,1),1,'last');                                    
                                    behEvents.startCueEndDist(numTr) = ...
                                        behEvents.wheel(indDist,3) * 0.04;
                                    numTr = numTr+1;
                                end
                                %%%%
                            else
                                behEvents.movieID(i) = 3; % black out
                            end
                        elseif(behEvents.movieOn(i,2) == 3)
                            behEvents.movieNumL = behEvents.movieNumL + 1;
                            indTr = find(behEvents.taskDescr(:,1) <= ...
                                behEvents.movieOn(i,1),1,'last');
                            if(~isempty(indTr))
                                behEvents.trStartCueLR(indTr) = 1; % left trial
                            end
                            behEvents.movieID(i) = 1;
                        elseif(behEvents.movieOn(i,2) == 4)
                            behEvents.movieNumR = behEvents.movieNumR + 1;
                            indTr = find(behEvents.taskDescr(:,1) <= ...
                                behEvents.movieOn(i,1),1,'last');
                            if(~isempty(indTr))
                                behEvents.trStartCueLR(indTr) = 0; % right trial
                            end
                            behEvents.movieID(i) = 0;
                        end                            
                  end
                  
%                   %% detect movie pulsed on the left screen
%                   [UpCrossings, DownCrossings] = ...
%                       getMovieOnTimeSmTr(baseFileName,nChannelsTot,sampleFreq,...
%                       behEvents.movieNumL);
%                   
%                   if(isempty(UpCrossings))
%                       disp(['Unmatching numbers of movie pulses on the left screen between the '...
%                           'behavior file and the recording. Please check.']);
%                       check = 1;
%                       return;
%                   end
%                   
%                   % store the movie pulse recording time into Arduino behavior file
%                   behEvents.TDTMovieOnLInd = DownCrossings;
%                   behEvents.TDTMovieOnL = DownCrossings / sampleFreq * 1000;
%                   behEvents.TDTMovieOffLInd = UpCrossings;
%                   behEvents.TDTMovieOffL = UpCrossings / sampleFreq * 1000;
%                   
%                   %% detect movie pulsed on the right screen
%                   [UpCrossings, DownCrossings] = ...
%                       getMovieOnTimeSmTr(baseFileName,nChannelsTot,sampleFreq,...
%                       behEvents.movieNumR);
%                   
%                   if(isempty(UpCrossings))
%                       disp(['Unmatching numbers of movie pulses on the right screen between the '...
%                           'behavior file and the recording. Please check.']);
%                       check = 1;
%                       return;
%                   end
%                   
%                   % store the movie pulse recording time into Arduino behavior file
%                   behEvents.TDTMovieOnRInd = DownCrossings;
%                   behEvents.TDTMovieOnR = DownCrossings / sampleFreq * 1000;
%                   behEvents.TDTMovieOffRInd = UpCrossings;
%                   behEvents.TDTMovieOffR = UpCrossings / sampleFreq * 1000;
              %%%%%%%
              
              else
                  numMovies = length(behEvents.movieOn(:,1));
                  [UpCrossings, DownCrossings] = ...
                      getMovieOnTimeSmTr(baseFileName,nChannelsTot,sampleFreq,...
                      numMovies);
                  
                  if(isempty(UpCrossings))
                      disp(['Unmatching numbers of movie pulses between the '...
                          'behavior file and the recording. Please check.']);
                      check = 1;
                      return;
                  end
                  
                  % store the movie pulse recording time into Arduino behavior file
                  behEvents.TDTMovieOnInd = DownCrossings;
                  behEvents.TDTMovieOn = DownCrossings / sampleFreq * 1000;
                  behEvents.TDTMovieOffInd = UpCrossings;
                  behEvents.TDTMovieOff = UpCrossings / sampleFreq * 1000;
              end
              save(fullNameB,'behEvents');
        end
        
        % convert arduino time to recording time
        if(~isempty(strfind(baseFileName,'A009-20190111')) || ...
                ~isempty(strfind(baseFileName,'A009-20190112')))
            Arduino2TDTtime_smTr_A009(baseFileName,sampleFreq,lfpFreq); 
        elseif(~isempty(strfind(baseFileName,'A001')) || ...
                ~isempty(strfind(baseFileName,'A002')))
            Arduino2TDTtime_smTr_A002(baseFileName,sampleFreq,lfpFreq); 
        elseif(~isempty(strfind(baseFileName,'A007'))) 
            Arduino2TDTtime_smTr_A007(baseFileName,sampleFreq,lfpFreq); 
        elseif(~isempty(strfind(baseFileName,'A004')))
            Arduino2TDTtime_smTr_A004(baseFileName,sampleFreq,lfpFreq);
        elseif(~isempty(strfind(baseFileName,'A011-20190215')) || ...
                ~isempty(strfind(baseFileName,'A011-20190218')) || ...
                ~isempty(strfind(baseFileName,'A011-20190219')) || ...
                ~isempty(strfind(baseFileName,'A012-20190220')))
            Arduino2TDTtime_smTr_A011_20190215(baseFileName,sampleFreq,lfpFreq); 
        elseif(~isempty(strfind(baseFileName,'A011-20190220')))
            Arduino2TDTtime_smTr_A011_20190220(baseFileName,sampleFreq,lfpFreq); 
        elseif(~isempty(strfind(baseFileName,'A012-20190225')))
            Arduino2TDTtime_smTr_A012_20190225(baseFileName,sampleFreq,lfpFreq);
        elseif(~isempty(strfind(baseFileName,'A011')) || ...
                ~isempty(strfind(baseFileName,'A012')))
            Arduino2TDTtime_smTr_A011(baseFileName,sampleFreq,lfpFreq);
        elseif(~isempty(strfind(baseFileName,'A015')) || ...
               ~isempty(strfind(baseFileName,'A014')) || ...
               ~isempty(strfind(baseFileName,'A013')) || ...
               ~isempty(strfind(baseFileName,'A016')))
            Arduino2TDTtime_smTr_A015(baseFileName,sampleFreq,lfpFreq);
        elseif(~isempty(strfind(baseFileName,'A022')) || ...
                ~isempty(strfind(baseFileName,'A023')) || ...
                ~isempty(strfind(baseFileName,'A024')) || ...
                ~isempty(strfind(baseFileName,'A025')) || ...
                ~isempty(strfind(baseFileName,'A029')) || ...
                ~isempty(strfind(baseFileName,'A028')))
            Arduino2TDTtime_smTr_opto(baseFileName,sampleFreq,lfpFreq);
         elseif(~isempty(strfind(baseFileName,'A027')) || ...
                ~isempty(strfind(baseFileName,'A026')))
            Arduino2TDTtime_smTr_passiveCue(baseFileName,sampleFreq,lfpFreq);
        else
            Arduino2TDTtime_smTr(baseFileName,sampleFreq,lfpFreq); 
        end
    end
