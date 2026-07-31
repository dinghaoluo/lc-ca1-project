function check = Arduino2RecordT_smTr_pix(baseFileName,sampleFreq,lfpFreq,nChannelsTot)
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
            % LoadBehMazeFile_smTr_pix, Dinghao, 16 Mar 2023
            % changed diodePar to report dummy values since pix does not 
            % report DS
            behEvents = LoadBehMazeFile_smTr_pix([baseFileName 'T.txt']);

            save([baseFileName 'B.mat'], 'behEvents');
        end
        
        if(~isfield(behEvents,'TDTsyncInd'))
            % get the recording time of the sync pulses            
            % changed 11 Mar 23, totChans = 3 (sync + mov + pulse), Dinghao
            if(baseFileName=="A069r-20230908-01" || baseFileName=="A069r-20230908-02")
                [UpCrossings] = getSyncPulseTimeSmTr_pix('analogin_unstable',3,...
                    sampleFreq);
            else
                [UpCrossings] = getSyncPulseTimeSmTr_pix('analogin',3,...
                    sampleFreq);
            end
            
            if(length(UpCrossings) ~= length(behEvents.ArdSyncMsec(:,1)))
                disp(['Unmatching numbers of sync pulses between the '...
                    'behavior file and the recording. Please check.']);
                check = 1;
            end
            
            % store the sync pulse recording time into Arduino behavior file
            behEvents.TDTsyncInd = UpCrossings;
            behEvents.TDTsyncMsec = UpCrossings / sampleFreq * 1000;
            
            save(fullNameB,'behEvents');
        end
        
        %% stimulation pulses 
        %%%% added by yingxue on 11/1/2019
        if(~isfield(behEvents,'TDTstimInd') && isfield(behEvents,'pulsePar'))
            % get the recording time of the sync pulses            
            % 11 Mar 23 same as above, totChans = 3, Dinghao 
            if(baseFileName=="A069r-20230908-01" || baseFileName=="A069r-20230908-02")
                [UpCrossings, DownCrossings] = getStimPulseTime_pix('analogin_unstable',3,...
                    sampleFreq);
            else
                [UpCrossings, DownCrossings] = getStimPulseTime_pix('analogin',3,...
                    sampleFreq); 
            end
          
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
        %%%%%% added by Yingxue on 03/30/2019
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
                      getMovieOnTimeSmTr('analogin',3,sampleFreq,...
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
        Arduino2TDTtime_smTr_opto(baseFileName,sampleFreq,lfpFreq);
         
    end
