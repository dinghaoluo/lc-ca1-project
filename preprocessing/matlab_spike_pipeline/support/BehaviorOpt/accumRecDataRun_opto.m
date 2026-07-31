function [recDataRunPre,recDataRunManip,recDataRunPost,...
        recDataRunPreOpt,recDataRunManipOpt,recDataRunPostOpt,...
        recDataRunPreOptCtrl,recDataRunManipOptCtrl,recDataRunPostOptCtrl] = ...
                accumRecDataRun_opto(path,recSess,manipSess,mazeSess,pulseW,stimCode)
    nRec = size(path,1);
    tmp.indRec = [];
    tmp.indSess = [];
    tmp.numSamplesMean = [];
    tmp.maxSpeedMean = [];
    tmp.meanSpeedMean = [];
    tmp.maxRunLenTMean = [];
    tmp.totRunLenTMean = [];
    tmp.numRunMean = [];
    tmp.maxAccMean = [];
    tmp.meanAccMean = [];
    tmp.totStopLenTMean = [];
    tmp.startCueToRunMean = [];
    tmp.numLicksBefRewMean = [];
    tmp.numLicksRewMean = [];
    tmp.med1stFiveLickDistMean = [];
    tmp.medLickDistMean = [];
    tmp.med1stFiveLickDistBefRewMean = [];
    tmp.medLickDistBefRewMean = [];
    tmp.percRewarded = [];
    tmp.percNonStop = [];
%     tmp.speedProfile = [];
%     tmp.lickProfile = [];
    tmp.pumpLfpIndMean = [];
    tmp.pumpMMMean = [];
    tmp.speedProfile = [];
    tmp.lickProfile = [];
    tmp.speedProfileRndSel = [];
    tmp.lickProfileRndSel = [];
    
    tmp.meanSpeedOverTimeRunBefRun = [];
    tmp.meanSpeedOverTimeRun0to1 = [];
    tmp.meanSpeedOverTimeRun3to5 = [];

    tmp.meanSpeedOverDistRun0to100 = [];
    tmp.meanSpeedOverDistRunAfter100 = [];

    tmp.meanRun30to100 = [];
    tmp.meanRun100to150 = [];
    tmp.meanRun150to180 = [];
    tmp.meanRun180to210 = [];
    
    %% added by Yingxue on 4/7/2021
    tmp.speedSimMean = [];
    tmp.lickSimMean = [];
    tmp.speedEucMean = [];
    tmp.lickEucMean = [];
    tmp.speedStdMean = [];
    tmp.lickStdMean = [];
    %%
    
    recDataRunPre = tmp;
    recDataRunManip = tmp;
    recDataRunPost = tmp;
    
    recDataRunPreOpt = tmp;
    recDataRunManipOpt = tmp;
    recDataRunPostOpt = tmp;
    
    recDataRunPreOptCtrl = tmp;
    recDataRunManipOptCtrl = tmp;
    recDataRunPostOptCtrl = tmp;
    
    for i = 1:nRec
        if(i == 9)
            a = 1;
        end
        ind = findstr(path(i,:),'\');
        recName = path(i,ind(end-1)+1:ind(end)-1);
        disp(recName);
        fullPath = [path(i,:) recName '_compSess.mat'];
        if(~exist(fullPath,'file'))
            disp([fullPath ' does not exist.']);
            return;
        end
        load(fullPath);
        
        if(stimCode == 5 || stimCode == 6) % stimulation across the running period or only around the reward
            indManip = find(manipSess{i} == stimCode);
        else
            indManip = find(manipSess{i} == stimCode); % & pulseW{i} == 3);
        end
        for j = 1:length(indManip)
            indPre = find(manipSess{i}(1:indManip(j)) == 0,1,'last');
%             if(~isempty(indPre))
%                 a= 1;
%             end
            if(isempty(indPre))
                indPre = NaN;
            end
            indPost = find(manipSess{i}(indManip(j)+1:end) == 0,1,'first')+indManip(j);
            if(isempty(indPost))
                indPost = NaN;
            end
            recDataRunPre = concatData(recDataRunPre,sessDataRun,...
                        sessDataLick,sessDataSpeed,indPre,i);
            recDataRunPost = concatData(recDataRunPost,sessDataRun,...
                        sessDataLick,sessDataSpeed,indPost,i);
            recDataRunManip = concatData(recDataRunManip,sessDataRun,...
                        sessDataLick,sessDataSpeed,indManip(j),i);
                                        
            recDataRunPreOpt = concatDataOpt(recDataRunPreOpt,sessDataRunOpt,...
                        sessDataLickOpt,sessDataSpeedOpt,indPre,i,1);
            recDataRunPostOpt = concatDataOpt(recDataRunPostOpt,sessDataRunOpt,...
                        sessDataLickOpt,sessDataSpeedOpt,indPost,i,1);
            recDataRunManipOpt = concatDataOpt(recDataRunManipOpt,sessDataRunOpt,...
                        sessDataLickOpt,sessDataSpeedOpt,indManip(j),i,1);
                    
            recDataRunPreOptCtrl = concatDataOpt(recDataRunPreOptCtrl,sessDataRunOpt,...
                        sessDataLickOpt,sessDataSpeedOpt,indPre,i,2);
            recDataRunPostOptCtrl = concatDataOpt(recDataRunPostOptCtrl,sessDataRunOpt,...
                        sessDataLickOpt,sessDataSpeedOpt,indPost,i,2);
            recDataRunManipOptCtrl = concatDataOpt(recDataRunManipOptCtrl,sessDataRunOpt,...
                        sessDataLickOpt,sessDataSpeedOpt,indManip(j),i,2);
                    
        end
    end

end

function recData = concatData(recData,sessData,sessDataLick,sessDataSpeed,indSess,indRec)
    randTrNo = 35;
    nSess = length(indSess);
    if(sum(isnan(indSess)) >= 1)
        recData.indRec = [recData.indRec NaN*ones(1,nSess)];
        recData.indSess = [recData.indSess NaN*ones(1,nSess)];
        recData.numSamplesMean = [recData.numSamplesMean NaN*ones(1,nSess)];
        recData.maxSpeedMean = [recData.maxSpeedMean NaN*ones(1,nSess)];
        recData.meanSpeedMean = [recData.meanSpeedMean NaN*ones(1,nSess)];
        recData.maxRunLenTMean = [recData.maxRunLenTMean NaN*ones(1,nSess)];
        recData.totRunLenTMean = [recData.totRunLenTMean NaN*ones(1,nSess)];
        recData.numRunMean = [recData.numRunMean NaN*ones(1,nSess)];
        recData.maxAccMean = [recData.maxAccMean NaN*ones(1,nSess)];
        recData.meanAccMean = [recData.meanAccMean NaN*ones(1,nSess)];
        recData.totStopLenTMean = [recData.totStopLenTMean NaN*ones(1,nSess)];
        recData.startCueToRunMean = [recData.startCueToRunMean NaN*ones(1,nSess)];
        recData.numLicksBefRewMean = [recData.numLicksBefRewMean NaN*ones(1,nSess)];
        recData.numLicksRewMean = [recData.numLicksRewMean NaN*ones(1,nSess)];
        recData.med1stFiveLickDistMean = [recData.med1stFiveLickDistMean NaN*ones(1,nSess)];
        recData.medLickDistMean = [recData.medLickDistMean NaN*ones(1,nSess)];
        recData.med1stFiveLickDistBefRewMean = [recData.med1stFiveLickDistBefRewMean NaN*ones(1,nSess)];
        recData.medLickDistBefRewMean = [recData.medLickDistBefRewMean NaN*ones(1,nSess)];
        recData.percRewarded = [recData.percRewarded NaN*ones(1,nSess)];
        recData.percNonStop = [recData.percNonStop NaN*ones(1,nSess)];
        recData.pumpLfpIndMean = [recData.pumpLfpIndMean NaN*ones(1,nSess)];
        recData.pumpMMMean = [recData.pumpMMMean NaN*ones(1,nSess)];
        
        %% added by Yingxue on 4/7/2021
        recData.speedSimMean = [recData.speedSimMean NaN*ones(1,nSess)];
        recData.lickSimMean = [recData.lickSimMean NaN*ones(1,nSess)];
        
        recData.speedEucMean = [recData.speedEucMean NaN*ones(1,nSess)];
        recData.lickEucMean = [recData.lickEucMean NaN*ones(1,nSess)];
        
        recData.speedStdMean = [recData.speedStdMean NaN*ones(1,nSess)];
        recData.lickStdMean = [recData.lickStdMean NaN*ones(1,nSess)];
    else
        recData.indRec = [recData.indRec indRec];
        recData.indSess = [recData.indSess indSess];
        recData.numSamplesMean = [recData.numSamplesMean sessData.numSamplesMean(indSess)];
        recData.maxSpeedMean = [recData.maxSpeedMean sessData.maxSpeedMean(indSess)];
        recData.meanSpeedMean = [recData.meanSpeedMean sessData.meanSpeedMean(indSess)];
        recData.maxRunLenTMean = [recData.maxRunLenTMean sessData.maxRunLenTMean(indSess)];
        recData.totRunLenTMean = [recData.totRunLenTMean sessData.totRunLenTMean(indSess)];
        recData.numRunMean = [recData.numRunMean sessData.numRunMean(indSess)];
        recData.maxAccMean = [recData.maxAccMean sessData.maxAccMean(indSess)];
        recData.meanAccMean = [recData.meanAccMean sessData.meanAccMean(indSess)];
        recData.totStopLenTMean = [recData.totStopLenTMean sessData.totStopLenTMean(indSess)];
        recData.startCueToRunMean = [recData.startCueToRunMean sessData.startCueToRunMean(indSess)];
        recData.numLicksBefRewMean = [recData.numLicksBefRewMean sessData.numLicksBefRewMean(indSess)];
        recData.numLicksRewMean = [recData.numLicksRewMean sessData.numLicksRewMean(indSess)];
        recData.med1stFiveLickDistMean = [recData.med1stFiveLickDistMean sessData.med1stFiveLickDistMean(indSess)];
        recData.medLickDistMean = [recData.medLickDistMean sessData.medLickDistMean(indSess)];
        recData.med1stFiveLickDistBefRewMean = [recData.med1stFiveLickDistBefRewMean sessData.med1stFiveLickDistBefRewMean(indSess)];
        recData.medLickDistBefRewMean = [recData.medLickDistBefRewMean sessData.medLickDistBefRewMean(indSess)];
        recData.percRewarded = [recData.percRewarded sessData.percRewarded(indSess)];
        recData.percNonStop = [recData.percNonStop sessData.percNonStop(indSess)];
        recData.spaceStepsLick = sessDataLick.spaceSteps;
        recData.spaceStepsSpeed = sessDataSpeed.spaceSteps;
        sess = unique(indSess);
        for i = 1:length(sess)
            recData.speedProfile = [recData.speedProfile; sessDataSpeed.meanRun{sess(i)}];
            recData.lickProfile = [recData.lickProfile; sessDataLick.meanRun{sess(i)}];
            
            indTrR = randperm(size(sessDataSpeed.Run{sess(i)},1));
            if(length(indTrR) > randTrNo)
                recData.speedProfileRndSel = [recData.speedProfileRndSel;...
                    sessDataSpeed.Run{sess(i)}(indTrR(1:randTrNo),:)];
                recData.lickProfileRndSel = [recData.lickProfileRndSel;...
                    sessDataLick.Run{sess(i)}(indTrR(1:randTrNo),:)];
            else
                recData.speedProfileRndSel = [recData.speedProfileRndSel;...
                    sessDataSpeed.Run{sess(i)}];
                recData.lickProfileRndSel = [recData.lickProfileRndSel;...
                    sessDataLick.Run{sess(i)}];
            end
        end
        recData.pumpLfpIndMean = [recData.pumpLfpIndMean sessData.pumpLfpIndMean(indSess)];
        recData.pumpMMMean = [recData.pumpMMMean sessData.pumpMMMean(indSess)];
        
        %% added by Yingxue on 4/7/2021
        recData.speedSimMean = [recData.speedSimMean sessData.speedSimMean(indSess)];
        recData.lickSimMean = [recData.lickSimMean sessData.lickSimMean(indSess)];
        
        recData.speedEucMean = [recData.speedEucMean sessData.speedEucMean(indSess)];
        recData.lickEucMean = [recData.lickEucMean sessData.lickEucMean(indSess)];
        
        recData.speedStdMean = [recData.speedStdMean mean(sessDataSpeed.stdCue{indSess})];
        recData.lickStdMean = [recData.lickStdMean mean(sessDataLick.stdCue{indSess})];
        %%
        
        recData.meanSpeedOverTimeRunBefRun = [recData.meanSpeedOverTimeRunBefRun ...
            sessDataSpeed.meanSpeedOverTimeRunBefRun{indSess}'];
        recData.meanSpeedOverTimeRun0to1 = [recData.meanSpeedOverTimeRun0to1 ...
            sessDataSpeed.meanSpeedOverTimeRun0to1{indSess}'];
        recData.meanSpeedOverTimeRun3to5 = [recData.meanSpeedOverTimeRun3to5 ...
            sessDataSpeed.meanSpeedOverTimeRun3to5{indSess}'];
        
        recData.meanSpeedOverDistRun0to100 = [recData.meanSpeedOverDistRun0to100 ...
            sessDataSpeed.meanSpeedOverDistRun0to100{indSess}'];
        recData.meanSpeedOverDistRunAfter100 = [recData.meanSpeedOverDistRunAfter100 ...
            sessDataSpeed.meanSpeedOverDistRunAfter100{indSess}'];
        
        recData.meanRun30to100 = [recData.meanRun30to100 ...
            sessDataLick.meanRun30to100{indSess}'];
        recData.meanRun100to150 = [recData.meanRun100to150 ...
            sessDataLick.meanRun100to150{indSess}'];
        recData.meanRun150to180 = [recData.meanRun150to180 ...
            sessDataLick.meanRun150to180{indSess}'];
        recData.meanRun180to210 = [recData.meanRun180to210 ...
            sessDataLick.meanRun180to210{indSess}'];
    end
end

function recData = concatDataOpt(recData,sessData,sessDataLick,sessDataSpeed,indSess,indRec,opt)
    randTrNo = 35;
    nSess = length(indSess);
    if(sum(isnan(indSess)) >= 1)
        recData.indRec = [recData.indRec NaN*ones(1,nSess)];
        recData.indSess = [recData.indSess NaN*ones(1,nSess)];
        recData.numSamplesMean = [recData.numSamplesMean NaN*ones(1,nSess)];
        recData.maxSpeedMean = [recData.maxSpeedMean NaN*ones(1,nSess)];
        recData.meanSpeedMean = [recData.meanSpeedMean NaN*ones(1,nSess)];
        recData.maxRunLenTMean = [recData.maxRunLenTMean NaN*ones(1,nSess)];
        recData.totRunLenTMean = [recData.totRunLenTMean NaN*ones(1,nSess)];
        recData.numRunMean = [recData.numRunMean NaN*ones(1,nSess)];
        recData.maxAccMean = [recData.maxAccMean NaN*ones(1,nSess)];
        recData.meanAccMean = [recData.meanAccMean NaN*ones(1,nSess)];
        recData.totStopLenTMean = [recData.totStopLenTMean NaN*ones(1,nSess)];
        recData.startCueToRunMean = [recData.startCueToRunMean NaN*ones(1,nSess)];
        recData.numLicksBefRewMean = [recData.numLicksBefRewMean NaN*ones(1,nSess)];
        recData.numLicksRewMean = [recData.numLicksRewMean NaN*ones(1,nSess)];
        recData.med1stFiveLickDistMean = [recData.med1stFiveLickDistMean NaN*ones(1,nSess)];
        recData.medLickDistMean = [recData.medLickDistMean NaN*ones(1,nSess)];
        recData.med1stFiveLickDistBefRewMean = [recData.med1stFiveLickDistBefRewMean NaN*ones(1,nSess)];
        recData.medLickDistBefRewMean = [recData.medLickDistBefRewMean NaN*ones(1,nSess)];
        recData.percRewarded = [recData.percRewarded NaN*ones(1,nSess)];
        recData.percNonStop = [recData.percNonStop NaN*ones(1,nSess)];
        recData.pumpLfpIndMean = [recData.pumpLfpIndMean NaN*ones(1,nSess)];
        recData.pumpMMMean = [recData.pumpMMMean NaN*ones(1,nSess)];
        
        %% added by Yingxue on 4/7/2021
        recData.speedSimMean = [recData.speedSimMean NaN*ones(1,nSess)];
        recData.lickSimMean = [recData.lickSimMean NaN*ones(1,nSess)];
        
        recData.speedEucMean = [recData.speedEucMean NaN*ones(1,nSess)];
        recData.lickEucMean = [recData.lickEucMean NaN*ones(1,nSess)];
        
        recData.speedStdMean = [recData.speedStdMean NaN*ones(1,nSess)];
        recData.lickStdMean = [recData.lickStdMean NaN*ones(1,nSess)];
    else
        recData.indRec = [recData.indRec indRec];
        recData.indSess = [recData.indSess indSess];
        recData.numSamplesMean = [recData.numSamplesMean sessData.numSamplesMean(opt,indSess)];
        recData.maxSpeedMean = [recData.maxSpeedMean sessData.maxSpeedMean(opt,indSess)];
        recData.meanSpeedMean = [recData.meanSpeedMean sessData.meanSpeedMean(opt,indSess)];
        recData.maxRunLenTMean = [recData.maxRunLenTMean sessData.maxRunLenTMean(opt,indSess)];
        recData.totRunLenTMean = [recData.totRunLenTMean sessData.totRunLenTMean(opt,indSess)];
        recData.numRunMean = [recData.numRunMean sessData.numRunMean(opt,indSess)];
        recData.maxAccMean = [recData.maxAccMean sessData.maxAccMean(opt,indSess)];
        recData.meanAccMean = [recData.meanAccMean sessData.meanAccMean(opt,indSess)];
        recData.totStopLenTMean = [recData.totStopLenTMean sessData.totStopLenTMean(opt,indSess)];
        recData.startCueToRunMean = [recData.startCueToRunMean sessData.startCueToRunMean(opt,indSess)];
        recData.numLicksBefRewMean = [recData.numLicksBefRewMean sessData.numLicksBefRewMean(opt,indSess)];
        recData.numLicksRewMean = [recData.numLicksRewMean sessData.numLicksRewMean(opt,indSess)];
        recData.med1stFiveLickDistMean = [recData.med1stFiveLickDistMean sessData.med1stFiveLickDistMean(opt,indSess)];
        recData.medLickDistMean = [recData.medLickDistMean sessData.medLickDistMean(opt,indSess)];
        recData.med1stFiveLickDistBefRewMean = [recData.med1stFiveLickDistBefRewMean sessData.med1stFiveLickDistBefRewMean(opt,indSess)];
        recData.medLickDistBefRewMean = [recData.medLickDistBefRewMean sessData.medLickDistBefRewMean(opt,indSess)];
        recData.percRewarded = [recData.percRewarded sessData.percRewarded(opt,indSess)];
        recData.percNonStop = [recData.percNonStop sessData.percNonStop(opt,indSess)];
        recData.spaceStepsLick = sessDataLick.spaceSteps;
        recData.spaceStepsSpeed = sessDataSpeed.spaceSteps;
        sess = unique(indSess);
        for i = 1:length(sess)
            recData.speedProfile = [recData.speedProfile; sessDataSpeed.meanRun{opt,sess(i)}];
            recData.lickProfile = [recData.lickProfile; sessDataLick.meanRun{opt,sess(i)}];
            
            indTrR = randperm(size(sessDataSpeed.Run{opt,sess(i)},1));
            if(length(indTrR) > randTrNo)
                recData.speedProfileRndSel = [recData.speedProfileRndSel;...
                    sessDataSpeed.Run{opt,sess(i)}(indTrR(1:randTrNo),:)];
                recData.lickProfileRndSel = [recData.lickProfileRndSel;...
                    sessDataLick.Run{opt,sess(i)}(indTrR(1:randTrNo),:)];
            else
                recData.speedProfileRndSel = [recData.speedProfileRndSel;...
                    sessDataSpeed.Run{opt,sess(i)}];
                recData.lickProfileRndSel = [recData.lickProfileRndSel;...
                    sessDataLick.Run{opt,sess(i)}];
            end
        end
        recData.pumpLfpIndMean = [recData.pumpLfpIndMean sessData.pumpLfpIndMean(opt,indSess)];
        recData.pumpMMMean = [recData.pumpMMMean sessData.pumpMMMean(opt,indSess)];
        
        %% added by Yingxue on 4/7/2021
        recData.speedSimMean = [recData.speedSimMean sessData.speedSimMean(opt,indSess)];
        recData.lickSimMean = [recData.lickSimMean sessData.lickSimMean(opt,indSess)];
        
        recData.speedEucMean = [recData.speedEucMean sessData.speedEucMean(opt,indSess)];
        recData.lickEucMean = [recData.lickEucMean sessData.lickEucMean(opt,indSess)];
        
        recData.speedStdMean = [recData.speedStdMean mean(sessDataSpeed.stdRun{opt,indSess})];
        recData.lickStdMean = [recData.lickStdMean mean(sessDataLick.stdRun{opt,indSess})];
        %%
        
        recData.meanSpeedOverTimeRunBefRun = [recData.meanSpeedOverTimeRunBefRun ...
            sessDataSpeed.meanSpeedOverTimeRunBefRun{opt,indSess}'];
        recData.meanSpeedOverTimeRun0to1 = [recData.meanSpeedOverTimeRun0to1 ...
            sessDataSpeed.meanSpeedOverTimeRun0to1{opt,indSess}'];
        recData.meanSpeedOverTimeRun3to5 = [recData.meanSpeedOverTimeRun3to5 ...
            sessDataSpeed.meanSpeedOverTimeRun3to5{opt,indSess}'];
        
        recData.meanSpeedOverDistRun0to100 = [recData.meanSpeedOverDistRun0to100 ...
            sessDataSpeed.meanSpeedOverDistRun0to100{opt,indSess}'];
        recData.meanSpeedOverDistRunAfter100 = [recData.meanSpeedOverDistRunAfter100 ...
            sessDataSpeed.meanSpeedOverDistRunAfter100{opt,indSess}'];
        
        recData.meanRun30to100 = [recData.meanRun30to100 ...
            sessDataLick.meanRun30to100{opt,indSess}'];
        recData.meanRun100to150 = [recData.meanRun100to150 ...
            sessDataLick.meanRun100to150{opt,indSess}'];
        recData.meanRun150to180 = [recData.meanRun150to180 ...
            sessDataLick.meanRun150to180{opt,indSess}'];
        recData.meanRun180to210 = [recData.meanRun180to210 ...
            sessDataLick.meanRun180to210{opt,indSess}'];
    end
end