function [recDataRewPre,recDataRewManip,recDataRewPost,...
    recDataRewPreOpt,recDataRewManipOpt,recDataRewPostOpt,...
        recDataRewPreOptCtrl,recDataRewManipOptCtrl,recDataRewPostOptCtrl] = ...
                accumRecDataRew_opto(path,recSess,manipSess,mazeSess,pulseW,stimCode)
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
    
    %% added by Yingxue on 4/7/2021
    tmp.speedSimMean = [];
    tmp.lickSimMean = [];
    tmp.speedEucMean = [];
    tmp.lickEucMean = [];
    %%
    
    tmp.meanSpeedOverTimeRew0to1 = [];
    tmp.meanSpeedOverTimeRew1to2 = [];
    tmp.meanSpeedOverTimeRew2to3 = [];
    tmp.meanSpeedOverTimeRew3to5 = [];
    
    recDataRewPre = tmp;
    recDataRewManip = tmp;
    recDataRewPost = tmp;
    
    recDataRewPreOpt = tmp;
    recDataRewManipOpt = tmp;
    recDataRewPostOpt = tmp;
    
    recDataRewPreOptCtrl = tmp;
    recDataRewManipOptCtrl = tmp;
    recDataRewPostOptCtrl = tmp;
    
    for i = 1:nRec
        ind = findstr(path(i,:),'\');
        recName = path(i,ind(end-1)+1:ind(end)-1);
        disp(recName);
        fullPath = [path(i,:) recName '_compSess.mat'];
        if(~exist(fullPath,'file'))
            disp([fullPath ' does not exist.']);
            return;
        end
        load(fullPath,'sessDataRew','sessDataSpeed','sessDataRewOpt');
        
        if(stimCode == 5 || stimCode == 6) % stimulation across the running period or only around the reward
            indManip = find(manipSess{i} == stimCode);
        else
            indManip = find(manipSess{i} == stimCode); % & pulseW{i} == 3);
        end
        for j = 1:length(indManip)
            indPre = find(manipSess{i}(1:indManip(j)) == 0,1,'last');
            if(isempty(indPre))
                indPre = NaN;
            end
            indPost = find(manipSess{i}(indManip(j)+1:end) == 0,1,'first')+indManip(j);
            if(isempty(indPost))
                indPost = NaN;
            end
            recDataRewPre = concatData(recDataRewPre,sessDataRew,sessDataSpeed,indPre,i);
            recDataRewPost = concatData(recDataRewPost,sessDataRew,sessDataSpeed,indPost,i);
            recDataRewManip = concatData(recDataRewManip,sessDataRew,sessDataSpeed,indManip(j),i);
            
            recDataRewPreOpt = concatDataOpt(recDataRewPreOpt,sessDataRewOpt,...
                        indPre,i,1);
            recDataRewPostOpt = concatDataOpt(recDataRewPostOpt,sessDataRewOpt,...
                        indPost,i,1);
            recDataRewManipOpt = concatDataOpt(recDataRewManipOpt,sessDataRewOpt,...
                        indManip(j),i,1);
                    
            recDataRewPreOptCtrl = concatDataOpt(recDataRewPreOptCtrl,sessDataRewOpt,...
                        indPre,i,2);
            recDataRewPostOptCtrl = concatDataOpt(recDataRewPostOptCtrl,sessDataRewOpt,...
                        indPost,i,2);
            recDataRewManipOptCtrl = concatDataOpt(recDataRewManipOptCtrl,sessDataRewOpt,...
                        indManip(j),i,2);
        end
    end

end

function recData = concatData(recData,sessData,sessDataSpeed,indSess,indRec)
    nSess = length(indSess);
    if(sum(isnan(indSess)) > 0)
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
        
        %% added by Yingxue on 4/7/2021
        recData.speedSimMean = [recData.speedSimMean NaN*ones(1,nSess)];
        recData.lickSimMean = [recData.lickSimMean NaN*ones(1,nSess)];
        
        recData.speedEucMean = [recData.speedEucMean NaN*ones(1,nSess)];
        recData.lickEucMean = [recData.lickEucMean NaN*ones(1,nSess)];
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
        
        %% added by Yingxue on 4/7/2021
        recData.speedSimMean = [recData.speedSimMean sessData.speedSimMean(indSess)];
        recData.lickSimMean = [recData.lickSimMean sessData.lickSimMean(indSess)];
        
        recData.speedEucMean = [recData.speedEucMean sessData.speedEucMean(indSess)];
        recData.lickEucMean = [recData.lickEucMean sessData.lickEucMean(indSess)];
        %%
        
        recData.meanSpeedOverTimeRew0to1 = [recData.meanSpeedOverTimeRew0to1 ...
            sessDataSpeed.meanSpeedOverTimeRew0to1{indSess}'];
        recData.meanSpeedOverTimeRew1to2 = [recData.meanSpeedOverTimeRew1to2 ...
            sessDataSpeed.meanSpeedOverTimeRew1to2{indSess}'];
        recData.meanSpeedOverTimeRew2to3 = [recData.meanSpeedOverTimeRew2to3 ...
            sessDataSpeed.meanSpeedOverTimeRew2to3{indSess}'];
        recData.meanSpeedOverTimeRew3to5 = [recData.meanSpeedOverTimeRew3to5 ...
            sessDataSpeed.meanSpeedOverTimeRew3to5{indSess}'];
    end
end

function recData = concatDataOpt(recData,sessData,indSess,indRec,opt)
    nSess = length(indSess);
    if(sum(isnan(indSess)) > 0)
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
        
        %% added by Yingxue on 4/7/2021
        recData.speedSimMean = [recData.speedSimMean NaN*ones(1,nSess)];
        recData.lickSimMean = [recData.lickSimMean NaN*ones(1,nSess)];
        
        recData.speedEucMean = [recData.speedEucMean NaN*ones(1,nSess)];
        recData.lickEucMean = [recData.lickEucMean NaN*ones(1,nSess)];
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
        
        %% added by Yingxue on 4/7/2021
        recData.speedSimMean = [recData.speedSimMean sessData.speedSimMean(opt,indSess)];
        recData.lickSimMean = [recData.lickSimMean sessData.lickSimMean(opt,indSess)];
        
        recData.speedEucMean = [recData.speedEucMean sessData.speedEucMean(opt,indSess)];
        recData.lickEucMean = [recData.lickEucMean sessData.lickEucMean(opt,indSess)];
        %%
    end
end
