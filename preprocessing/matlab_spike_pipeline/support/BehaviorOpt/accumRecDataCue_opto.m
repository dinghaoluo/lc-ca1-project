function [recDataCuePre,recDataCueManip,recDataCuePost,...
    recDataCuePreOpt,recDataCueManipOpt,recDataCuePostOpt,...
        recDataCuePreOptCtrl,recDataCueManipOptCtrl,recDataCuePostOptCtrl] = ...
                accumRecDataCue_opto(path,recSess,manipSess,mazeSess,pulseW,stimCode)
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
    tmp.med1stFiveLickDistMean = [];
    tmp.medLickDistMean = [];
    
    %% added by Yingxue on 4/7/2021
    tmp.speedSimMean = [];
    tmp.lickSimMean = [];
    tmp.speedEucMean = [];
    tmp.lickEucMean = [];
    tmp.speedStdMean = [];
    tmp.lickStdMean = [];
    %%
    
    recDataCuePre = tmp;
    recDataCueManip = tmp;
    recDataCuePost = tmp;
    
    recDataCuePreOpt = tmp;
    recDataCueManipOpt = tmp;
    recDataCuePostOpt = tmp;
    
    recDataCuePreOptCtrl = tmp;
    recDataCueManipOptCtrl = tmp;
    recDataCuePostOptCtrl = tmp;
    
    for i = 1:nRec
        ind = findstr(path(i,:),'\');
        recName = path(i,ind(end-1)+1:ind(end)-1);
        disp(recName);
        fullPath = [path(i,:) recName '_compSess.mat'];
        if(~exist(fullPath,'file'))
            disp([fullPath ' does not exist.']);
            return;
        end
        load(fullPath,'sessDataCue','sessDataSpeed','sessDataLick',...
            'sessDataCueOpt','sessDataSpeedOpt','sessDataLickOpt');
        
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
            recDataCuePre = concatData(recDataCuePre,sessDataCue,sessDataSpeed,sessDataLick,indPre,i);
            recDataCuePost = concatData(recDataCuePost,sessDataCue,sessDataSpeed,sessDataLick,indPost,i);
            recDataCueManip = concatData(recDataCueManip,sessDataCue,sessDataSpeed,sessDataLick,indManip(j),i);
            
            recDataCuePreOpt = concatDataOpt(recDataCuePreOpt,sessDataCueOpt,sessDataSpeedOpt,sessDataLickOpt,...
                        indPre,i,1);
            recDataCuePostOpt = concatDataOpt(recDataCuePostOpt,sessDataCueOpt,sessDataSpeedOpt,sessDataLickOpt,...
                        indPost,i,1);
            recDataCueManipOpt = concatDataOpt(recDataCueManipOpt,sessDataCueOpt,sessDataSpeedOpt,sessDataLickOpt,...
                        indManip(j),i,1);
                    
            recDataCuePreOptCtrl = concatDataOpt(recDataCuePreOptCtrl,sessDataCueOpt,sessDataSpeedOpt,sessDataLickOpt,...
                        indPre,i,2);
            recDataCuePostOptCtrl = concatDataOpt(recDataCuePostOptCtrl,sessDataCueOpt,sessDataSpeedOpt,sessDataLickOpt,...
                        indPost,i,2);
            recDataCueManipOptCtrl = concatDataOpt(recDataCueManipOptCtrl,sessDataCueOpt,sessDataSpeedOpt,sessDataLickOpt,...
                        indManip(j),i,2);
        end
    end

end

function recData = concatData(recData,sessData,sessDataSpeed,sessDataLick,indSess,indRec)
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
        recData.med1stFiveLickDistMean = [recData.med1stFiveLickDistMean NaN*ones(1,nSess)];
        recData.medLickDistMean = [recData.medLickDistMean NaN*ones(1,nSess)];
        
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
        recData.med1stFiveLickDistMean = [recData.med1stFiveLickDistMean sessData.med1stFiveLickDistMean(indSess)];
        recData.medLickDistMean = [recData.medLickDistMean sessData.medLickDistMean(indSess)];
        
        %% added by Yingxue on 4/7/2021
        recData.speedSimMean = [recData.speedSimMean sessData.speedSimMean(indSess)];
        recData.lickSimMean = [recData.lickSimMean sessData.lickSimMean(indSess)];
        
        recData.speedEucMean = [recData.speedEucMean sessData.speedEucMean(indSess)];
        recData.lickEucMean = [recData.lickEucMean sessData.lickEucMean(indSess)];
        
        recData.speedStdMean = [recData.speedStdMean mean(sessDataSpeed.stdCue{indSess})];
        recData.lickStdMean = [recData.lickStdMean mean(sessDataLick.stdCue{indSess})];
        %%
    end
end

function recData = concatDataOpt(recData,sessData,sessDataSpeed,sessDataLick,indSess,indRec,opt)
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
        recData.med1stFiveLickDistMean = [recData.med1stFiveLickDistMean NaN*ones(1,nSess)];
        recData.medLickDistMean = [recData.medLickDistMean NaN*ones(1,nSess)];
        
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
        recData.med1stFiveLickDistMean = [recData.med1stFiveLickDistMean sessData.med1stFiveLickDistMean(opt,indSess)];
        recData.medLickDistMean = [recData.medLickDistMean sessData.medLickDistMean(opt,indSess)];
        
        %% added by Yingxue on 4/7/2021
        recData.speedSimMean = [recData.speedSimMean sessData.speedSimMean(opt,indSess)];
        recData.lickSimMean = [recData.lickSimMean sessData.lickSimMean(opt,indSess)];
        
        recData.speedEucMean = [recData.speedEucMean sessData.speedEucMean(opt,indSess)];
        recData.lickEucMean = [recData.lickEucMean sessData.lickEucMean(opt,indSess)];
        
        %%
    end
end
