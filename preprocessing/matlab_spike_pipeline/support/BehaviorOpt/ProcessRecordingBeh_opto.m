function ProcessRecordingBeh_opto(path,recName,sessNo,mazeSessArr,mazeTypeArr)
% process individual files in a recording
% e.g.: ProcessRecordingBeh_opto('Z:\Yingxue\DataAnalysisXiaoliang\A517-20191114\','A517-20191114',[1 3 4 5],[])
% e.g.:
% ProcessRecordingBeh_opto(activeLickOptoSSTPath(1,:),activeLickOptoSSTPath(1,end-13:end-1),ALRecSessionsOptoSST{1},ALOptoMazeSessionSST{1},ALOptoMazeTypeSST{1})

    onlyRun = 0;
    
    if(isempty(mazeSessArr))
        mazeSessArr = ones(1,length(sessNo));
        mazeTypeArr = zeros(1,length(sessNo));
    end
    sessNoTmp = unique(sessNo);
    for i = 1:length(sessNoTmp)
        fileName = [recName '-0' num2str(sessNoTmp(i))];
        pathTmp = [path fileName '\'];
        disp(['Session ' num2str(sessNoTmp(i))]);
        indMazeSess = mazeSessArr(sessNo == sessNoTmp(i));
        ProcessingAlignedBeh_opto(pathTmp,fileName,onlyRun,indMazeSess);
    end

    % compare between different sessions
    sessDataRun = compBehRun(path,recName,sessNo,mazeSessArr,mazeTypeArr,onlyRun);
    close all;    
    sessDataCue = compBehCue(path,recName,sessNo,mazeSessArr,onlyRun);
    close all;
    sessDataRew = compBehRew(path,recName,sessNo,mazeSessArr,onlyRun);
    
    sessDataLick = compBehLick(path,recName,sessNo,mazeSessArr);
    close all;
    sessDataSpeed = compBehSpeed(path,recName,sessNo,mazeSessArr,onlyRun);
    close all;   
    
    %added by Yingxue on 10/01/2020
    sessDataRunOpt = compBehRunOpt(path,recName,sessNo,mazeSessArr,mazeTypeArr,onlyRun);
    close all;
    sessDataCueOpt = compBehCueOpt(path,recName,sessNo,mazeSessArr,onlyRun);
    close all;
    sessDataRewOpt = compBehRewOpt(path,recName,sessNo,mazeSessArr,onlyRun);
    
    sessDataLickOpt = compBehLickOpt(path,recName,sessNo,mazeSessArr);
    close all;
    sessDataSpeedOpt = compBehSpeedOpt(path,recName,sessNo,mazeSessArr,onlyRun);
    close all;  

    fullPath = [path '\' recName '_compSess.mat'];
    save(fullPath, 'sessDataRun','sessDataLick','sessDataSpeed');
    save(fullPath,'sessDataCue','sessDataRew','-append');
    
    % added by Yingxue on 10/01/2020
    save(fullPath, 'sessDataLickOpt','sessDataSpeedOpt','-append');
    save(fullPath,'sessDataRunOpt','sessDataCueOpt','sessDataRewOpt','-append');
    
end

function sessData = compBehRun(path,recName,sessNo,mazeSess,mazeType,onlyRun)
    nSess = length(sessNo);
    sessData.sessNo = sessNo;
    sessData.numSamples = [];
    sessData.totStopLenT = [];
    sessData.numSamplesMean = zeros(1,nSess);
    sessData.numSamplesSEM = zeros(1,nSess);
    sessData.maxSpeedMean = zeros(1,nSess);
    sessData.maxSpeedSEM = zeros(1,nSess);
    sessData.meanSpeedMean = zeros(1,nSess);
    sessData.meanSpeedSEM = zeros(1,nSess);
    sessData.maxSpeedRunMean = zeros(1,nSess);
    sessData.maxSpeedRunSEM = zeros(1,nSess);
    sessData.meanSpeedRunMean = zeros(1,nSess);
    sessData.meanSpeedRunSEM = zeros(1,nSess);
    sessData.maxRunLenTMean = zeros(1,nSess);
    sessData.maxRunLenTSEM = zeros(1,nSess);
    sessData.totRunLenTMean = zeros(1,nSess);
    sessData.totRunLenTSEM = zeros(1,nSess);
    sessData.numRunMean = zeros(1,nSess);
    sessData.numRunSEM = zeros(1,nSess);
    sessData.maxAccMean = zeros(1,nSess);
    sessData.maxAccSEM = zeros(1,nSess);    
    sessData.meanAccMean = zeros(1,nSess);
    sessData.meanAccSEM = zeros(1,nSess);
    sessData.totStopLenTMean = zeros(1,nSess);
    sessData.totStopLenTSEM = zeros(1,nSess);
    sessData.percStopTr =  zeros(1,nSess);
    sessData.startCueToRunMean = zeros(1,nSess);
    sessData.startCueToRunSEM = zeros(1,nSess);
    sessData.med1stFiveLickDistMean = zeros(1,nSess);
    sessData.med1stFiveLickDistSEM = zeros(1,nSess);
    sessData.medLickDistMean = zeros(1,nSess);
    sessData.medLickDistSEM = zeros(1,nSess);
    sessData.med1stFiveLickDistBefRewMean = zeros(1,nSess);
    sessData.med1stFiveLickDistBefRewSEM = zeros(1,nSess);
    sessData.medLickDistBefRewMean = zeros(1,nSess);
    sessData.medLickDistBefRewSEM = zeros(1,nSess);
    sessData.percNonStop = zeros(1,nSess);
    sessData.meanSpeedProfile = [];
    sessData.rewToRunMean = zeros(1,nSess);
    sessData.rewToRunSEM = zeros(1,nSess);
    
    sessData.speedSimMean = zeros(1,nSess);
    sessData.speedSimSem = zeros(1,nSess);
    sessData.lickSimMean = zeros(1,nSess);
    sessData.lickSimSem = zeros(1,nSess);
    
    sessData.speedEucMean = zeros(1,nSess);
    sessData.speedEucSem = zeros(1,nSess);
    sessData.lickEucMean = zeros(1,nSess);
    sessData.lickEucSem = zeros(1,nSess);
        
    for i = 1:nSess
        fileName = [recName '-0' num2str(sessNo(i))];
        pathTmp = [path fileName '\'];

        fullpath = [pathTmp fileName '_behPar_msess' num2str(mazeSess(i)) '.mat'];
        if(~exist(fullpath,'file'))
            disp([fullpath ' does not exist.']);
            return;
        end
        load(fullpath,'behParRun','param');
        
        fullPath = [pathTmp fileName '_alignRun_msess' num2str(mazeSess(i)) '.mat']; 
        if(exist(fullPath) == 0)
            disp('The aligned to run file does not exist');
            return;
        end
        load(fullPath,'trialsRun');
        
        fullPath = [pathTmp fileName '_lickDistSim_msess' num2str(mazeSess(i)) '.mat'];
        if(exist(fullPath) == 0)
            disp('The lick over distance similarity file does not exist');
            return;
        end
        load(fullPath,'lickOverDistSim');
        
        fullPath = [pathTmp fileName '_speedDistSim_msess' num2str(mazeSess(i)) '_Run' num2str(onlyRun) '.mat'];
        if(exist(fullPath) == 0)
            disp('The speed over distance similarity file does not exist');
            return;
        end
        load(fullPath,'speedOverDistSim');
        
        trInd = lickOverDistSim.trInd; % param.startTr:param.endTr;
        nTr = length(trInd);
            
        sessData.mazeType(i) = mazeType(i);
        % added by Yingxue on 7/6/2022
        sessData.numSamples{i} = behParRun.numSamples(trInd);
        sessData.totStopLenT{i} = behParRun.totStopLenT(trInd);
        %
        sessData.numSamplesMean(i) = mean(behParRun.numSamples(trInd));
        sessData.numSamplesSEM(i) = std(behParRun.numSamples(trInd))/sqrt(nTr);
        sessData.maxSpeedMean(i) = mean(behParRun.maxSpeed(trInd));
        sessData.maxSpeedSEM(i) = std(behParRun.maxSpeed(trInd))/sqrt(nTr);
        sessData.meanSpeedMean(i) = mean(behParRun.meanSpeed(trInd));
        sessData.meanSpeedSEM(i) = std(behParRun.meanSpeed(trInd))/sqrt(nTr);
        sessData.maxSpeedRunMean(i) = mean(behParRun.maxSpeedRun(trInd));
        sessData.maxSpeedRunSEM(i) = std(behParRun.maxSpeedRun(trInd))/sqrt(nTr);
        sessData.meanSpeedRunMean(i) = mean(behParRun.meanSpeedRun(trInd));
        sessData.meanSpeedRunSEM(i) = std(behParRun.meanSpeedRun(trInd))/sqrt(nTr);
        sessData.maxRunLenTMean(i) = mean(behParRun.maxRunLenT(trInd));
        sessData.maxRunLenTSEM(i) = std(behParRun.maxRunLenT(trInd))/sqrt(nTr);
        sessData.totRunLenTMean(i) = mean(behParRun.totRunLenT(trInd));
        sessData.totRunLenTSEM(i) = std(behParRun.totRunLenT(trInd))/sqrt(nTr);  
        sessData.numRunMean(i) = mean(behParRun.numRun(trInd));
        sessData.numRunSEM(i) = std(behParRun.numRun(trInd))/sqrt(nTr);
        sessData.maxAccMean(i) = mean(behParRun.maxAcc(trInd));
        sessData.maxAccSEM(i) = std(behParRun.maxAcc(trInd))/sqrt(nTr);
        sessData.meanAccMean(i) = mean(behParRun.meanAcc(trInd));
        sessData.meanAccSEM(i) = std(behParRun.meanAcc(trInd))/sqrt(nTr);
        sessData.totStopLenTMean(i) = mean(behParRun.totStopLenT(trInd));
        sessData.totStopLenTSEM(i) = std(behParRun.totStopLenT(trInd))/sqrt(nTr);
        sessData.percStopTr(i) = sum(behParRun.totStopLenT(trInd) > 0)/nTr; % stop longer than 0 ms, added on 8/13/2022
        sessData.startCueToRunMean(i) = mean(behParRun.startCueToRun(trInd));
        sessData.startCueToRunSEM(i) = std(behParRun.startCueToRun(trInd))/sqrt(nTr);
        sessData.numLicksBefRewMean(i) = mean(behParRun.numLicksBefRew(trInd));
        sessData.numLicksBefRewSEM(i) = std(behParRun.numLicksBefRew(trInd))/sqrt(nTr);
        sessData.numLicksRewMean(i) = mean(behParRun.numLicksRew(trInd));
        sessData.numLicksRewSEM(i) = std(behParRun.numLicksRew(trInd))/sqrt(nTr);
        sessData.rewToRunMean(i) = mean(behParRun.rewardToRun(trInd));
        sessData.rewToRunSEM(i) = mean(behParRun.rewardToRun(trInd))/sqrt(nTr);
        
        %% added by Yingxue on 4/19/2021
        sessData.lickSimMean(i) = lickOverDistSim.meanRun;
        sessData.lickSimSem(i) = lickOverDistSim.semRun;
        
        sessData.speedSimMean(i) = speedOverDistSim.meanRun;
        sessData.speedSimSem(i) = speedOverDistSim.semRun;
        
        sessData.lickEucMean(i) = lickOverDistSim.meanEucRun;
        sessData.lickEucSem(i) = lickOverDistSim.semEucRun;
        
        sessData.speedEucMean(i) = speedOverDistSim.meanEucRun;
        sessData.speedEucSem(i) = speedOverDistSim.semEucRun;
        %%
        
        indN = ~isnan(behParRun.med1stFiveLickDist(trInd));
        med1stFiveLickDist = behParRun.med1stFiveLickDist(trInd);
        sessData.med1stFiveLickDistMean(i) = mean(med1stFiveLickDist(indN));
        sessData.med1stFiveLickDistSEM(i) = std(med1stFiveLickDist(indN))/sqrt(sum(indN));
        
        indN = ~isnan(behParRun.medLickDist(trInd));
        medLickDist = behParRun.medLickDist(trInd);
        sessData.medLickDistMean(i) = mean(medLickDist(indN));
        sessData.medLickDistSEM(i) = std(medLickDist(indN))/sqrt(sum(indN));  
        
        indN = ~isnan(behParRun.med1stFiveLickDistBefRew(trInd));
        med1stFiveLickDistBefRew = behParRun.med1stFiveLickDistBefRew(trInd);
        sessData.med1stFiveLickDistBefRewMean(i) = mean(med1stFiveLickDistBefRew(indN));
        sessData.med1stFiveLickDistBefRewSEM(i) = std(med1stFiveLickDistBefRew(indN))/sqrt(sum(indN));
        
        indN = ~isnan(behParRun.medLickDistBefRew(trInd));
        medLickDistBefRew = behParRun.medLickDistBefRew(trInd);
        sessData.medLickDistBefRewMean(i) = mean(medLickDistBefRew(indN));
        sessData.medLickDistBefRewSEM(i) = std(medLickDistBefRew(indN))/sqrt(sum(indN));
        
        sessData.percRewarded(i) = sum(behParRun.rewarded(trInd)==1)/nTr;
        sessData.percNonStop(i) = sum(behParRun.nonStop(trInd))/nTr;
        sessData.speedProfile{i} = behParRun.speedProfile(trInd,:)/10;
        sessData.lickProfile{i} = behParRun.lickProfile(trInd,:);
        
        indN = behParRun.pumpLfpInd(trInd) > 0;
        sessData.pumpLfpIndMean(i) = mean(behParRun.pumpLfpInd(trInd(indN)));
        sessData.pumpLfpIndSEM(i) = std(behParRun.pumpLfpInd(trInd(indN)))/sqrt(sum(indN));
        sessData.pumpMMMean(i) = mean(behParRun.pumpMM(trInd(indN)));
        sessData.pumpMMSEM(i) = std(behParRun.pumpMM(trInd(indN)))/sqrt(sum(indN));
        
    end
    
%     barPlot(1:nSess,sessData.numSamplesMean,sessData.numSamplesSEM,'Session','No. samples / Trial','Run onset');
%     print ('-painters', '-dpdf', [path recName '_numSamplesRun'], '-r600')
    barPlot(1:nSess,sessData.meanSpeedMean,sessData.meanSpeedSEM,'Session','Mean speed (mm/s)','Run onset');
    print ('-painters', '-dpdf', [path recName '_meanSpeedRun'], '-r600')
    barPlot(1:nSess,sessData.meanSpeedRunMean,sessData.meanSpeedRunSEM,'Session','Mean speed run only (mm/s)','Run onset');
    print ('-painters', '-dpdf', [path recName '_meanSpeedRunRun'], '-r600')
%     barPlot(1:nSess,sessData.maxSpeedMean,sessData.maxSpeedSEM,'Session','Max speed (mm/s)','Run onset');
%     print ('-painters', '-dpdf', [path recName '_maxSpeedRun'], '-r600')
%     barPlot(1:nSess,sessData.maxRunLenTMean,sessData.maxRunLenTSEM,'Session','Max run length (s)','Run onset');
%     print ('-painters', '-dpdf', [path recName '_maxRunLenTRun'], '-r600')
    barPlot(1:nSess,sessData.totRunLenTMean,sessData.totRunLenTSEM,'Session','Total run length (s)','Run onset');
    print ('-painters', '-dpdf', [path recName '_totRunLenTRun'], '-r600')
    barPlot(1:nSess,sessData.numRunMean,sessData.numRunSEM,'Session','No. run segments','Run onset');
    print ('-painters', '-dpdf', [path recName '_numRunRun'], '-r600')
%     barPlot(1:nSess,sessData.meanAccMean,sessData.meanAccSEM,'Session','Mean acceleration (mm/s^2)','Run onset');
%     print ('-painters', '-dpdf', [path recName '_meanAccRun'], '-r600')
    barPlot(1:nSess,sessData.totStopLenTMean,sessData.totStopLenTSEM,'Session','Total stop time (s)','Run onset');
    print ('-painters', '-dpdf', [path recName '_totStopLenTRun'], '-r600')
%     barPlot(1:nSess,sessData.startCueToRunMean,sessData.startCueToRunSEM,'Session','Start cue to run time (s)','Run onset');
%     print ('-painters', '-dpdf', [path recName '_startCueToRun'], '-r600')
    
    barPlot(1:nSess,sessData.numLicksBefRewMean,sessData.numLicksBefRewSEM,'Session','Num licks bef. rew.','Run onset');
    print ('-painters', '-dpdf', [path recName '_numLicksBefRew'], '-r600')
    barPlot(1:nSess,sessData.numLicksRewMean,sessData.numLicksRewSEM,'Session','Num licks aft. rew.','Run onset');
    print ('-painters', '-dpdf', [path recName '_numLicksRew'], '-r600')
    
%     barPlot(1:nSess,sessData.med1stFiveLickDistMean,sessData.med1stFiveLickDistSEM,'Session','Med. dist. first five licks (mm)','Run onset');
%     print ('-painters', '-dpdf', [path recName '_med1stFiveLickDistRun'], '-r600')
%     barPlot(1:nSess,sessData.medLickDistMean,sessData.medLickDistSEM,'Session','Med. dist. licks (mm)','Run onset');
%     print ('-painters', '-dpdf', [path recName '_medLickDistRun'], '-r600')
    
    barPlot(1:nSess,sessData.med1stFiveLickDistBefRewMean,sessData.med1stFiveLickDistBefRewSEM,'Session','Med. dist. first five licks before reward (mm)','Run onset');
    print ('-painters', '-dpdf', [path recName '_med1stFiveLickDistBefRew'], '-r600')
%     barPlot(1:nSess,sessData.medLickDistBefRewMean,sessData.medLickDistBefRewSEM,'Session','Med. dist. licks before reward (mm)','Run onset');
%     print ('-painters', '-dpdf', [path recName '_medLickDistBefRewRun'], '-r600')

    barPlot(1:nSess,sessData.percRewarded,zeros(1,nSess),'Session','Performance (ratio correct)','Run onset');
    print ('-painters', '-dpdf', [path recName '_percRewardedRun'], '-r600')
%     barPlot(1:nSess,sessData.percNonStop,zeros(1,nSess),'Session','Perc. non-stop','Run onset');
%     print ('-painters', '-dpdf', [path recName '_percNonstopRun'], '-r600')
%     runPlot(sessData.speedProfile,(1:size(sessData.speedProfile{1},2))/1250,'Time (s)','Speed (cm/s)','Run onset');

    barPlot(1:nSess,sessData.pumpMMMean,sessData.pumpMMSEM,'Session','pump location (mm)','Run onset');
    print ('-painters', '-dpdf', [path recName '_pumpMM'], '-r600')
    
    barPlot(1:nSess,sessData.pumpLfpIndMean/1250,sessData.pumpLfpIndSEM/1250,'Session','pump time (s)','Run onset');
    print ('-painters', '-dpdf', [path recName '_pumpTime'], '-r600')
        
    barPlot(1:nSess,sessData.speedSimMean,sessData.speedSimSem,'Session','speed over distance similarity','Run onset');
    print ('-painters', '-dpdf', [path recName '_speedSimRun'], '-r600')
    
    barPlot(1:nSess,sessData.lickSimMean,sessData.lickSimSem,'Session','lick over distance similarity','Run onset');
    print ('-painters', '-dpdf', [path recName '_lickSimRun'], '-r600')
    
    barPlot(1:nSess,sessData.speedEucMean,sessData.speedEucSem,'Session','speed over distance euclidean','Run onset');
    print ('-painters', '-dpdf', [path recName '_speedEucRun'], '-r600')
    
    barPlot(1:nSess,sessData.lickEucMean,sessData.lickEucSem,'Session','lick over distance euclidean','Run onset');
    print ('-painters', '-dpdf', [path recName '_lickEucRun'], '-r600')
    
%     lickPlot(sessData.lickProfile,(1:size(sessData.speedProfile{1},2))/1250,'Time (s)','No. licks','Run onset');
end

%% added by Yingxue on 10/01/2020
function sessData = compBehRunOpt(path,recName,sessNo,mazeSess,mazeType,onlyRun)
    nSess = length(sessNo);
    sessData.sessNo = sessNo;
    sessData.mazeType = mazeType;
    sessData.numSamples = cell(1,nSess);
    sessData.totStopLenT = cell(1,nSess);
    sessData.numSamplesMean = zeros(3,nSess);
    sessData.numSamplesSEM = zeros(3,nSess);
    sessData.maxSpeedMean = zeros(3,nSess);
    sessData.maxSpeedSEM = zeros(3,nSess);
    sessData.meanSpeedMean = zeros(3,nSess);
    sessData.meanSpeedSEM = zeros(3,nSess);
    sessData.maxSpeedRunMean = zeros(3,nSess);
    sessData.maxSpeedRunSEM = zeros(3,nSess);
    sessData.meanSpeedRunMean = zeros(3,nSess);
    sessData.meanSpeedRunSEM = zeros(3,nSess);
    sessData.maxRunLenTMean = zeros(3,nSess);
    sessData.maxRunLenTSEM = zeros(3,nSess);
    sessData.totRunLenTMean = zeros(3,nSess);
    sessData.totRunLenTSEM = zeros(3,nSess);
    sessData.numRunMean = zeros(3,nSess);
    sessData.numRunSEM = zeros(3,nSess);
    sessData.maxAccMean = zeros(3,nSess);
    sessData.maxAccSEM = zeros(3,nSess);    
    sessData.meanAccMean = zeros(3,nSess);
    sessData.meanAccSEM = zeros(3,nSess);
    sessData.totStopLenTMean = zeros(3,nSess);
    sessData.totStopLenTSEM = zeros(3,nSess);
    sessData.percStopTr =  zeros(3,nSess);
    sessData.startCueToRunMean = zeros(3,nSess);
    sessData.startCueToRunSEM = zeros(3,nSess);
    sessData.numLicksBefRewMean = zeros(3,nSess);
    sessData.numLicksBefRewSEM = zeros(3,nSess);
    sessData.numLicksRewMean = zeros(3,nSess);
    sessData.numLicksRewSEM = zeros(3,nSess);
    sessData.med1stFiveLickDistMean = zeros(3,nSess);
    sessData.med1stFiveLickDistSEM = zeros(3,nSess);
    sessData.medLickDistMean = zeros(3,nSess);
    sessData.medLickDistSEM = zeros(3,nSess);
    sessData.med1stFiveLickDistBefRewMean = zeros(3,nSess);
    sessData.med1stFiveLickDistBefRewSEM = zeros(3,nSess);
    sessData.medLickDistBefRewMean = zeros(3,nSess);
    sessData.medLickDistBefRewSEM = zeros(3,nSess); 
    sessData.rewToRunMean = zeros(3,nSess);
    sessData.rewToRunSEM = zeros(3,nSess);
    
    sessData.speedSimMean = zeros(2,nSess);
    sessData.speedSimSem = zeros(2,nSess);
    sessData.lickSimMean = zeros(2,nSess);
    sessData.lickSimSem = zeros(2,nSess);
    
    sessData.speedEucMean = zeros(2,nSess);
    sessData.speedEucSem = zeros(2,nSess);
    sessData.lickEucMean = zeros(2,nSess);
    sessData.lickEucSem = zeros(2,nSess);
    
    sessData.percRewarded = zeros(3,nSess);
    sessData.percNonStop = zeros(3,nSess);

    sessData.pumpLfpIndMean = zeros(3,nSess);
    sessData.pumpLfpIndSEM = zeros(3,nSess);
    sessData.pumpMMMean = zeros(3,nSess);
    sessData.pumpMMSEM = zeros(3,nSess);
    
    for i = 1:nSess
        fileName = [recName '-0' num2str(sessNo(i))];
        pathTmp = [path fileName '\'];

        fullpath = [pathTmp fileName '_behPar_msess' num2str(mazeSess(i)) '.mat'];
        if(~exist(fullpath,'file'))
            disp([fullpath ' does not exist.']);
            return;
        end
        load(fullpath,'behParRun','param');
        
        fullPath = [pathTmp fileName '_alignRun_msess' num2str(mazeSess(i)) '.mat']; 
        if(exist(fullPath) == 0)
            disp('The aligned to run file does not exist');
            return;
        end
        load(fullPath,'trialsRun');
        
        fullPath = [pathTmp fileName '_lickDist_msess' num2str(mazeSess(i)) '.mat']; 
        if(exist(fullPath) == 0)
            disp('The _lickDist file does not exist');
            return;
        end
        load(fullPath,'lickOverDistOpt','lickOverDistOptCtrl','lickOverDistOptCtrl1');
        
        fullPath = [pathTmp fileName '_lickDistSim_msess' num2str(mazeSess(i)) '.mat'];
        if(exist(fullPath) == 0)
            disp('The lick over distance similarity file does not exist');
            return;
        end
        load(fullPath,'lickOverDistSim','lickOverDistSimOpt','lickOverDistSimOptCtrl');
        
        fullPath = [pathTmp fileName '_speedDistSim_msess' num2str(mazeSess(i)) '_Run' num2str(onlyRun) '.mat'];
        if(exist(fullPath) == 0)
            disp('The speed over distance similarity file does not exist');
            return;
        end
        load(fullPath,'speedOverDistSim','speedOverDistSimOpt','speedOverDistSimOptCtrl');
        
        if(~isempty(lickOverDistOpt.trInd))
            trIndOpt = lickOverDistOpt.trInd;
            trIndCtrl = lickOverDistOptCtrl.trInd;
            trIndCtrl1 = lickOverDistOptCtrl1.trInd; % added by Yingxue on 1/25/2022
            nTrOpt = length(trIndOpt);
            nTrCtrl = length(trIndCtrl);
            nTrCtrl1 = length(trIndCtrl1); % added by Yingxue on 1/25/2022
            
            %% added by Yingxue on 1/25/2022
            mazeT = unique(behParRun.mazeType);
            mazeT = mazeT(mazeT ~= 0);
            if(mazeT == 4)
                trIndOptTmp = trIndOpt-1;
            else
                trIndOptTmp = trIndOpt;
            end
            sessData.stimOnLfpInd{i} =  behParRun.stimOnLfpInd(trIndOptTmp);
            sessData.stimOffLfpInd{i} =  behParRun.stimOffLfpInd(trIndOptTmp);
            sessData.stimAcrossTrialRun{i} =  behParRun.stimAcrossTrialRun(trIndOptTmp);
            sessData.stimAcrossTrialCue{i} =  behParRun.stimAcrossTrialCue(trIndOptTmp);
            sessData.stimAcrossTrialRew{i} =  behParRun.stimAcrossTrialRew(trIndOptTmp);
            
            % added by Yingxue on 7/6/2022
            sessData.numSamples{i} = [{behParRun.numSamples(trIndOpt)} ...
                {behParRun.numSamples(trIndCtrl)} ...
                {behParRun.numSamples(trIndCtrl1)}];
            sessData.totStopLenT{i} = [{behParRun.totStopLenT(trIndOpt)} ...
                {behParRun.totStopLenT(trIndCtrl)} ...
                {behParRun.totStopLenT(trIndCtrl1)}];
            %
            
            %% changed by Yingxue on 1/25/2022, added the features for control trials immediately 
            % after the stim trials
            sessData.numSamplesMean(:,i) = [mean(behParRun.numSamples(trIndOpt));...
                mean(behParRun.numSamples(trIndCtrl));...
                mean(behParRun.numSamples(trIndCtrl1))];
            sessData.numSamplesSEM(:,i) = [std(behParRun.numSamples(trIndOpt))/sqrt(nTrOpt);...
                std(behParRun.numSamples(trIndCtrl))/sqrt(nTrCtrl);...
                std(behParRun.numSamples(trIndCtrl1))/sqrt(nTrCtrl1)];
            sessData.maxSpeedMean(:,i) = [mean(behParRun.maxSpeed(trIndOpt));...
                mean(behParRun.maxSpeed(trIndCtrl));...
                mean(behParRun.maxSpeed(trIndCtrl1))];
            sessData.maxSpeedSEM(:,i) = [std(behParRun.maxSpeed(trIndOpt))/sqrt(nTrOpt);...
                std(behParRun.maxSpeed(trIndCtrl))/sqrt(nTrCtrl);...
                std(behParRun.maxSpeed(trIndCtrl1))/sqrt(nTrCtrl1)];
            sessData.meanSpeedMean(:,i) = [mean(behParRun.meanSpeed(trIndOpt));...
                mean(behParRun.meanSpeed(trIndCtrl));...
                mean(behParRun.meanSpeed(trIndCtrl1))];
            sessData.meanSpeedSEM(:,i) = [std(behParRun.meanSpeed(trIndOpt))/sqrt(nTrOpt);...
                std(behParRun.meanSpeed(trIndCtrl))/sqrt(nTrCtrl);...
                std(behParRun.meanSpeed(trIndCtrl1))/sqrt(nTrCtrl1)];
            sessData.maxSpeedRunMean(:,i) = [mean(behParRun.maxSpeedRun(trIndOpt));...
                mean(behParRun.maxSpeedRun(trIndCtrl));...
                mean(behParRun.maxSpeedRun(trIndCtrl1))];
            sessData.maxSpeedRunSEM(:,i) = [std(behParRun.maxSpeedRun(trIndOpt))/sqrt(nTrOpt);...
                std(behParRun.maxSpeedRun(trIndCtrl))/sqrt(nTrCtrl);...
                std(behParRun.maxSpeedRun(trIndCtrl1))/sqrt(nTrCtrl1)];
            sessData.meanSpeedRunMean(:,i) = [mean(behParRun.meanSpeedRun(trIndOpt));...
                mean(behParRun.meanSpeedRun(trIndCtrl));...
                mean(behParRun.meanSpeedRun(trIndCtrl1))];
            sessData.meanSpeedRunSEM(:,i) = [std(behParRun.meanSpeedRun(trIndOpt))/sqrt(nTrOpt);...
                std(behParRun.meanSpeedRun(trIndCtrl))/sqrt(nTrCtrl);...
                std(behParRun.meanSpeedRun(trIndCtrl1))/sqrt(nTrCtrl1)];
            sessData.maxRunLenTMean(:,i) = [mean(behParRun.maxRunLenT(trIndOpt));...
                mean(behParRun.maxRunLenT(trIndCtrl));...
                mean(behParRun.maxRunLenT(trIndCtrl1))];
            sessData.maxRunLenTSEM(:,i) = [std(behParRun.maxRunLenT(trIndOpt))/sqrt(nTrOpt);...
                std(behParRun.maxRunLenT(trIndCtrl))/sqrt(nTrCtrl);...
                std(behParRun.maxRunLenT(trIndCtrl1))/sqrt(nTrCtrl1)];
            sessData.totRunLenTMean(:,i) = [mean(behParRun.totRunLenT(trIndOpt));...
                mean(behParRun.totRunLenT(trIndCtrl));...
                mean(behParRun.totRunLenT(trIndCtrl1))];
            sessData.totRunLenTSEM(:,i) = [std(behParRun.totRunLenT(trIndOpt))/sqrt(nTrOpt);...
                std(behParRun.totRunLenT(trIndCtrl))/sqrt(nTrCtrl);...
                std(behParRun.totRunLenT(trIndCtrl1))/sqrt(nTrCtrl1)];      
            sessData.numRunMean(:,i) = [mean(behParRun.numRun(trIndOpt));...
                mean(behParRun.numRun(trIndCtrl));...
                mean(behParRun.numRun(trIndCtrl1))];
            sessData.numRunSEM(:,i) = [std(behParRun.numRun(trIndOpt))/sqrt(nTrOpt);...
                std(behParRun.numRun(trIndCtrl))/sqrt(nTrCtrl);...
                std(behParRun.numRun(trIndCtrl1))/sqrt(nTrCtrl1)];
            sessData.maxAccMean(:,i) = [mean(behParRun.maxAcc(trIndOpt));...
                mean(behParRun.maxAcc(trIndCtrl));...
                mean(behParRun.maxAcc(trIndCtrl1))];
            sessData.maxAccSEM(:,i) = [std(behParRun.maxAcc(trIndOpt))/sqrt(nTrOpt);...
                std(behParRun.maxAcc(trIndCtrl))/sqrt(nTrCtrl);...
                std(behParRun.maxAcc(trIndCtrl1))/sqrt(nTrCtrl1)];
            sessData.meanAccMean(:,i) = [mean(behParRun.meanAcc(trIndOpt));...
                mean(behParRun.meanAcc(trIndCtrl));...
                mean(behParRun.meanAcc(trIndCtrl1))];
            sessData.meanAccSEM(:,i) = [std(behParRun.meanAcc(trIndOpt))/sqrt(nTrOpt);...
                std(behParRun.meanAcc(trIndCtrl))/sqrt(nTrCtrl);...
                std(behParRun.meanAcc(trIndCtrl1))/sqrt(nTrCtrl1)];
            sessData.totStopLenTMean(:,i) = [mean(behParRun.totStopLenT(trIndOpt));...
                mean(behParRun.totStopLenT(trIndCtrl));...
                mean(behParRun.totStopLenT(trIndCtrl1))];
            sessData.totStopLenTSEM(:,i) = [std(behParRun.totStopLenT(trIndOpt))/sqrt(nTrOpt);...
                std(behParRun.totStopLenT(trIndCtrl))/sqrt(nTrCtrl);...
                std(behParRun.totStopLenT(trIndCtrl1))/sqrt(nTrCtrl1)];
            sessData.percStopTr(:,i) = [sum(behParRun.totStopLenT(trIndOpt) > 0)/nTrOpt;...
                sum(behParRun.totStopLenT(trIndCtrl) > 0)/nTrCtrl;...
                sum(behParRun.totStopLenT(trIndCtrl1) > 0)/nTrCtrl1]; % stop longer than 0 ms, added on 8/13/2022
            sessData.startCueToRunMean(:,i) = [mean(behParRun.startCueToRun(trIndOpt));...
                mean(behParRun.startCueToRun(trIndCtrl));...
                mean(behParRun.startCueToRun(trIndCtrl1))];
            sessData.startCueToRunSEM(:,i) = [std(behParRun.startCueToRun(trIndOpt))/sqrt(nTrOpt);...
                std(behParRun.startCueToRun(trIndCtrl))/sqrt(nTrCtrl);...
                std(behParRun.startCueToRun(trIndCtrl1))/sqrt(nTrCtrl1)];
            sessData.numLicksBefRewMean(:,i) = [mean(behParRun.numLicksBefRew(trIndOpt));...
                mean(behParRun.numLicksBefRew(trIndCtrl));...
                mean(behParRun.numLicksBefRew(trIndCtrl1))];
            sessData.numLicksBefRewSEM(:,i) = [std(behParRun.numLicksBefRew(trIndOpt))/sqrt(nTrOpt);...
                std(behParRun.numLicksBefRew(trIndCtrl))/sqrt(nTrCtrl);...
                std(behParRun.numLicksBefRew(trIndCtrl1))/sqrt(nTrCtrl1)];
            sessData.numLicksRewMean(:,i) = [mean(behParRun.numLicksRew(trIndOpt));...
                mean(behParRun.numLicksRew(trIndCtrl));...
                mean(behParRun.numLicksRew(trIndCtrl1))];
            sessData.numLicksRewSEM(:,i) = [std(behParRun.numLicksRew(trIndOpt))/sqrt(nTrOpt);...
                std(behParRun.numLicksRew(trIndCtrl))/sqrt(nTrCtrl);...
                std(behParRun.numLicksRew(trIndCtrl1))/sqrt(nTrCtrl1)];
            sessData.rewToRunMean(:,i) = [mean(behParRun.rewardToRun(trIndOpt));...
                mean(behParRun.rewardToRun(trIndCtrl));...
                mean(behParRun.rewardToRun(trIndCtrl1))];
            sessData.rewToRunSEM(:,i) = [mean(behParRun.rewardToRun(trIndOpt))/sqrt(nTrOpt);...
                mean(behParRun.rewardToRun(trIndCtrl))/sqrt(nTrCtrl);...
                mean(behParRun.rewardToRun(trIndCtrl1))/sqrt(nTrCtrl1)];
            
            %% added by Yingxue on 4/19/2021
            sessData.lickSimMean(:,i) = [lickOverDistSimOpt.meanRun;lickOverDistSimOptCtrl.meanRun];
            sessData.lickSimSem(:,i) = [lickOverDistSimOpt.semRun;lickOverDistSimOptCtrl.semRun];

            sessData.speedSimMean(:,i) = [speedOverDistSimOpt.meanRun;speedOverDistSimOptCtrl.meanRun];
            sessData.speedSimSem(:,i) = [speedOverDistSimOpt.semRun;speedOverDistSimOptCtrl.semRun];
            
            sessData.lickEucMean(:,i) = [lickOverDistSimOpt.meanEucRun;lickOverDistSimOptCtrl.meanEucRun];
            sessData.lickEucSem(:,i) = [lickOverDistSimOpt.semEucRun;lickOverDistSimOptCtrl.semEucRun];

            sessData.speedEucMean(:,i) = [speedOverDistSimOpt.meanEucRun;speedOverDistSimOptCtrl.meanEucRun];
            sessData.speedEucSem(:,i) = [speedOverDistSimOpt.semEucRun;speedOverDistSimOptCtrl.semEucRun];
            %%
            
            indN = ~isnan(behParRun.med1stFiveLickDist(trIndOpt));
            med1stFiveLickDist = behParRun.med1stFiveLickDist(trIndOpt);
            indNCtrl = ~isnan(behParRun.med1stFiveLickDist(trIndCtrl));
            med1stFiveLickDistCtrl = behParRun.med1stFiveLickDist(trIndCtrl);
            indNCtrl1 = ~isnan(behParRun.med1stFiveLickDist(trIndCtrl1));
            med1stFiveLickDistCtrl1 = behParRun.med1stFiveLickDist(trIndCtrl1);
            sessData.med1stFiveLickDistMean(:,i) = [mean(med1stFiveLickDist(indN));...
                mean(med1stFiveLickDistCtrl(indNCtrl));...
                mean(med1stFiveLickDistCtrl1(indNCtrl1))];
            sessData.med1stFiveLickDistSEM(:,i) = [std(med1stFiveLickDist(indN))/sqrt(sum(indN));...
                std(med1stFiveLickDistCtrl(indNCtrl))/sqrt(sum(indNCtrl));...
                std(med1stFiveLickDistCtrl1(indNCtrl1))/sqrt(sum(indNCtrl1))];

            indN = ~isnan(behParRun.medLickDist(trIndOpt));
            medLickDist = behParRun.medLickDist(trIndOpt);
            indNCtrl = ~isnan(behParRun.medLickDist(trIndCtrl));
            medLickDistCtrl = behParRun.medLickDist(trIndCtrl);
            indNCtrl1 = ~isnan(behParRun.medLickDist(trIndCtrl1));
            medLickDistCtrl1 = behParRun.medLickDist(trIndCtrl1);
            sessData.medLickDistMean(:,i) = [mean(medLickDist(indN));...
                mean(medLickDistCtrl(indNCtrl));...
                mean(medLickDistCtrl1(indNCtrl1))];
            sessData.medLickDistSEM(:,i) = [std(medLickDist(indN))/sqrt(sum(indN));...
                std(medLickDistCtrl(indNCtrl))/sqrt(sum(indNCtrl));...
                std(medLickDistCtrl1(indNCtrl1))/sqrt(sum(indNCtrl1))];  

            indN = ~isnan(behParRun.med1stFiveLickDistBefRew(trIndOpt));
            med1stFiveLickDistBefRew = behParRun.med1stFiveLickDistBefRew(trIndOpt);
            indNCtrl = ~isnan(behParRun.med1stFiveLickDistBefRew(trIndCtrl));
            med1stFiveLickDistBefRewCtrl = behParRun.med1stFiveLickDistBefRew(trIndCtrl);
            indNCtrl1 = ~isnan(behParRun.med1stFiveLickDistBefRew(trIndCtrl1));
            med1stFiveLickDistBefRewCtrl1 = behParRun.med1stFiveLickDistBefRew(trIndCtrl1);
            sessData.med1stFiveLickDistBefRewMean(:,i) = [mean(med1stFiveLickDistBefRew(indN));...
                mean(med1stFiveLickDistBefRewCtrl(indNCtrl));...
                mean(med1stFiveLickDistBefRewCtrl1(indNCtrl1))];
            sessData.med1stFiveLickDistBefRewSEM(:,i) = [std(med1stFiveLickDistBefRew(indN))/sqrt(sum(indN));...
                std(med1stFiveLickDistBefRewCtrl(indNCtrl))/sqrt(sum(indNCtrl));...
                std(med1stFiveLickDistBefRewCtrl1(indNCtrl1))/sqrt(sum(indNCtrl1))];

            indN = ~isnan(behParRun.medLickDistBefRew(trIndOpt));
            medLickDistBefRew = behParRun.medLickDistBefRew(trIndOpt);
            indNCtrl = ~isnan(behParRun.medLickDistBefRew(trIndCtrl));
            medLickDistBefRewCtrl = behParRun.medLickDistBefRew(trIndCtrl);
            indNCtrl1 = ~isnan(behParRun.medLickDistBefRew(trIndCtrl1));
            medLickDistBefRewCtrl1 = behParRun.medLickDistBefRew(trIndCtrl1);
            sessData.medLickDistBefRewMean(:,i) = [mean(medLickDistBefRew(indN));...
                mean(medLickDistBefRewCtrl(indNCtrl));...
                mean(medLickDistBefRewCtrl1(indNCtrl1))];
            sessData.medLickDistBefRewSEM(:,i) = [std(medLickDistBefRew(indN))/sqrt(sum(indN));...
                std(medLickDistBefRewCtrl(indNCtrl))/sqrt(sum(indNCtrl));...
                std(medLickDistBefRewCtrl1(indNCtrl1))/sqrt(sum(indNCtrl1))];

            sessData.percRewarded(:,i) = [sum(behParRun.rewarded(trIndOpt)==1)/nTrOpt;...
                sum(behParRun.rewarded(trIndCtrl)==1)/nTrCtrl;...
                sum(behParRun.rewarded(trIndCtrl1)==1)/nTrCtrl1];
            sessData.percNonStop(:,i) = [sum(behParRun.nonStop(trIndOpt))/nTrOpt;...
                sum(behParRun.nonStop(trIndCtrl))/nTrCtrl;...
                sum(behParRun.nonStop(trIndCtrl1))/nTrCtrl1];
            sessData.speedProfile(:,i) = [{behParRun.speedProfile(trIndOpt,:)/10};...
                {behParRun.speedProfile(trIndCtrl,:)/10};...
                {behParRun.speedProfile(trIndCtrl1,:)/10}];
            sessData.lickProfile(:,i) = [{behParRun.lickProfile(trIndOpt,:)};...
                {behParRun.lickProfile(trIndCtrl,:)};...
                {behParRun.lickProfile(trIndCtrl1,:)}];

            indN = behParRun.pumpLfpInd(trIndOpt) > 0;
            indNCtrl = behParRun.pumpLfpInd(trIndCtrl) > 0;
            indNCtrl1 = behParRun.pumpLfpInd(trIndCtrl1) > 0;
            sessData.pumpLfpIndMean(:,i) = [mean(behParRun.pumpLfpInd(trIndOpt(indN)));...
                mean(behParRun.pumpLfpInd(trIndCtrl(indNCtrl)));...
                mean(behParRun.pumpLfpInd(trIndCtrl1(indNCtrl1)))];
            sessData.pumpLfpIndSEM(:,i) = [std(behParRun.pumpLfpInd(trIndOpt(indN)))/sqrt(sum(indN));...
                std(behParRun.pumpLfpInd(trIndCtrl(indNCtrl)))/sqrt(sum(indNCtrl));...
                std(behParRun.pumpLfpInd(trIndCtrl1(indNCtrl1)))/sqrt(sum(indNCtrl1))];
            sessData.pumpMMMean(:,i) = [mean(behParRun.pumpMM(trIndOpt(indN)));...
                mean(behParRun.pumpMM(trIndCtrl(indNCtrl)));...
                mean(behParRun.pumpMM(trIndCtrl1(indNCtrl1)))];
            sessData.pumpMMSEM(:,i) = [std(behParRun.pumpMM(trIndOpt(indN)))/sqrt(sum(indN));...
                std(behParRun.pumpMM(trIndCtrl(indNCtrl)))/sqrt(sum(indNCtrl));...
                std(behParRun.pumpMM(trIndCtrl1(indNCtrl1)))/sqrt(sum(indNCtrl1))];
        else
            trInd = param.startTr:param.endTr;
            nTr = length(trInd);
            
            %% added by Yingxue on 2/10/2022
            sessData.stimOnLfpInd{i} =  [];
            sessData.stimOffLfpInd{i} =  [];
            sessData.stimAcrossTrialRun{i} =  [];
            sessData.stimAcrossTrialCue{i} =  [];
            sessData.stimAcrossTrialRew{i} =  [];
            
            % added by Yingxue on 7/6/2022
            sessData.numSamples{i} = behParRun.numSamples(trInd);
            sessData.totStopLenT{i} = behParRun.totStopLenT(trInd);
            %
        
            sessData.numSamplesMean(:,i) = mean(behParRun.numSamples(trInd));
            sessData.numSamplesSEM(:,i) = std(behParRun.numSamples(trInd))/sqrt(nTr);
            sessData.maxSpeedMean(:,i) = mean(behParRun.maxSpeed(trInd));
            sessData.maxSpeedSEM(:,i) = std(behParRun.maxSpeed(trInd))/sqrt(nTr);
            sessData.meanSpeedMean(:,i) = mean(behParRun.meanSpeed(trInd));
            sessData.meanSpeedSEM(:,i) = std(behParRun.meanSpeed(trInd))/sqrt(nTr);
            sessData.maxSpeedRunMean(:,i) = mean(behParRun.maxSpeedRun(trInd));
            sessData.maxSpeedRunSEM(:,i) = std(behParRun.maxSpeedRun(trInd))/sqrt(nTr);
            sessData.meanSpeedRunMean(:,i) = mean(behParRun.meanSpeedRun(trInd));
            sessData.meanSpeedRunSEM(:,i) = std(behParRun.meanSpeedRun(trInd))/sqrt(nTr);
            sessData.maxRunLenTMean(:,i) = mean(behParRun.maxRunLenT(trInd));
            sessData.maxRunLenTSEM(:,i) = std(behParRun.maxRunLenT(trInd))/sqrt(nTr);
            sessData.totRunLenTMean(:,i) = mean(behParRun.totRunLenT(trInd));
            sessData.totRunLenTSEM(:,i) = std(behParRun.totRunLenT(trInd))/sqrt(nTr);
            sessData.numRunMean(:,i) = mean(behParRun.numRun(trInd));
            sessData.numRunSEM(:,i) = std(behParRun.numRun(trInd))/sqrt(nTr);
            sessData.maxAccMean(:,i) = mean(behParRun.maxAcc(trInd));
            sessData.maxAccSEM(:,i) = std(behParRun.maxAcc(trInd))/sqrt(nTr);
            sessData.meanAccMean(:,i) = mean(behParRun.meanAcc(trInd));
            sessData.meanAccSEM(:,i) = std(behParRun.meanAcc(trInd))/sqrt(nTr);
            sessData.totStopLenTMean(:,i) = mean(behParRun.totStopLenT(trInd));
            sessData.percStopTr(:,i) =  sum(behParRun.totStopLenT(trInd) > 0)/nTr;
            sessData.totStopLenTSEM(:,i) = std(behParRun.totStopLenT(trInd))/sqrt(nTr);
            sessData.startCueToRunMean(:,i) = mean(behParRun.startCueToRun(trInd));
            sessData.startCueToRunSEM(:,i) = std(behParRun.startCueToRun(trInd))/sqrt(nTr);
            sessData.numLicksBefRewMean(:,i) = mean(behParRun.numLicksBefRew(trInd));
            sessData.numLicksBefRewSEM(:,i) = std(behParRun.numLicksBefRew(trInd))/sqrt(nTr);
            sessData.numLicksRewMean(:,i) = mean(behParRun.numLicksRew(trInd));
            sessData.numLicksRewSEM(:,i) = std(behParRun.numLicksRew(trInd))/sqrt(nTr);
            sessData.rewToRunMean(:,i) = mean(behParRun.rewardToRun(trInd));
            sessData.rewToRunSEM(:,i) = mean(behParRun.rewardToRun(trInd))/sqrt(nTr);

            %% added by Yingxue on 4/19/2021
            sessData.lickSimMean(:,i) = lickOverDistSim.meanRun;
            sessData.lickSimSem(:,i) = lickOverDistSim.semRun;

            sessData.speedSimMean(:,i) = speedOverDistSim.meanRun;
            sessData.speedSimSem(:,i) = speedOverDistSim.semRun;
            
            sessData.lickEucMean(:,i) = lickOverDistSim.meanEucRun;
            sessData.lickEucSem(:,i) = lickOverDistSim.semEucRun;

            sessData.speedEucMean(:,i) = speedOverDistSim.meanEucRun;
            sessData.speedEucSem(:,i) = speedOverDistSim.semEucRun;
            %%
            
            indN = ~isnan(behParRun.med1stFiveLickDist(trInd));
            med1stFiveLickDist = behParRun.med1stFiveLickDist(trInd);
            sessData.med1stFiveLickDistMean(:,i) = mean(med1stFiveLickDist(indN));
            sessData.med1stFiveLickDistSEM(:,i) = std(med1stFiveLickDist(indN))/sqrt(sum(indN));

            indN = ~isnan(behParRun.medLickDist(trInd));
            medLickDist = behParRun.medLickDist(trInd);
            sessData.medLickDistMean(:,i) = mean(medLickDist(indN));
            sessData.medLickDistSEM(:,i) = std(medLickDist(indN))/sqrt(sum(indN));  

            indN = ~isnan(behParRun.med1stFiveLickDistBefRew(trInd));
            med1stFiveLickDistBefRew = behParRun.med1stFiveLickDistBefRew(trInd);
            sessData.med1stFiveLickDistBefRewMean(:,i) = mean(med1stFiveLickDistBefRew(indN));
            sessData.med1stFiveLickDistBefRewSEM(:,i) = std(med1stFiveLickDistBefRew(indN))/sqrt(sum(indN));

            indN = ~isnan(behParRun.medLickDistBefRew(trInd));
            medLickDistBefRew = behParRun.medLickDistBefRew(trInd);
            sessData.medLickDistBefRewMean(:,i) = mean(medLickDistBefRew(indN));
            sessData.medLickDistBefRewSEM(:,i) = std(medLickDistBefRew(indN))/sqrt(sum(indN));

            sessData.percRewarded(:,i) = sum(behParRun.rewarded(trInd)==1)/nTr;
            sessData.percNonStop(:,i) = sum(behParRun.nonStop(trInd))/nTr;
            sessData.speedProfile(:,i) = [{behParRun.speedProfile(trInd,:)/10};...
                {behParRun.speedProfile(trInd,:)/10};...
                {behParRun.speedProfile(trInd,:)/10}];
            sessData.lickProfile(:,i) = [{behParRun.lickProfile(trInd,:)};...
                {behParRun.lickProfile(trInd,:)};...
                {behParRun.lickProfile(trInd,:)}];

            indN = behParRun.pumpLfpInd(trInd) > 0;
            sessData.pumpLfpIndMean(:,i) = mean(behParRun.pumpLfpInd(trInd(indN)));
            sessData.pumpLfpIndSEM(:,i) = std(behParRun.pumpLfpInd(trInd(indN)))/sqrt(sum(indN));
            sessData.pumpMMMean(:,i) = mean(behParRun.pumpMM(trInd(indN)));
            sessData.pumpMMSEM(:,i) = std(behParRun.pumpMM(trInd(indN)))/sqrt(sum(indN));
        end
    end
    
    %% opto
    barPlot(1:nSess,sessData.meanSpeedMean(1,:),sessData.meanSpeedSEM(1,:),'Session','Mean speed (mm/s)','Run onset Opt');
    print ('-painters', '-dpdf', [path recName '_meanSpeedRunOpt'], '-r600')
    
    barPlot(1:nSess,sessData.numRunMean(1,:),sessData.numRunSEM(1,:),'Session','No. run segments','Run onset Opt');
    print ('-painters', '-dpdf', [path recName '_numRunRunOpt'], '-r600')
    
    barPlot(1:nSess,sessData.numLicksBefRewMean(1,:),sessData.numLicksBefRewSEM(1,:),'Session','Num licks bef. rew.','Run onset Opt');
    print ('-painters', '-dpdf', [path recName '_numLicksBefRewOpt'], '-r600')
   
    barPlot(1:nSess,sessData.med1stFiveLickDistBefRewMean(1,:),sessData.med1stFiveLickDistBefRewSEM(1,:),'Session','Med. dist. first five licks before reward (mm)','Run onset Opt');
    print ('-painters', '-dpdf', [path recName '_med1stFiveLickDistBefRewOpt'], '-r600')

    barPlot(1:nSess,sessData.percRewarded(1,:),zeros(1,nSess),'Session','Performance (ratio correct)','Run onset Opt');
    print ('-painters', '-dpdf', [path recName '_percRewardedRunOpt'], '-r600')

    barPlot(1:nSess,sessData.pumpMMMean(1,:),sessData.pumpMMSEM(1,:),'Session','pump location (mm)','Run onset Opt');
    print ('-painters', '-dpdf', [path recName '_pumpMMOpt'], '-r600')
    
    barPlot(1:nSess,sessData.pumpLfpIndMean(1,:)/1250,sessData.pumpLfpIndSEM(1,:)/1250,'Session','pump time (s)','Run onset Opt');
    print ('-painters', '-dpdf', [path recName '_pumpTimeOpt'], '-r600')
    
    barPlot(1:nSess,sessData.speedSimMean(1,:),sessData.speedSimSem(1,:),'Session','speed over distance similarity','Run onset Opt');
    print ('-painters', '-dpdf', [path recName '_speedSimRunOpt'], '-r600')
    
    barPlot(1:nSess,sessData.lickSimMean(1,:),sessData.lickSimSem(1,:),'Session','lick over distance similarity','Run onset Opt');
    print ('-painters', '-dpdf', [path recName '_lickSimRunOpt'], '-r600')
    
    barPlot(1:nSess,sessData.speedEucMean(1,:),sessData.speedEucSem(1,:),'Session','speed over distance euclidean','Run onset Opt');
    print ('-painters', '-dpdf', [path recName '_speedEucRunOpt'], '-r600')
    
    barPlot(1:nSess,sessData.lickEucMean(1,:),sessData.lickEucSem(1,:),'Session','lick over distance euclidean','Run onset Opt');
    print ('-painters', '-dpdf', [path recName '_lickEucRunOpt'], '-r600')
    
    barPlot(1:nSess,sessData.totStopLenTMean(1,:),sessData.totStopLenTSEM(1,:),'Session','Total stop time (s)','Run onset Opt');
    print ('-painters', '-dpdf', [path recName '_totStopLenTRunOpt'], '-r600')
        
    %% opto ctrl
    barPlot(1:nSess,sessData.meanSpeedMean(2,:),sessData.meanSpeedSEM(2,:),'Session','Mean speed (mm/s)','Run onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_meanSpeedRunOptCtrl'], '-r600')

    barPlot(1:nSess,sessData.numRunMean(2,:),sessData.numRunSEM(2,:),'Session','No. run segments','Run onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_numRunRunOptCtrl'], '-r600')

    barPlot(1:nSess,sessData.numLicksBefRewMean(2,:),sessData.numLicksBefRewSEM(2,:),'Session','Num licks bef. rew.','Run onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_numLicksBefRewOptCtrl'], '-r600')

    barPlot(1:nSess,sessData.med1stFiveLickDistBefRewMean(2,:),sessData.med1stFiveLickDistBefRewSEM(2,:),'Session','Med. dist. first five licks before reward (mm)','Run onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_med1stFiveLickDistBefRewOptCtrl'], '-r600')

    barPlot(1:nSess,sessData.percRewarded(2,:),zeros(1,nSess),'Session','Performance (ratio correct)','Run onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_percRewardedRunOptCtrl'], '-r600')

    barPlot(1:nSess,sessData.pumpMMMean(2,:),sessData.pumpMMSEM(2,:),'Session','pump location (mm)','Run onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_pumpMMOptCtrl'], '-r600')

    barPlot(1:nSess,sessData.pumpLfpIndMean(2,:)/1250,sessData.pumpLfpIndSEM(2,:)/1250,'Session','pump time (s)','Run onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_pumpTimeOptCtrl'], '-r600')

    barPlot(1:nSess,sessData.speedSimMean(2,:),sessData.speedSimSem(2,:),'Session','speed over distance similarity','Run onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_speedSimRunOptCtrl'], '-r600')

    barPlot(1:nSess,sessData.lickSimMean(2,:),sessData.lickSimSem(2,:),'Session','lick over distance similarity','Run onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_lickSimRunOptCtrl'], '-r600')
    
    barPlot(1:nSess,sessData.speedEucMean(2,:),sessData.speedEucSem(2,:),'Session','speed over distance euclidean','Run onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_speedEucRunOptCtrl'], '-r600')
    
    barPlot(1:nSess,sessData.lickEucMean(2,:),sessData.lickEucSem(2,:),'Session','lick over distance euclidean','Run onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_lickEucRunOptCtrl'], '-r600')
    
    barPlot(1:nSess,sessData.totStopLenTMean(2,:),sessData.totStopLenTSEM(2,:),'Session','Total stop time (s)','Run onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_totStopLenTRunOptCtrl'], '-r600')
    
    %% opto ctrl1, added by Yingxue on 1/25/2022
    barPlot(1:nSess,sessData.meanSpeedMean(3,:),sessData.meanSpeedSEM(3,:),'Session','Mean speed (mm/s)','Run onset Ctrl aft Opt');
    print ('-painters', '-dpdf', [path recName '_meanSpeedRunOptCtrl1'], '-r600')

    barPlot(1:nSess,sessData.numRunMean(3,:),sessData.numRunSEM(3,:),'Session','No. run segments','Run onset Ctrl aft Opt');
    print ('-painters', '-dpdf', [path recName '_numRunRunOptCtrl1'], '-r600')

    barPlot(1:nSess,sessData.numLicksBefRewMean(3,:),sessData.numLicksBefRewSEM(3,:),'Session','Num licks bef. rew.','Run onset Ctrl aft Opt');
    print ('-painters', '-dpdf', [path recName '_numLicksBefRewOptCtrl1'], '-r600')

    barPlot(1:nSess,sessData.med1stFiveLickDistBefRewMean(3,:),sessData.med1stFiveLickDistBefRewSEM(3,:),'Session','Med. dist. first five licks before reward (mm)','Run onset Ctrl aft Opt');
    print ('-painters', '-dpdf', [path recName '_med1stFiveLickDistBefRewOptCtrl1'], '-r600')

    barPlot(1:nSess,sessData.percRewarded(3,:),zeros(1,nSess),'Session','Performance (ratio correct)','Run onset Ctrl aft Opt');
    print ('-painters', '-dpdf', [path recName '_percRewardedRunOptCtrl1'], '-r600')

    barPlot(1:nSess,sessData.pumpMMMean(3,:),sessData.pumpMMSEM(3,:),'Session','pump location (mm)','Run onset Ctrl aft Opt');
    print ('-painters', '-dpdf', [path recName '_pumpMMOptCtrl1'], '-r600')

    barPlot(1:nSess,sessData.pumpLfpIndMean(3,:)/1250,sessData.pumpLfpIndSEM(3,:)/1250,'Session','pump time (s)','Run onset Ctrl aft Opt');
    print ('-painters', '-dpdf', [path recName '_pumpTimeOptCtrl1'], '-r600')

    barPlot(1:nSess,sessData.totStopLenTMean(3,:),sessData.totStopLenTSEM(3,:),'Session','Total stop time (s)','Run onset Ctrl aft Opt');
    print ('-painters', '-dpdf', [path recName '_totStopLenTRunOptCtrl1'], '-r600')
end

function sessData = compBehRew(path,recName,sessNo,mazeSess,onlyRun)
    nSess = length(sessNo);
    sessData.sessNo = sessNo;
    sessData.numSamplesMean = zeros(1,nSess);
    sessData.numSamplesSEM = zeros(1,nSess);
    sessData.maxSpeedMean = zeros(1,nSess);
    sessData.maxSpeedSEM = zeros(1,nSess);
    sessData.meanSpeedMean = zeros(1,nSess);
    sessData.meanSpeedSEM = zeros(1,nSess);
    sessData.maxSpeedRunMean = zeros(1,nSess);
    sessData.maxSpeedRunSEM = zeros(1,nSess);
    sessData.meanSpeedRunMean = zeros(1,nSess);
    sessData.meanSpeedRunSEM = zeros(1,nSess);
    sessData.maxRunLenTMean = zeros(1,nSess);
    sessData.maxRunLenTSEM = zeros(1,nSess);
    sessData.totRunLenTMean = zeros(1,nSess);
    sessData.totRunLenTSEM = zeros(1,nSess);
    sessData.numRunMean = zeros(1,nSess);
    sessData.numRunSEM = zeros(1,nSess);
    sessData.maxAccMean = zeros(1,nSess);
    sessData.maxAccSEM = zeros(1,nSess);
    sessData.meanAccMean = zeros(1,nSess);
    sessData.meanAccSEM = zeros(1,nSess);
    sessData.totStopLenTMean = zeros(1,nSess);
    sessData.totStopLenTSEM = zeros(1,nSess); 
    
    sessData.speedSimMean = zeros(1,nSess);
    sessData.speedSimSem = zeros(1,nSess);
    sessData.lickSimMean = zeros(1,nSess);
    sessData.lickSimSem = zeros(1,nSess);
    
    sessData.speedEucMean = zeros(1,nSess);
    sessData.speedEucSem = zeros(1,nSess);
    sessData.lickEucMean = zeros(1,nSess);
    sessData.lickEucSem = zeros(1,nSess);
    
    for i = 1:nSess
        fileName = [recName '-0' num2str(sessNo(i))];
        pathTmp = [path fileName '\'];
        
        fullpath = [pathTmp fileName '_behPar_msess' num2str(mazeSess(i)) '.mat'];
        if(~exist(fullpath,'file'))
            disp([fullpath ' does not exist.']);
            return;
        end
        load(fullpath,'behParRew','param');
        
        fullPath = [pathTmp fileName '_lickDistSim_msess' num2str(mazeSess(i)) '.mat'];
        if(exist(fullPath) == 0)
            disp('The lick over distance similarity file does not exist');
            return;
        end
        load(fullPath,'lickOverDistSim');
        
        fullPath = [pathTmp fileName '_speedDistSim_msess' num2str(mazeSess(i)) '_Run' num2str(onlyRun) '.mat'];
        if(exist(fullPath) == 0)
            disp('The speed over distance similarity file does not exist');
            return;
        end
        load(fullPath,'speedOverDistSim');
        
        trInd = lickOverDistSim.trInd; %param.startTr:param.endTr;
        nTr = length(trInd);
        
        sessData.numSamplesMean(i) = mean(behParRew.numSamples(trInd));
        sessData.numSamplesSEM(i) = std(behParRew.numSamples(trInd))/sqrt(nTr);
        sessData.maxSpeedMean(i) = mean(behParRew.maxSpeed(trInd));
        sessData.maxSpeedSEM(i) = std(behParRew.maxSpeed(trInd))/sqrt(nTr);
        sessData.meanSpeedMean(i) = mean(behParRew.meanSpeed(trInd));
        sessData.meanSpeedSEM(i) = std(behParRew.meanSpeed(trInd))/sqrt(nTr);
        sessData.maxSpeedRunMean(i) = mean(behParRew.maxSpeedRun(trInd));
        sessData.maxSpeedRunSEM(i) = std(behParRew.maxSpeedRun(trInd))/sqrt(nTr);
        sessData.meanSpeedRunMean(i) = mean(behParRew.meanSpeedRun(trInd));
        sessData.meanSpeedRunSEM(i) = std(behParRew.meanSpeedRun(trInd))/sqrt(nTr);
        sessData.maxRunLenTMean(i) = mean(behParRew.maxRunLenT(trInd));
        sessData.maxRunLenTSEM(i) = std(behParRew.maxRunLenT(trInd))/sqrt(nTr);
        sessData.totRunLenTMean(i) = mean(behParRew.totRunLenT(trInd));
        sessData.totRunLenTSEM(i) = std(behParRew.totRunLenT(trInd))/sqrt(nTr);
        sessData.numRunMean(i) = mean(behParRew.numRun(trInd));
        sessData.numRunSEM(i) = std(behParRew.numRun(trInd))/sqrt(nTr);
        sessData.maxAccMean(i) = mean(behParRew.maxAcc(trInd));
        sessData.maxAccSEM(i) = std(behParRew.maxAcc(trInd))/sqrt(nTr);
        sessData.meanAccMean(i) = mean(behParRew.meanAcc(trInd));
        sessData.meanAccSEM(i) = std(behParRew.meanAcc(trInd))/sqrt(nTr);
        sessData.totStopLenTMean(i) = mean(behParRew.totStopLenTRew(trInd));
        sessData.totStopLenTSEM(i) = std(behParRew.totStopLenTRew(trInd))/sqrt(nTr);
        
        %% added by Yingxue on 4/19/2021
        sessData.lickSimMean(i) = lickOverDistSim.meanRew;
        sessData.lickSimSem(i) = lickOverDistSim.semRew;
        
        sessData.speedSimMean(i) = speedOverDistSim.meanRew;
        sessData.speedSimSem(i) = speedOverDistSim.semRew;  
        
        sessData.lickEucMean(i) = lickOverDistSim.meanEucRew;
        sessData.lickEucSem(i) = lickOverDistSim.semEucRew;
        
        sessData.speedEucMean(i) = speedOverDistSim.meanEucRew;
        sessData.speedEucSem(i) = speedOverDistSim.semEucRew;  
    end
    
    barPlot(1:nSess,sessData.numSamplesMean,sessData.numSamplesSEM,'Session','No. samples / Trial','Reward');
    print ('-painters', '-dpdf', [path recName '_numSamplesRew'], '-r600')
    barPlot(1:nSess,sessData.meanSpeedMean,sessData.meanSpeedSEM,'Session','Mean speed (mm/s)','Reward');
    print ('-painters', '-dpdf', [path recName '_meanSpeedRew'], '-r600')
    barPlot(1:nSess,sessData.meanSpeedRunMean,sessData.meanSpeedRunSEM,'Session','Mean speed run only (mm/s)','Reward');
    print ('-painters', '-dpdf', [path recName '_meanSpeedRunRew'], '-r600')
    barPlot(1:nSess,sessData.maxRunLenTMean,sessData.maxRunLenTSEM,'Session','Max run length (s)','Reward');
    print ('-painters', '-dpdf', [path recName '_maxRunLenTRew'], '-r600')
    barPlot(1:nSess,sessData.totRunLenTMean,sessData.totRunLenTSEM,'Session','Total run length (s)','Reward');
    print ('-painters', '-dpdf', [path recName '_totRunLenTRew'], '-r600')
    barPlot(1:nSess,sessData.numRunMean,sessData.numRunSEM,'Session','No. run segments','Reward');
    print ('-painters', '-dpdf', [path recName '_numRunRew'], '-r600')
    barPlot(1:nSess,sessData.meanAccMean,sessData.meanAccSEM,'Session','Mean acceleration (mm/s^2)','Reward');
    print ('-painters', '-dpdf', [path recName '_meanAccRew'], '-r600')
    barPlot(1:nSess,sessData.totStopLenTMean,sessData.totStopLenTSEM,'Session','Total stop time (s)','Reward');
    print ('-painters', '-dpdf', [path recName '_totStopLenTRew'], '-r600')
    barPlot(1:nSess,sessData.speedSimMean,sessData.speedSimSem,'Session','speed over distance similarity','Reward');
    print ('-painters', '-dpdf', [path recName '_speedSimRew'], '-r600')    
    barPlot(1:nSess,sessData.lickSimMean,sessData.lickSimSem,'Session','lick over distance similarity','Reward');
    print ('-painters', '-dpdf', [path recName '_lickSimRew'], '-r600')
    barPlot(1:nSess,sessData.speedEucMean,sessData.speedEucSem,'Session','speed over distance euclidean','Reward');
    print ('-painters', '-dpdf', [path recName '_speedEucRew'], '-r600')    
    barPlot(1:nSess,sessData.lickEucMean,sessData.lickEucSem,'Session','lick over distance euclidean','Reward');
    print ('-painters', '-dpdf', [path recName '_lickEucRew'], '-r600')
end

%% added by Yingxue on 10/01/2020
function sessData = compBehRewOpt(path,recName,sessNo,mazeSess,onlyRun)
    nSess = length(sessNo);
    sessData.sessNo = sessNo;
    sessData.numSamplesMean = zeros(2,nSess);
    sessData.numSamplesSEM = zeros(2,nSess);
    sessData.maxSpeedMean = zeros(2,nSess);
    sessData.maxSpeedSEM = zeros(2,nSess);
    sessData.meanSpeedMean = zeros(2,nSess);
    sessData.meanSpeedSEM = zeros(2,nSess);
    sessData.maxSpeedRunMean = zeros(2,nSess);
    sessData.maxSpeedRunSEM = zeros(2,nSess);
    sessData.meanSpeedRunMean = zeros(2,nSess);
    sessData.meanSpeedRunSEM = zeros(2,nSess);
    sessData.maxRunLenTMean = zeros(2,nSess);
    sessData.maxRunLenTSEM = zeros(2,nSess);
    sessData.totRunLenTMean = zeros(2,nSess);
    sessData.totRunLenTSEM = zeros(2,nSess);
    sessData.numRunMean = zeros(2,nSess);
    sessData.numRunSEM = zeros(2,nSess);
    sessData.maxAccMean = zeros(2,nSess);
    sessData.maxAccSEM = zeros(2,nSess);
    sessData.meanAccMean = zeros(2,nSess);
    sessData.meanAccSEM = zeros(2,nSess);
    sessData.totStopLenTMean = zeros(2,nSess);
    sessData.totStopLenTSEM = zeros(2,nSess); 
    
    sessData.speedSimMean = zeros(2,nSess);
    sessData.speedSimSem = zeros(2,nSess);
    sessData.lickSimMean = zeros(2,nSess);
    sessData.lickSimSem = zeros(2,nSess);
    
    sessData.speedEucMean = zeros(2,nSess);
    sessData.speedEucSem = zeros(2,nSess);
    sessData.lickEucMean = zeros(2,nSess);
    sessData.lickEucSem = zeros(2,nSess);
    
    for i = 1:nSess
        fileName = [recName '-0' num2str(sessNo(i))];
        pathTmp = [path fileName '\'];
        
        fullpath = [pathTmp fileName '_behPar_msess' num2str(mazeSess(i)) '.mat'];
        if(~exist(fullpath,'file'))
            disp([fullpath ' does not exist.']);
            return;
        end
        load(fullpath,'behParRew','param');
        
        fullPath = [pathTmp fileName '_lickDist_msess' num2str(mazeSess(i)) '.mat']; 
        if(exist(fullPath) == 0)
            disp('The _lickDist file does not exist');
            return;
        end
        load(fullPath,'lickOverDistOpt','lickOverDistOptCtrl');
        
        fullPath = [pathTmp fileName '_lickDistSim_msess' num2str(mazeSess(i)) '.mat'];
        if(exist(fullPath) == 0)
            disp('The lick over distance similarity file does not exist');
            return;
        end
        load(fullPath,'lickOverDistSim','lickOverDistSimOpt','lickOverDistSimOptCtrl');
        
        fullPath = [pathTmp fileName '_speedDistSim_msess' num2str(mazeSess(i)) '_Run' num2str(onlyRun) '.mat'];
        if(exist(fullPath) == 0)
            disp('The speed over distance similarity file does not exist');
            return;
        end
        load(fullPath,'speedOverDistSim','speedOverDistSimOpt','speedOverDistSimOptCtrl');
        
        if(~isempty(lickOverDistOpt.trInd))
            trIndOpt = lickOverDistOpt.trInd;
            trIndCtrl = lickOverDistOptCtrl.trInd;
            nTrOpt = length(trIndOpt);
            nTrCtrl = length(trIndCtrl);
            
            sessData.numSamplesMean(:,i) = [mean(behParRew.numSamples(trIndOpt));...
                mean(behParRew.numSamples(trIndCtrl))];
            sessData.numSamplesSEM(:,i) = [std(behParRew.numSamples(trIndOpt))/sqrt(nTrOpt);...
                std(behParRew.numSamples(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.maxSpeedMean(:,i) = [mean(behParRew.maxSpeed(trIndOpt));...
                mean(behParRew.maxSpeed(trIndCtrl))];
            sessData.maxSpeedSEM(:,i) = [std(behParRew.maxSpeed(trIndOpt))/sqrt(nTrOpt);...
                std(behParRew.maxSpeed(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.meanSpeedMean(:,i) = [mean(behParRew.meanSpeed(trIndOpt));...
                mean(behParRew.meanSpeed(trIndCtrl))];
            sessData.meanSpeedSEM(:,i) = [std(behParRew.meanSpeed(trIndOpt))/sqrt(nTrOpt);...
                std(behParRew.meanSpeed(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.maxSpeedRunMean(:,i) = [mean(behParRew.maxSpeedRun(trIndOpt));...
                mean(behParRew.maxSpeedRun(trIndCtrl))];
            sessData.maxSpeedRunSEM(:,i) = [std(behParRew.maxSpeedRun(trIndOpt))/sqrt(nTrOpt);...
                mean(behParRew.maxSpeedRun(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.meanSpeedRunMean(:,i) = [mean(behParRew.meanSpeedRun(trIndOpt));...
                mean(behParRew.meanSpeedRun(trIndCtrl))];
            sessData.meanSpeedRunSEM(:,i) = [std(behParRew.meanSpeedRun(trIndOpt))/sqrt(nTrOpt);...
                std(behParRew.meanSpeedRun(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.maxRunLenTMean(:,i) = [mean(behParRew.maxRunLenT(trIndOpt));...
                mean(behParRew.maxRunLenT(trIndCtrl))];
            sessData.maxRunLenTSEM(:,i) = [std(behParRew.maxRunLenT(trIndOpt))/sqrt(nTrOpt);...
                std(behParRew.maxRunLenT(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.totRunLenTMean(:,i) = [mean(behParRew.totRunLenT(trIndOpt));...
                mean(behParRew.totRunLenT(trIndCtrl))];
            sessData.totRunLenTSEM(:,i) = [std(behParRew.totRunLenT(trIndOpt))/sqrt(nTrOpt);...
                std(behParRew.totRunLenT(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.numRunMean(:,i) = [mean(behParRew.numRun(trIndOpt));...
                mean(behParRew.numRun(trIndCtrl))];
            sessData.numRunSEM(:,i) = [std(behParRew.numRun(trIndOpt))/sqrt(nTrOpt);...
                std(behParRew.numRun(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.maxAccMean(:,i) = [mean(behParRew.maxAcc(trIndOpt));...
                mean(behParRew.maxAcc(trIndCtrl))];
            sessData.maxAccSEM(:,i) = [std(behParRew.maxAcc(trIndOpt))/sqrt(nTrOpt);...
                std(behParRew.maxAcc(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.meanAccMean(:,i) = [mean(behParRew.meanAcc(trIndOpt));...
                mean(behParRew.meanAcc(trIndCtrl))];
            sessData.meanAccSEM(:,i) = [std(behParRew.meanAcc(trIndOpt))/sqrt(nTrOpt);...
                std(behParRew.meanAcc(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.totStopLenTMean(:,i) = [mean(behParRew.totStopLenTRew(trIndOpt));...
                mean(behParRew.totStopLenTRew(trIndCtrl))];
            sessData.totStopLenTSEM(:,i) = [std(behParRew.totStopLenTRew(trIndOpt))/sqrt(nTrOpt);...
                std(behParRew.totStopLenTRew(trIndCtrl))/sqrt(nTrCtrl)];
            
            %% added by Yingxue on 4/19/2021
            sessData.lickSimMean(:,i) = [lickOverDistSimOpt.meanRew;lickOverDistSimOptCtrl.meanRew];
            sessData.lickSimSem(:,i) = [lickOverDistSimOpt.semRew;lickOverDistSimOptCtrl.semRew];

            sessData.speedSimMean(:,i) = [speedOverDistSimOpt.meanRew;speedOverDistSimOptCtrl.meanRew];
            sessData.speedSimSem(:,i) = [speedOverDistSimOpt.semRew;speedOverDistSimOptCtrl.semRew];  
            
            sessData.lickEucMean(:,i) = [lickOverDistSimOpt.meanEucRew;lickOverDistSimOptCtrl.meanEucRew];
            sessData.lickEucSem(:,i) = [lickOverDistSimOpt.semEucRew;lickOverDistSimOptCtrl.semEucRew];

            sessData.speedEucMean(:,i) = [speedOverDistSimOpt.meanEucRew;speedOverDistSimOptCtrl.meanEucRew];
            sessData.speedEucSem(:,i) = [speedOverDistSimOpt.semEucRew;speedOverDistSimOptCtrl.semEucRew]; 
        else
            trInd = param.startTr:param.endTr;
            nTr = length(trInd);

            sessData.numSamplesMean(:,i) = mean(behParRew.numSamples(trInd));
            sessData.numSamplesSEM(:,i) = std(behParRew.numSamples(trInd))/sqrt(nTr);
            sessData.maxSpeedMean(:,i) = mean(behParRew.maxSpeed(trInd));
            sessData.maxSpeedSEM(:,i) = std(behParRew.maxSpeed(trInd))/sqrt(nTr);
            sessData.meanSpeedMean(:,i) = mean(behParRew.meanSpeed(trInd));
            sessData.meanSpeedSEM(:,i) = std(behParRew.meanSpeed(trInd))/sqrt(nTr);
            sessData.maxSpeedRunMean(:,i) = mean(behParRew.maxSpeedRun(trInd));
            sessData.maxSpeedRunSEM(:,i) = std(behParRew.maxSpeedRun(trInd))/sqrt(nTr);
            sessData.meanSpeedRunMean(:,i) = mean(behParRew.meanSpeedRun(trInd));
            sessData.meanSpeedRunSEM(:,i) = std(behParRew.meanSpeedRun(trInd))/sqrt(nTr);
            sessData.maxRunLenTMean(:,i) = mean(behParRew.maxRunLenT(trInd));
            sessData.maxRunLenTSEM(:,i) = std(behParRew.maxRunLenT(trInd))/sqrt(nTr);
            sessData.totRunLenTMean(:,i) = mean(behParRew.totRunLenT(trInd));
            sessData.totRunLenTSEM(:,i) = std(behParRew.totRunLenT(trInd))/sqrt(nTr);
            sessData.numRunMean(:,i) = mean(behParRew.numRun(trInd));
            sessData.numRunSEM(:,i) = std(behParRew.numRun(trInd))/sqrt(nTr);
            sessData.maxAccMean(:,i) = mean(behParRew.maxAcc(trInd));
            sessData.maxAccSEM(:,i) = std(behParRew.maxAcc(trInd))/sqrt(nTr);
            sessData.meanAccMean(:,i) = mean(behParRew.meanAcc(trInd));
            sessData.meanAccSEM(:,i) = std(behParRew.meanAcc(trInd))/sqrt(nTr);
            sessData.totStopLenTMean(:,i) = mean(behParRew.totStopLenTRew(trInd));
            sessData.totStopLenTSEM(:,i) = std(behParRew.totStopLenTRew(trInd))/sqrt(nTr);
            
            %% added by Yingxue on 4/19/2021
            sessData.lickSimMean(:,i) = lickOverDistSim.meanRew;
            sessData.lickSimSem(:,i) = lickOverDistSim.semRew;

            sessData.speedSimMean(:,i) = speedOverDistSim.meanRew;
            sessData.speedSimSem(:,i) = speedOverDistSim.semRew;
            
            sessData.lickEucMean(:,i) = lickOverDistSim.meanEucRew;
            sessData.lickEucSem(:,i) = lickOverDistSim.semEucRew;

            sessData.speedEucMean(:,i) = speedOverDistSim.meanEucRew;
            sessData.speedEucSem(:,i) = speedOverDistSim.semEucRew;
        end
    end
    
    %% opto
    barPlot(1:nSess,sessData.numSamplesMean(1,:),sessData.numSamplesSEM(1,:),'Session','No. samples / Trial','Reward Opt');
    print ('-painters', '-dpdf', [path recName '_numSamplesRewOpt'], '-r600')
    barPlot(1:nSess,sessData.meanSpeedMean(1,:),sessData.meanSpeedSEM(1,:),'Session','Mean speed (mm/s)','Reward Opt');
    print ('-painters', '-dpdf', [path recName '_meanSpeedRewOpt'], '-r600')
    barPlot(1:nSess,sessData.meanSpeedRunMean(1,:),sessData.meanSpeedRunSEM(1,:),'Session','Mean speed run only (mm/s)','Reward Opt');
    print ('-painters', '-dpdf', [path recName '_meanSpeedRunRewOpt'], '-r600')
    barPlot(1:nSess,sessData.maxRunLenTMean(1,:),sessData.maxRunLenTSEM(1,:),'Session','Max run length (s)','Reward Opt');
    print ('-painters', '-dpdf', [path recName '_maxRunLenTRewOpt'], '-r600')
    barPlot(1:nSess,sessData.totRunLenTMean(1,:),sessData.totRunLenTSEM(1,:),'Session','Total run length (s)','Reward Opt');
    print ('-painters', '-dpdf', [path recName '_totRunLenTRewOpt'], '-r600')
    barPlot(1:nSess,sessData.numRunMean(1,:),sessData.numRunSEM(1,:),'Session','No. run segments','Reward Opt');
    print ('-painters', '-dpdf', [path recName '_numRunRewOpt'], '-r600')
    barPlot(1:nSess,sessData.meanAccMean(1,:),sessData.meanAccSEM(1,:),'Session','Mean acceleration (mm/s^2)','Reward Opt');
    print ('-painters', '-dpdf', [path recName '_meanAccRewOpt'], '-r600')
    barPlot(1:nSess,sessData.totStopLenTMean(1,:),sessData.totStopLenTSEM(1,:),'Session','Total stop time (s)','Reward Opt');
    print ('-painters', '-dpdf', [path recName '_totStopLenTRewOpt'], '-r600')
    barPlot(1:nSess,sessData.speedSimMean(1,:),sessData.speedSimSem(1,:),'Session','speed over distance similarity','Reward Opt');
    print ('-painters', '-dpdf', [path recName '_speedSimRewOpt'], '-r600')    
    barPlot(1:nSess,sessData.lickSimMean(1,:),sessData.lickSimSem(1,:),'Session','lick over distance similarity','Reward Opt');
    print ('-painters', '-dpdf', [path recName '_lickSimRewOpt'], '-r600')
    barPlot(1:nSess,sessData.speedEucMean(1,:),sessData.speedEucSem(1,:),'Session','speed over distance euclidean','Reward Opt');
    print ('-painters', '-dpdf', [path recName '_speedEucRewOpt'], '-r600')    
    barPlot(1:nSess,sessData.lickEucMean(1,:),sessData.lickEucSem(1,:),'Session','lick over distance euclidean','Reward Opt');
    print ('-painters', '-dpdf', [path recName '_lickEucRewOpt'], '-r600')
    
    %% opto ctrl
    barPlot(1:nSess,sessData.numSamplesMean(2,:),sessData.numSamplesSEM(2,:),'Session','No. samples / Trial','Reward Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_numSamplesRewOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.meanSpeedMean(2,:),sessData.meanSpeedSEM(2,:),'Session','Mean speed (mm/s)','Reward Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_meanSpeedRewOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.meanSpeedRunMean(2,:),sessData.meanSpeedRunSEM(2,:),'Session','Mean speed run only (mm/s)','Reward Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_meanSpeedRunRewOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.maxRunLenTMean(2,:),sessData.maxRunLenTSEM(2,:),'Session','Max run length (s)','Reward Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_maxRunLenTRewOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.totRunLenTMean(2,:),sessData.totRunLenTSEM(2,:),'Session','Total run length (s)','Reward Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_totRunLenTRewOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.numRunMean(2,:),sessData.numRunSEM(2,:),'Session','No. run segments','Reward Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_numRunRewOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.meanAccMean(2,:),sessData.meanAccSEM(2,:),'Session','Mean acceleration (mm/s^2)','Reward Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_meanAccRewOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.totStopLenTMean(2,:),sessData.totStopLenTSEM(2,:),'Session','Total stop time (s)','Reward Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_totStopLenTRewOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.speedSimMean(2,:),sessData.speedSimSem(2,:),'Session','speed over distance similarity','Reward Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_speedSimRewOptCtrl'], '-r600')    
    barPlot(1:nSess,sessData.lickSimMean(2,:),sessData.lickSimSem(2,:),'Session','lick over distance similarity','Reward Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_lickSimRewOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.speedEucMean(2,:),sessData.speedEucSem(2,:),'Session','speed over distance eulicdean','Reward Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_speedEucRewOptCtrl'], '-r600')    
    barPlot(1:nSess,sessData.lickEucMean(2,:),sessData.lickEucSem(2,:),'Session','lick over distance euclidean','Reward Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_lickEucRewOptCtrl'], '-r600')
end

function sessData = compBehCue(path,recName,sessNo,mazeSess,onlyRun)
    nSess = length(sessNo);
    sessData.sessNo = sessNo;
    sessData.numSamplesMean = zeros(1,nSess);
    sessData.numSamplesSEM = zeros(1,nSess);
    sessData.maxSpeedMean = zeros(1,nSess);
    sessData.maxSpeedSEM = zeros(1,nSess);
    sessData.meanSpeedMean = zeros(1,nSess);
    sessData.meanSpeedSEM = zeros(1,nSess);
    sessData.maxSpeedRunMean = zeros(1,nSess);
    sessData.maxSpeedRunSEM = zeros(1,nSess);
    sessData.meanSpeedRunMean = zeros(1,nSess);
    sessData.meanSpeedRunSEM = zeros(1,nSess);
    sessData.maxRunLenTMean = zeros(1,nSess);
    sessData.maxRunLenTSEM = zeros(1,nSess);
    sessData.totRunLenTMean = zeros(1,nSess);
    sessData.totRunLenTSEM = zeros(1,nSess);
    sessData.numRunMean = zeros(1,nSess);
    sessData.numRunSEM = zeros(1,nSess);
    sessData.maxAccMean = zeros(1,nSess);
    sessData.maxAccSEM = zeros(1,nSess);
    sessData.meanAccMean = zeros(1,nSess);
    sessData.meanAccSEM = zeros(1,nSess);
    sessData.totStopLenTMean = zeros(1,nSess);
    sessData.totStopLenTSEM = zeros(1,nSess);
    sessData.med1stFiveLickDistMean = zeros(1,nSess);
    sessData.med1stFiveLickDistSEM = zeros(1,nSess);
    sessData.medLickDistMean = zeros(1,nSess);
    sessData.medLickDistSEM = zeros(1,nSess);
    
    sessData.speedSimMean = zeros(1,nSess);
    sessData.speedSimSem = zeros(1,nSess);
    sessData.lickSimMean = zeros(1,nSess);
    sessData.lickSimSem = zeros(1,nSess);
    
    sessData.speedEucMean = zeros(1,nSess);
    sessData.speedEucSem = zeros(1,nSess);
    sessData.lickEucMean = zeros(1,nSess);
    sessData.lickEucSem = zeros(1,nSess);
    
    for i = 1:nSess
        fileName = [recName '-0' num2str(sessNo(i))];
        pathTmp = [path fileName '\'];
        
        fullpath = [pathTmp fileName '_behPar_msess' num2str(mazeSess(i)) '.mat'];
        if(~exist(fullpath,'file'))
            disp([fullpath ' does not exist.']);
            return;
        end
        load(fullpath,'behParCue','param');
        
        fullPath = [pathTmp fileName '_lickDistSim_msess' num2str(mazeSess(i)) '.mat'];
        if(exist(fullPath) == 0)
            disp('The lick over distance similarity file does not exist');
            return;
        end
        load(fullPath,'lickOverDistSim');
        
        fullPath = [pathTmp fileName '_speedDistSim_msess' num2str(mazeSess(i)) '_Run' num2str(onlyRun) '.mat'];
        if(exist(fullPath) == 0)
            disp('The speed over distance similarity file does not exist');
            return;
        end
        load(fullPath,'speedOverDistSim');
        
        trInd = lickOverDistSim.trInd; %param.startTr:param.endTr;
        nTr = length(trInd);
        
        sessData.numSamplesMean(i) = mean(behParCue.numSamples(trInd));
        sessData.numSamplesSEM(i) = std(behParCue.numSamples(trInd))/sqrt(nTr);
        sessData.maxSpeedMean(i) = mean(behParCue.maxSpeed(trInd));
        sessData.maxSpeedSEM(i) = std(behParCue.maxSpeed(trInd))/sqrt(nTr);
        sessData.meanSpeedMean(i) = mean(behParCue.meanSpeed(trInd));
        sessData.meanSpeedSEM(i) = std(behParCue.meanSpeed(trInd))/sqrt(nTr);
        sessData.maxSpeedRunMean(i) = mean(behParCue.maxSpeedRun(trInd));
        sessData.maxSpeedRunSEM(i) = std(behParCue.maxSpeedRun(trInd))/sqrt(nTr);
        sessData.meanSpeedRunMean(i) = mean(behParCue.meanSpeedRun(trInd));
        sessData.meanSpeedRunSEM(i) = std(behParCue.meanSpeedRun(trInd))/sqrt(nTr);
        sessData.maxRunLenTMean(i) = mean(behParCue.maxRunLenT(trInd));
        sessData.maxRunLenTSEM(i) = std(behParCue.maxRunLenT(trInd))/sqrt(nTr);
        sessData.totRunLenTMean(i) = mean(behParCue.totRunLenT(trInd));
        sessData.totRunLenTSEM(i) = std(behParCue.totRunLenT(trInd))/sqrt(nTr);
        sessData.numRunMean(i) = mean(behParCue.numRun(trInd));
        sessData.numRunSEM(i) = std(behParCue.numRun(trInd))/sqrt(nTr);
        sessData.maxAccMean(i) = mean(behParCue.maxAcc(trInd));
        sessData.maxAccSEM(i) = std(behParCue.maxAcc(trInd))/sqrt(nTr);
        sessData.meanAccMean(i) = mean(behParCue.meanAcc(trInd));
        sessData.meanAccSEM(i) = std(behParCue.meanAcc(trInd))/sqrt(nTr);
        sessData.totStopLenTMean(i) = mean(behParCue.totStopLenT(trInd));
        sessData.totStopLenTSEM(i) = std(behParCue.totStopLenT(trInd))/sqrt(nTr);
        sessData.numLicksBefRewMean(i) = mean(behParCue.numLicksBefRew(trInd));
        sessData.numLicksBefRewSEM(i) = std(behParCue.numLicksBefRew(trInd))/sqrt(nTr);
        sessData.numLicksRewMean(i) = mean(behParCue.numLicksRew(trInd));
        sessData.numLicksRewSEM(i) = std(behParCue.numLicksRew(trInd))/sqrt(nTr);
        
        %% added by Yingxue on 4/19/2021
        sessData.lickSimMean(i) = lickOverDistSim.meanCue;
        sessData.lickSimSem(i) = lickOverDistSim.semCue;
        
        sessData.speedSimMean(i) = speedOverDistSim.meanCue;
        sessData.speedSimSem(i) = speedOverDistSim.semCue;
        
        sessData.lickEucMean(i) = lickOverDistSim.meanEucCue;
        sessData.lickEucSem(i) = lickOverDistSim.semEucCue;
        
        sessData.speedEucMean(i) = speedOverDistSim.meanEucCue;
        sessData.speedEucSem(i) = speedOverDistSim.semEucCue;
        %%
        
        indN = ~isnan(behParCue.med1stFiveLickDist(trInd));
        med1stFiveLickDist = behParCue.med1stFiveLickDist(trInd);
        sessData.med1stFiveLickDistMean(i) = mean(med1stFiveLickDist(indN));
        sessData.med1stFiveLickDistSEM(i) = std(med1stFiveLickDist(indN))/sqrt(sum(indN));
        
        indN = ~isnan(behParCue.medLickDist(trInd));
        medLickDist = behParCue.medLickDist(trInd);
        sessData.medLickDistMean(i) = mean(medLickDist(indN));
        sessData.medLickDistSEM(i) = std(medLickDist(indN))/sqrt(sum(indN));
        
        indN = ~isnan(behParCue.med1stFiveLickDistBefRew(trInd));
        med1stFiveLickDist = behParCue.med1stFiveLickDist(trInd);
        sessData.med1stFiveLickDistBefRewMean(i) = mean(med1stFiveLickDist(indN));
        sessData.med1stFiveLickDistBefRewSEM(i) = std(med1stFiveLickDist(indN))/sqrt(sum(indN));
        
        indN = ~isnan(behParCue.medLickDistBefRew(trInd));
        medLickDistBefRew = behParCue.medLickDistBefRew(trInd);
        sessData.medLickDistBefRewMean(i) = mean(medLickDistBefRew(indN));
        sessData.medLickDistBefRewSEM(i) = std(medLickDistBefRew(indN))/sqrt(sum(indN));
    end
    
    barPlot(1:nSess,sessData.numSamplesMean,sessData.numSamplesSEM,'Session','No. samples / Trial','Cue onset');
    print ('-painters', '-dpdf', [path recName '_numSamplesCue'], '-r600')
    barPlot(1:nSess,sessData.meanSpeedMean,sessData.meanSpeedSEM,'Session','Mean speed (mm/s)','Cue onset');
    print ('-painters', '-dpdf', [path recName '_meanSpeedCue'], '-r600')
    barPlot(1:nSess,sessData.meanSpeedRunMean,sessData.meanSpeedRunSEM,'Session','Mean speed run only (mm/s)','Cue onset');
    print ('-painters', '-dpdf', [path recName '_meanSpeedRunCue'], '-r600')
    barPlot(1:nSess,sessData.maxRunLenTMean,sessData.maxRunLenTSEM,'Session','Max run length (s)','Cue onset');
    print ('-painters', '-dpdf', [path recName '_maxRunLenTCue'], '-r600')
    barPlot(1:nSess,sessData.totRunLenTMean,sessData.totRunLenTSEM,'Session','Total run length (s)','Cue onset');
    print ('-painters', '-dpdf', [path recName '_totRunLenTCue'], '-r600')
    barPlot(1:nSess,sessData.numRunMean,sessData.numRunSEM,'Session','No. run segments','Cue onset');
    print ('-painters', '-dpdf', [path recName '_numRunCue'], '-r600')
    barPlot(1:nSess,sessData.meanAccMean,sessData.meanAccSEM,'Session','Mean acceleration (mm/s^2)','Cue onset');
    print ('-painters', '-dpdf', [path recName '_meanAccCue'], '-r600')
    barPlot(1:nSess,sessData.totStopLenTMean,sessData.totStopLenTSEM,'Session','Total stop time (s)','Cue onset');
    print ('-painters', '-dpdf', [path recName '_totStopLenTCue'], '-r600')
    barPlot(1:nSess,sessData.med1stFiveLickDistMean,sessData.med1stFiveLickDistSEM,'Session','Med. dist. first five licks (mm)','Cue onset');
    print ('-painters', '-dpdf', [path recName '_med1stFiveLickDistCue'], '-r600')
    barPlot(1:nSess,sessData.medLickDistMean,sessData.medLickDistSEM,'Session','Med. dist. licks (mm)','Cue onset');
    print ('-painters', '-dpdf', [path recName '_medLickDistCue'], '-r600')
    barPlot(1:nSess,sessData.speedSimMean,sessData.speedSimSem,'Session','speed over distance similarity','Cue onset');
    print ('-painters', '-dpdf', [path recName '_speedSimCue'], '-r600')    
    barPlot(1:nSess,sessData.lickSimMean,sessData.lickSimSem,'Session','lick over distance similarity','Cue onset');
    print ('-painters', '-dpdf', [path recName '_lickSimCue'], '-r600')
    barPlot(1:nSess,sessData.speedEucMean,sessData.speedEucSem,'Session','speed over distance euclidean','Cue onset');
    print ('-painters', '-dpdf', [path recName '_speedEucCue'], '-r600')    
    barPlot(1:nSess,sessData.lickEucMean,sessData.lickEucSem,'Session','lick over distance euclidean','Cue onset');
    print ('-painters', '-dpdf', [path recName '_lickEucCue'], '-r600')
end

%% added by Yingxue on 10/01/2020
function sessData = compBehCueOpt(path,recName,sessNo,mazeSess,onlyRun)
    nSess = length(sessNo);
    sessData.sessNo = sessNo;
    sessData.numSamplesMean = zeros(2,nSess);
    sessData.numSamplesSEM = zeros(2,nSess);
    sessData.maxSpeedMean = zeros(2,nSess);
    sessData.maxSpeedSEM = zeros(2,nSess);
    sessData.meanSpeedMean = zeros(2,nSess);
    sessData.meanSpeedSEM = zeros(2,nSess);
    sessData.maxSpeedRunMean = zeros(2,nSess);
    sessData.maxSpeedRunSEM = zeros(2,nSess);
    sessData.meanSpeedRunMean = zeros(2,nSess);
    sessData.meanSpeedRunSEM = zeros(2,nSess);
    sessData.maxRunLenTMean = zeros(2,nSess);
    sessData.maxRunLenTSEM = zeros(2,nSess);
    sessData.totRunLenTMean = zeros(2,nSess);
    sessData.totRunLenTSEM = zeros(2,nSess);
    sessData.numRunMean = zeros(2,nSess);
    sessData.numRunSEM = zeros(2,nSess);
    sessData.maxAccMean = zeros(2,nSess);
    sessData.maxAccSEM = zeros(2,nSess);
    sessData.meanAccMean = zeros(2,nSess);
    sessData.meanAccSEM = zeros(2,nSess);
    sessData.totStopLenTMean = zeros(2,nSess);
    sessData.totStopLenTSEM = zeros(2,nSess);
    sessData.numLicksBefRewMean = zeros(2,nSess);
    sessData.numLicksBefRewSEM = zeros(2,nSess);
    sessData.numLicksRewMean = zeros(2,nSess);
    sessData.numLicksRewSEM = zeros(2,nSess);
    sessData.med1stFiveLickDistMean = zeros(2,nSess);
    sessData.med1stFiveLickDistSEM = zeros(2,nSess);
    sessData.medLickDistMean = zeros(2,nSess);
    sessData.medLickDistSEM = zeros(2,nSess); 
    sessData.med1stFiveLickDistBefRewMean = zeros(2,nSess);
    sessData.med1stFiveLickDistBefRewSEM = zeros(2,nSess);
    sessData.medLickDistBefRewMean = zeros(2,nSess);
    sessData.medLickDistBefRewSEM = zeros(2,nSess); 
    sessData.meanRunSim = zeros(2,nSess);
    sessData.meanRewSim = zeros(2,nSess);
    sessData.meanCueSim = zeros(2,nSess);
    
    sessData.speedSimMean = zeros(2,nSess);
    sessData.speedSimSem = zeros(2,nSess);
    sessData.lickSimMean = zeros(2,nSess);
    sessData.lickSimSem = zeros(2,nSess);
    
    sessData.speedEucMean = zeros(2,nSess);
    sessData.speedEucSem = zeros(2,nSess);
    sessData.lickEucMean = zeros(2,nSess);
    sessData.lickEucSem = zeros(2,nSess);
    
    for i = 1:nSess
        fileName = [recName '-0' num2str(sessNo(i))];
        pathTmp = [path fileName '\'];
        
        fullpath = [pathTmp fileName '_behPar_msess' num2str(mazeSess(i)) '.mat'];
        if(~exist(fullpath,'file'))
            disp([fullpath ' does not exist.']);
            return;
        end
        load(fullpath,'behParCue','param');
        
        fullPath = [pathTmp fileName '_lickDist_msess' num2str(mazeSess(i)) '.mat']; 
        if(exist(fullPath) == 0)
            disp('The _lickDist file does not exist');
            return;
        end
        load(fullPath,'lickOverDistOpt','lickOverDistOptCtrl');
        
        fullPath = [pathTmp fileName '_lickDistSim_msess' num2str(mazeSess(i)) '.mat'];
        if(exist(fullPath) == 0)
            disp('The lick over distance similarity file does not exist');
            return;
        end
        load(fullPath,'lickOverDistSim','lickOverDistSimOpt','lickOverDistSimOptCtrl');
        
        fullPath = [pathTmp fileName '_speedDistSim_msess' num2str(mazeSess(i)) '_Run' num2str(onlyRun) '.mat'];
        if(exist(fullPath) == 0)
            disp('The speed over distance similarity file does not exist');
            return;
        end
        load(fullPath,'speedOverDistSim','speedOverDistSimOpt','speedOverDistSimOptCtrl');
        
        if(~isempty(lickOverDistOpt.trInd))
            trIndOpt = lickOverDistOpt.trInd;
            trIndCtrl = lickOverDistOptCtrl.trInd;
            nTrOpt = length(trIndOpt);
            nTrCtrl = length(trIndCtrl);
            
            sessData.numSamplesMean(:,i) = [mean(behParCue.numSamples(trIndOpt));...
                mean(behParCue.numSamples(trIndCtrl))];
            sessData.numSamplesSEM(:,i) = [std(behParCue.numSamples(trIndOpt))/sqrt(nTrOpt);...
                std(behParCue.numSamples(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.maxSpeedMean(:,i) = [mean(behParCue.maxSpeed(trIndOpt));...
                mean(behParCue.maxSpeed(trIndCtrl))];
            sessData.maxSpeedSEM(:,i) = [std(behParCue.maxSpeed(trIndOpt))/sqrt(nTrOpt);...
                std(behParCue.maxSpeed(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.meanSpeedMean(:,i) = [mean(behParCue.meanSpeed(trIndOpt));...
                mean(behParCue.meanSpeed(trIndCtrl))];
            sessData.meanSpeedSEM(:,i) = [std(behParCue.meanSpeed(trIndOpt))/sqrt(nTrOpt);...
                std(behParCue.meanSpeed(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.maxSpeedRunMean(:,i) = [mean(behParCue.maxSpeedRun(trIndOpt));...
                mean(behParCue.maxSpeedRun(trIndCtrl))];
            sessData.maxSpeedRunSEM(:,i) = [std(behParCue.maxSpeedRun(trIndOpt))/sqrt(nTrOpt);...
                std(behParCue.maxSpeedRun(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.meanSpeedRunMean(:,i) = [mean(behParCue.meanSpeedRun(trIndOpt));...
                mean(behParCue.meanSpeedRun(trIndCtrl))];
            sessData.meanSpeedRunSEM(:,i) = [std(behParCue.meanSpeedRun(trIndOpt))/sqrt(nTrOpt);...
                std(behParCue.meanSpeedRun(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.maxRunLenTMean(:,i) = [mean(behParCue.maxRunLenT(trIndOpt));...
                mean(behParCue.maxRunLenT(trIndCtrl))];
            sessData.maxRunLenTSEM(:,i) = [std(behParCue.maxRunLenT(trIndOpt))/sqrt(nTrOpt);...
                std(behParCue.maxRunLenT(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.totRunLenTMean(:,i) = [mean(behParCue.totRunLenT(trIndOpt));...
                mean(behParCue.totRunLenT(trIndCtrl))];
            sessData.totRunLenTSEM(:,i) = [std(behParCue.totRunLenT(trIndOpt))/sqrt(nTrOpt);...
                std(behParCue.totRunLenT(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.numRunMean(:,i) = [mean(behParCue.numRun(trIndOpt));...
                mean(behParCue.numRun(trIndCtrl))];
            sessData.numRunSEM(:,i) = [std(behParCue.numRun(trIndOpt))/sqrt(nTrOpt);...
                std(behParCue.numRun(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.maxAccMean(:,i) = [mean(behParCue.maxAcc(trIndOpt));...
                mean(behParCue.maxAcc(trIndCtrl))];
            sessData.maxAccSEM(:,i) = [std(behParCue.maxAcc(trIndOpt))/sqrt(nTrOpt);...
                std(behParCue.maxAcc(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.meanAccMean(:,i) = [mean(behParCue.meanAcc(trIndOpt));...
                mean(behParCue.meanAcc(trIndCtrl))];
            sessData.meanAccSEM(:,i) = [std(behParCue.meanAcc(trIndOpt))/sqrt(nTrOpt);...
                std(behParCue.meanAcc(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.totStopLenTMean(:,i) = [mean(behParCue.totStopLenT(trIndOpt));...
                mean(behParCue.totStopLenT(trIndCtrl))];
            sessData.totStopLenTSEM(:,i) = [std(behParCue.totStopLenT(trIndOpt))/sqrt(nTrOpt);...
                std(behParCue.totStopLenT(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.numLicksBefRewMean(:,i) = [mean(behParCue.numLicksBefRew(trIndOpt));...
                mean(behParCue.numLicksBefRew(trIndCtrl))];
            sessData.numLicksBefRewSEM(:,i) = [std(behParCue.numLicksBefRew(trIndOpt))/sqrt(nTrOpt);...
                std(behParCue.numLicksBefRew(trIndCtrl))/sqrt(nTrCtrl)];
            sessData.numLicksRewMean(:,i) = [mean(behParCue.numLicksRew(trIndOpt));...
                mean(behParCue.numLicksRew(trIndCtrl))];
            sessData.numLicksRewSEM(:,i) = [std(behParCue.numLicksRew(trIndOpt))/sqrt(nTrOpt);...
                std(behParCue.numLicksRew(trIndCtrl))/sqrt(nTrCtrl)];
            
            %% added by Yingxue 0n 4/19/2021
            sessData.lickSimMean(:,i) = [lickOverDistSimOpt.meanCue;lickOverDistSimOptCtrl.meanCue];
            sessData.lickSimSem(:,i) = [lickOverDistSimOpt.semCue;lickOverDistSimOptCtrl.semCue];

            sessData.speedSimMean(:,i) = [speedOverDistSimOpt.meanCue;speedOverDistSimOptCtrl.meanCue];
            sessData.speedSimSem(:,i) = [speedOverDistSimOpt.semCue;speedOverDistSimOptCtrl.semCue];
            
            sessData.lickEucMean(:,i) = [lickOverDistSimOpt.meanEucCue;lickOverDistSimOptCtrl.meanEucCue];
            sessData.lickEucSem(:,i) = [lickOverDistSimOpt.semEucCue;lickOverDistSimOptCtrl.semEucCue];

            sessData.speedEucMean(:,i) = [speedOverDistSimOpt.meanEucCue;speedOverDistSimOptCtrl.meanEucCue];
            sessData.speedEucSem(:,i) = [speedOverDistSimOpt.semEucCue;speedOverDistSimOptCtrl.semEucCue];
            %%
            
            indN = ~isnan(behParCue.med1stFiveLickDist(trIndOpt));
            med1stFiveLickDist = behParCue.med1stFiveLickDist(trIndOpt);
            indNCtrl = ~isnan(behParCue.med1stFiveLickDist(trIndCtrl));
            med1stFiveLickDistCtrl = behParCue.med1stFiveLickDist(trIndCtrl);
            sessData.med1stFiveLickDistMean(:,i) = [mean(med1stFiveLickDist(indN));...
                mean(med1stFiveLickDistCtrl(indNCtrl))];
            sessData.med1stFiveLickDistSEM(:,i) = [std(med1stFiveLickDist(indN))/sqrt(sum(indN));...
                std(med1stFiveLickDistCtrl(indNCtrl))/sqrt(sum(indNCtrl))];

            indN = ~isnan(behParCue.medLickDist(trIndOpt));
            medLickDist = behParCue.medLickDist(trIndOpt);
            indNCtrl = ~isnan(behParCue.medLickDist(trIndCtrl));
            medLickDistCtrl = behParCue.medLickDist(trIndCtrl);
            sessData.medLickDistMean(:,i) = [mean(medLickDist(indN));...
                mean(medLickDistCtrl(indNCtrl))];
            sessData.medLickDistSEM(:,i) = [std(medLickDist(indN))/sqrt(sum(indN));...
                std(medLickDistCtrl(indNCtrl))/sqrt(sum(indNCtrl))];

            indN = ~isnan(behParCue.med1stFiveLickDistBefRew(trIndOpt));
            med1stFiveLickDist = behParCue.med1stFiveLickDist(trIndOpt);
            indNCtrl = ~isnan(behParCue.med1stFiveLickDistBefRew(trIndCtrl));
            med1stFiveLickDistCtrl = behParCue.med1stFiveLickDist(trIndCtrl);
            sessData.med1stFiveLickDistBefRewMean(:,i) = [mean(med1stFiveLickDist(indN));...
                mean(med1stFiveLickDistCtrl(indNCtrl))];
            sessData.med1stFiveLickDistBefRewSEM(:,i) = [std(med1stFiveLickDist(indN))/sqrt(sum(indN));...
                std(med1stFiveLickDistCtrl(indNCtrl))/sqrt(sum(indNCtrl))];

            indN = ~isnan(behParCue.medLickDistBefRew(trIndOpt));
            medLickDistBefRew = behParCue.medLickDistBefRew(trIndOpt);
            indNCtrl = ~isnan(behParCue.medLickDistBefRew(trIndCtrl));
            medLickDistBefRewCtrl = behParCue.medLickDistBefRew(trIndCtrl);
            sessData.medLickDistBefRewMean(:,i) = [mean(medLickDistBefRew(indN));...
                mean(medLickDistBefRewCtrl(indNCtrl))];
            sessData.medLickDistBefRewSEM(:,i) = [std(medLickDistBefRew(indN))/sqrt(sum(indN));...
                std(medLickDistBefRewCtrl(indNCtrl))/sqrt(sum(indNCtrl))];
        else
            trInd = param.startTr:param.endTr;
            nTr = length(trInd);

            sessData.numSamplesMean(:,i) = mean(behParCue.numSamples(trInd));
            sessData.numSamplesSEM(:,i) = std(behParCue.numSamples(trInd))/sqrt(nTr);
            sessData.maxSpeedMean(:,i) = mean(behParCue.maxSpeed(trInd));
            sessData.maxSpeedSEM(:,i) = std(behParCue.maxSpeed(trInd))/sqrt(nTr);
            sessData.meanSpeedMean(:,i) = mean(behParCue.meanSpeed(trInd));
            sessData.meanSpeedSEM(:,i) = std(behParCue.meanSpeed(trInd))/sqrt(nTr);
            sessData.maxSpeedRunMean(:,i) = mean(behParCue.maxSpeedRun(trInd));
            sessData.maxSpeedRunSEM(:,i) = std(behParCue.maxSpeedRun(trInd))/sqrt(nTr);
            sessData.meanSpeedRunMean(:,i) = mean(behParCue.meanSpeedRun(trInd));
            sessData.meanSpeedRunSEM(:,i) = std(behParCue.meanSpeedRun(trInd))/sqrt(nTr);
            sessData.maxRunLenTMean(:,i) = mean(behParCue.maxRunLenT(trInd));
            sessData.maxRunLenTSEM(:,i) = std(behParCue.maxRunLenT(trInd))/sqrt(nTr);
            sessData.totRunLenTMean(:,i) = mean(behParCue.totRunLenT(trInd));
            sessData.totRunLenTSEM(:,i) = std(behParCue.totRunLenT(trInd))/sqrt(nTr);
            sessData.numRunMean(:,i) = mean(behParCue.numRun(trInd));
            sessData.numRunSEM(:,i) = std(behParCue.numRun(trInd))/sqrt(nTr);
            sessData.maxAccMean(:,i) = mean(behParCue.maxAcc(trInd));
            sessData.maxAccSEM(:,i) = std(behParCue.maxAcc(trInd))/sqrt(nTr);
            sessData.meanAccMean(:,i) = mean(behParCue.meanAcc(trInd));
            sessData.meanAccSEM(:,i) = std(behParCue.meanAcc(trInd))/sqrt(nTr);
            sessData.totStopLenTMean(:,i) = mean(behParCue.totStopLenT(trInd));
            sessData.totStopLenTSEM(:,i) = std(behParCue.totStopLenT(trInd))/sqrt(nTr);
            sessData.numLicksBefRewMean(:,i) = mean(behParCue.numLicksBefRew(trInd));
            sessData.numLicksBefRewSEM(:,i) = std(behParCue.numLicksBefRew(trInd))/sqrt(nTr);
            sessData.numLicksRewMean(:,i) = mean(behParCue.numLicksRew(trInd));
            sessData.numLicksRewSEM(:,i) = std(behParCue.numLicksRew(trInd))/sqrt(nTr);
            
            sessData.lickSimMean(:,i) = lickOverDistSim.meanCue;
            sessData.lickSimSem(:,i) = lickOverDistSim.semCue;

            sessData.speedSimMean(:,i) = speedOverDistSim.meanCue;
            sessData.speedSimSem(:,i) = speedOverDistSim.semCue;
            
            %% added by Yingxue 0n 4/19/2021
            sessData.lickEucMean(:,i) = lickOverDistSim.meanEucCue;
            sessData.lickEucSem(:,i) = lickOverDistSim.semEucCue;

            sessData.speedEucMean(:,i) = speedOverDistSim.meanEucCue;
            sessData.speedEucSem(:,i) = speedOverDistSim.semEucCue;
            %%

            indN = ~isnan(behParCue.med1stFiveLickDist(trInd));
            med1stFiveLickDist = behParCue.med1stFiveLickDist(trInd);
            sessData.med1stFiveLickDistMean(:,i) = mean(med1stFiveLickDist(indN));
            sessData.med1stFiveLickDistSEM(:,i) = std(med1stFiveLickDist(indN))/sqrt(sum(indN));

            indN = ~isnan(behParCue.medLickDist(trInd));
            medLickDist = behParCue.medLickDist(trInd);
            sessData.medLickDistMean(:,i) = mean(medLickDist(indN));
            sessData.medLickDistSEM(:,i) = std(medLickDist(indN))/sqrt(sum(indN));

            indN = ~isnan(behParCue.med1stFiveLickDistBefRew(trInd));
            med1stFiveLickDist = behParCue.med1stFiveLickDist(trInd);
            sessData.med1stFiveLickDistBefRewMean(:,i) = mean(med1stFiveLickDist(indN));
            sessData.med1stFiveLickDistBefRewSEM(:,i) = std(med1stFiveLickDist(indN))/sqrt(sum(indN));

            indN = ~isnan(behParCue.medLickDistBefRew(trInd));
            medLickDistBefRew = behParCue.medLickDistBefRew(trInd);
            sessData.medLickDistBefRewMean(:,i) = mean(medLickDistBefRew(indN));
            sessData.medLickDistBefRewSEM(:,i) = std(medLickDistBefRew(indN))/sqrt(sum(indN));
        end
    end
    
    %% opto
    barPlot(1:nSess,sessData.numSamplesMean(1,:),sessData.numSamplesSEM(1,:),'Session','No. samples / Trial','Cue onset Opt');
    print ('-painters', '-dpdf', [path recName '_numSamplesCueOpt'], '-r600')
    barPlot(1:nSess,sessData.meanSpeedMean(1,:),sessData.meanSpeedSEM(1,:),'Session','Mean speed (mm/s)','Cue onset Opt');
    print ('-painters', '-dpdf', [path recName '_meanSpeedCueOpt'], '-r600')
    barPlot(1:nSess,sessData.meanSpeedRunMean(1,:),sessData.meanSpeedRunSEM(1,:),'Session','Mean speed run only (mm/s)','Cue onset Opt');
    print ('-painters', '-dpdf', [path recName '_meanSpeedRunCueOpt'], '-r600')
    barPlot(1:nSess,sessData.maxRunLenTMean(1,:),sessData.maxRunLenTSEM(1,:),'Session','Max run length (s)','Cue onset Opt');
    print ('-painters', '-dpdf', [path recName '_maxRunLenTCueOpt'], '-r600')
    barPlot(1:nSess,sessData.totRunLenTMean(1,:),sessData.totRunLenTSEM(1,:),'Session','Total run length (s)','Cue onset Opt');
    print ('-painters', '-dpdf', [path recName '_totRunLenTCueOpt'], '-r600')
    barPlot(1:nSess,sessData.numRunMean(1,:),sessData.numRunSEM(1,:),'Session','No. run segments','Cue onset Opt');
    print ('-painters', '-dpdf', [path recName '_numRunCueOpt'], '-r600')
    barPlot(1:nSess,sessData.meanAccMean(1,:),sessData.meanAccSEM(1,:),'Session','Mean acceleration (mm/s^2)','Cue onset Opt');
    print ('-painters', '-dpdf', [path recName '_meanAccCueOpt'], '-r600')
    barPlot(1:nSess,sessData.totStopLenTMean(1,:),sessData.totStopLenTSEM(1,:),'Session','Total stop time (s)','Cue onset Opt');
    print ('-painters', '-dpdf', [path recName '_totStopLenTCueOpt'], '-r600')
    barPlot(1:nSess,sessData.med1stFiveLickDistMean(1,:),sessData.med1stFiveLickDistSEM(1,:),'Session','Med. dist. first five licks (mm)','Cue onset Opt');
    print ('-painters', '-dpdf', [path recName '_med1stFiveLickDistCueOpt'], '-r600')
    barPlot(1:nSess,sessData.medLickDistMean(1,:),sessData.medLickDistSEM(1,:),'Session','Med. dist. licks (mm)','Cue onset Opt');
    print ('-painters', '-dpdf', [path recName '_medLickDistCueOpt'], '-r600')
    barPlot(1:nSess,sessData.speedSimMean(1,:),sessData.speedSimSem(1,:),'Session','speed over distance similarity','Cue onset Opt');
    print ('-painters', '-dpdf', [path recName '_speedSimCueOpt'], '-r600')
    barPlot(1:nSess,sessData.lickSimMean(1,:),sessData.lickSimSem(1,:),'Session','lick over distance similarity','Cue onset Opt');
    print ('-painters', '-dpdf', [path recName '_lickSimCueOpt'], '-r600')
    barPlot(1:nSess,sessData.speedEucMean(1,:),sessData.speedEucSem(1,:),'Session','speed over distance euclidean','Cue onset Opt');
    print ('-painters', '-dpdf', [path recName '_speedEucCueOpt'], '-r600')
    barPlot(1:nSess,sessData.lickEucMean(1,:),sessData.lickEucSem(1,:),'Session','lick over distance euclidean','Cue onset Opt');
    print ('-painters', '-dpdf', [path recName '_lickEucCueOpt'], '-r600')
    
    %% opto ctrl
    barPlot(1:nSess,sessData.numSamplesMean(2,:),sessData.numSamplesSEM(2,:),'Session','No. samples / Trial','Cue onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_numSamplesCueOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.meanSpeedMean(2,:),sessData.meanSpeedSEM(2,:),'Session','Mean speed (mm/s)','Cue onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_meanSpeedCueOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.meanSpeedRunMean(2,:),sessData.meanSpeedRunSEM(2,:),'Session','Mean speed run only (mm/s)','Cue onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_meanSpeedRunCueOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.maxRunLenTMean(2,:),sessData.maxRunLenTSEM(2,:),'Session','Max run length (s)','Cue onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_maxRunLenTCueOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.totRunLenTMean(2,:),sessData.totRunLenTSEM(2,:),'Session','Total run length (s)','Cue onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_totRunLenTCueOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.numRunMean(2,:),sessData.numRunSEM(2,:),'Session','No. run segments','Cue onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_numRunCueOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.meanAccMean(2,:),sessData.meanAccSEM(2,:),'Session','Mean acceleration (mm/s^2)','Cue onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_meanAccCueOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.totStopLenTMean(2,:),sessData.totStopLenTSEM(2,:),'Session','Total stop time (s)','Cue onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_totStopLenTCueOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.med1stFiveLickDistMean(2,:),sessData.med1stFiveLickDistSEM(2,:),'Session','Med. dist. first five licks (mm)','Cue onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_med1stFiveLickDistCueOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.medLickDistMean(2,:),sessData.medLickDistSEM(2,:),'Session','Med. dist. licks (mm)','Cue onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_medLickDistCueOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.speedSimMean(2,:),sessData.speedSimSem(2,:),'Session','speed over distance similarity','Cue onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_speedSimCueOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.lickSimMean(2,:),sessData.lickSimSem(2,:),'Session','lick over distance similarity','Cue onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_lickSimCueOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.speedEucMean(2,:),sessData.speedEucSem(2,:),'Session','speed over distance euclidean','Cue onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_speedEucCueOptCtrl'], '-r600')
    barPlot(1:nSess,sessData.lickEucMean(2,:),sessData.lickEucSem(2,:),'Session','lick over distance euclidean','Cue onset Opt Ctrl');
    print ('-painters', '-dpdf', [path recName '_lickEucCueOptCtrl'], '-r600')
end

function sessData = compBehLick(path,recName,sessNo,mazeSess)
    nSess = length(sessNo);
    sessData.sessNo = sessNo;
    sessData.Run = [];
    sessData.meanRun = [];
    sessData.stdRun = [];
    sessData.SEMRun = [];
    sessData.Cue = [];
    sessData.meanCue = [];
    sessData.stdCue = [];
    sessData.SEMCue = [];
        
    sessData.meanRun30to100 = [];
    sessData.meanRun30to120 = [];
    sessData.meanRun100to150 = [];
    sessData.meanRun150to180 = [];
    sessData.meanRun180to210 = [];
    for i = 1:nSess
        fileName = [recName '-0' num2str(sessNo(i))];
        pathTmp = [path fileName '\'];
        
        fullpath = [pathTmp fileName '_lickDist_msess' num2str(mazeSess(i)) '.mat'];
        if(~exist(fullpath,'file'))
            disp([fullpath ' does not exist.']);
            return;
        end
        load(fullpath,'lickOverDist','param');
        trInd = lickOverDist.trInd; %param.startTr:param.endTr;
        
        sessData.spaceSteps = param.spaceSteps;
        sessData.Run{i} = lickOverDist.Run(trInd,:);
        sessData.meanRun{i} = lickOverDist.meanRun;
        sessData.stdRun{i} = lickOverDist.stdRun;
        sessData.SEMRun{i} = lickOverDist.SEMRun;
        
        sessData.Cue{i} = lickOverDist.Cue(trInd,:);
        sessData.meanCue{i} = lickOverDist.meanCue;
        sessData.stdCue{i} = lickOverDist.stdCue;
        sessData.SEMCue{i} = lickOverDist.SEMCue;
        
        sessData.meanRun30to100{i} = lickOverDist.meanRun30to100;
        sessData.meanRun30to120{i} = lickOverDist.meanRun30to120;
        sessData.meanRun100to150{i} = lickOverDist.meanRun100to150;
        sessData.meanRun150to180{i} = lickOverDist.meanRun150to180;
        sessData.meanRun180to210{i} = lickOverDist.meanRun180to210;
                
%         %% added by Yingxue on 9/12/2020
%         Run30to100{i} = sessData.Run{i}(:,lickOverDist.ind30to100);
%         Run100to150{i} = sessData.Run{i}(:,lickOverDist.ind100to150);
%         Run150to180{i} = sessData.Run{i}(:,lickOverDist.ind150to180);
%         Run180to210{i} = sessData.Run{i}(:,lickOverDist.ind180to210);
    end
    
    lickPlot(sessData.Run,param.spaceSteps/10,'Dist (cm)','No. licks','Run onset');
    saveas(gcf,[path recName '_LickDistRun']);
    print ('-painters', '-dpdf', [path recName '_LickDistRun'], '-r600')    
    
    plotBoxPlot(sessData.meanRun30to100,'Lick 30-100 cm',[recName '_LickDistRun30to100'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun30to120,'Lick 30-120 cm',[recName '_LickDistRun30to120'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun100to150,'Lick 100-150 cm',[recName '_LickDistRun100to150'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun150to180,'Lick 150-180 cm',[recName '_LickDistRun150to180'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun180to210,'Lick 180-210 cm',[recName '_LickDistRun180to210'],path,[-1 10]);
    
    lickPlot(sessData.Cue,param.spaceSteps/10,'Dist (cm)','No. licks','Cue onset');
    saveas(gcf,[path recName '_LickDistCue']);
    print ('-painters', '-dpdf', [path recName '_LickDistCue'], '-r600')    
end

%% added by Yingxue on 10/01/2020
function sessData = compBehLickOpt(path,recName,sessNo,mazeSess)
    nSess = length(sessNo);
    sessData.sessNo = sessNo;
    
    for i = 1:nSess
        fileName = [recName '-0' num2str(sessNo(i))];
        pathTmp = [path fileName '\'];
        
        fullpath = [pathTmp fileName '_lickDist_msess' num2str(mazeSess(i)) '.mat'];
        if(~exist(fullpath,'file'))
            disp([fullpath ' does not exist.']);
            return;
        end
        load(fullpath,'lickOverDist','lickOverDistOpt','lickOverDistOptCtrl',...
            'lickOverDistOptCtrl1','param');
                
        if(~isempty(lickOverDistOpt.trInd))
            trIndOpt = lickOverDistOpt.trInd;
            trIndCtrl = lickOverDistOptCtrl.trInd;
            trIndCtrl1 = lickOverDistOptCtrl1.trInd; % added by Yingxue on 1/25/2022

            sessData.spaceSteps = param.spaceSteps;
            sessData.Run(:,i) = [{lickOverDist.Run(trIndOpt,:)};...
                {lickOverDist.Run(trIndCtrl,:)};...
                {lickOverDist.Run(trIndCtrl1,:)}];
            sessData.meanRun(:,i) = [{lickOverDistOpt.meanRun};...
                {lickOverDistOptCtrl.meanRun};...
                {lickOverDistOptCtrl1.meanRun}];
            sessData.stdRun(:,i) = [{lickOverDistOpt.stdRun};...
                {lickOverDistOptCtrl.stdRun};...
                {lickOverDistOptCtrl1.stdRun}];
            sessData.SEMRun(:,i) = [{lickOverDistOpt.SEMRun};...
                {lickOverDistOptCtrl.SEMRun};...
                {lickOverDistOptCtrl1.SEMRun}];

            sessData.meanRun30to100(:,i) = [{lickOverDistOpt.meanRun30to100};...
                {lickOverDistOptCtrl.meanRun30to100};...
                {lickOverDistOptCtrl1.meanRun30to100}];
            sessData.meanRun30to120(:,i) = [{lickOverDistOpt.meanRun30to120};...
                {lickOverDistOptCtrl.meanRun30to120};...
                {lickOverDistOptCtrl1.meanRun30to120}];
            sessData.meanRun100to150(:,i) = [{lickOverDistOpt.meanRun100to150};...
                {lickOverDistOptCtrl.meanRun100to150};...
                {lickOverDistOptCtrl1.meanRun100to150}];
            sessData.meanRun150to180(:,i) = [{lickOverDistOpt.meanRun150to180};...
                {lickOverDistOptCtrl.meanRun150to180};...
                {lickOverDistOptCtrl1.meanRun150to180}];
            sessData.meanRun180to210(:,i) = [{lickOverDistOpt.meanRun180to210};
                {lickOverDistOptCtrl.meanRun180to210};...
                {lickOverDistOptCtrl1.meanRun180to210}];
            
%         %% added by Yingxue on 9/12/2020
%         Run30to100{i} = sessData.Run{i}(:,lickOverDist.ind30to100);
%         Run100to150{i} = sessData.Run{i}(:,lickOverDist.ind100to150);
%         Run150to180{i} = sessData.Run{i}(:,lickOverDist.ind150to180);
%         Run180to210{i} = sessData.Run{i}(:,lickOverDist.ind180to210);
        else
            trInd = param.startTr:param.endTr;

            sessData.spaceSteps = param.spaceSteps;
            sessData.Run(:,i) = [{lickOverDist.Run(trInd,:)};...
                {lickOverDist.Run(trInd,:)};...
                {lickOverDist.Run(trInd,:)}];
            sessData.meanRun(:,i) = [{lickOverDist.meanRun};...
                {lickOverDist.meanRun};...
                {lickOverDist.meanRun}];
            sessData.stdRun(:,i) = [{lickOverDist.stdRun};...
                {lickOverDist.stdRun};...
                {lickOverDist.stdRun}];
            sessData.SEMRun(:,i) = [{lickOverDist.SEMRun};...
                {lickOverDist.SEMRun};...
                {lickOverDist.SEMRun}];

            sessData.meanRun30to100(:,i) = [{lickOverDist.meanRun30to100};...
                {lickOverDist.meanRun30to100};...
                {lickOverDist.meanRun30to100}];
            sessData.meanRun30to120(:,i) = [{lickOverDist.meanRun30to120};...
                {lickOverDist.meanRun30to120};...
                {lickOverDist.meanRun30to120}];
            sessData.meanRun100to150(:,i) = [{lickOverDist.meanRun100to150};...
                {lickOverDist.meanRun100to150};...
                {lickOverDist.meanRun100to150}];
            sessData.meanRun150to180(:,i) = [{lickOverDist.meanRun150to180};...
                {lickOverDist.meanRun150to180};...
                {lickOverDist.meanRun150to180}];
            sessData.meanRun180to210(:,i) = [{lickOverDist.meanRun180to210};
                {lickOverDist.meanRun180to210};...
                {lickOverDist.meanRun180to210}];
        end
    end
    
    lickPlot(sessData.Run(1,:),param.spaceSteps/10,'Dist (cm)','No. licks','Run onset Opt');
    saveas(gcf,[path recName '_LickDistRunOpt']);
    print ('-painters', '-dpdf', [path recName '_LickDistRunOpt'], '-r600')    
    
    plotBoxPlot(sessData.meanRun30to100(1,:),'Lick 30-100 cm',[recName '_LickDistRun30to100Opt'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun30to120(1,:),'Lick 30-120 cm',[recName '_LickDistRun30to120Opt'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun100to150(1,:),'Lick 100-150 cm',[recName '_LickDistRun100to150Opt'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun150to180(1,:),'Lick 150-180 cm',[recName '_LickDistRun150to180Opt'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun180to210(1,:),'Lick 180-210 cm',[recName '_LickDistRun180to210Opt'],path,[-1 10]);
    
    lickPlot(sessData.Run(2,:),param.spaceSteps/10,'Dist (cm)','No. licks','Run onset Opt Ctrl');
    saveas(gcf,[path recName '_LickDistRunOptCtrl']);
    print ('-painters', '-dpdf', [path recName '_LickDistRunOptCtrl'], '-r600')    

    plotBoxPlot(sessData.meanRun30to100(2,:),'Lick 30-100 cm',[recName '_LickDistRun30to100OptCtrl'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun30to120(2,:),'Lick 30-120 cm',[recName '_LickDistRun30to120OptCtrl'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun100to150(2,:),'Lick 100-150 cm',[recName '_LickDistRun100to150OptCtrl'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun150to180(2,:),'Lick 150-180 cm',[recName '_LickDistRun150to180OptCtrl'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun180to210(2,:),'Lick 180-210 cm',[recName '_LickDistRun180to210OptCtrl'],path,[-1 10]);
    
    lickPlot(sessData.Run(3,:),param.spaceSteps/10,'Dist (cm)','No. licks','Run onset Ctrl after opt');
    saveas(gcf,[path recName '_LickDistRunOptCtrl1']);
    print ('-painters', '-dpdf', [path recName '_LickDistRunOptCtrl1'], '-r600')    

    plotBoxPlot(sessData.meanRun30to100(3,:),'Lick 30-100 cm',[recName '_LickDistRun30to100OptCtrl1'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun30to120(2,:),'Lick 30-120 cm',[recName '_LickDistRun30to120OptCtrl1'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun100to150(3,:),'Lick 100-150 cm',[recName '_LickDistRun100to150OptCtrl1'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun150to180(3,:),'Lick 150-180 cm',[recName '_LickDistRun150to180OptCtrl1'],path,[-1 10]);
    plotBoxPlot(sessData.meanRun180to210(3,:),'Lick 180-210 cm',[recName '_LickDistRun180to210OptCtrl1'],path,[-1 10]);
end

function plotBoxPlot(x,yl,fn,pathAnal,ylim)
    figure;
    colorArr = [...
                127 127 229;...
                125 125 168]/255;
    xa = [];
    g = [];
    for i = 1:length(x)
        xa = [xa;x{i}(:)];
        g = [g;repmat({num2str(i)},length(x{i}(:)),1)];
    end
    boxplot(xa,g,'Notch','on','Widths',0.3,'Symbol','');
    h = findobj(gca,'Tag','Box');
    for j = 1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colorArr(1,:),'FaceAlpha',0.5);
    end
    
    ylabel(yl);
    set(gca,'YLim',ylim)
    
    fileName1 = [pathAnal fn];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function lickPlot(data,x,xl,yl,ti)
    figure
    nSess = length(data);
    for i = 1:nSess
        options.handle     = subplot(1,nSess,i);
        options.color_area = [224 191 232]./255; %[128 193 219]./255;    % Blue theme
        options.color_line = [129 98 137]./255; %[ 52 148 186]./255;
        options.alpha      = 0.5;
        options.line_width = 2;
        options.error      = 'sem';
        options.x_axis = x;
        plot_areaerrorbarSub(data{i}, options);
        set(gca,'XLim',[0 max(x)],'YLim',[-2 10])
        xlabel(xl)
        ylabel(yl)
    end
end

function sessData = compBehSpeed(path,recName,sessNo,mazeSess,onlyRun)
    nSess = length(sessNo);
    sessData.sessNo = sessNo;
    sessData.Run = [];
    sessData.meanRun = [];
    sessData.stdRun = [];
    sessData.SEMRun = [];
    sessData.Cue = [];
    sessData.meanCue = [];
    sessData.stdCue = [];
    sessData.SEMCue = [];
    sessData.meanSpeedOverDistRun0to100 = [];
    sessData.meanSpeedOverDistRunAfter100 = [];
    sessData.meanSpeedOverTimeRunBL = [];
    sessData.meanSpeedOverTimeRunBefRun = [];
    sessData.meanSpeedOverTimeRun0to1 = [];
    sessData.meanSpeedOverTimeRun3to5 = [];
    
    for i = 1:nSess
        fileName = [recName '-0' num2str(sessNo(i))];
        pathTmp = [path fileName '\'];
        
        fullpath = [pathTmp fileName '_runSpeedDist_Run' num2str(onlyRun) '_msess' num2str(mazeSess(i)) '.mat'];
        if(~exist(fullpath,'file'))
            disp([fullpath ' does not exist.']);
            return;
        end
        load(fullpath,'speedOverDist','param');
        
        fullpath = [pathTmp fileName '_runSpeedTime_Run' num2str(onlyRun) '_msess' num2str(mazeSess(i)) '.mat'];
        if(~exist(fullpath,'file'))
            disp([fullpath ' does not exist.']);
            return;
        end
        load(fullpath,'speedOverTime');
        trInd = speedOverTime.trInd; %param.startTr:param.endTr;
       
        sessData.spaceSteps = param.spaceSteps;
        sessData.Run{i} = speedOverDist.Run(trInd,:)/10;
        sessData.meanRun{i} = speedOverDist.meanRun/10;
        sessData.stdRun{i} = speedOverDist.stdRun/10;
        sessData.SEMRun{i} = speedOverDist.SEMRun/10;
        
        %% added on 8/13/2022
        sessData.RunLowSpeed{i} = speedOverDist.lowSpeed(trInd,:);
        sessData.meanRunLowSpeed{i} = speedOverDist.meanRunLowSpeed;
        sessData.stdRunLowSpeed{i} = speedOverDist.stdRunLowSpeed;
        sessData.SEMRunLowSpeed{i} = speedOverDist.SEMRunLowSpeed;
        %%
        
        sessData.Cue{i} = speedOverDist.Cue(trInd,:)/10;
        sessData.meanCue{i} = speedOverDist.meanCue/10;
        sessData.stdCue{i} = speedOverDist.stdCue/10;
        sessData.SEMCue{i} = speedOverDist.SEMCue/10;
                
        sessData.meanSpeedOverDistRun0to100{i} = speedOverDist.meanRun30to100/10;
        sessData.meanSpeedOverDistRunAfter100{i} = speedOverDist.meanRunAfter100/10;
        
        sessData.meanSpeedOverTimeRunBL{i} = speedOverTime.meanRunSpeedBL/10;
        sessData.meanSpeedOverTimeRunBefRun{i} = speedOverTime.meanRunSpeedBefRun/10;
        sessData.meanSpeedOverTimeRun0to1{i} = speedOverTime.meanRunSpeed0to1/10;
        sessData.meanSpeedOverTimeRun3to5{i} = speedOverTime.meanRunSpeed3to5/10;
        
        sessData.meanSpeedOverTimeRew0to1{i} = speedOverTime.meanRewSpeed0to1/10;
        sessData.meanSpeedOverTimeRew1to2{i} = speedOverTime.meanRewSpeed1to2/10;
        sessData.meanSpeedOverTimeRew2to3{i} = speedOverTime.meanRewSpeed2to3/10;
        sessData.meanSpeedOverTimeRew3to5{i} = speedOverTime.meanRewSpeed3to5/10;
        
        %% added by Yingxue on 9/12/2020
        SpeedOverDistRun0to100{i} = sessData.Run{i}(:,speedOverDist.ind30to100);
        SpeedOverDistRunAfter100{i} = sessData.Run{i}(:,speedOverDist.indAfter100);
        
        SpeedOverTimeRunBL{i} = speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.indBaseline)/10;
        SpeedOverTimeRunBefRun{i} = speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.indBefRun)/10;
        SpeedOverTimeRun0to1{i} = speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.ind0to1)/10;
        SpeedOverTimeRun3to5{i} = speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.ind3to5)/10;
        
        SpeedOverTimeRew0to1{i} = speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind0to1Rew)/10;
        SpeedOverTimeRew1to2{i} = speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind1to2Rew)/10;
        SpeedOverTimeRew2to3{i} = speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind2to3Rew)/10;
        SpeedOverTimeRew3to5{i} = speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind3to5Rew)/10;
    end
    
    runPlot(sessData.Run,param.spaceSteps/10,'Dist (cm)','Speed (cm/s)','Run onset');
    saveas(gcf,[path recName '_SpeedDistRun']);
    print ('-painters', '-dpdf', [path recName '_SpeedDistRun'], '-r600')    
    runPlot(sessData.Cue,param.spaceSteps/10,'Dist (cm)','Speed (cm/s)','Cue onset');
    saveas(gcf,[path recName '_SpeedDistCue']);
    print ('-painters', '-dpdf', [path recName '_SpeedDistCue'], '-r600')    
    
    % changed by Yingxue 9/12/2020, using SpeedOverDistRun0to100 instead of
    % sessData.meanSpeedOverDistRun0to100 
    plotBoxPlot(sessData.meanSpeedOverDistRun0to100,'Speed (cm/s) RunOn 0-100cm',[recName '_SpeedDistRun0to100'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverDistRunAfter100,'Speed (cm/s) RunOn After100cm',[recName '_SpeedDistRunAfter100'],path,[-1 100]);
    
    plotBoxPlot(sessData.meanSpeedOverTimeRunBL,'Speed (cm/s) RunOn BL',[recName '_SpeedTimeRunBL'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverTimeRunBefRun,'Speed (cm/s) RunOn BefRun',[recName '_SpeedTimeRunBefRun'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverTimeRun0to1,'Speed (cm/s) RunOn 0-1s',[recName '_SpeedTimeRun0to1s'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverTimeRun3to5,'Speed (cm/s) RunOn 3-5s',[recName '_SpeedTimeRun3to5s'],path,[-1 100]);
    
    plotBoxPlot(sessData.meanSpeedOverTimeRew0to1,'Speed (cm/s) Rew 0-1s',[recName '_SpeedTimeRew0to1s'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverTimeRew1to2,'Speed (cm/s) Rew 1-2s',[recName '_SpeedTimeRew1to2s'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverTimeRew2to3,'Speed (cm/s) Rew 2-3s',[recName '_SpeedTimeRew2to3s'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverTimeRew3to5,'Speed (cm/s) Rew 3-5s',[recName '_SpeedTimeRew3to5s'],path,[-1 100]);
end

%% added by Yingxue on 10/01/2020
function sessData = compBehSpeedOpt(path,recName,sessNo,mazeSess,onlyRun)
    nSess = length(sessNo);
    sessData.sessNo = sessNo;
    
    for i = 1:nSess
        fileName = [recName '-0' num2str(sessNo(i))];
        pathTmp = [path fileName '\'];
        
        fullpath = [pathTmp fileName '_runSpeedDist_Run' num2str(onlyRun) '_msess' num2str(mazeSess(i)) '.mat'];
        if(~exist(fullpath,'file'))
            disp([fullpath ' does not exist.']);
            return;
        end
        load(fullpath,'speedOverDist','speedOverDistOpt','speedOverDistOptCtrl',...
            'speedOverDistOptCtrl1','param');
        
        fullpath = [pathTmp fileName '_runSpeedTime_Run' num2str(onlyRun) '_msess' num2str(mazeSess(i)) '.mat'];
        if(~exist(fullpath,'file'))
            disp([fullpath ' does not exist.']);
            return;
        end
        load(fullpath,'speedOverTime','speedOverTimeOpt','speedOverTimeOptCtrl',...
            'speedOverTimeOptCtrl1');
        
        if(~isempty(speedOverDistOpt.trInd))
            trIndOpt = speedOverDistOpt.trInd;
            trIndCtrl = speedOverDistOptCtrl.trInd;
            trIndCtrl1 = speedOverDistOptCtrl1.trInd; 
                % changed by Yingxue on 1/25/2022 to add the ctrl trial after opt
       
            sessData.spaceSteps = param.spaceSteps;
            sessData.Run(:,i) = [{speedOverDist.Run(trIndOpt,:)/10};...
                {speedOverDist.Run(trIndCtrl,:)/10};...
                {speedOverDist.Run(trIndCtrl1,:)/10}];
            sessData.meanRun(:,i) = [{speedOverDistOpt.meanRun/10};...
                {speedOverDistOptCtrl.meanRun/10};...
                {speedOverDistOptCtrl1.meanRun/10}];
            sessData.stdRun(:,i) = [{speedOverDistOpt.stdRun/10};...
                {speedOverDistOptCtrl.stdRun/10};...
                {speedOverDistOptCtrl1.stdRun/10}];
            sessData.SEMRun(:,i) = [{speedOverDistOpt.SEMRun/10};...
                {speedOverDistOptCtrl.SEMRun/10};...
                {speedOverDistOptCtrl1.SEMRun/10}];
            
            %% added on 8/13/2022
            sessData.RunLowSpeed(:,i) = [{speedOverDist.lowSpeed(trIndOpt,:)};...
                {speedOverDist.lowSpeed(trIndCtrl,:)};...
                {speedOverDist.lowSpeed(trIndCtrl1,:)}];
            sessData.meanRunLowSpeed(:,i) = [{speedOverDistOpt.meanRunLowSpeed};...
                {speedOverDistOptCtrl.meanRunLowSpeed};...
                {speedOverDistOptCtrl1.meanRunLowSpeed}];
            sessData.stdRunLowSpeed(:,i) = [{speedOverDistOpt.stdRunLowSpeed};...
                {speedOverDistOptCtrl.stdRunLowSpeed};...
                {speedOverDistOptCtrl1.stdRunLowSpeed}];
            sessData.SEMRunLowSpeed(:,i) = [{speedOverDistOpt.SEMRunLowSpeed};...
                {speedOverDistOptCtrl.SEMRunLowSpeed};...
                {speedOverDistOptCtrl1.SEMRunLowSpeed}];
            %%

            sessData.Cue(:,i) = [{speedOverDist.Cue(trIndOpt,:)/10};...
                {speedOverDist.Cue(trIndCtrl,:)/10};...
                {speedOverDist.Cue(trIndCtrl1,:)/10}];

            sessData.meanSpeedOverDistRun0to100(:,i) = [{speedOverDistOpt.meanRun30to100/10};...
                {speedOverDistOptCtrl.meanRun30to100/10};...
                {speedOverDistOptCtrl1.meanRun30to100/10}];
            sessData.meanSpeedOverDistRunAfter100(:,i) = [{speedOverDistOpt.meanRunAfter100/10};...
                {speedOverDistOptCtrl.meanRunAfter100/10};...
                {speedOverDistOptCtrl1.meanRunAfter100/10}];

            sessData.meanSpeedOverTimeRunBL(:,i) = [{speedOverTimeOpt.meanRunSpeedBL/10};...
                {speedOverTimeOptCtrl.meanRunSpeedBL/10};...
                {speedOverTimeOptCtrl1.meanRunSpeedBL/10}];
            sessData.meanSpeedOverTimeRunBefRun(:,i) = [{speedOverTimeOpt.meanRunSpeedBefRun/10};...
                {speedOverTimeOptCtrl.meanRunSpeedBefRun/10};...
                {speedOverTimeOptCtrl1.meanRunSpeedBefRun/10}];
            sessData.meanSpeedOverTimeRun0to1(:,i) = [{speedOverTimeOpt.meanRunSpeed0to1/10};...
                {speedOverTimeOptCtrl.meanRunSpeed0to1/10};...
                {speedOverTimeOptCtrl1.meanRunSpeed0to1/10}];
            sessData.meanSpeedOverTimeRun3to5(:,i) = [{speedOverTimeOpt.meanRunSpeed3to5/10};...
                {speedOverTimeOptCtrl.meanRunSpeed3to5/10};...
                {speedOverTimeOptCtrl1.meanRunSpeed3to5/10}];


            %% added by Yingxue on 9/12/2020
            SpeedOverDistRun0to100(:,i) = [{sessData.Run{1,i}(:,speedOverDist.ind30to100)};...
                {sessData.Run{2,i}(:,speedOverDist.ind30to100)};...
                {sessData.Run{3,i}(:,speedOverDist.ind30to100)}];
            SpeedOverDistRunAfter100(:,i) = [{sessData.Run{1,i}(:,speedOverDist.indAfter100)};...
                {sessData.Run{2,i}(:,speedOverDist.indAfter100)};...
                {sessData.Run{3,i}(:,speedOverDist.indAfter100)}];

            SpeedOverTimeRunBL(:,i) = [{speedOverTime.runSpeedAlignedRun(trIndOpt,speedOverTime.indBaseline)/10};...
                {speedOverTime.runSpeedAlignedRun(trIndCtrl,speedOverTime.indBaseline)/10};...
                {speedOverTime.runSpeedAlignedRun(trIndCtrl1,speedOverTime.indBaseline)/10}];
            SpeedOverTimeRunBefRun(:,i) = [{speedOverTime.runSpeedAlignedRun(trIndOpt,speedOverTime.indBefRun)/10};...
                {speedOverTime.runSpeedAlignedRun(trIndCtrl,speedOverTime.indBefRun)/10};...
                {speedOverTime.runSpeedAlignedRun(trIndCtrl1,speedOverTime.indBefRun)/10}];
            SpeedOverTimeRun0to1(:,i) = [{speedOverTime.runSpeedAlignedRun(trIndOpt,speedOverTime.ind0to1)/10};...
                {speedOverTime.runSpeedAlignedRun(trIndCtrl,speedOverTime.ind0to1)/10};...
                {speedOverTime.runSpeedAlignedRun(trIndCtrl1,speedOverTime.ind0to1)/10}];
            SpeedOverTimeRun3to5(:,i) = [{speedOverTime.runSpeedAlignedRun(trIndOpt,speedOverTime.ind3to5)/10};...
                {speedOverTime.runSpeedAlignedRun(trIndCtrl,speedOverTime.ind3to5)/10};...
                {speedOverTime.runSpeedAlignedRun(trIndCtrl1,speedOverTime.ind3to5)/10}];

            SpeedOverTimeRew0to1(:,i) = [{speedOverTime.runSpeedAlignedRew(trIndOpt,speedOverTime.ind0to1Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trIndCtrl,speedOverTime.ind0to1Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trIndCtrl1,speedOverTime.ind0to1Rew)/10}];
            SpeedOverTimeRew1to2(:,i) = [{speedOverTime.runSpeedAlignedRew(trIndOpt,speedOverTime.ind1to2Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trIndCtrl,speedOverTime.ind1to2Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trIndCtrl1,speedOverTime.ind1to2Rew)/10}];
            SpeedOverTimeRew2to3(:,i) = [{speedOverTime.runSpeedAlignedRew(trIndOpt,speedOverTime.ind2to3Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trIndCtrl,speedOverTime.ind2to3Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trIndCtrl1,speedOverTime.ind2to3Rew)/10}];
            SpeedOverTimeRew3to5(:,i) = [{speedOverTime.runSpeedAlignedRew(trIndOpt,speedOverTime.ind3to5Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trIndCtrl,speedOverTime.ind3to5Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trIndCtrl1,speedOverTime.ind3to5Rew)/10}];
        else
            trInd = param.startTr:param.endTr;

            sessData.spaceSteps = param.spaceSteps;
            sessData.Run(:,i) = [{speedOverDist.Run(trInd,:)/10};...
                {speedOverDist.Run(trInd,:)/10}; ...
                {speedOverDist.Run(trInd,:)/10}];
            sessData.meanRun(:,i) = [{speedOverDist.meanRun/10};...
                {speedOverDist.meanRun/10};...
                {speedOverDist.meanRun/10}];
            sessData.stdRun(:,i) = [{speedOverDist.stdRun/10};...
                {speedOverDist.stdRun/10};...
                {speedOverDist.stdRun/10}];
            sessData.SEMRun(:,i) = [{speedOverDist.SEMRun/10};...
                {speedOverDist.SEMRun/10};...
                {speedOverDist.SEMRun/10}];
            
            %% added on 8/13/2022
            sessData.RunLowSpeed(:,i) = [{speedOverDist.lowSpeed(trInd,:)};...
                {speedOverDist.lowSpeed(trInd,:)};...
                {speedOverDist.lowSpeed(trInd,:)}];
            sessData.meanRunLowSpeed(:,i) = [{speedOverDist.meanRunLowSpeed};...
                {speedOverDist.meanRunLowSpeed};...
                {speedOverDist.meanRunLowSpeed}];
            sessData.stdRunLowSpeed(:,i) = [{speedOverDist.stdRunLowSpeed};...
                {speedOverDist.stdRunLowSpeed};...
                {speedOverDist.stdRunLowSpeed}];
            sessData.SEMRunLowSpeed(:,i) = [{speedOverDist.SEMRunLowSpeed};...
                {speedOverDist.SEMRunLowSpeed};...
                {speedOverDist.SEMRunLowSpeed}];
            %%

            sessData.Cue(:,i) = [{speedOverDist.Cue(trInd,:)/10};...
                {speedOverDist.Cue(trInd,:)/10};...
                {speedOverDist.Cue(trInd,:)/10}];

            sessData.meanSpeedOverDistRun0to100(:,i) = [{speedOverDist.meanRun30to100/10};...
                {speedOverDist.meanRun30to100/10};...
                {speedOverDist.meanRun30to100/10}];
            sessData.meanSpeedOverDistRunAfter100(:,i) = [{speedOverDist.meanRunAfter100/10};...
                {speedOverDist.meanRunAfter100/10};...
                {speedOverDist.meanRunAfter100/10}];

            sessData.meanSpeedOverTimeRunBL(:,i) = [{speedOverTime.meanRunSpeedBL/10};...
                {speedOverTime.meanRunSpeedBL/10};...
                {speedOverTime.meanRunSpeedBL/10}];
            sessData.meanSpeedOverTimeRunBefRun(:,i) = [{speedOverTime.meanRunSpeedBefRun/10};...
                {speedOverTime.meanRunSpeedBefRun/10};...
                {speedOverTime.meanRunSpeedBefRun/10}];
            sessData.meanSpeedOverTimeRun0to1(:,i) = [{speedOverTime.meanRunSpeed0to1/10};...
                {speedOverTime.meanRunSpeed0to1/10};...
                {speedOverTime.meanRunSpeed0to1/10}];
            sessData.meanSpeedOverTimeRun3to5(:,i) = [{speedOverTime.meanRunSpeed3to5/10};...
                {speedOverTime.meanRunSpeed3to5/10};...
                {speedOverTime.meanRunSpeed3to5/10}];

            %% added by Yingxue on 9/12/2020
            SpeedOverDistRun0to100(:,i) = [{sessData.Run{1,i}(:,speedOverDist.ind30to100)};...
                {sessData.Run{2,i}(:,speedOverDist.ind30to100)};...
                {sessData.Run{3,i}(:,speedOverDist.ind30to100)}];
            SpeedOverDistRunAfter100(:,i) = [{sessData.Run{1,i}(:,speedOverDist.indAfter100)};...
                {sessData.Run{2,i}(:,speedOverDist.indAfter100)};...
                {sessData.Run{3,i}(:,speedOverDist.indAfter100)}];

            SpeedOverTimeRunBL(:,i) = [{speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.indBaseline)/10};...
                {speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.indBaseline)/10};...
                {speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.indBaseline)/10}];
            SpeedOverTimeRunBefRun(:,i) = [{speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.indBefRun)/10};...
                {speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.indBefRun)/10};...
                {speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.indBefRun)/10}];
            SpeedOverTimeRun0to1(:,i) = [{speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.ind0to1)/10};...
                {speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.ind0to1)/10};...
                {speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.ind0to1)/10}];
            SpeedOverTimeRun3to5(:,i) = [{speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.ind3to5)/10};...
                {speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.ind3to5)/10};...
                {speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.ind3to5)/10}];

            SpeedOverTimeRew0to1(:,i) = [{speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind0to1Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind0to1Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind0to1Rew)/10}];
            SpeedOverTimeRew1to2(:,i) = [{speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind1to2Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind1to2Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind1to2Rew)/10}];
            SpeedOverTimeRew2to3(:,i) = [{speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind2to3Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind2to3Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind2to3Rew)/10}];
            SpeedOverTimeRew3to5(:,i) = [{speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind3to5Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind3to5Rew)/10};...
                {speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind3to5Rew)/10}];
        end
    end
    
    runPlot(sessData.Run(1,:),param.spaceSteps/10,'Dist (cm)','Speed (cm/s)','Run onset');
    saveas(gcf,[path recName '_SpeedDistRunOpt']);
    print ('-painters', '-dpdf', [path recName '_SpeedDistRunOpt'], '-r600')    
    runPlot(sessData.Cue(1,:),param.spaceSteps/10,'Dist (cm)','Speed (cm/s)','Cue onset');
    saveas(gcf,[path recName '_SpeedDistCueOpt']);
    print ('-painters', '-dpdf', [path recName '_SpeedDistCueOpt'], '-r600')    
    
    % changed by Yingxue 9/12/2020, using SpeedOverDistRun0to100 instead of
    % sessData.meanSpeedOverDistRun0to100 
    plotBoxPlot(sessData.meanSpeedOverDistRun0to100(1,:),'Speed (cm/s) RunOn 0-100cm',[recName '_SpeedDistRun0to100Opt'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverDistRunAfter100(1,:),'Speed (cm/s) RunOn After100cm',[recName '_SpeedDistRunAfter100Opt'],path,[-1 100]);
    
    plotBoxPlot(sessData.meanSpeedOverTimeRunBL(1,:),'Speed (cm/s) RunOn BL',[recName '_SpeedTimeRunBLOpt'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverTimeRunBefRun(1,:),'Speed (cm/s) RunOn BefRun',[recName '_SpeedTimeRunBefRunOpt'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverTimeRun0to1(1,:),'Speed (cm/s) RunOn 0-1s',[recName '_SpeedTimeRun0to1sOpt'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverTimeRun3to5(1,:),'Speed (cm/s) RunOn 3-5s',[recName '_SpeedTimeRun3to5sOpt'],path,[-1 100]);
    
    runPlot(sessData.Run(2,:),param.spaceSteps/10,'Dist (cm)','Speed (cm/s)','Run onset (opto ctrl)');
    saveas(gcf,[path recName '_SpeedDistRunOptCtrl']);
    print ('-painters', '-dpdf', [path recName '_SpeedDistRunOptCtrl'], '-r600')    
    runPlot(sessData.Cue(2,:),param.spaceSteps/10,'Dist (cm)','Speed (cm/s)','Cue onset (opto ctrl)');
    saveas(gcf,[path recName '_SpeedDistCueOptCtrl']);
    print ('-painters', '-dpdf', [path recName '_SpeedDistCueOptCtrl'], '-r600')    

    % changed by Yingxue 9/12/2020, using SpeedOverDistRun0to100 instead of
    % sessData.meanSpeedOverDistRun0to100 
    plotBoxPlot(sessData.meanSpeedOverDistRun0to100(2,:),'Speed (cm/s) RunOn 0-100cm',[recName '_SpeedDistRun0to100OptCtrl'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverDistRunAfter100(2,:),'Speed (cm/s) RunOn After100cm',[recName '_SpeedDistRunAfter100OptCtrl'],path,[-1 100]);

    plotBoxPlot(sessData.meanSpeedOverTimeRunBL(2,:),'Speed (cm/s) RunOn BL',[recName '_SpeedTimeRunBLOptCtrl'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverTimeRunBefRun(2,:),'Speed (cm/s) RunOn BefRun',[recName '_SpeedTimeRunBefRunOptCtrl'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverTimeRun0to1(2,:),'Speed (cm/s) RunOn 0-1s',[recName '_SpeedTimeRun0to1sOptCtrl'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverTimeRun3to5(2,:),'Speed (cm/s) RunOn 3-5s',[recName '_SpeedTimeRun3to5sOptCtrl'],path,[-1 100]);
    
    %% added by Yingxue 1/25/2022 
    runPlot(sessData.Run(3,:),param.spaceSteps/10,'Dist (cm)','Speed (cm/s)','Run onset (ctrl after opto)');
    saveas(gcf,[path recName '_SpeedDistRunOptCtrl1']);
    print ('-painters', '-dpdf', [path recName '_SpeedDistRunOptCtrl1'], '-r600')    
    runPlot(sessData.Cue(3,:),param.spaceSteps/10,'Dist (cm)','Speed (cm/s)','Cue onset (ctrl after opto)');
    saveas(gcf,[path recName '_SpeedDistCueOptCtrl1']);
    print ('-painters', '-dpdf', [path recName '_SpeedDistCueOptCtrl1'], '-r600')    

    plotBoxPlot(sessData.meanSpeedOverDistRun0to100(3,:),'Speed (cm/s) RunOn 0-100cm',[recName '_SpeedDistRun0to100OptCtrl1'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverDistRunAfter100(3,:),'Speed (cm/s) RunOn After100cm',[recName '_SpeedDistRunAfter100OptCtrl1'],path,[-1 100]);

    plotBoxPlot(sessData.meanSpeedOverTimeRunBL(3,:),'Speed (cm/s) RunOn BL',[recName '_SpeedTimeRunBLOptCtrl1'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverTimeRunBefRun(3,:),'Speed (cm/s) RunOn BefRun',[recName '_SpeedTimeRunBefRunOptCtrl1'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverTimeRun0to1(3,:),'Speed (cm/s) RunOn 0-1s',[recName '_SpeedTimeRun0to1sOptCtrl1'],path,[-1 100]);
    plotBoxPlot(sessData.meanSpeedOverTimeRun3to5(3,:),'Speed (cm/s) RunOn 3-5s',[recName '_SpeedTimeRun3to5sOptCtrl1'],path,[-1 100]);
end
 
function runPlot(data,x,xl,yl,ti)
    figure
    nSess = length(data);
    for i = 1:nSess
        options.handle     = subplot(1,nSess,i);
        options.color_area = [224 191 232]./255; %[128 193 219]./255;    % Blue theme
        options.color_line = [129 98 137]./255; %[ 52 148 186]./255;
        options.alpha      = 0.5;
        options.line_width = 2;
        options.error      = 'sem';
        options.x_axis = x;
        plot_areaerrorbarSub(data{i}, options);
        set(gca,'XLim',[0 max(x)],'YLim',[-2 100])
        xlabel(xl)
        ylabel(yl)
    end
end

function sessData = compRecTime(path,recName,sessNo,injectionSess)
    nSess = length(sessNo);
    sessData.sessNo = sessNo;
    sessData.timeEnd = [];
    sessData.timeStart = [];
    sessData.duration = zeros(1,nSess);
    sessData.tDiffInjStart = zeros(1,nSess);
    sessData.tDiffInjEnd = zeros(1,nSess);
    
    fullpath = [path recName '_recTime.mat'];
    if(~exist(fullpath,'file'))
        disp([fullpath ' does not exist.']);
        return;
    end
    load(fullpath,'fInfo');
    ind = fInfo.sessNo == injectionSess;
    injTime = datetime(fInfo.time{ind});
    
    for i = 1:nSess
        fileName = [recName '-0' num2str(sessNo(i))];
        pathTmp = [path fileName '\'];
        
        fullpath = [pathTmp fileName 'B.mat'];
        if(~exist(fullpath,'file'))
            disp([fullpath ' does not exist.']);
            return;
        end
        load(fullpath,'behEvents');
        timeStart = behEvents.trialT(1,1);
        fnames = fieldnames(behEvents);
        timeEnd = timeStart;
        for n = 1:length(fnames)
            if(strcmp(fnames{n},'movieTDescr'))
                continue;
            end
            timeEnd = max(timeEnd, behEvents.(fnames{n})(end,1));
        end
        sessData.duration(i) = (timeEnd - timeStart)/1000; % sec
        
        ind = fInfo.sessNo == sessNo(i);
        sessData.timeEnd{i} = datetime(fInfo.time{ind});
        sessData.timeStart{i} = datetime(fInfo.time{ind}) - ...
            seconds(sessData.duration(i));
        
        sessData.tDiffInjStart(i) = seconds(sessData.timeStart{i} - injTime);
        sessData.tDiffInjEnd(i) = seconds(sessData.timeEnd{i} - injTime);
    end
end