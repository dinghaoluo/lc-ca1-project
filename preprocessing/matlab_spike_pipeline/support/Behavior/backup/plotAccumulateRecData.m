function plotAccumulateRecData()
    RecordingListNT;
    
    fullPath = [folderPath 'allRecData.mat'];
    if(~exist(fullPath,'file'))
        disp([fullPath ' does not exist.']);
        return;
    end
    load(fullPath,'recDataRunPre','recDataRunManip','recDataRunPost');
    
    fullPath = [folderPath 'allRecDataStats.mat'];
    load(fullPath,'meanRecDataRunPre','meanRecDataRunManip','meanRecDataRunPost',...
        'semRecDataRunPre','semRecDataRunManip','semRecDataRunPost',...
        'rankRecDataRunPrePost','rankRecDataRunPreM','rankRecDataRunPostM',...
        'anovaRecDataRun');
    
    nCond = length(meanRecDataRunPre);
    for i = 10 %:nCond
        path = [folderPath '_numSamples_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.numSamples;anovaRecDataRun{i}.numSamplesPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.numSamples;anovaRecDataRun{i}.numSamplesPrePost];
            rankYZ = [rankRecDataRunPostM{i}.numSamples;anovaRecDataRun{i}.numSamplesPostM];
        end
        plotStat(recDataRunPre{i}.numSamplesMean,recDataRunManip{i}.numSamplesMean,...
            recDataRunPost{i}.numSamplesMean,meanRecDataRunPre{i}.numSamples,...
            meanRecDataRunManip{i}.numSamples,meanRecDataRunPost{i}.numSamples,...
            semRecDataRunPre{i}.numSamples,semRecDataRunManip{i}.numSamples,...
            semRecDataRunPost{i}.numSamples,...
            rankXY,rankXZ,rankYZ,'Session','No. samples',['Cond' num2str(i)],path);
        
        path = [folderPath '_maxSpeed_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.maxSpeed;anovaRecDataRun{i}.maxSpeedPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.maxSpeed;anovaRecDataRun{i}.maxSpeedPrePost];
            rankYZ = [rankRecDataRunPostM{i}.maxSpeed;anovaRecDataRun{i}.maxSpeedPostM];
        end
        plotStat(recDataRunPre{i}.maxSpeedMean,recDataRunManip{i}.maxSpeedMean,...
            recDataRunPost{i}.maxSpeedMean,meanRecDataRunPre{i}.maxSpeed,...
            meanRecDataRunManip{i}.maxSpeed,meanRecDataRunPost{i}.maxSpeed,...
            semRecDataRunPre{i}.maxSpeed,semRecDataRunManip{i}.maxSpeed,...
            semRecDataRunPost{i}.maxSpeed,...
            rankXY,rankXZ,rankYZ,'Session','Max speed (mm/s)',['Cond' num2str(i)],path);
        
        path = [folderPath '_meanSpeed_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.meanSpeed;anovaRecDataRun{i}.meanSpeedPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.meanSpeed;anovaRecDataRun{i}.meanSpeedPrePost];
            rankYZ = [rankRecDataRunPostM{i}.meanSpeed;anovaRecDataRun{i}.meanSpeedPostM];
        end
        plotStat(recDataRunPre{i}.meanSpeedMean,recDataRunManip{i}.meanSpeedMean,...
            recDataRunPost{i}.meanSpeedMean,meanRecDataRunPre{i}.meanSpeed,...
            meanRecDataRunManip{i}.meanSpeed,meanRecDataRunPost{i}.meanSpeed,...
            semRecDataRunPre{i}.meanSpeed,semRecDataRunManip{i}.meanSpeed,...
            semRecDataRunPost{i}.meanSpeed,...
            rankXY,rankXZ,rankYZ,'Session','Mean speed (mm/s)',['Cond' num2str(i)],path);
            
        path = [folderPath '_maxRunLenT_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.maxRunLenT;anovaRecDataRun{i}.maxRunLenTPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.maxRunLenT;anovaRecDataRun{i}.maxRunLenTPrePost];
            rankYZ = [rankRecDataRunPostM{i}.maxRunLenT;anovaRecDataRun{i}.maxRunLenTPostM];
        end
        plotStat(recDataRunPre{i}.maxRunLenTMean,recDataRunManip{i}.maxRunLenTMean,...
            recDataRunPost{i}.maxRunLenTMean,meanRecDataRunPre{i}.maxRunLenT,...
            meanRecDataRunManip{i}.maxRunLenT,meanRecDataRunPost{i}.maxRunLenT,...
            semRecDataRunPre{i}.maxRunLenT,semRecDataRunManip{i}.maxRunLenT,...
            semRecDataRunPost{i}.maxRunLenT,...
            rankXY,rankXZ,rankYZ,'Session','Longest running bout (s)',['Cond' num2str(i)],path);
        
        path = [folderPath '_totRunLenT_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.totRunLenT;anovaRecDataRun{i}.totRunLenTPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.totRunLenT;anovaRecDataRun{i}.totRunLenTPrePost];
            rankYZ = [rankRecDataRunPostM{i}.totRunLenT;anovaRecDataRun{i}.totRunLenTPostM];
        end
        plotStat(recDataRunPre{i}.totRunLenTMean,recDataRunManip{i}.totRunLenTMean,...
            recDataRunPost{i}.totRunLenTMean,meanRecDataRunPre{i}.totRunLenT,...
            meanRecDataRunManip{i}.totRunLenT,meanRecDataRunPost{i}.totRunLenT,...
            semRecDataRunPre{i}.totRunLenT,semRecDataRunManip{i}.totRunLenT,...
            semRecDataRunPost{i}.totRunLenT,...
            rankXY,rankXZ,rankYZ,'Session','Total run time (s)',['Cond' num2str(i)],path);
        
        path = [folderPath '_numRun_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.numRun;anovaRecDataRun{i}.numRunPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.numRun;anovaRecDataRun{i}.numRunPrePost];
            rankYZ = [rankRecDataRunPostM{i}.numRun;anovaRecDataRun{i}.numRunPostM];
        end
        plotStat(recDataRunPre{i}.numRunMean,recDataRunManip{i}.numRunMean,...
            recDataRunPost{i}.numRunMean,meanRecDataRunPre{i}.numRun,...
            meanRecDataRunManip{i}.numRun,meanRecDataRunPost{i}.numRun,...
            semRecDataRunPre{i}.numRun,semRecDataRunManip{i}.numRun,...
            semRecDataRunPost{i}.numRun,...
            rankXY,rankXZ,rankYZ,'Session','No. running bouts',['Cond' num2str(i)],path);
        
        path = [folderPath '_maxAcc_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.maxAcc;anovaRecDataRun{i}.maxAccPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.maxAcc;anovaRecDataRun{i}.maxAccPrePost];
            rankYZ = [rankRecDataRunPostM{i}.maxAcc;anovaRecDataRun{i}.maxAccPostM];
        end
        plotStat(recDataRunPre{i}.maxAccMean,recDataRunManip{i}.maxAccMean,...
            recDataRunPost{i}.maxAccMean,meanRecDataRunPre{i}.maxAcc,...
            meanRecDataRunManip{i}.maxAcc,meanRecDataRunPost{i}.maxAcc,...
            semRecDataRunPre{i}.maxAcc,semRecDataRunManip{i}.maxAcc,...
            semRecDataRunPost{i}.maxAcc,...
            rankXY,rankXZ,rankYZ,'Session','Max acceleration (mm/s^2)',['Cond' num2str(i)],path);
        
        path = [folderPath '_meanAcc_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.meanAcc;anovaRecDataRun{i}.meanAccPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.meanAcc;anovaRecDataRun{i}.meanAccPrePost];
            rankYZ = [rankRecDataRunPostM{i}.meanAcc;anovaRecDataRun{i}.meanAccPostM];
        end
        plotStat(recDataRunPre{i}.meanAccMean,recDataRunManip{i}.meanAccMean,...
            recDataRunPost{i}.meanAccMean,meanRecDataRunPre{i}.meanAcc,...
            meanRecDataRunManip{i}.meanAcc,meanRecDataRunPost{i}.meanAcc,...
            semRecDataRunPre{i}.meanAcc,semRecDataRunManip{i}.meanAcc,...
            semRecDataRunPost{i}.meanAcc,...
            rankXY,rankXZ,rankYZ,'Session','Mean acceleration (mm/s^2)',['Cond' num2str(i)],path);
        
        path = [folderPath '_totStopLenT_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.totStopLenT;anovaRecDataRun{i}.totStopLenTPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.totStopLenT;anovaRecDataRun{i}.totStopLenTPrePost];
            rankYZ = [rankRecDataRunPostM{i}.totStopLenT;anovaRecDataRun{i}.totStopLenTPostM];
        end
        plotStat(recDataRunPre{i}.totStopLenTMean,recDataRunManip{i}.totStopLenTMean,...
            recDataRunPost{i}.totStopLenTMean,meanRecDataRunPre{i}.totStopLenT,...
            meanRecDataRunManip{i}.totStopLenT,meanRecDataRunPost{i}.totStopLenT,...
            semRecDataRunPre{i}.totStopLenT,semRecDataRunManip{i}.totStopLenT,...
            semRecDataRunPost{i}.totStopLenT,...
            rankXY,rankXZ,rankYZ,'Session','Total stop time (s)',['Cond' num2str(i)],path);
        
        path = [folderPath '_startCueToRun_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.startCueToRun;anovaRecDataRun{i}.startCueToRunPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.startCueToRun;anovaRecDataRun{i}.startCueToRunPrePost];
            rankYZ = [rankRecDataRunPostM{i}.startCueToRun;anovaRecDataRun{i}.startCueToRunPostM];
        end
        plotStat(recDataRunPre{i}.startCueToRunMean,recDataRunManip{i}.startCueToRunMean,...
            recDataRunPost{i}.startCueToRunMean,meanRecDataRunPre{i}.startCueToRun,...
            meanRecDataRunManip{i}.startCueToRun,meanRecDataRunPost{i}.startCueToRun,...
            semRecDataRunPre{i}.startCueToRun,semRecDataRunManip{i}.startCueToRun,...
            semRecDataRunPost{i}.startCueToRun,...
            rankXY,rankXZ,rankYZ,'Session','Start cue to run onset (s)',['Cond' num2str(i)],path);
        
        path = [folderPath '_numLicksBefRew_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.numLicksBefRew;anovaRecDataRun{i}.numLicksBefRewPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.numLicksBefRew;anovaRecDataRun{i}.numLicksBefRewPrePost];
            rankYZ = [rankRecDataRunPostM{i}.numLicksBefRew;anovaRecDataRun{i}.numLicksBefRewPostM];
        end
        plotStat(recDataRunPre{i}.numLicksBefRewMean,recDataRunManip{i}.numLicksBefRewMean,...
            recDataRunPost{i}.numLicksBefRewMean,meanRecDataRunPre{i}.numLicksBefRew,...
            meanRecDataRunManip{i}.numLicksBefRew,meanRecDataRunPost{i}.numLicksBefRew,...
            semRecDataRunPre{i}.numLicksBefRew,semRecDataRunManip{i}.numLicksBefRew,...
            semRecDataRunPost{i}.numLicksBefRew,...
            rankXY,rankXZ,rankYZ,'Session','No. licks before reward',['Cond' num2str(i)],path);
        
        path = [folderPath '_numLicksRew_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.numLicksRew;anovaRecDataRun{i}.numLicksRewPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.numLicksRew;anovaRecDataRun{i}.numLicksRewPrePost];
            rankYZ = [rankRecDataRunPostM{i}.numLicksRew;anovaRecDataRun{i}.numLicksRewPostM];
        end
        plotStat(recDataRunPre{i}.numLicksRewMean,recDataRunManip{i}.numLicksRewMean,...
            recDataRunPost{i}.numLicksRewMean,meanRecDataRunPre{i}.numLicksRew,...
            meanRecDataRunManip{i}.numLicksRew,meanRecDataRunPost{i}.numLicksRew,...
            semRecDataRunPre{i}.numLicksRew,semRecDataRunManip{i}.numLicksRew,...
            semRecDataRunPost{i}.numLicksRew,...
            rankXY,rankXZ,rankYZ,'Session','No. licks after reward',['Cond' num2str(i)],path);
        
        path = [folderPath '_med1stFiveLickDist_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.med1stFiveLickDist;anovaRecDataRun{i}.med1stFiveLickDistPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.med1stFiveLickDist;anovaRecDataRun{i}.med1stFiveLickDistPrePost];
            rankYZ = [rankRecDataRunPostM{i}.med1stFiveLickDist;anovaRecDataRun{i}.med1stFiveLickDistPostM];
        end
        plotStat(recDataRunPre{i}.med1stFiveLickDistMean,recDataRunManip{i}.med1stFiveLickDistMean,...
            recDataRunPost{i}.med1stFiveLickDistMean,meanRecDataRunPre{i}.med1stFiveLickDist,...
            meanRecDataRunManip{i}.med1stFiveLickDist,meanRecDataRunPost{i}.med1stFiveLickDist,...
            semRecDataRunPre{i}.med1stFiveLickDist,semRecDataRunManip{i}.med1stFiveLickDist,...
            semRecDataRunPost{i}.med1stFiveLickDist,...
            rankXY,rankXZ,rankYZ,'Session','Median first-5-lick distance (mm)',['Cond' num2str(i)],path);
        
        path = [folderPath '_medLickDist_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.medLickDist;anovaRecDataRun{i}.medLickDistPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.medLickDist;anovaRecDataRun{i}.medLickDistPrePost];
            rankYZ = [rankRecDataRunPostM{i}.medLickDist;anovaRecDataRun{i}.medLickDistPostM];
        end
        plotStat(recDataRunPre{i}.medLickDistMean,recDataRunManip{i}.medLickDistMean,...
            recDataRunPost{i}.medLickDistMean,meanRecDataRunPre{i}.medLickDist,...
            meanRecDataRunManip{i}.medLickDist,meanRecDataRunPost{i}.medLickDist,...
            semRecDataRunPre{i}.medLickDist,semRecDataRunManip{i}.medLickDist,...
            semRecDataRunPost{i}.medLickDist,...
            rankXY,rankXZ,rankYZ,'Session','Median lick distance (mm)',['Cond' num2str(i)],path);
        
        path = [folderPath '_med1stFiveLickDistBefRew_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.med1stFiveLickDistBefRew;anovaRecDataRun{i}.med1stFiveLickDistBefRewPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.med1stFiveLickDistBefRew;anovaRecDataRun{i}.med1stFiveLickDistBefRewPrePost];
            rankYZ = [rankRecDataRunPostM{i}.med1stFiveLickDistBefRew;anovaRecDataRun{i}.med1stFiveLickDistBefRewPostM];
        end
        plotStat(recDataRunPre{i}.med1stFiveLickDistBefRewMean,recDataRunManip{i}.med1stFiveLickDistBefRewMean,...
            recDataRunPost{i}.med1stFiveLickDistBefRewMean,meanRecDataRunPre{i}.med1stFiveLickDistBefRew,...
            meanRecDataRunManip{i}.med1stFiveLickDistBefRew,meanRecDataRunPost{i}.med1stFiveLickDistBefRew,...
            semRecDataRunPre{i}.med1stFiveLickDistBefRew,semRecDataRunManip{i}.med1stFiveLickDistBefRew,...
            semRecDataRunPost{i}.med1stFiveLickDistBefRew,...
            rankXY,rankXZ,rankYZ,'Session','Median first-5-lick distance befrore reward (mm)',['Cond' num2str(i)],path);
        
        path = [folderPath '_medLickDistBefRew_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.medLickDistBefRew;anovaRecDataRun{i}.medLickDistBefRewPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.medLickDistBefRew;anovaRecDataRun{i}.medLickDistBefRewPrePost];
            rankYZ = [rankRecDataRunPostM{i}.medLickDistBefRew;anovaRecDataRun{i}.medLickDistBefRewPostM];
        end
        plotStat(recDataRunPre{i}.medLickDistBefRewMean,recDataRunManip{i}.medLickDistBefRewMean,...
            recDataRunPost{i}.medLickDistBefRewMean,meanRecDataRunPre{i}.medLickDistBefRew,...
            meanRecDataRunManip{i}.medLickDistBefRew,meanRecDataRunPost{i}.medLickDistBefRew,...
            semRecDataRunPre{i}.medLickDistBefRew,semRecDataRunManip{i}.medLickDistBefRew,...
            semRecDataRunPost{i}.medLickDistBefRew,...
            rankXY,rankXZ,rankYZ,'Session','Median lick distance befrore reward (mm)',['Cond' num2str(i)],path);
        
        path = [folderPath '_percRewarded_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.percRewarded;anovaRecDataRun{i}.percRewardedPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.percRewarded;anovaRecDataRun{i}.percRewardedPrePost];
            rankYZ = [rankRecDataRunPostM{i}.percRewarded;anovaRecDataRun{i}.percRewardedPostM];
        end
        plotStat(recDataRunPre{i}.percRewarded,recDataRunManip{i}.percRewarded,...
            recDataRunPost{i}.percRewarded,meanRecDataRunPre{i}.percRewarded,...
            meanRecDataRunManip{i}.percRewarded,meanRecDataRunPost{i}.percRewarded,...
            semRecDataRunPre{i}.percRewarded,semRecDataRunManip{i}.percRewarded,...
            semRecDataRunPost{i}.percRewarded,...
            rankXY,rankXZ,rankYZ,'Session','Percent. rewarded',['Cond' num2str(i)],path);
        
        path = [folderPath '_percNonStop_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.percNonStop;anovaRecDataRun{i}.percNonStopPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.percNonStop;anovaRecDataRun{i}.percNonStopPrePost];
            rankYZ = [rankRecDataRunPostM{i}.percNonStop;anovaRecDataRun{i}.percNonStopPostM];
        end
        plotStat(recDataRunPre{i}.percNonStop,recDataRunManip{i}.percNonStop,...
            recDataRunPost{i}.percNonStop,meanRecDataRunPre{i}.percNonStop,...
            meanRecDataRunManip{i}.percNonStop,meanRecDataRunPost{i}.percNonStop,...
            semRecDataRunPre{i}.percNonStop,semRecDataRunManip{i}.percNonStop,...
            semRecDataRunPost{i}.percNonStop,...
            rankXY,rankXZ,rankYZ,'Session','Percent. non-stop trials',['Cond' num2str(i)],path);
        
        path = [folderPath '_pumpLfpInd_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.pumpLfpInd;anovaRecDataRun{i}.pumpLfpIndPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.pumpLfpInd;anovaRecDataRun{i}.pumpLfpIndPrePost];
            rankYZ = [rankRecDataRunPostM{i}.pumpLfpInd;anovaRecDataRun{i}.pumpLfpIndPostM];
        end
        plotStat(recDataRunPre{i}.pumpLfpIndMean,recDataRunManip{i}.pumpLfpIndMean,...
            recDataRunPost{i}.pumpLfpIndMean,meanRecDataRunPre{i}.pumpLfpInd,...
            meanRecDataRunManip{i}.pumpLfpInd,meanRecDataRunPost{i}.pumpLfpInd,...
            semRecDataRunPre{i}.pumpLfpInd,semRecDataRunManip{i}.pumpLfpInd,...
            semRecDataRunPost{i}.pumpLfpInd,...
            rankXY,rankXZ,rankYZ,'Session','Pump on lfp index',['Cond' num2str(i)],path);
        
        path = [folderPath '_pumpMM_cond' num2str(i)];
        rankXY = [rankRecDataRunPreM{i}.pumpMM;anovaRecDataRun{i}.pumpMMPreM];
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRunPrePost{i}.pumpMM;anovaRecDataRun{i}.pumpMMPrePost];
            rankYZ = [rankRecDataRunPostM{i}.pumpMM;anovaRecDataRun{i}.pumpMMPostM];
        end
        plotStat(recDataRunPre{i}.pumpMMMean,recDataRunManip{i}.pumpMMMean,...
            recDataRunPost{i}.pumpMMMean,meanRecDataRunPre{i}.pumpMM,...
            meanRecDataRunManip{i}.pumpMM,meanRecDataRunPost{i}.pumpMM,...
            semRecDataRunPre{i}.pumpMM,semRecDataRunManip{i}.pumpMM,...
            semRecDataRunPost{i}.pumpMM,...
            rankXY,rankXZ,rankYZ,'Session','Pump on distance (mm)',['Cond' num2str(i)],path);
        
%         pause;
        close all;
    end
end

function plotStat(dataX,dataY,dataZ,meanX,meanY,meanZ,semX,semY,semZ,...
                rankXY,rankXZ,rankYZ,xlab,ylab,ti,path)    
    if(~isempty(rankXZ))  
        dataArr = [dataX' dataY' dataZ'];
        meanArr = [meanX meanY meanZ];
        errBar = [semX semY semZ];
        
        barPlotWithStat(1:3,meanArr,errBar,dataArr,xlab,ylab,ti,rankXY,rankXZ,rankYZ);
    else
        dataArr = [dataX' dataY'];
        meanArr = [meanX meanY];
        errBar = [semX semY];
        
        barPlotWithStat(1:2,meanArr,errBar,dataArr,xlab,ylab,ti,rankXY,rankXZ,rankYZ);       
    end
    print('-painters', '-dpdf', path, '-r600')
    savefig([path '.fig']);
end