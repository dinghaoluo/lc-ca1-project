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
    for i = 3 %1:nCond
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
            rankXY,rankXZ,rankYZ,'Session','Percent. rewarded',['Cond' num2str(i)],path,[0 1]);
        
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
%         
        colorSel = 0;
        rankXY = rankRecDataRunPreM{i}.pRSMeanLickRndSelMeanPerRec;
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
            meanLickPost = [];
        else
            rankXZ = rankRecDataRunPrePost{i}.pRSMeanLickRndSelMeanPerRec;
            rankYZ = rankRecDataRunPostM{i}.pRSMeanLickRndSelMeanPerRec;
            meanLickPost = mean(recDataRunPost{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLickAfter30),2);
        end
        plotBoxPlot(mean(recDataRunPre{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLickAfter30),2),...
            meanLickPost,...
            mean(recDataRunManip{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLickAfter30),2),...
            'No. licks',['lickProfileRndSelBox' num2str(i)],...
            folderPath,[-0.1 1.0],rankXY,rankXZ,rankYZ,colorSel);
        
        rankXY = rankRecDataRunPreM{i}.pRSMeanLickRndSelMeanPerRec30to150;
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
            meanLickPost = [];
        else
            rankXZ = rankRecDataRunPrePost{i}.pRSMeanLickRndSelMeanPerRec30to150;
            rankYZ = rankRecDataRunPostM{i}.pRSMeanLickRndSelMeanPerRec30to150;
            meanLickPost = mean(recDataRunPost{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLick30to150),2);
        end
        plotBoxPlot(mean(recDataRunPre{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLick30to150),2),...
            meanLickPost,...
            mean(recDataRunManip{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLick30to150),2),...
            'No. licks (30-150 cm)',['lickProfileRndSel30to150Box' num2str(i)],...
            folderPath,[-0.1 1.0],rankXY,rankXZ,rankYZ,colorSel);
        
        rankXY = rankRecDataRunPreM{i}.pRSMeanLickRndSelMeanPerRec30to100;
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
            meanLickPost = [];
        else
            rankXZ = rankRecDataRunPrePost{i}.pRSMeanLickRndSelMeanPerRec30to100;
            rankYZ = rankRecDataRunPostM{i}.pRSMeanLickRndSelMeanPerRec30to100;
            meanLickPost = mean(recDataRunPost{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLick30to100),2);
        end
        plotBoxPlot(mean(recDataRunPre{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLick30to100),2),...
            meanLickPost,...
            mean(recDataRunManip{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLick30to100),2),...
            'No. licks (30-100 cm)',['lickProfileRndSel30to100Box' num2str(i)],...
            folderPath,[-0.1 0.1],rankXY,rankXZ,rankYZ,colorSel);
        
        rankXY = rankRecDataRunPreM{i}.pRSMeanLickRndSelMeanPerRec100to150;
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
            meanLickPost = [];
        else
            rankXZ = rankRecDataRunPrePost{i}.pRSMeanLickRndSelMeanPerRec100to150;
            rankYZ = rankRecDataRunPostM{i}.pRSMeanLickRndSelMeanPerRec100to150;
            meanLickPost = mean(recDataRunPost{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLick100to150),2);
        end
        plotBoxPlot(mean(recDataRunPre{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLick100to150),2),...
            meanLickPost,...
            mean(recDataRunManip{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLick100to150),2),...
            'No. licks (100-150 cm)',['lickProfileRndSel100to150Box' num2str(i)],...
            folderPath,[-0.1 1.5],rankXY,rankXZ,rankYZ,colorSel);
        
        rankXY = rankRecDataRunPreM{i}.pRSMeanLickRndSelMeanPerRec150to180;
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
            meanLickPost = [];
        else
            rankXZ = rankRecDataRunPrePost{i}.pRSMeanLickRndSelMeanPerRec150to180;
            rankYZ = rankRecDataRunPostM{i}.pRSMeanLickRndSelMeanPerRec150to180;
            meanLickPost = mean(recDataRunPost{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLick150to180),2);
        end
        plotBoxPlot(mean(recDataRunPre{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLick150to180),2),...
            meanLickPost,...
            mean(recDataRunManip{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLick150to180),2),...
            'No. licks (150-180 cm)',['lickProfileRndSel150to180Box' num2str(i)],...
            folderPath,[-0.1 2.5],rankXY,rankXZ,rankYZ,colorSel);
        
        rankXY = rankRecDataRunPreM{i}.pRSMeanLickRndSelMeanPerRecAfter180;
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
            meanLickPost = [];
        else
            rankXZ = rankRecDataRunPrePost{i}.pRSMeanLickRndSelMeanPerRecAfter180;
            rankYZ = rankRecDataRunPostM{i}.pRSMeanLickRndSelMeanPerRecAfter180;
            meanLickPost = mean(recDataRunPost{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLickAfter180),2);
        end
        plotBoxPlot(mean(recDataRunPre{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLickAfter180),2),...
            meanLickPost,...
            mean(recDataRunManip{i}.lickProfileRndSel(:,rankRecDataRunPreM{i}.indLickAfter180),2),...
            'No. licks (>180 cm)',['lickProfileRndSelAfter180Box' num2str(i)],...
            folderPath,[-0.1 3.5],rankXY,rankXZ,rankYZ,colorSel);
    
        rankXY = rankRecDataRunPreM{i}.pRSMeanSpeedRndSelMeanPerRec;
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
            meanSpeedPost = [];
        else
            rankXZ = rankRecDataRunPrePost{i}.pRSMeanSpeedRndSelMeanPerRec;
            rankYZ = rankRecDataRunPostM{i}.pRSMeanSpeedRndSelMeanPerRec;
            meanSpeedPost = mean(recDataRunPost{i}.speedProfileRndSel,2);
        end
        plotBoxPlot(mean(recDataRunPre{i}.speedProfileRndSel,2),...
            meanSpeedPost,...
            mean(recDataRunManip{i}.speedProfileRndSel,2),...
            'Speed (cm/s)',['speedProfileRndSelBox' num2str(i)],...
            folderPath,[10 125],rankXY,rankXZ,rankYZ,colorSel);
        
        rankXY = rankRecDataRunPreM{i}.pRSMeanSpeedRndSelMeanPerRec30to100;
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
            meanSpeedPost = [];
        else
            rankXZ = rankRecDataRunPrePost{i}.pRSMeanSpeedRndSelMeanPerRec30to100;
            rankYZ = rankRecDataRunPostM{i}.pRSMeanSpeedRndSelMeanPerRec30to100;
            meanSpeedPost = mean(recDataRunPost{i}.speedProfileRndSel(:,rankRecDataRunPreM{i}.indSpeed30to100),2);
        end
        plotBoxPlot(mean(recDataRunPre{i}.speedProfileRndSel(:,rankRecDataRunPreM{i}.indSpeed30to100),2),...
            meanSpeedPost,...
            mean(recDataRunManip{i}.speedProfileRndSel(:,rankRecDataRunPreM{i}.indSpeed30to100),2),...
            'Speed (cm/s) 30-100 cm',['speedProfileRndSel30to100Box' num2str(i)],...
            folderPath,[10 120],rankXY,rankXZ,rankYZ,colorSel);
        
        rankXY = rankRecDataRunPreM{i}.pRSMeanSpeedRndSelMeanPerRecAfter100;
        if(isempty(rankRecDataRunPrePost{i}))
            rankXZ = [];
            rankYZ = [];
            meanSpeedPost = [];
        else
            rankXZ = rankRecDataRunPrePost{i}.pRSMeanSpeedRndSelMeanPerRecAfter100;
            rankYZ = rankRecDataRunPostM{i}.pRSMeanSpeedRndSelMeanPerRecAfter100;
            meanSpeedPost = mean(recDataRunPost{i}.speedProfileRndSel(:,rankRecDataRunPreM{i}.indSpeedAfter100),2);
        end
        plotBoxPlot(mean(recDataRunPre{i}.speedProfileRndSel(:,rankRecDataRunPreM{i}.indSpeedAfter100),2),...
            meanSpeedPost,...
            mean(recDataRunManip{i}.speedProfileRndSel(:,rankRecDataRunPreM{i}.indSpeedAfter100),2),...
            'Speed (cm/s) >100 cm',['speedProfileRndSelAfter100Box' num2str(i)],...
            folderPath,[10 100],rankXY,rankXZ,rankYZ,colorSel);
        
        plotTrace(folderPath,recDataRunManip{i}.spaceStepsLick/10,...
            recDataRunManip{i}.lickProfileRndSel,...
            recDataRunPre{i}.lickProfileRndSel,...
            [30,210],[0 2],'No. lick',['lickProfileRndSel' num2str(i)])
        
        plotTrace(folderPath,recDataRunManip{i}.spaceStepsSpeed/10,...
            recDataRunManip{i}.speedProfileRndSel,...
            recDataRunPre{i}.speedProfileRndSel,...
            [0,180],[0 100],'Speed (cm/sec)',['speedProfileRndSel' num2str(i)])
        
%         pause;
        close all;
    end
end

function plotTrace(pathAnal,timeSteps,mani,pre,xli,yli,yl,filename)
    options.handle     = figure;
    set(options.handle,'OuterPosition',...
        [500 500 280 280])
    options.color_areaX = [27 115 187]./255;    % Blue theme
    options.color_lineX = [ 39 169 225]./255;
    options.color_areaY = [187 189 192]./255;    % Orange theme
    options.color_lineY = [167 169  151]./255;
    options.alpha      = 0.5;
    options.line_width = 0.5;
    options.error      = 'sem';
    options.x_axisX = timeSteps;
    options.x_axisY = timeSteps;
    plot_areaerrorbarXY(mani, pre,...
        options);
    set(gca,'XLim',xli);
    set(gca,'YLim',yli);
    xlabel('Dist (cm)')
    ylabel(yl)
    legend('','M','','Pre')
    
    fileName1 = [pathAnal filename...
        '_Mani_Pre'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotBoxPlot(x1,x2,x3,yl,fn,pathAnal,ylimit,p12,p13,p23,colorSel)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 200 400])
    if(colorSel == 0)
        colorArr = [163 207 98;...
                234 131 114;...
                163 207 98]/255;
    elseif(colorSel == 1)            
        colorArr = [234 131 114;...
                116 53 61;...
                234 131 114]/255;
    else        
        colorArr = [163 207 98;... 
            63 79 37;...
            163 207 98]/255;
    end
    x = [x1;x3;x2];
    g = [repmat(1,length(x1),1);...
        repmat(2,length(x3),1);...
        repmat(3,length(x2),1)];
    boxplot(x,g,'Notch','on','Widths',0.3,'Symbol','');
    h = findobj(gca,'Tag','Box');
    for j = 1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colorArr(j,:),'FaceAlpha',0.5);
    end
    if(~isempty(ylimit))
        set(gca,'YLim',ylimit);
    end
    ylabel(yl);
    text(1,ylimit(2)*0.75,num2str(p12,'p = %f'));
    if(~isempty(p23))
        text(2.2,ylimit(2)*0.8,num2str(p23,'p = %f'));
    end
    if(~isempty(p13))
        text(1.5,ylimit(2)*0.85,num2str(p13,'p = %f'));
    end
    
    fileName1 = [pathAnal fn];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotStat(dataX,dataY,dataZ,meanX,meanY,meanZ,semX,semY,semZ,...
                rankXY,rankXZ,rankYZ,xlab,ylab,ti,path,ylimit) 
    if(nargin == 16)
        ylimit = [];
    end
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
    if(~isempty(ylimit))
        set(gca,'YLim',ylimit);
    end
    print('-painters', '-dpdf', path, '-r600')
    print('-painters', '-dpdf', path, '-r600')
    savefig([path '.fig']);
end