function plotAccumulateRecDataCue()
    RecordingListNT;
    
    fullPath = [folderPath 'allRecData.mat'];
    if(~exist(fullPath,'file'))
        disp([fullPath ' does not exist.']);
        return;
    end
    load(fullPath,'recDataCuePre','recDataCueManip','recDataCuePost');
    
    fullPath = [folderPath 'allRecDataStats.mat'];
    load(fullPath,'meanRecDataCuePre','meanRecDataCueManip','meanRecDataCuePost',...
        'semRecDataCuePre','semRecDataCueManip','semRecDataCuePost',...
        'rankRecDataCuePrePost','rankRecDataCuePreM','rankRecDataCuePostM',...
        'anovaRecDataCue');
    
    nCond = length(meanRecDataCuePre);
    for i = 17
%         path = [folderPath '_numSamples_cond' num2str(i) '_Cue'];
%         rankXY = [rankRecDataCuePreM{i}.numSamples;anovaRecDataCue{i}.numSamplesPreM];
%         if(isempty(rankRecDataCuePrePost{i}))
%             rankXZ = [];
%             rankYZ = [];
%         else
%             rankXZ = [rankRecDataCuePrePost{i}.numSamples;anovaRecDataCue{i}.numSamplesPrePost];
%             rankYZ = [rankRecDataCuePostM{i}.numSamples;anovaRecDataCue{i}.numSamplesPostM];
%         end
%         plotStat(recDataCuePre{i}.numSamplesMean,recDataCueManip{i}.numSamplesMean,...
%             recDataCuePost{i}.numSamplesMean,meanRecDataCuePre{i}.numSamples,...
%             meanRecDataCueManip{i}.numSamples,meanRecDataCuePost{i}.numSamples,...
%             semRecDataCuePre{i}.numSamples,semRecDataCueManip{i}.numSamples,...
%             semRecDataCuePost{i}.numSamples,...
%             rankXY,rankXZ,rankYZ,'Session','No. samples',['Cond' num2str(i)],path);
%         
%         path = [folderPath '_maxSpeed_cond' num2str(i) '_Cue'];
%         rankXY = [rankRecDataCuePreM{i}.maxSpeed;anovaRecDataCue{i}.maxSpeedPreM];
%         if(isempty(rankRecDataCuePrePost{i}))
%             rankXZ = [];
%             rankYZ = [];
%         else
%             rankXZ = [rankRecDataCuePrePost{i}.maxSpeed;anovaRecDataCue{i}.maxSpeedPrePost];
%             rankYZ = [rankRecDataCuePostM{i}.maxSpeed;anovaRecDataCue{i}.maxSpeedPostM];
%         end
%         plotStat(recDataCuePre{i}.maxSpeedMean,recDataCueManip{i}.maxSpeedMean,...
%             recDataCuePost{i}.maxSpeedMean,meanRecDataCuePre{i}.maxSpeed,...
%             meanRecDataCueManip{i}.maxSpeed,meanRecDataCuePost{i}.maxSpeed,...
%             semRecDataCuePre{i}.maxSpeed,semRecDataCueManip{i}.maxSpeed,...
%             semRecDataCuePost{i}.maxSpeed,...
%             rankXY,rankXZ,rankYZ,'Session','Max speed (mm/s)',['Cond' num2str(i)],path);
%         
%         path = [folderPath '_meanSpeed_cond' num2str(i) '_Cue'];
%         rankXY = [rankRecDataCuePreM{i}.meanSpeed;anovaRecDataCue{i}.meanSpeedPreM];
%         if(isempty(rankRecDataCuePrePost{i}))
%             rankXZ = [];
%             rankYZ = [];
%         else
%             rankXZ = [rankRecDataCuePrePost{i}.meanSpeed;anovaRecDataCue{i}.meanSpeedPrePost];
%             rankYZ = [rankRecDataCuePostM{i}.meanSpeed;anovaRecDataCue{i}.meanSpeedPostM];
%         end
%         plotStat(recDataCuePre{i}.meanSpeedMean,recDataCueManip{i}.meanSpeedMean,...
%             recDataCuePost{i}.meanSpeedMean,meanRecDataCuePre{i}.meanSpeed,...
%             meanRecDataCueManip{i}.meanSpeed,meanRecDataCuePost{i}.meanSpeed,...
%             semRecDataCuePre{i}.meanSpeed,semRecDataCueManip{i}.meanSpeed,...
%             semRecDataCuePost{i}.meanSpeed,...
%             rankXY,rankXZ,rankYZ,'Session','Mean speed (mm/s)',['Cond' num2str(i)],path);
%             
%         path = [folderPath '_maxRunLenT_cond' num2str(i) '_Cue'];
%         rankXY = [rankRecDataCuePreM{i}.maxRunLenT;anovaRecDataCue{i}.maxRunLenTPreM];
%         if(isempty(rankRecDataCuePrePost{i}))
%             rankXZ = [];
%             rankYZ = [];
%         else
%             rankXZ = [rankRecDataCuePrePost{i}.maxRunLenT;anovaRecDataCue{i}.maxRunLenTPrePost];
%             rankYZ = [rankRecDataCuePostM{i}.maxRunLenT;anovaRecDataCue{i}.maxRunLenTPostM];
%         end
%         plotStat(recDataCuePre{i}.maxRunLenTMean,recDataCueManip{i}.maxRunLenTMean,...
%             recDataCuePost{i}.maxRunLenTMean,meanRecDataCuePre{i}.maxRunLenT,...
%             meanRecDataCueManip{i}.maxRunLenT,meanRecDataCuePost{i}.maxRunLenT,...
%             semRecDataCuePre{i}.maxRunLenT,semRecDataCueManip{i}.maxRunLenT,...
%             semRecDataCuePost{i}.maxRunLenT,...
%             rankXY,rankXZ,rankYZ,'Session','Longest running bout (s)',['Cond' num2str(i)],path);
%         
%         path = [folderPath '_totRunLenT_cond' num2str(i) '_Cue'];
%         rankXY = [rankRecDataCuePreM{i}.totRunLenT;anovaRecDataCue{i}.totRunLenTPreM];
%         if(isempty(rankRecDataCuePrePost{i}))
%             rankXZ = [];
%             rankYZ = [];
%         else
%             rankXZ = [rankRecDataCuePrePost{i}.totRunLenT;anovaRecDataCue{i}.totRunLenTPrePost];
%             rankYZ = [rankRecDataCuePostM{i}.totRunLenT;anovaRecDataCue{i}.totRunLenTPostM];
%         end
%         plotStat(recDataCuePre{i}.totRunLenTMean,recDataCueManip{i}.totRunLenTMean,...
%             recDataCuePost{i}.totRunLenTMean,meanRecDataCuePre{i}.totRunLenT,...
%             meanRecDataCueManip{i}.totRunLenT,meanRecDataCuePost{i}.totRunLenT,...
%             semRecDataCuePre{i}.totRunLenT,semRecDataCueManip{i}.totRunLenT,...
%             semRecDataCuePost{i}.totRunLenT,...
%             rankXY,rankXZ,rankYZ,'Session','Total run time (s)',['Cond' num2str(i)],path);
%         
%         path = [folderPath '_numRun_cond' num2str(i) '_Cue'];
%         rankXY = [rankRecDataCuePreM{i}.numRun;anovaRecDataCue{i}.numRunPreM];
%         if(isempty(rankRecDataCuePrePost{i}))
%             rankXZ = [];
%             rankYZ = [];
%         else
%             rankXZ = [rankRecDataCuePrePost{i}.numRun;anovaRecDataCue{i}.numRunPrePost];
%             rankYZ = [rankRecDataCuePostM{i}.numRun;anovaRecDataCue{i}.numRunPostM];
%         end
%         plotStat(recDataCuePre{i}.numRunMean,recDataCueManip{i}.numRunMean,...
%             recDataCuePost{i}.numRunMean,meanRecDataCuePre{i}.numRun,...
%             meanRecDataCueManip{i}.numRun,meanRecDataCuePost{i}.numRun,...
%             semRecDataCuePre{i}.numRun,semRecDataCueManip{i}.numRun,...
%             semRecDataCuePost{i}.numRun,...
%             rankXY,rankXZ,rankYZ,'Session','No. running bouts',['Cond' num2str(i)],path);
%         
%         path = [folderPath '_maxAcc_cond' num2str(i) '_Cue'];
%         rankXY = [rankRecDataCuePreM{i}.maxAcc;anovaRecDataCue{i}.maxAccPreM];
%         if(isempty(rankRecDataCuePrePost{i}))
%             rankXZ = [];
%             rankYZ = [];
%         else
%             rankXZ = [rankRecDataCuePrePost{i}.maxAcc;anovaRecDataCue{i}.maxAccPrePost];
%             rankYZ = [rankRecDataCuePostM{i}.maxAcc;anovaRecDataCue{i}.maxAccPostM];
%         end
%         plotStat(recDataCuePre{i}.maxAccMean,recDataCueManip{i}.maxAccMean,...
%             recDataCuePost{i}.maxAccMean,meanRecDataCuePre{i}.maxAcc,...
%             meanRecDataCueManip{i}.maxAcc,meanRecDataCuePost{i}.maxAcc,...
%             semRecDataCuePre{i}.maxAcc,semRecDataCueManip{i}.maxAcc,...
%             semRecDataCuePost{i}.maxAcc,...
%             rankXY,rankXZ,rankYZ,'Session','Max acceleration (mm/s^2)',['Cond' num2str(i)],path);
%         
%         path = [folderPath '_meanAcc_cond' num2str(i) '_Cue'];
%         rankXY = [rankRecDataCuePreM{i}.meanAcc;anovaRecDataCue{i}.meanAccPreM];
%         if(isempty(rankRecDataCuePrePost{i}))
%             rankXZ = [];
%             rankYZ = [];
%         else
%             rankXZ = [rankRecDataCuePrePost{i}.meanAcc;anovaRecDataCue{i}.meanAccPrePost];
%             rankYZ = [rankRecDataCuePostM{i}.meanAcc;anovaRecDataCue{i}.meanAccPostM];
%         end
%         plotStat(recDataCuePre{i}.meanAccMean,recDataCueManip{i}.meanAccMean,...
%             recDataCuePost{i}.meanAccMean,meanRecDataCuePre{i}.meanAcc,...
%             meanRecDataCueManip{i}.meanAcc,meanRecDataCuePost{i}.meanAcc,...
%             semRecDataCuePre{i}.meanAcc,semRecDataCueManip{i}.meanAcc,...
%             semRecDataCuePost{i}.meanAcc,...
%             rankXY,rankXZ,rankYZ,'Session','Mean acceleration (mm/s^2)',['Cond' num2str(i)],path);
        
        path = [folderPath '_totStopLenT_cond' num2str(i) '_Cue'];
        rankXY = [rankRecDataCuePreM{i}.totStopLenT;anovaRecDataCue{i}.totStopLenTPreM];
        if(isempty(rankRecDataCuePrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataCuePrePost{i}.totStopLenT;anovaRecDataCue{i}.totStopLenTPrePost];
            rankYZ = [rankRecDataCuePostM{i}.totStopLenT;anovaRecDataCue{i}.totStopLenTPostM];
        end
        plotStat(recDataCuePre{i}.totStopLenTMean,recDataCueManip{i}.totStopLenTMean,...
            recDataCuePost{i}.totStopLenTMean,meanRecDataCuePre{i}.totStopLenT,...
            meanRecDataCueManip{i}.totStopLenT,meanRecDataCuePost{i}.totStopLenT,...
            semRecDataCuePre{i}.totStopLenT,semRecDataCueManip{i}.totStopLenT,...
            semRecDataCuePost{i}.totStopLenT,...
            rankXY,rankXZ,rankYZ,'Session','Total stop time (s)',['Cond' num2str(i)],path);
        
%         path = [folderPath '_med1stFiveLickDist_cond' num2str(i) '_Cue'];
%         rankXY = [rankRecDataCuePreM{i}.med1stFiveLickDist;anovaRecDataCue{i}.med1stFiveLickDistPreM];
%         if(isempty(rankRecDataCuePrePost{i}))
%             rankXZ = [];
%             rankYZ = [];
%         else
%             rankXZ = [rankRecDataCuePrePost{i}.med1stFiveLickDist;anovaRecDataCue{i}.med1stFiveLickDistPrePost];
%             rankYZ = [rankRecDataCuePostM{i}.med1stFiveLickDist;anovaRecDataCue{i}.med1stFiveLickDistPostM];
%         end
%         plotStat(recDataCuePre{i}.med1stFiveLickDistMean,recDataCueManip{i}.med1stFiveLickDistMean,...
%             recDataCuePost{i}.med1stFiveLickDistMean,meanRecDataCuePre{i}.med1stFiveLickDist,...
%             meanRecDataCueManip{i}.med1stFiveLickDist,meanRecDataCuePost{i}.med1stFiveLickDist,...
%             semRecDataCuePre{i}.med1stFiveLickDist,semRecDataCueManip{i}.med1stFiveLickDist,...
%             semRecDataCuePost{i}.med1stFiveLickDist,...
%             rankXY,rankXZ,rankYZ,'Session','Median first-5-lick distance (mm)',['Cond' num2str(i)],path);
%         
%         path = [folderPath '_medLickDist_cond' num2str(i) '_Cue'];
%         rankXY = [rankRecDataCuePreM{i}.medLickDist;anovaRecDataCue{i}.medLickDistPreM];
%         if(isempty(rankRecDataCuePrePost{i}))
%             rankXZ = [];
%             rankYZ = [];
%         else
%             rankXZ = [rankRecDataCuePrePost{i}.medLickDist;anovaRecDataCue{i}.medLickDistPrePost];
%             rankYZ = [rankRecDataCuePostM{i}.medLickDist;anovaRecDataCue{i}.medLickDistPostM];
%         end
%         plotStat(recDataCuePre{i}.medLickDistMean,recDataCueManip{i}.medLickDistMean,...
%             recDataCuePost{i}.medLickDistMean,meanRecDataCuePre{i}.medLickDist,...
%             meanRecDataCueManip{i}.medLickDist,meanRecDataCuePost{i}.medLickDist,...
%             semRecDataCuePre{i}.medLickDist,semRecDataCueManip{i}.medLickDist,...
%             semRecDataCuePost{i}.medLickDist,...
%             rankXY,rankXZ,rankYZ,'Session','Median lick distance (mm)',['Cond' num2str(i)],path);
% 
% %         pause;
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