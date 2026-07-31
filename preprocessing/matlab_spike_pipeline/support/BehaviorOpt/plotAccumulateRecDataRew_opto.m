function plotAccumulateRecDataRew_opto(mode,opt)
    if(mode == 1)
        RecordingListNT_opto1; % including the continous stimulation
    elseif(mode == 2)
        RecordingListNT_opto2; % without the continous stimulation
    else
        RecordingListNT_opto;
    end
    
    fullPath = [folderPath 'allRecData.mat'];
    if(~exist(fullPath,'file'))
        disp([fullPath ' does not exist.']);
        return;
    end
    load(fullPath,'recDataRewPre','recDataRewManip','recDataRewPost');
    if(opt == 1)
        load(fullPath,'recDataRewManipOpt');
        recDataRewManip = recDataRewManipOpt;
    elseif(opt == 2)
        load(fullPath,'recDataRewManipOptCtrl');
        recDataRewManip = recDataRewManipOptCtrl;
    end
    
    fullPath = [folderPath 'allRecDataStats.mat'];
    if(opt == 0)        
        load(fullPath,'meanRecDataRewPre','meanRecDataRewManip','meanRecDataRewPost',...
            'semRecDataRewPre','semRecDataRewManip','semRecDataRewPost',...
            'rankRecDataRewPrePost','rankRecDataRewPreM','rankRecDataRewPostM',...
            'anovaRecDataRew');
    elseif(opt == 1)
        load(fullPath,'meanRecDataRewPreOpt','meanRecDataRewManipOpt','meanRecDataRewPostOpt',...
            'semRecDataRewPreOpt','semRecDataRewManipOpt','semRecDataRewPostOpt',...
            'rankRecDataRewPrePostOpt','rankRecDataRewPreMOpt','rankRecDataRewPostMOpt',...
            'rankRecDataRewOptCtrl','anovaRecDataRewOpt');
        meanRecDataRewPre = meanRecDataRewPreOpt;
        meanRecDataRewManip = meanRecDataRewManipOpt;
        meanRecDataRewPost = meanRecDataRewPostOpt;
        semRecDataRewPre = semRecDataRewPreOpt;
        semRecDataRewManip = semRecDataRewManipOpt;
        semRecDataRewPost = semRecDataRewPostOpt;
        rankRecDataRewPrePost = rankRecDataRewPrePostOpt;
        rankRecDataRewPreM = rankRecDataRewPreMOpt;
        rankRecDataRewPostM = rankRecDataRewPostMOpt;
        anovaRecDataRew = anovaRecDataRewOpt;
    else
        load(fullPath,'meanRecDataRewPreOptCtrl','meanRecDataRewManipOptCtrl','meanRecDataRewPostOptCtrl',...
            'semRecDataRewPreOptCtrl','semRecDataRewManipOptCtrl','semRecDataRewPostOptCtrl',...
            'rankRecDataRewPrePostOptCtrl','rankRecDataRewPreMOptCtrl','rankRecDataRewPostMOptCtrl',...
            'anovaRecDataRewOptCtrl');
        meanRecDataRewPre = meanRecDataRewPreOptCtrl;
        meanRecDataRewManip = meanRecDataRewManipOptCtrl;
        meanRecDataRewPost = meanRecDataRewPostOptCtrl;
        semRecDataRewPre = semRecDataRewPreOptCtrl;
        semRecDataRewManip = semRecDataRewManipOptCtrl;
        semRecDataRewPost = semRecDataRewPostOptCtrl;
        rankRecDataRewPrePost = rankRecDataRewPrePostOptCtrl;
        rankRecDataRewPreM = rankRecDataRewPreMOptCtrl;
        rankRecDataRewPostM = rankRecDataRewPostMOptCtrl;
        anovaRecDataRew = anovaRecDataRewOptCtrl;
    end
    
    nCond = length(meanRecDataRewPre);
    for i = 1:nCond
        disp(['Condition ' num2str(i)])
        if(length(recDataRewPre{i}.indRec)<=1)
            continue;
        end
        path = [folderPath 'Figure\_numSamplesRew_cond' num2str(i) '_opt' num2str(opt)];
        rankXY = [rankRecDataRewPreM{i}.numSamples;anovaRecDataRew{i}.numSamplesPreM];
        if(isempty(rankRecDataRewPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRewPrePost{i}.numSamples;anovaRecDataRew{i}.numSamplesPrePost];
            rankYZ = [rankRecDataRewPostM{i}.numSamples;anovaRecDataRew{i}.numSamplesPostM];
        end
        plotStat(recDataRewPre{i}.numSamplesMean,recDataRewManip{i}.numSamplesMean,...
            recDataRewPost{i}.numSamplesMean,meanRecDataRewPre{i}.numSamples,...
            meanRecDataRewManip{i}.numSamples,meanRecDataRewPost{i}.numSamples,...
            semRecDataRewPre{i}.numSamples,semRecDataRewManip{i}.numSamples,...
            semRecDataRewPost{i}.numSamples,...
            rankXY,rankXZ,rankYZ,'Session','No. samples (Rew)',['Cond' num2str(i)],path);
        
        path = [folderPath 'Figure\_maxSpeedRew_cond' num2str(i) '_opt' num2str(opt)];
        rankXY = [rankRecDataRewPreM{i}.maxSpeed;anovaRecDataRew{i}.maxSpeedPreM];
        if(isempty(rankRecDataRewPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRewPrePost{i}.maxSpeed;anovaRecDataRew{i}.maxSpeedPrePost];
            rankYZ = [rankRecDataRewPostM{i}.maxSpeed;anovaRecDataRew{i}.maxSpeedPostM];
        end
        plotStat(recDataRewPre{i}.maxSpeedMean,recDataRewManip{i}.maxSpeedMean,...
            recDataRewPost{i}.maxSpeedMean,meanRecDataRewPre{i}.maxSpeed,...
            meanRecDataRewManip{i}.maxSpeed,meanRecDataRewPost{i}.maxSpeed,...
            semRecDataRewPre{i}.maxSpeed,semRecDataRewManip{i}.maxSpeed,...
            semRecDataRewPost{i}.maxSpeed,...
            rankXY,rankXZ,rankYZ,'Session','Max speed (mm/s) (Rew)',['Cond' num2str(i)],path);
        
        path = [folderPath 'Figure\_meanSpeedRew_cond' num2str(i) '_opt' num2str(opt)];
        rankXY = [rankRecDataRewPreM{i}.meanSpeed;anovaRecDataRew{i}.meanSpeedPreM];
        if(isempty(rankRecDataRewPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRewPrePost{i}.meanSpeed;anovaRecDataRew{i}.meanSpeedPrePost];
            rankYZ = [rankRecDataRewPostM{i}.meanSpeed;anovaRecDataRew{i}.meanSpeedPostM];
        end
        plotStat(recDataRewPre{i}.meanSpeedMean,recDataRewManip{i}.meanSpeedMean,...
            recDataRewPost{i}.meanSpeedMean,meanRecDataRewPre{i}.meanSpeed,...
            meanRecDataRewManip{i}.meanSpeed,meanRecDataRewPost{i}.meanSpeed,...
            semRecDataRewPre{i}.meanSpeed,semRecDataRewManip{i}.meanSpeed,...
            semRecDataRewPost{i}.meanSpeed,...
            rankXY,rankXZ,rankYZ,'Session','Mean speed (mm/s) (Rew)',['Cond' num2str(i)],path);
            
        path = [folderPath 'Figure\_maxRunLenTRew_cond' num2str(i) '_opt' num2str(opt)];
        rankXY = [rankRecDataRewPreM{i}.maxRunLenT;anovaRecDataRew{i}.maxRunLenTPreM];
        if(isempty(rankRecDataRewPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRewPrePost{i}.maxRunLenT;anovaRecDataRew{i}.maxRunLenTPrePost];
            rankYZ = [rankRecDataRewPostM{i}.maxRunLenT;anovaRecDataRew{i}.maxRunLenTPostM];
        end
        plotStat(recDataRewPre{i}.maxRunLenTMean,recDataRewManip{i}.maxRunLenTMean,...
            recDataRewPost{i}.maxRunLenTMean,meanRecDataRewPre{i}.maxRunLenT,...
            meanRecDataRewManip{i}.maxRunLenT,meanRecDataRewPost{i}.maxRunLenT,...
            semRecDataRewPre{i}.maxRunLenT,semRecDataRewManip{i}.maxRunLenT,...
            semRecDataRewPost{i}.maxRunLenT,...
            rankXY,rankXZ,rankYZ,'Session','Longest running bout (s) (Rew)',['Cond' num2str(i)],path);
        
        path = [folderPath 'Figure\_totRunLenTRew_cond' num2str(i) '_opt' num2str(opt)];
        rankXY = [rankRecDataRewPreM{i}.totRunLenT;anovaRecDataRew{i}.totRunLenTPreM];
        if(isempty(rankRecDataRewPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRewPrePost{i}.totRunLenT;anovaRecDataRew{i}.totRunLenTPrePost];
            rankYZ = [rankRecDataRewPostM{i}.totRunLenT;anovaRecDataRew{i}.totRunLenTPostM];
        end
        plotStat(recDataRewPre{i}.totRunLenTMean,recDataRewManip{i}.totRunLenTMean,...
            recDataRewPost{i}.totRunLenTMean,meanRecDataRewPre{i}.totRunLenT,...
            meanRecDataRewManip{i}.totRunLenT,meanRecDataRewPost{i}.totRunLenT,...
            semRecDataRewPre{i}.totRunLenT,semRecDataRewManip{i}.totRunLenT,...
            semRecDataRewPost{i}.totRunLenT,...
            rankXY,rankXZ,rankYZ,'Session','Total run time (s) (Rew)',['Cond' num2str(i)],path);
        
        path = [folderPath 'Figure\_numRunRew_cond' num2str(i) '_opt' num2str(opt)];
        rankXY = [rankRecDataRewPreM{i}.numRun;anovaRecDataRew{i}.numRunPreM];
        if(isempty(rankRecDataRewPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRewPrePost{i}.numRun;anovaRecDataRew{i}.numRunPrePost];
            rankYZ = [rankRecDataRewPostM{i}.numRun;anovaRecDataRew{i}.numRunPostM];
        end
        plotStat(recDataRewPre{i}.numRunMean,recDataRewManip{i}.numRunMean,...
            recDataRewPost{i}.numRunMean,meanRecDataRewPre{i}.numRun,...
            meanRecDataRewManip{i}.numRun,meanRecDataRewPost{i}.numRun,...
            semRecDataRewPre{i}.numRun,semRecDataRewManip{i}.numRun,...
            semRecDataRewPost{i}.numRun,...
            rankXY,rankXZ,rankYZ,'Session','No. running bouts (Rew)',['Cond' num2str(i)],path);
        
        path = [folderPath 'Figure\_maxAccRew_cond' num2str(i) '_opt' num2str(opt)];
        rankXY = [rankRecDataRewPreM{i}.maxAcc;anovaRecDataRew{i}.maxAccPreM];
        if(isempty(rankRecDataRewPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRewPrePost{i}.maxAcc;anovaRecDataRew{i}.maxAccPrePost];
            rankYZ = [rankRecDataRewPostM{i}.maxAcc;anovaRecDataRew{i}.maxAccPostM];
        end
        plotStat(recDataRewPre{i}.maxAccMean,recDataRewManip{i}.maxAccMean,...
            recDataRewPost{i}.maxAccMean,meanRecDataRewPre{i}.maxAcc,...
            meanRecDataRewManip{i}.maxAcc,meanRecDataRewPost{i}.maxAcc,...
            semRecDataRewPre{i}.maxAcc,semRecDataRewManip{i}.maxAcc,...
            semRecDataRewPost{i}.maxAcc,...
            rankXY,rankXZ,rankYZ,'Session','Max acceleration (mm/s^2) (Rew)',['Cond' num2str(i)],path);
        
        path = [folderPath 'Figure\_meanAccRew_cond' num2str(i) '_opt' num2str(opt)];
        rankXY = [rankRecDataRewPreM{i}.meanAcc;anovaRecDataRew{i}.meanAccPreM];
        if(isempty(rankRecDataRewPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRewPrePost{i}.meanAcc;anovaRecDataRew{i}.meanAccPrePost];
            rankYZ = [rankRecDataRewPostM{i}.meanAcc;anovaRecDataRew{i}.meanAccPostM];
        end
        plotStat(recDataRewPre{i}.meanAccMean,recDataRewManip{i}.meanAccMean,...
            recDataRewPost{i}.meanAccMean,meanRecDataRewPre{i}.meanAcc,...
            meanRecDataRewManip{i}.meanAcc,meanRecDataRewPost{i}.meanAcc,...
            semRecDataRewPre{i}.meanAcc,semRecDataRewManip{i}.meanAcc,...
            semRecDataRewPost{i}.meanAcc,...
            rankXY,rankXZ,rankYZ,'Session','Mean acceleration (mm/s^2) (Rew)',['Cond' num2str(i)],path);
        
        path = [folderPath 'Figure\_totStopLenTRew_cond' num2str(i) '_opt' num2str(opt)];
        rankXY = [rankRecDataRewPreM{i}.totStopLenT;anovaRecDataRew{i}.totStopLenTPreM];
        if(isempty(rankRecDataRewPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRewPrePost{i}.totStopLenT;anovaRecDataRew{i}.totStopLenTPrePost];
            rankYZ = [rankRecDataRewPostM{i}.totStopLenT;anovaRecDataRew{i}.totStopLenTPostM];
        end
        plotStat(recDataRewPre{i}.totStopLenTMean,recDataRewManip{i}.totStopLenTMean,...
            recDataRewPost{i}.totStopLenTMean,meanRecDataRewPre{i}.totStopLenT,...
            meanRecDataRewManip{i}.totStopLenT,meanRecDataRewPost{i}.totStopLenT,...
            semRecDataRewPre{i}.totStopLenT,semRecDataRewManip{i}.totStopLenT,...
            semRecDataRewPost{i}.totStopLenT,...
            rankXY,rankXZ,rankYZ,'Session','Total stop time (s) (Rew)',['Cond' num2str(i)],path);
        
       path = [folderPath 'Figure\_speedSimMeanRew_cond' num2str(i) '_opt' num2str(opt)];
        rankXY = [rankRecDataRewPreM{i}.speedSimMean;anovaRecDataRew{i}.speedSimMeanPreM];
        if(isempty(rankRecDataRewPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRewPrePost{i}.speedSimMean;anovaRecDataRew{i}.speedSimMeanPrePost];
            rankYZ = [rankRecDataRewPostM{i}.speedSimMean;anovaRecDataRew{i}.speedSimMeanPostM];
        end        
        plotStat(recDataRewPre{i}.speedSimMean,recDataRewManip{i}.speedSimMean,...
            recDataRewPost{i}.speedSimMean,meanRecDataRewPre{i}.speedSimMean,...
            meanRecDataRewManip{i}.speedSimMean,meanRecDataRewPost{i}.speedSimMean,...
            semRecDataRewPre{i}.speedSimMean,semRecDataRewManip{i}.speedSimMean,...
            semRecDataRewPost{i}.speedSimMean,...
            rankXY,rankXZ,rankYZ,'Session','speed Sim Mean',['Cond' num2str(i)],path);
        
        path = [folderPath 'Figure\_lickSimMeanRew_cond' num2str(i) '_opt' num2str(opt)];
        rankXY = [rankRecDataRewPreM{i}.lickSimMean;anovaRecDataRew{i}.lickSimMeanPreM];
        if(isempty(rankRecDataRewPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRewPrePost{i}.lickSimMean;anovaRecDataRew{i}.lickSimMeanPrePost];
            rankYZ = [rankRecDataRewPostM{i}.lickSimMean;anovaRecDataRew{i}.lickSimMeanPostM];
        end        
        plotStat(recDataRewPre{i}.lickSimMean,recDataRewManip{i}.lickSimMean,...
            recDataRewPost{i}.lickSimMean,meanRecDataRewPre{i}.lickSimMean,...
            meanRecDataRewManip{i}.lickSimMean,meanRecDataRewPost{i}.lickSimMean,...
            semRecDataRewPre{i}.lickSimMean,semRecDataRewManip{i}.lickSimMean,...
            semRecDataRewPost{i}.lickSimMean,...
            rankXY,rankXZ,rankYZ,'Session','lick Sim Mean',['Cond' num2str(i)],path);

        path = [folderPath 'Figure\_speedEucMeanRew_cond' num2str(i) '_opt' num2str(opt)];
        rankXY = [rankRecDataRewPreM{i}.speedEucMean;anovaRecDataRew{i}.speedEucMeanPreM];
        if(isempty(rankRecDataRewPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRewPrePost{i}.speedEucMean;anovaRecDataRew{i}.speedEucMeanPrePost];
            rankYZ = [rankRecDataRewPostM{i}.speedEucMean;anovaRecDataRew{i}.speedEucMeanPostM];
        end        
        plotStat(recDataRewPre{i}.speedEucMean,recDataRewManip{i}.speedEucMean,...
            recDataRewPost{i}.speedEucMean,meanRecDataRewPre{i}.speedEucMean,...
            meanRecDataRewManip{i}.speedEucMean,meanRecDataRewPost{i}.speedEucMean,...
            semRecDataRewPre{i}.speedEucMean,semRecDataRewManip{i}.speedEucMean,...
            semRecDataRewPost{i}.speedEucMean,...
            rankXY,rankXZ,rankYZ,'Session','speed Euc Mean',['Cond' num2str(i)],path);
        
        path = [folderPath 'Figure\_lickEucMeanRew_cond' num2str(i) '_opt' num2str(opt)];
        rankXY = [rankRecDataRewPreM{i}.lickEucMean;anovaRecDataRew{i}.lickEucMeanPreM];
        if(isempty(rankRecDataRewPrePost{i}))
            rankXZ = [];
            rankYZ = [];
        else
            rankXZ = [rankRecDataRewPrePost{i}.lickEucMean;anovaRecDataRew{i}.lickEucMeanPrePost];
            rankYZ = [rankRecDataRewPostM{i}.lickEucMean;anovaRecDataRew{i}.lickEucMeanPostM];
        end        
        plotStat(recDataRewPre{i}.lickEucMean,recDataRewManip{i}.lickEucMean,...
            recDataRewPost{i}.lickEucMean,meanRecDataRewPre{i}.lickEucMean,...
            meanRecDataRewManip{i}.lickEucMean,meanRecDataRewPost{i}.lickEucMean,...
            semRecDataRewPre{i}.lickEucMean,semRecDataRewManip{i}.lickEucMean,...
            semRecDataRewPost{i}.lickEucMean,...
            rankXY,rankXZ,rankYZ,'Session','lick Euc Mean',['Cond' num2str(i)],path);
        
%         colorSel = 0;
%         
%         rankXY = rankRecDataRewPreM{i}.pRSMeanSpeedOverTimeRew0to1;
%         if(isempty(rankRecDataRewPrePost{i}))
%             rankXZ = [];
%             rankYZ = [];
%             meanSpeedPost = [];
%         else
%             rankXZ = rankRecDataRewPrePost{i}.pRSMeanSpeedOverTimeRew0to1;
%             rankYZ = rankRecDataRewPostM{i}.pRSMeanSpeedOverTimeRew0to1;
%             meanSpeedPost = recDataRewPost{i}.meanSpeedOverTimeRew0to1';
%         end
%         plotBoxPlot(recDataRewPre{i}.meanSpeedOverTimeRew0to1',...
%             meanSpeedPost,...
%             recDataRewManip{i}.meanSpeedOverTimeRew0to1',...
%             'Speed (cm/s) 0to1s (Rew)',['speedProfileRew0to1sBox' num2str(i) '_opt' num2str(opt)],...
%             [folderPath 'Figure\'],[-10 100],rankXY,rankXZ,rankYZ,colorSel);
%         
%        rankXY = rankRecDataRewPreM{i}.pRSMeanSpeedOverTimeRew1to2;
%         if(isempty(rankRecDataRewPrePost{i}))
%             rankXZ = [];
%             rankYZ = [];
%             meanSpeedPost = [];
%         else
%             rankXZ = rankRecDataRewPrePost{i}.pRSMeanSpeedOverTimeRew1to2;
%             rankYZ = rankRecDataRewPostM{i}.pRSMeanSpeedOverTimeRew1to2;
%             meanSpeedPost = recDataRewPost{i}.meanSpeedOverTimeRew1to2';
%         end
%         plotBoxPlot(recDataRewPre{i}.meanSpeedOverTimeRew1to2',...
%             meanSpeedPost,...
%             recDataRewManip{i}.meanSpeedOverTimeRew1to2',...
%             'Speed (cm/s) 1to2s (Rew)',['speedProfileRew1to2sBox' num2str(i) '_opt' num2str(opt)],...
%             [folderPath 'Figure\'],[-10 100],rankXY,rankXZ,rankYZ,colorSel); 
%         
%         rankXY = rankRecDataRewPreM{i}.pRSMeanSpeedOverTimeRew2to3;
%         if(isempty(rankRecDataRewPrePost{i}))
%             rankXZ = [];
%             rankYZ = [];
%             meanSpeedPost = [];
%         else
%             rankXZ = rankRecDataRewPrePost{i}.pRSMeanSpeedOverTimeRew2to3;
%             rankYZ = rankRecDataRewPostM{i}.pRSMeanSpeedOverTimeRew2to3;
%             meanSpeedPost = recDataRewPost{i}.meanSpeedOverTimeRew2to3';
%         end
%         plotBoxPlot(recDataRewPre{i}.meanSpeedOverTimeRew2to3',...
%             meanSpeedPost,...
%             recDataRewManip{i}.meanSpeedOverTimeRew2to3',...
%             'Speed (cm/s) 2to3s (Rew)',['speedProfileRew2to3sBox' num2str(i) '_opt' num2str(opt)],...
%             [folderPath 'Figure\'],[-10 100],rankXY,rankXZ,rankYZ,colorSel); 
%         
%         rankXY = rankRecDataRewPreM{i}.pRSMeanSpeedOverTimeRew3to5;
%         if(isempty(rankRecDataRewPrePost{i}))
%             rankXZ = [];
%             rankYZ = [];
%             meanSpeedPost = [];
%         else
%             rankXZ = rankRecDataRewPrePost{i}.pRSMeanSpeedOverTimeRew3to5;
%             rankYZ = rankRecDataRewPostM{i}.pRSMeanSpeedOverTimeRew3to5;
%             meanSpeedPost = recDataRewPost{i}.meanSpeedOverTimeRew3to5';
%         end
%         plotBoxPlot(recDataRewPre{i}.meanSpeedOverTimeRew3to5',...
%             meanSpeedPost,...
%             recDataRewManip{i}.meanSpeedOverTimeRew3to5',...
%             'Speed (cm/s) 3to5s (Rew)',['speedProfileRew3to5sBox' num2str(i) '_opt' num2str(opt)],...
%             [folderPath 'Figure\'],[10 100],rankXY,rankXZ,rankYZ,colorSel); 
        
%         pause;
        close all;
    end
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