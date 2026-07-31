function plotAccumulateRecDataRew_optovsCtrl(mode)
% plot the comparison between optogenetic trials and the control trials
% during the optogenetic session
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
    load(fullPath,'recDataRewManipOpt','recDataRewManipOptCtrl');
    
    fullPath = [folderPath 'allRecDataStats.mat'];
    load(fullPath,'meanRecDataRewManipOpt','semRecDataRewManipOpt',...
            'rankRecDataRewOptCtrl','meanRecDataRewManipOptCtrl','semRecDataRewManipOptCtrl');
            
    nCond = length(meanRecDataRewManipOpt);
    for i = 1:nCond
        disp(['Condition ' num2str(i)])
        if(length(recDataRewManipOpt{i}.indRec)<=1)
            continue;
        end
        path = [folderPath 'Figure\_numSamplesRew_cond' num2str(i) '_OptVsOptCtrl'];
        plotStat(recDataRewManipOpt{i}.numSamplesMean,recDataRewManipOptCtrl{i}.numSamplesMean,[],...
            meanRecDataRewManipOpt{i}.numSamples,meanRecDataRewManipOptCtrl{i}.numSamples,[],...
            semRecDataRewManipOpt{i}.numSamples,semRecDataRewManipOptCtrl{i}.numSamples,[],...
            rankRecDataRewOptCtrl{i}.numSamples,[],[],'Session','No. samples',['Cond' num2str(i)],path);
        
        path = [folderPath 'Figure\_maxSpeedRew_cond' num2str(i) '_OptVsOptCtrl'];
        plotStat(recDataRewManipOpt{i}.maxSpeedMean,recDataRewManipOptCtrl{i}.maxSpeedMean,[],...
            meanRecDataRewManipOpt{i}.maxSpeed,meanRecDataRewManipOptCtrl{i}.maxSpeed,[],...
            semRecDataRewManipOpt{i}.maxSpeed,semRecDataRewManipOptCtrl{i}.maxSpeed,[],...
            rankRecDataRewOptCtrl{i}.maxSpeed,[],[],'Session','Max speed (mm/s)',['Cond' num2str(i)],path);
                
        path = [folderPath 'Figure\_meanSpeedRew_cond' num2str(i) '_OptVsOptCtrl'];
        plotStat(recDataRewManipOpt{i}.meanSpeedMean,recDataRewManipOptCtrl{i}.meanSpeedMean,[],...
            meanRecDataRewManipOpt{i}.meanSpeed,meanRecDataRewManipOptCtrl{i}.meanSpeed,[],...
            semRecDataRewManipOpt{i}.meanSpeed,semRecDataRewManipOptCtrl{i}.meanSpeed,[],...
            rankRecDataRewOptCtrl{i}.meanSpeed,[],[],'Session','Mean speed (mm/s)',['Cond' num2str(i)],path);
                    
        path = [folderPath 'Figure\_maxRunLenTRew_cond' num2str(i) '_OptVsOptCtrl'];
        plotStat(recDataRewManipOpt{i}.maxRunLenTMean,recDataRewManipOptCtrl{i}.maxRunLenTMean,[],...
            meanRecDataRewManipOpt{i}.maxRunLenT,meanRecDataRewManipOptCtrl{i}.maxRunLenT,[],...
            semRecDataRewManipOpt{i}.maxRunLenT,semRecDataRewManipOptCtrl{i}.maxRunLenT,[],...
            rankRecDataRewOptCtrl{i}.maxRunLenT,[],[],'Session','Longest running bout (s)',['Cond' num2str(i)],path);
                
        path = [folderPath 'Figure\_totRunLenTRew_cond' num2str(i) '_OptVsOptCtrl'];
        plotStat(recDataRewManipOpt{i}.totRunLenTMean,recDataRewManipOptCtrl{i}.totRunLenTMean,[],...
            meanRecDataRewManipOpt{i}.totRunLenT,meanRecDataRewManipOptCtrl{i}.totRunLenT,[],...
            semRecDataRewManipOpt{i}.totRunLenT,semRecDataRewManipOptCtrl{i}.totRunLenT,[],...
            rankRecDataRewOptCtrl{i}.totRunLenT,[],[],'Session','Total run time (s)',['Cond' num2str(i)],path);
                
        path = [folderPath 'Figure\_numRunRew_cond' num2str(i) '_OptVsOptCtrl'];
        plotStat(recDataRewManipOpt{i}.numRunMean,recDataRewManipOptCtrl{i}.numRunMean,[],...
            meanRecDataRewManipOpt{i}.numRun,meanRecDataRewManipOptCtrl{i}.numRun,[],...
            semRecDataRewManipOpt{i}.numRun,semRecDataRewManipOptCtrl{i}.numRun,[],...
            rankRecDataRewOptCtrl{i}.numRun,[],[],'Session','No. running bouts',['Cond' num2str(i)],path);
                
        path = [folderPath 'Figure\_maxAccRew_cond' num2str(i) '_OptVsOptCtrl'];
        plotStat(recDataRewManipOpt{i}.maxAccMean,recDataRewManipOptCtrl{i}.maxAccMean,[],...
            meanRecDataRewManipOpt{i}.maxAcc,meanRecDataRewManipOptCtrl{i}.maxAcc,[],...
            semRecDataRewManipOpt{i}.maxAcc,semRecDataRewManipOptCtrl{i}.maxAcc,[],...
            rankRecDataRewOptCtrl{i}.maxAcc,[],[],'Session','Max acceleration (mm/s^2)',['Cond' num2str(i)],path);
        
        path = [folderPath 'Figure\_meanAccRew_cond' num2str(i) '_OptVsOptCtrl'];
        plotStat(recDataRewManipOpt{i}.meanAccMean,recDataRewManipOptCtrl{i}.meanAccMean,[],...
            meanRecDataRewManipOpt{i}.meanAcc,meanRecDataRewManipOptCtrl{i}.meanAcc,[],...
            semRecDataRewManipOpt{i}.meanAcc,semRecDataRewManipOptCtrl{i}.meanAcc,[],...
            rankRecDataRewOptCtrl{i}.meanAcc,[],[],'Session','Mean acceleration (mm/s^2)',['Cond' num2str(i)],path);
        
        path = [folderPath 'Figure\_totStopLenTRew_cond' num2str(i) '_OptVsOptCtrl'];
        plotStat(recDataRewManipOpt{i}.totStopLenTMean,recDataRewManipOptCtrl{i}.totStopLenTMean,[],...
            meanRecDataRewManipOpt{i}.totStopLenT,meanRecDataRewManipOptCtrl{i}.totStopLenT,[],...
            semRecDataRewManipOpt{i}.totStopLenT,semRecDataRewManipOptCtrl{i}.totStopLenT,[],...
            rankRecDataRewOptCtrl{i}.totStopLenT,[],[],'Session','Total stop time (s)',['Cond' num2str(i)],path);
                       
        path = [folderPath 'Figure\_speedSimMeanRew_cond' num2str(i) '_OptVsOptCtrl'];
        plotStat(recDataRewManipOpt{i}.speedSimMean,recDataRewManipOptCtrl{i}.speedSimMean,[],...
            meanRecDataRewManipOpt{i}.speedSimMean,meanRecDataRewManipOptCtrl{i}.speedSimMean,[],...
            semRecDataRewManipOpt{i}.speedSimMean,semRecDataRewManipOptCtrl{i}.speedSimMean,[],...
            rankRecDataRewOptCtrl{i}.speedSimMean,[],[],'Session','speed Sim Mean',['Cond' num2str(i)],path);
        
        path = [folderPath 'Figure\_lickSimMeanRew_cond' num2str(i) '_OptVsOptCtrl'];
        plotStat(recDataRewManipOpt{i}.lickSimMean,recDataRewManipOptCtrl{i}.lickSimMean,[],...
            meanRecDataRewManipOpt{i}.lickSimMean,meanRecDataRewManipOptCtrl{i}.lickSimMean,[],...
            semRecDataRewManipOpt{i}.lickSimMean,semRecDataRewManipOptCtrl{i}.lickSimMean,[],...
            rankRecDataRewOptCtrl{i}.lickSimMean,[],[],'Session','lick Sim Mean',['Cond' num2str(i)],path);
                
        path = [folderPath 'Figure\_speedEucMeanRew_cond' num2str(i) '_OptVsOptCtrl'];
        plotStat(recDataRewManipOpt{i}.speedEucMean,recDataRewManipOptCtrl{i}.speedEucMean,[],...
            meanRecDataRewManipOpt{i}.speedEucMean,meanRecDataRewManipOptCtrl{i}.speedEucMean,[],...
            semRecDataRewManipOpt{i}.speedEucMean,semRecDataRewManipOptCtrl{i}.speedEucMean,[],...
            rankRecDataRewOptCtrl{i}.speedEucMean,[],[],'Session','speed Euc Mean',['Cond' num2str(i)],path);
        
        path = [folderPath 'Figure\_lickEucMeanRew_cond' num2str(i) '_OptVsOptCtrl'];
        plotStat(recDataRewManipOpt{i}.lickEucMean,recDataRewManipOptCtrl{i}.lickEucMean,[],...
            meanRecDataRewManipOpt{i}.lickEucMean,meanRecDataRewManipOptCtrl{i}.lickEucMean,[],...
            semRecDataRewManipOpt{i}.lickEucMean,semRecDataRewManipOptCtrl{i}.lickEucMean,[],...
            rankRecDataRewOptCtrl{i}.lickEucMean,[],[],'Session','lick Euc Mean',['Cond' num2str(i)],path);
        
%         pause;
        close all;
    end
end

function plotBoxPlot(x1,x2,yl,fn,pathAnal,ylimit,p12,colorSel)
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
    x = [x1;x2];
    g = [repmat(1,length(x1),1);...
        repmat(2,length(x2),1)];
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
        
        barPlotWithStat(1:2,meanArr,errBar,dataArr,xlab,ylab,ti,rankXY,[],[]);       
    end
    print('-painters', '-dpdf', path, '-r600')
    savefig([path '.fig']);
end