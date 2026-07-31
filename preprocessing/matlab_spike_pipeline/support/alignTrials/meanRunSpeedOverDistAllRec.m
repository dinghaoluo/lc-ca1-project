function meanRunSpeedOverDistAllRec(onlyRun)
% calculate mean speed over all the recordings

    RecordingList;
    
    speedOverDistMeanAllArr = [];
    speedOverDistStdAllArr = [];
    speedOverDistAllArr = [];
    
    indTr = 5:55;
    nTrials = length(indTr);
    
    %% active licking
    numRec = size(listRecordingsActiveLickPath,1);
    totNTr = 0;
    for i = 1:numRec
        path = listRecordingsActiveLickPath(i,:);
        fileName = listRecordingsActiveLickFileName(i,:);
        
        fullPath = [path,fileName, '_runSpeedDist_Run' num2str(onlyRun) '.mat'];
        load(fullPath,'speedOverDistMean','speedOverDistStd','speedOverDist','param');
        
        speedOverDistMeanAllArr = [speedOverDistMeanAllArr; speedOverDistMean];
        speedOverDistStdAllArr = [speedOverDistStdAllArr; speedOverDistStd];
        if(size(speedOverDist,1) < max(indTr))
            speedOverDistAllArr = [speedOverDistAllArr; speedOverDist];
            totNTr = totNTr+size(speedOverDist,1);
        else
            speedOverDistAllArr = [speedOverDistAllArr; speedOverDist(indTr,:)]; 
            totNTr = totNTr+nTrials;
        end
    end
    
    ind = speedOverDistAllArr(:,1) > 3000 | speedOverDistAllArr(:,5) > 3000;
    speedOverDistMeanAll = mean(speedOverDistMeanAllArr);
    speedOverDistSEMAll = sqrt(sum(speedOverDistStdAllArr.^2)/numRec);
    speedOverDistAllArr = speedOverDistAllArr(ind == 0,:);
    speedOverDistMeanAll1 = mean(speedOverDistAllArr);
    speedOverDistSEMAll1 = std(speedOverDistAllArr)/sqrt(sum(ind == 0));
    
    fullPath = ['Z:\Yingxue\DataAnalysisRaphi\Behavior\ALRunSpeedDistMeanAllRec_Run' num2str(onlyRun) '.mat'];
    save(fullPath, 'speedOverDistMeanAllArr','speedOverDistStdAllArr',...
        'speedOverDistMeanAll','speedOverDistSEMAll','speedOverDistMeanAll1',...
        'speedOverDistSEMAll1');
    speedOverDistAllArrAL = speedOverDistAllArr;
    
    %% passive licking
    numRec = size(listRecordingsPassiveLickPath,1);
    totNTr = 0;
    speedOverDistAllArr = [];
    for i = 1:numRec
        path = listRecordingsPassiveLickPath(i,:);
        fileName = listRecordingsPassiveLickFileName(i,:);
        
        fullPath = [path,fileName, '_runSpeedDist_Run' num2str(onlyRun) '.mat'];
        load(fullPath,'speedOverDistMean','speedOverDistStd','speedOverDist');
        
        speedOverDistMeanAllArr = [speedOverDistMeanAllArr; speedOverDistMean];
        speedOverDistStdAllArr = [speedOverDistStdAllArr; speedOverDistStd];
        if(size(speedOverDist,1) < max(indTr))
            speedOverDistAllArr = [speedOverDistAllArr; speedOverDist];
            totNTr = totNTr+size(speedOverDist,1);
        else
            speedOverDistAllArr = [speedOverDistAllArr; speedOverDist(indTr,:)]; 
            totNTr = totNTr+nTrials;
        end
    end
    
    speedOverDistMeanAll = mean(speedOverDistMeanAllArr);
    speedOverDistSEMAll = sqrt(sum(speedOverDistStdAllArr.^2)/numRec);
    speedOverDistMeanAll1 = mean(speedOverDistAllArr);
    speedOverDistSEMAll1 = std(speedOverDistAllArr)/sqrt(totNTr);
    
    fullPath = ['Z:\Yingxue\DataAnalysisRaphi\Behavior\PLRunSpeedDistMeanAllRec_Run' num2str(onlyRun) '.mat'];
    save(fullPath, 'speedOverDistMeanAllArr','speedOverDistStdAllArr',...
        'speedOverDistMeanAll','speedOverDistSEMAll','speedOverDistMeanAll1',...
        'speedOverDistSEMAll1');
    
    options.handle     = figure;
    set(options.handle,'OuterPosition',...
        [500 500 280 280])
    options.color_areaX = [27 117 187]./255;    % Blue theme
    options.color_lineX = [ 39 169 225]./255;
    options.color_areaY = [187 189 192]./255;    % Orange theme
    options.color_lineY = [167 169  171]./255;
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'sem';
    options.x_axisX = param.spaceSteps/10;
    options.x_axisY = param.spaceSteps/10;
    plot_areaerrorbarXY(speedOverDistAllArrAL/10, speedOverDistAllArr/10,...
        options);
    set(gca,'YLim',[0 60],'XLim',[0 180])
    xlabel('Dist (cm)')
    ylabel('Speed (cm/s)')
    legend('','Active','','Passive')

    fileName1 = ['Z:\Yingxue\DataAnalysisRaphi\Behavior\' ...
        'SpeedVsDist-AL-PL'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    %% no cue
    numRec = size(listRecordingsNoCuePath,1);
    totNTr = 0;
    speedOverDistAllArr = [];
    for i = 1:numRec
        path = listRecordingsNoCuePath(i,:);
        fileName = listRecordingsNoCueFileName(i,:);
        
        fullPath = [path,fileName, '_runSpeedDist_Run' num2str(onlyRun) '.mat'];
        load(fullPath,'speedOverDistMean','speedOverDistStd','speedOverDist');
        
        speedOverDistMeanAllArr = [speedOverDistMeanAllArr; speedOverDistMean];
        speedOverDistStdAllArr = [speedOverDistStdAllArr; speedOverDistStd];
        if(size(speedOverDist,1) < max(indTr))
            speedOverDistAllArr = [speedOverDistAllArr; speedOverDist];
            totNTr = totNTr+size(speedOverDist,1);
        else
            speedOverDistAllArr = [speedOverDistAllArr; speedOverDist(indTr,:)]; 
            totNTr = totNTr+nTrials;
        end
    end
    
    ind = speedOverDistAllArr(:,5) > 3000 | speedOverDistAllArr(:,6) > 1000 ...
                | speedOverDistAllArr(:,7) > 3000;
    speedOverDistMeanAll = mean(speedOverDistMeanAllArr);
    speedOverDistSEMAll = sqrt(sum(speedOverDistStdAllArr.^2)/numRec);
    speedOverDistAllArr = speedOverDistAllArr(ind == 0,:);
    speedOverDistMeanAll1 = mean(speedOverDistAllArr);
    speedOverDistSEMAll1 = std(speedOverDistAllArr)/sqrt(sum(ind == 0));
    
    fullPath = ['Z:\Yingxue\DataAnalysisRaphi\Behavior\NCRunSpeedDistMeanAllRec_Run' num2str(onlyRun) '.mat'];
    save(fullPath, 'speedOverDistMeanAllArr','speedOverDistStdAllArr',...
        'speedOverDistMeanAll','speedOverDistSEMAll','speedOverDistMeanAll1',...
        'speedOverDistSEMAll1');
    
    options.handle     = figure;
    set(options.handle,'OuterPosition',...
        [500 500 280 280])
    options.color_areaX = [27 117 187]./255;    % Blue theme
    options.color_lineX = [ 39 169 225]./255;
    options.color_areaY = [187 189 192]./255;    % Orange theme
    options.color_lineY = [167 169  171]./255;
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'sem';
    options.x_axisX = param.spaceSteps/10;
    options.x_axisY = param.spaceSteps/10;
    plot_areaerrorbarXY(speedOverDistAllArrAL/10, speedOverDistAllArr/10,...
        options);
    set(gca,'YLim',[0 60],'XLim',[0 180])
    xlabel('Dist (cm)')
    ylabel('Speed (cm/s)')
    legend('','Active','','Passive no cue')
    
    fileName1 = ['Z:\Yingxue\DataAnalysisRaphi\Behavior\' ...
        'SpeedVsDist-AL-PNoCue'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end
