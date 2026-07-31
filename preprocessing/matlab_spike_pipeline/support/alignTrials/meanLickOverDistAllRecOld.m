function meanLickOverDistAllRec
% calculate mean lick over all the recordings

    RecordingList;
    
    indTr = 5:55;
    nTrials = length(indTr);
    
    %% active licking
    numRec = size(listRecordingsActiveLickPath,1);
    totNTr = 0;
    lickOverDistMeanAllArr = [];
    lickOverDistStdAllArr = [];
    lickOverDistAllArr = [];
    for i = 1:numRec
        path = listRecordingsActiveLickPath(i,:);
        fileName = listRecordingsActiveLickFileName(i,:);
        
        fullPath = [path,fileName, '_lickDist.mat'];
        load(fullPath,'lickOverDistMean','lickOverDistStd','lickOverDist','param');
        
        lickOverDistMeanAllArr = [lickOverDistMeanAllArr; lickOverDistMean];
        lickOverDistStdAllArr = [lickOverDistStdAllArr; lickOverDistStd];
        if(size(lickOverDist,1) < max(indTr))
            lickOverDistAllArr = [lickOverDistAllArr; lickOverDist];
            totNTr = totNTr+size(lickOverDist,1);
        else
            lickOverDistAllArr = [lickOverDistAllArr; lickOverDist(indTr,:)]; 
            totNTr = totNTr+nTrials;
        end
    end
    
    lickOverDistMeanAll = mean(lickOverDistMeanAllArr);
    lickOverDistSEMAll = sqrt(sum(lickOverDistStdAllArr.^2)/numRec);
    lickOverDistMeanAll1 = mean(lickOverDistAllArr);
    lickOverDistSEMAll1 = std(lickOverDistAllArr)/sqrt(nTrials*numRec);
    
    fullPath = ['Z:\Yingxue\DataAnalysisRaphi\ALLickDistMeanAllRec.mat'];
    save(fullPath, 'lickOverDistMeanAllArr','lickOverDistStdAllArr',...
        'lickOverDistMeanAll','lickOverDistSEMAll','lickOverDistMeanAll1',...
        'lickOverDistSEMAll1','param');
    
    lickOverDistAllArrAL = lickOverDistAllArr;
    
    %% passive licking
    numRec = size(listRecordingsPassiveLickPath,1);
    totNTr = 0;
    lickOverDistMeanAllArr = [];
    lickOverDistStdAllArr = [];
    lickOverDistAllArr = [];
    for i = 1:numRec
        path = listRecordingsPassiveLickPath(i,:);
        fileName = listRecordingsPassiveLickFileName(i,:);
        
        fullPath = [path,fileName, '_lickDist.mat'];
        load(fullPath,'lickOverDistMean','lickOverDistStd','lickOverDist');
        
        lickOverDistMeanAllArr = [lickOverDistMeanAllArr; lickOverDistMean];
        lickOverDistStdAllArr = [lickOverDistStdAllArr; lickOverDistStd];
        if(size(lickOverDist,1) < max(indTr))
            lickOverDistAllArr = [lickOverDistAllArr; lickOverDist];
            totNTr = totNTr+size(lickOverDist,1);
        else
            lickOverDistAllArr = [lickOverDistAllArr; lickOverDist(indTr,:)]; 
            totNTr = totNTr+nTrials;
        end
    end
    
    lickOverDistMeanAll = mean(lickOverDistMeanAllArr);
    lickOverDistSEMAll = sqrt(sum(lickOverDistStdAllArr.^2)/numRec);
    lickOverDistMeanAll1 = mean(lickOverDistAllArr);
    lickOverDistSEMAll1 = std(lickOverDistAllArr)/sqrt(nTrials*numRec);
    
    fullPath = ['Z:\Yingxue\DataAnalysisRaphi\PLLickDistMeanAllRec.mat'];
    save(fullPath, 'lickOverDistMeanAllArr','lickOverDistStdAllArr',...
        'lickOverDistMeanAll','lickOverDistSEMAll','lickOverDistMeanAll1',...
        'lickOverDistSEMAll1');
    
    options.handle     = figure;
    options.color_areaX = [27 117 187]./255;    % Blue theme
    options.color_lineX = [ 39 169 225]./255;
    options.color_areaY = [187 189 192]./255;    % Orange theme
    options.color_lineY = [167 169  171]./255;
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'sem';
    options.x_axisX = param.spaceSteps/10;
    options.x_axisY = param.spaceSteps/10;
    plot_areaerrorbarXY(lickOverDistAllArrAL, lickOverDistAllArr,...
        options);
    xlabel('Dist (cm)')
    ylabel('No. licks ')
    legend('','Active','','Passive')
    
    %% no cue
    %% passive licking
    numRec = size(listRecordingsNoCuePath,1);
    totNTr = 0;
    lickOverDistMeanAllArr = [];
    lickOverDistStdAllArr = [];
    lickOverDistAllArr = [];
    for i = 1:numRec
        path = listRecordingsNoCuePath(i,:);
        fileName = listRecordingsNoCueFileName(i,:);
        
        fullPath = [path,fileName, '_lickDist.mat'];
        load(fullPath,'lickOverDistMean','lickOverDistStd','lickOverDist');
        
        lickOverDistMeanAllArr = [lickOverDistMeanAllArr; lickOverDistMean];
        lickOverDistStdAllArr = [lickOverDistStdAllArr; lickOverDistStd];
        if(size(lickOverDist,1) < max(indTr))
            lickOverDistAllArr = [lickOverDistAllArr; lickOverDist];
            totNTr = totNTr+size(lickOverDist,1);
        else
            lickOverDistAllArr = [lickOverDistAllArr; lickOverDist(indTr,:)]; 
            totNTr = totNTr+nTrials;
        end
    end
    
    lickOverDistMeanAll = mean(lickOverDistMeanAllArr);
    lickOverDistSEMAll = sqrt(sum(lickOverDistStdAllArr.^2)/numRec);
    lickOverDistMeanAll1 = mean(lickOverDistAllArr);
    lickOverDistSEMAll1 = std(lickOverDistAllArr)/sqrt(nTrials*numRec);
    
    fullPath = ['Z:\Yingxue\DataAnalysisRaphi\NCLickDistMeanAllRec.mat'];
    save(fullPath, 'lickOverDistMeanAllArr','lickOverDistStdAllArr',...
        'lickOverDistMeanAll','lickOverDistSEMAll','lickOverDistMeanAll1',...
        'lickOverDistSEMAll1');
    
    options.handle     = figure;
    options.color_areaX = [27 117 187]./255;    % Blue theme
    options.color_lineX = [ 39 169 225]./255;
    options.color_areaY = [187 189 192]./255;    % Orange theme
    options.color_lineY = [167 169  171]./255;
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'sem';
    options.x_axisX = param.spaceSteps/10;
    options.x_axisY = param.spaceSteps/10;
    plot_areaerrorbarXY(lickOverDistAllArrAL, lickOverDistAllArr,...
        options);
    xlabel('Dist (cm)')
    ylabel('No. licks ')
    legend('','Active','','Passive no cue')
end
