function meanLickAlignedOverDistAllRec
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
    totNTr1 = 0;
    lickAlignedOverDistMeanAllArr = [];
    lickAlignedOverDistStdAllArr = [];
    lickAlignedOverDistAllArr = [];
    for i = 1:numRec  %% only taking the first 10 for now
        path = listRecordingsActiveLickPath(i,:);
        fileName = listRecordingsActiveLickFileName(i,:);
        
        fullPath = [path,fileName, '_lickDist.mat'];
        load(fullPath,'lickAlignedOverDistMean','lickAlignedOverDistStd','lickAlignedOverDist',...
            'lickOverDistMean','lickOverDistStd','lickOverDist','param');
        
        lickOverDistMeanAllArr = [lickOverDistMeanAllArr;lickOverDistMean];
        lickOverDistStdAllArr = [lickOverDistStdAllArr;lickOverDistStd];
        if(size(lickOverDist,1) < max(indTr))
            lickOverDistAllArr = [lickOverDistAllArr; lickOverDist];
            totNTr = totNTr+size(lickOverDist,1);
        else
            lickOverDistAllArr = [lickOverDistAllArr; lickOverDist(indTr,:)]; 
            totNTr = totNTr+nTrials;
        end
        
        lickAlignedOverDistMeanAllArr = [lickAlignedOverDistMeanAllArr; lickAlignedOverDistMean];
        lickAlignedOverDistStdAllArr = [lickAlignedOverDistStdAllArr; lickAlignedOverDistStd];
        if(size(lickAlignedOverDist,1) < max(indTr))
            lickAlignedOverDistAllArr = [lickAlignedOverDistAllArr; lickAlignedOverDist];
            totNTr1 = totNTr1+size(lickAlignedOverDist,1);
        else
            lickAlignedOverDistAllArr = [lickAlignedOverDistAllArr; lickAlignedOverDist(indTr,:)]; 
            totNTr1 = totNTr1+nTrials;
        end
    end
    
    lickOverDistMeanAll = mean(lickOverDistMeanAllArr);
    lickOverDistSEMAll = sqrt(sum(lickOverDistStdAllArr.^2)/numRec);
    lickOverDistMeanAll1 = mean(lickOverDistAllArr);
    lickOverDistSEMAll1 = std(lickOverDistAllArr)/sqrt(totNTr);
    
    lickAlignedOverDistMeanAll = mean(lickAlignedOverDistMeanAllArr);
    lickAlignedOverDistSEMAll = sqrt(sum(lickAlignedOverDistStdAllArr.^2)/numRec);
    lickAlignedOverDistMeanAll1 = mean(lickAlignedOverDistAllArr);
    lickAlignedOverDistSEMAll1 = std(lickAlignedOverDistAllArr)/sqrt(totNTr1);
    
    fullPath = ['Z:\Yingxue\DataAnalysisRaphi\Behavior\ALLickAlignedDistMeanAllRec.mat'];
    save(fullPath, 'lickOverDistAllArr','lickOverDistMeanAllArr','lickOverDistStdAllArr',...
        'lickOverDistMeanAll','lickOverDistSEMAll','lickOverDistMeanAll1',...
        'lickOverDistSEMAll1','lickAlignedOverDistAllArr',...
        'lickAlignedOverDistMeanAllArr','lickAlignedOverDistStdAllArr',...
        'lickAlignedOverDistMeanAll','lickAlignedOverDistSEMAll','lickAlignedOverDistMeanAll1',...
        'lickAlignedOverDistSEMAll1','param');
    
    lickOverDistAllArrAL = lickOverDistAllArr;
    lickAlignedOverDistAllArrAL = lickAlignedOverDistAllArr;
    
    %% passive licking
    numRec = size(listRecordingsPassiveLickPath,1);
    totNTr = 0;
    lickOverDistMeanAllArr = [];
    lickOverDistStdAllArr = [];
    lickOverDistAllArr = [];
    totNTr1 = 0;
    lickAlignedOverDistMeanAllArr = [];
    lickAlignedOverDistStdAllArr = [];
    lickAlignedOverDistAllArr = [];
    for i = 1:numRec
        path = listRecordingsPassiveLickPath(i,:);
        fileName = listRecordingsPassiveLickFileName(i,:);
        
        fullPath = [path,fileName, '_lickDist.mat'];
        load(fullPath,'lickAlignedOverDistMean','lickAlignedOverDistStd','lickAlignedOverDist',...
            'lickOverDistMean','lickOverDistStd','lickOverDist');
        
        lickOverDistMeanAllArr = [lickOverDistMeanAllArr;lickOverDistMean];
        lickOverDistStdAllArr = [lickOverDistStdAllArr;lickOverDistStd];
        if(size(lickOverDist,1) < max(indTr))
            lickOverDistAllArr = [lickOverDistAllArr; lickOverDist];
            totNTr = totNTr+size(lickOverDist,1);
        else
            lickOverDistAllArr = [lickOverDistAllArr; lickOverDist(indTr,:)]; 
            totNTr = totNTr+nTrials;
        end
        
        lickAlignedOverDistMeanAllArr = [lickAlignedOverDistMeanAllArr; lickAlignedOverDistMean];
        lickAlignedOverDistStdAllArr = [lickAlignedOverDistStdAllArr; lickAlignedOverDistStd];
        if(size(lickAlignedOverDist,1) < max(indTr))
            lickAlignedOverDistAllArr = [lickAlignedOverDistAllArr; lickAlignedOverDist];
            totNTr1 = totNTr1+size(lickAlignedOverDist,1);
        else
            lickAlignedOverDistAllArr = [lickAlignedOverDistAllArr; lickAlignedOverDist(indTr,:)]; 
            totNTr1 = totNTr1+nTrials;
        end
    end
    
    lickOverDistMeanAll = mean(lickOverDistMeanAllArr);
    lickOverDistSEMAll = sqrt(sum(lickOverDistStdAllArr.^2)/numRec);
    lickOverDistMeanAll1 = mean(lickOverDistAllArr);
    lickOverDistSEMAll1 = std(lickOverDistAllArr)/sqrt(totNTr);
    
    lickAlignedOverDistMeanAll = mean(lickAlignedOverDistMeanAllArr);
    lickAlignedOverDistSEMAll = sqrt(sum(lickAlignedOverDistStdAllArr.^2)/numRec);
    lickAlignedOverDistMeanAll1 = mean(lickAlignedOverDistAllArr);
    lickAlignedOverDistSEMAll1 = std(lickAlignedOverDistAllArr)/sqrt(totNTr1);
    
    fullPath = ['Z:\Yingxue\DataAnalysisRaphi\Behavior\PLLickAlignedDistMeanAllRec.mat'];
    save(fullPath, 'lickOverDistAllArr','lickOverDistMeanAllArr','lickOverDistStdAllArr',...
        'lickOverDistMeanAll','lickOverDistSEMAll','lickOverDistMeanAll1',...
        'lickOverDistSEMAll1','lickAlignedOverDistAllArr',...
        'lickAlignedOverDistMeanAllArr','lickAlignedOverDistStdAllArr',...
        'lickAlignedOverDistMeanAll','lickAlignedOverDistSEMAll','lickAlignedOverDistMeanAll1',...
        'lickAlignedOverDistSEMAll1');
    
    options.handle     = figure;
    set(options.handle,'OuterPosition',...
        [500 500 280 280])
    options.color_areaX = [27 117 187]./255;    % Blue theme
    options.color_lineX = [ 39 169 225]./255;
    options.color_areaY = [187 189 192]./255;    % Orange theme
    options.color_lineY = [167 169  171]./255;
    options.alpha      = 0.5;
    options.line_width = 0.5;
    options.error      = 'sem';
    options.x_axisX = param.spaceStepsAligned/10;
    options.x_axisY = param.spaceStepsAligned/10;
    plot_areaerrorbarXY(lickAlignedOverDistAllArrAL, lickAlignedOverDistAllArr,...
        options);
    xlabel('Dist (cm)')
    ylabel('No. licks ')
    legend('','Active','','Passive')
    
    fileName1 = ['Z:\Yingxue\DataAnalysisRaphi\Behavior\' ...
        'LickAlignedToPumpVsDist-AL-PL'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    options.handle     = figure;
    set(options.handle,'OuterPosition',...
        [500 500 280 280])
    options.color_areaX = [27 117 187]./255;    % Blue theme
    options.color_lineX = [ 39 169 225]./255;
    options.color_areaY = [187 189 192]./255;    % Orange theme
    options.color_lineY = [167 169  171]./255;
    options.alpha      = 0.5;
    options.line_width = 0.5;
    options.error      = 'sem';
    options.x_axisX = param.spaceSteps/10;
    options.x_axisY = param.spaceSteps/10;
    plot_areaerrorbarXY(lickOverDistAllArrAL, lickOverDistAllArr,...
        options);
    set(gca,'YLim',[0 0.5],'XLim',[30 180])
    xlabel('Dist (cm)')
    ylabel('No. licks ')
    legend('','Active','','Passive')
    
    fileName1 = ['Z:\Yingxue\DataAnalysisRaphi\Behavior\' ...
        'LickVsDist-AL-PL'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    %% no cue
    %% passive licking
    numRec = size(listRecordingsNoCuePath,1);
    totNTr = 0;    
    lickOverDistMeanAllArr = [];
    lickOverDistStdAllArr = [];
    lickOverDistAllArr = [];
    totNTr1 = 0;
    lickAlignedOverDistMeanAllArr = [];
    lickAlignedOverDistStdAllArr = [];
    lickAlignedOverDistAllArr = [];
    for i = 1:numRec
        path = listRecordingsNoCuePath(i,:);
        fileName = listRecordingsNoCueFileName(i,:);
        
        fullPath = [path,fileName, '_lickDist.mat'];
        load(fullPath,'lickAlignedOverDistMean','lickAlignedOverDistStd','lickAlignedOverDist',...
            'lickOverDistMean','lickOverDistStd','lickOverDist');
        
        lickOverDistMeanAllArr = [lickOverDistMeanAllArr;lickOverDistMean];
        lickOverDistStdAllArr = [lickOverDistStdAllArr;lickOverDistStd];
        if(size(lickOverDist,1) < max(indTr))
            lickOverDistAllArr = [lickOverDistAllArr; lickOverDist];
            totNTr = totNTr+size(lickOverDist,1);
        else
            lickOverDistAllArr = [lickOverDistAllArr; lickOverDist(indTr,:)]; 
            totNTr = totNTr+nTrials;
        end
        
        lickAlignedOverDistMeanAllArr = [lickAlignedOverDistMeanAllArr; lickAlignedOverDistMean];
        lickAlignedOverDistStdAllArr = [lickAlignedOverDistStdAllArr; lickAlignedOverDistStd];
        if(size(lickAlignedOverDist,1) < max(indTr))
            lickAlignedOverDistAllArr = [lickAlignedOverDistAllArr; lickAlignedOverDist];
            totNTr1 = totNTr1+size(lickAlignedOverDist,1);
        else
            lickAlignedOverDistAllArr = [lickAlignedOverDistAllArr; lickAlignedOverDist(indTr,:)]; 
            totNTr1 = totNTr1+nTrials;
        end
    end
    
    lickOverDistMeanAll = mean(lickOverDistMeanAllArr);
    lickOverDistSEMAll = sqrt(sum(lickOverDistStdAllArr.^2)/numRec);
    lickOverDistMeanAll1 = mean(lickOverDistAllArr);
    lickOverDistSEMAll1 = std(lickOverDistAllArr)/sqrt(totNTr);
    
    lickAlignedOverDistMeanAll = mean(lickAlignedOverDistMeanAllArr);
    lickAlignedOverDistSEMAll = sqrt(sum(lickAlignedOverDistStdAllArr.^2)/numRec);
    lickAlignedOverDistMeanAll1 = mean(lickAlignedOverDistAllArr);
    lickAlignedOverDistSEMAll1 = std(lickAlignedOverDistAllArr)/sqrt(totNTr1);
    
    fullPath = ['Z:\Yingxue\DataAnalysisRaphi\Behavior\NCLickAlignedDistMeanAllRec.mat'];
    save(fullPath, 'lickOverDistAllArr','lickOverDistMeanAllArr','lickOverDistStdAllArr',...
        'lickOverDistMeanAll','lickOverDistSEMAll','lickOverDistMeanAll1',...
        'lickOverDistSEMAll1','lickAlignedOverDistAllArr',...
        'lickAlignedOverDistMeanAllArr','lickAlignedOverDistStdAllArr',...
        'lickAlignedOverDistMeanAll','lickAlignedOverDistSEMAll','lickAlignedOverDistMeanAll1',...
        'lickAlignedOverDistSEMAll1');
    
    options.handle     = figure;
    set(options.handle,'OuterPosition',...
        [500 500 280 280])
    options.color_areaX = [27 117 187]./255;    % Blue theme
    options.color_lineX = [ 39 169 225]./255;
    options.color_areaY = [187 189 192]./255;    % Orange theme
    options.color_lineY = [167 169  171]./255;
    options.alpha      = 0.5;
    options.line_width = 0.5;
    options.error      = 'sem';
    options.x_axisX = param.spaceStepsAligned/10;
    options.x_axisY = param.spaceStepsAligned/10;
    plot_areaerrorbarXY(lickAlignedOverDistAllArrAL, lickAlignedOverDistAllArr,...
        options);
    xlabel('Dist (cm)')
    ylabel('No. licks ')
    legend('','Active','','Passive no cue')
    
    fileName1 = ['Z:\Yingxue\DataAnalysisRaphi\Behavior\' ...
        'LickAlignedToPumpVsDist-AL-PNoCue'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    options.handle     = figure;
    set(options.handle,'OuterPosition',...
        [500 500 280 280])
    options.color_areaX = [27 117 187]./255;    % Blue theme
    options.color_lineX = [ 39 169 225]./255;
    options.color_areaY = [187 189 192]./255;    % Orange theme
    options.color_lineY = [167 169  171]./255;
    options.alpha      = 0.5;
    options.line_width = 0.5;
    options.error      = 'sem';
    options.x_axisX = param.spaceSteps/10;
    options.x_axisY = param.spaceSteps/10;
    plot_areaerrorbarXY(lickOverDistAllArrAL, lickOverDistAllArr,...
        options);
    set(gca,'YLim',[0 0.5],'XLim',[30 180])
    xlabel('Dist (cm)')
    ylabel('No. licks ')
    legend('','Active','','Passive no cue')
    
    fileName1 = ['Z:\Yingxue\DataAnalysisRaphi\Behavior\' ...
        'LickVsDist-AL-PNoCue'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end
