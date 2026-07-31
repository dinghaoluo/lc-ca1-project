function plotLickSpeedProfileSingleRec(path,fileNameBase,sessions)
% plot the licking and speed profile for a single recording
% e.g.: plotLickSpeedProfileSingleRec('Z:\Yingxue\DataAnalysisXiaoliang\ANMZ529\A529-20200710\','A529-20200710',[1 2])
    fullPath = [path '\' fileNameBase '_compSess.mat']; 
    load(fullPath, 'sessDataRun','sessDataLick','sessDataSpeed');
    folderPath = 'Z:\Yingxue\DataAnalysisXiaoliang\';
    plotTrace(folderPath,sessDataLick.spaceSteps/10,...
            sessDataLick.Run{sessions(2)},...
            sessDataLick.Run{sessions(1)},...
            [30,210],[0 2],'No. lick',['lickProfile' fileNameBase])
        
    plotTrace(folderPath,sessDataSpeed.spaceSteps/10,...
            sessDataSpeed.Run{sessions(2)},...
            sessDataSpeed.Run{sessions(1)},...
            [0,180],[0 80],'Speed (cm/sec)',['speedProfile' fileNameBase])
end

function plotTrace(pathAnal,timeSteps,mani,pre,xli,yli,yl,filename)
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
    options.x_axisX = timeSteps;
    options.x_axisY = timeSteps;
    plot_areaerrorbarXY(mani, pre,...
        options);
    set(gca,'XLim',xli);
    set(gca,'YLim',yli);
    xlabel('Time (s)')
    ylabel(yl)
    legend('','M','','Pre')
    
    fileName1 = [pathAnal filename...
        '_Mani_Pre'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end
