function PyrIntInitPeakSegment(taskSel)
% compare Pyramidal neurons and PV interneurons on their initial peak
    
    pathAnal0 = 'Z:\Yingxue\Draft\SST\Pyramidal\';
    pathAnalInt0 = 'Z:\Yingxue\Draft\SST\Interneuron\';
    
    if(exist([pathAnal0 'initPeakPyrAllRec.mat']))
        load([pathAnal0 'initPeakPyrAllRec.mat']);
    end
    if(exist([pathAnalInt0 'initPeakIntAllRec.mat']))
        load([pathAnalInt0 'initPeakIntAllRec.mat']);
    end
    pathAnalPeak0 = 'Z:\Yingxue\Draft\SST\PyramidalInitPeak\2\';
    if(exist([pathAnalPeak0 'initPeakPyrAllRec_km2.mat']))
        load([pathAnalPeak0 'initPeakPyrAllRec_km2.mat'],'FRProfileMean');
    end
    
    if(taskSel == 1) % including all the neurons
        pathAnal = 'Z:\Yingxue\Draft\SST\PyramidalIntInitPeak\';
    elseif(taskSel == 2) % including AL and PL neurons
        pathAnal = 'Z:\Yingxue\Draft\SST\PyramidalIntInitPeakALPL\';
    else
        pathAnal = 'Z:\Yingxue\Draft\SST\PyramidalIntInitPeakAL\';
    end
    if(exist([pathAnal 'initPeakPyrIntAllRec.mat']))
        load([pathAnal 'initPeakPyrIntAllRec.mat']);
    end
    if(exist([pathAnal 'initPeakPyrIntAllRecSeg.mat']))
        load([pathAnal 'initPeakPyrIntAllRecSeg.mat']);
    end
    
    if(exist(pathAnal) == 0)
        mkdir(pathAnal);
    end
                    
    %% pyramidal rise and down
    indTmp = modPyr1NoCue.timeStepRun > 0;
    [~,indMax] = max(InitAll.avgFRProfile(:,indTmp)');
    [indMaxOrd,indOrd] = sort(indMax);
    avgFRProfileOrdered = InitAll.avgFRProfile(indOrd,:);
    InitAll1.indOrdered = indOrd;
    InitAll1.indMaxOrdered = indMaxOrd;
    
    mean0to1 = mean(avgFRProfileOrdered(:,FRProfileMean.indFR0to1),2);
    meanBefRun = mean(avgFRProfileOrdered(:,FRProfileMean.indFRBefRun),2);
    InitAll1.ratio0to1BefRun = mean0to1./meanBefRun;
    
    InitAll = PyrIntInitPeakByType(taskSel,modPyr1NoCue,modPyr1AL,modPyr1PL,...
                                    modInt1NoCue,modInt1AL,modInt1PL);

    InitField.avgFRProfile = InitAll.avgFRProfile(InitAll.isNeuWithFieldAligned==1,:);
    [~,indMax] = max(InitField.avgFRProfile(:,indTmp)');
    [indMaxOrd,indOrd] = sort(indMax);
    avgFRProfileOrdered = InitField.avgFRProfile(indOrd,:);
    InitField.indOrdered = indOrd;
    InitField.indMaxOrdered = indMaxOrd;

    mean0to1 = mean(avgFRProfileOrdered(:,FRProfileMean.indFR0to1),2);
    meanBefRun = mean(avgFRProfileOrdered(:,FRProfileMean.indFRBefRun),2);
    InitField.ratio0to1BefRun = mean0to1./meanBefRun;

    idx = find(InitField.ratio0to1BefRun >= 1.5);
    idx1 = find(InitField.ratio0to1BefRun <= 2/3);
        
    save([pathAnal 'initPeakPyrIntAllRecSeg.mat'],'InitAll1','InitField'); 
    
    %% plot depth of rise and down pyramidal neurons
    % plot change around the run onset
    timeBin = modInt1NoCue.timeStepRun(2) - modInt1NoCue.timeStepRun(1);
    plotChangeRunOnset(pathAnal,InitAll1.indMaxOrdered*timeBin,log(InitAll1.ratio0to1BefRun),...
            'Time (s)','log (FR 0to1s / BefRun)','FRRatio0to1vsBefRunPyrAllOrdered');
    
    plotChangeRunOnset(pathAnal,InitField.indMaxOrdered*timeBin,log(InitField.ratio0to1BefRun),...
            'Time (s)','log (FR 0to1s / BefRun) Field','FRRatio0to1vsBefRunPyrFieldOrdered');
        
    plotChangeRunOnsetRD(pathAnal,InitField.indMaxOrdered*timeBin,log(InitField.ratio0to1BefRun),...
            'Time (s)','log (FR 0to1s / BefRun) Field','FRRatio0to1vsBefRunPyrFieldOrderedRD',idx,idx1);
end

function plotChangeRunOnset(pathAnal,x,y,xl,yl,fileName)
    handle1 = figure;
    set(handle1,'OuterPosition',[500 500 280 280]);
    h = plot(x,y,'.');
    xlabel(xl);
    ylabel(yl);
    set(h,'MarkerSize',7,'MarkerFaceColor',[0.4 0.4 0.4],...
        'MarkerEdgeColor',[0.4 0.4 0.4]);
    set(gca,'FontSize',10,'XLim',[0 4],'YLim',[-5 5]);
    ind = y > 5 | y < -5;
    ind = sum(ind)/length(y);
    title(['Percent not displayed ' num2str(ind*100) '%'])
    fileName1 = [pathAnal fileName];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotChangeRunOnsetRD(pathAnal,x,y,xl,yl,fileName,idx,idx1)
    handle1 = figure;
    set(handle1,'OuterPosition',[500 500 280 280]);
    h = plot(x,y,'.');
    xlabel(xl);
    ylabel(yl);
    set(h,'MarkerSize',7,'MarkerFaceColor',[188 212 238]/255,...
        'MarkerEdgeColor',[188 212 238]/255);
    set(gca,'FontSize',10,'XLim',[0 4],'YLim',[-5 5]);
    
    hold on;
    h = plot(x(idx),y(idx),'r.');
    set(h,'MarkerSize',7,'MarkerFaceColor',[251 204 171]/255,...
        'MarkerEdgeColor',[251 204 171]/255);
    h = plot(x(idx1),y(idx1),'g.');
    set(h,'MarkerSize',7,'MarkerFaceColor',[214 170 207]/255,...
        'MarkerEdgeColor',[214 170 207]/255);
    
    ind = y > 5 | y < -5;
    ind = sum(ind)/length(y);
    title(['Percent not displayed ' num2str(ind*100) '%'])
    fileName1 = [pathAnal fileName];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end
