function plotPhasePolarEg()
    
    onlyRun = 1;
    pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\Pyramidal\';
    
    load([pathAnal 'autoCorrPyrAllRec.mat']);
    mod.task = [autoCorrPyrNoCue.task autoCorrPyrAL.task autoCorrPyrPL.task];
    mod.indRec = [autoCorrPyrNoCue.indRec autoCorrPyrAL.indRec autoCorrPyrPL.indRec];
    mod.indNeu = [autoCorrPyrNoCue.indNeu autoCorrPyrAL.indNeu autoCorrPyrPL.indNeu];
    mod.isNeuWithField = [modPyrNoCue.isNeuWithField modPyrAL.isNeuWithField modPyrPL.isNeuWithField];
    mod.idxC = [autoCorrPyrAll.idxC1(autoCorrPyrAll.task == autoCorrPyrNoCue.task(1))' ...
        autoCorrPyrAll.idxC1(autoCorrPyrAll.task == autoCorrPyrAL.task(1))' ...
        autoCorrPyrAll.idxC1(autoCorrPyrAll.task == autoCorrPyrPL.task(1))'];
    mod.phaseMeanResultantLen = [modPyrNoCue.phaseMeanResultantLen modPyrAL.phaseMeanResultantLen modPyrPL.phaseMeanResultantLen];
    mod.phaseMeanDire = [modPyrNoCue.phaseMeanDire modPyrAL.phaseMeanDire modPyrPL.phaseMeanDire];
    mod.burstMeanDire = [modPyrNoCue.burstMeanDire modPyrAL.burstMeanDire modPyrPL.burstMeanDire];
    mod.fractBurst = [modPyrNoCue.fractBurst modPyrAL.fractBurst modPyrPL.fractBurst];
    
%     histPhaseMeanDire(mod.phaseMeanDire,'Mean theta phase (deg)','Num. neurons',...
%         'Mean theta phase - all cells',pathAnal,'Pyr_PhaseMeanDireAll');
%     
%     histPhaseMeanDire(mod.burstMeanDire(mod.fractBurst>0),'Mean burst theta phase (deg)','Num. neurons',...
%         'Mean burst theta phase - all cells',pathAnal,'Pyr_BurstMeanDireAll');
%     
%     histPhaseDireClusters(mod.phaseMeanDire(mod.idxC == 2),...
%         mod.phaseMeanDire(mod.idxC == 1),'Mean theta phase (deg)','Num. neurons',...
%         'Mean theta phase - clusters','Pyr_PhaseMeanDireC1C2Hist',pathAnal,modPyrStatsC.pWWPhaseMeanDireC);
%     
%     histPhaseDireClusters(mod.burstMeanDire(mod.idxC == 2 & mod.fractBurst>0),...
%         mod.burstMeanDire(mod.idxC == 1 & mod.fractBurst>0),'Mean burst theta phase (deg)','Num. neurons',...
%         'Mean burst theta phase - clusters','Pyr_BurstMeanDireC1C2Hist',pathAnal,modPyrStatsC.pWWburstMeanDireC);
    
%     clu1 = 2;
%     ind = find(mod.idxC == clu1 & mod.isNeuWithField == 1);
%     [phaseLenField,indP] = sort(mod.phaseMeanResultantLen(ind));
%     for i = 1:10
%         indPhaseLenFieldMax = ind(indP(end-i+1));
%         taskField = mod.task(indPhaseLenFieldMax);
%         indRecField = mod.indRec(indPhaseLenFieldMax);
%         indNeuField = mod.indNeu(indPhaseLenFieldMax);
% %         plotPolarPlotNeu(taskField,indRecField,indNeuField,onlyRun,'Phase histogram',...
% %             ['Pyr_EgC' num2str(clu1) 'FieldHighModPolar'],pathAnal,...
% %             mod.phaseMeanResultantLen(indPhaseLenFieldMax),...
% %             mod.phaseMeanDire(indPhaseLenFieldMax));
% %         plotPolarBurstNeu(taskField,indRecField,indNeuField,onlyRun,'Burst phase histogram',...
% %             ['Pyr_EgC' num2str(clu1) 'FieldHighBurstModPolar'],pathAnal);
%         plotPhaseHistNeu(taskField,indRecField,indNeuField,onlyRun,'Phase histogram',...
%             ['Pyr_EgC' num2str(clu1) 'FieldHighModHist'],pathAnal,...
%             mod.phaseMeanDire(indPhaseLenFieldMax))
%         plotBurstHistNeu(taskField,indRecField,indNeuField,onlyRun,'Burst phase histogram',...
%             ['Pyr_EgC' num2str(clu1) 'FieldHighBurstModHist'],pathAnal,...
%             mod.burstMeanDire(indPhaseLenFieldMax))
%     end
    
    clu2 = 1;
    ind = find(mod.idxC == clu2 & mod.isNeuWithField == 1);
    [phaseLenField,indP] = sort(mod.phaseMeanResultantLen(ind));
    for i = 1:10
        indPhaseLenFieldMax = ind(indP(end-i+1));
        taskField = mod.task(indPhaseLenFieldMax);
        indRecField = mod.indRec(indPhaseLenFieldMax);
        indNeuField = mod.indNeu(indPhaseLenFieldMax);
%         plotPolarPlotNeu(taskField,indRecField,indNeuField,onlyRun,'Phase histogram',...
%             ['Pyr_EgC' num2str(clu1) 'FieldHighModPolar'],pathAnal);
%         plotPolarBurstNeu(taskField,indRecField,indNeuField,onlyRun,'Burst phase histogram',...
%             ['Pyr_EgC' num2str(clu1) 'FieldHighBurstModPolar'],pathAnal);
        plotPhaseHistNeu(taskField,indRecField,indNeuField,onlyRun,'Phase histogram',...
            ['Pyr_EgC' num2str(clu2) 'FieldHighModHist'],pathAnal,...
            mod.phaseMeanDire(indPhaseLenFieldMax))
        plotBurstHistNeu(taskField,indRecField,indNeuField,onlyRun,'Burst phase histogram',...
            ['Pyr_EgC' num2str(clu2) 'FieldHighBurstModHist'],pathAnal,...
            mod.burstMeanDire(indPhaseLenFieldMax))
    end

end

function histPhaseDireClusters(x1,x2,xl,yl,ti,fn,pathAnal,p)
    [figNew,pos] = CreateFig();
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    colorArr = [163 207 98;...
                234 131 114]/255;
    
    x1 = [x1 x1+2*pi]/pi*180;
    x2 = [x2 x2+2*pi]/pi*180;
    steps = 0:10:720;
    hist(x1,steps);
    h = findobj(gca,'Type','patch');
    h.FaceColor = colorArr(2,:);
    h.EdgeColor = [0.5 0.5 0.5];
    h.FaceAlpha = 0.5;
    hold on;
    if(~isempty(x2))
        hist(x2,steps);
        h = findobj(gca,'Type','patch');
        h(1).FaceColor = colorArr(1,:);
        h(1).EdgeColor = [0.5 0.5 0.5];
        h(1).FaceAlpha = 0.5;
        title([ti ' p = ' num2str(p)])  
    else
        title(ti)
    end
    set(gca,'XLim',[0 720],'XTick',[0 180 360 540 720]);
    plot([360 360],[0 max(yticks)],':','Color',[0.4 0.4 0.4],'LineWidth',2);
    xlabel(xl);
    ylabel(yl);
    
    fileName1 = [pathAnal fn];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function histPhaseMeanDire(phaseMeanDire,xl,yl,ti,pathAnal,fn)
    x = [phaseMeanDire, 2*pi+phaseMeanDire];
    x = x/pi*180;
    steps = 0:10:720;
    [figNew,pos] = CreateFig();
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    hist(x,steps);
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0.7 0.7 0.7];
    h.EdgeColor = [0.5 0.5 0.5];
    hold on;
    plot([360 360],[0 max(yticks)],':','Color',[1 0 0],'LineWidth',2);
    set(gca,'XLim',[0 720],'XTick',[0 180 360 540 720]);
    xlabel(xl);
    ylabel(yl);
    title(ti);
    fileName1 = [pathAnal fn];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotPhaseHistNeu(task,indRec,indNeu,onlyRun,ti,fn,pathAnal,meanDire)
    RecordingList;
    if(task == 1)
        path = listRecordingsNoCuePath(indRec,:);
        fileName = listRecordingsNoCueFileName(indRec,:);
        mazeSess = mazeSessionNoCue(indRec);
    elseif(task == 2)
        path = listRecordingsActiveLickPath(indRec,:);
        fileName = listRecordingsActiveLickFileName(indRec,:);
        mazeSess = mazeSessionActiveLick(indRec);
    else
        path = listRecordingsPassiveLickPath(indRec,:);
        fileName = listRecordingsPassiveLickFileName(indRec,:);
        mazeSess = mazeSessionPassiveLick(indRec);
    end
    
    fileNameInfo = [fileName '_Info.mat'];
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        disp('_Info.mat file does not exist.');
        return;
    end
    load(fullPath,'beh'); 
   
	fileNameThetaPhase = [fileName '_ThetaPhaseL_Run' num2str(onlyRun) '.mat'];
    fullPath = [path fileNameThetaPhase];
    if(exist(fullPath) == 0)
        disp('_ThetaPhaseL file does not exist.');
        return;
    end
    load(fullPath,'spikeThetaPhaseStruct','spikeThetaPhaseStructSess');
    if(length(beh.mazeSessAll) > 1)
        spikeThetaPhase = spikeThetaPhaseStructSess{mazeSess};
    else
        spikeThetaPhase = spikeThetaPhaseStruct;
    end
    [figNew,pos] = CreateFig();
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    ind = find(spikeThetaPhase.posPhase >= 0 & spikeThetaPhase.posPhase <= 720);  
    hold on;
    plot([360 360],[0 0.04],':','Color',[1 0 0],'LineWidth',2);
    plot(spikeThetaPhase.posPhase(ind),...
        spikeThetaPhase.histPhasePerNeuron{indNeu}(ind)/sum(spikeThetaPhase.histPhasePerNeuron{indNeu}(ind)),...
        'Color',[0.6 0.6 0.6],'LineWidth',2);
    hold on;
    h = plot([meanDire/pi*180 meanDire/pi*180],[0 0.02],'k-');
    set(h,'LineWidth',3);
    h = plot([meanDire/pi*180+360 meanDire/pi*180+360],[0 0.02],'k-');
    set(h,'LineWidth',3);
    
    set(gca,'XLim',[0 720],'XTick',[0 180 360 540 720])
    xlabel('Theta phase (deg)')
    ylabel('Probability')
    title(ti);    
    
    fileName1 = [pathAnal fn '-task' num2str(task) '-rec' num2str(indRec) '-neu' num2str(indNeu)];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotBurstHistNeu(task,indRec,indNeu,onlyRun,ti,fn,pathAnal,meanDire)
    RecordingList;
    if(task == 1)
        path = listRecordingsNoCuePath(indRec,:);
        fileName = listRecordingsNoCueFileName(indRec,:);
        mazeSess = mazeSessionNoCue(indRec);
    elseif(task == 2)
        path = listRecordingsActiveLickPath(indRec,:);
        fileName = listRecordingsActiveLickFileName(indRec,:);
        mazeSess = mazeSessionActiveLick(indRec);
    else
        path = listRecordingsPassiveLickPath(indRec,:);
        fileName = listRecordingsPassiveLickFileName(indRec,:);
        mazeSess = mazeSessionPassiveLick(indRec);
    end
    
    fileNameInfo = [fileName '_Info.mat'];
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        disp('_Info.mat file does not exist.');
        return;
    end
    load(fullPath,'beh'); 
   
	fileNameBurst = [fileName '_burstAll_THL_Run' num2str(onlyRun) ...
                     '.mat'];
    fullPath = [path fileNameBurst];
    if(exist(fullPath) == 0)
        disp('_bustAll file does not exist.');
        return;
    end
    load(fullPath,'burstIsiPerNeuron','burstIsiPerNeuronSess');
    if(length(beh.mazeSessAll) > 1)
        burstIsi = burstIsiPerNeuronSess{mazeSess};
    else
        burstIsi = burstIsiPerNeuron;
    end
        
    [figNew,pos] = CreateFig();
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
%     colorArr = [163 207 98;...
%                 234 131 114]/255;
    
    stepPhase = 5;
    thetaSpikesMulti = getMultiCycles(burstIsi.phaseBurst{indNeu}', 2);
    [histPhasePerNeuron,posPhase] =...
            hist(thetaSpikesMulti,[stepPhase/2:stepPhase:720-stepPhase/2]); 
    hold on;
    plot([360 360],[0 0.04],':','Color',[1 0 0],'LineWidth',2);
    plot(posPhase,...
        histPhasePerNeuron/sum(histPhasePerNeuron),...
        'Color',[0.6 0.6 0.6],'LineWidth',2);
    hold on;
    h = plot([meanDire/pi*180 meanDire/pi*180],[0 0.02],'k-');
    set(h,'LineWidth',3);
    h = plot([meanDire/pi*180+360 meanDire/pi*180+360],[0 0.02],'k-');
    set(h,'LineWidth',3);
    
    set(gca,'XLim',[0 720],'XTick',[0 180 360 540 720])
    xlabel('Theta phase (deg)')
    ylabel('Probability')
    title(ti);    
    
    fileName1 = [pathAnal fn '-task' num2str(task) '-rec' num2str(indRec) '-neu' num2str(indNeu)];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotPolarPlotNeu(task,indRec,indNeu,onlyRun,ti,fn,pathAnal,phaseLen,phaseMean)
    
    RecordingList;
    if(task == 1)
        path = listRecordingsNoCuePath(indRec,:);
        fileName = listRecordingsNoCueFileName(indRec,:);
        mazeSess = mazeSessionNoCue(indRec);
    elseif(task == 2)
        path = listRecordingsActiveLickPath(indRec,:);
        fileName = listRecordingsActiveLickFileName(indRec,:);
        mazeSess = mazeSessionActiveLick(indRec);
    else
        path = listRecordingsPassiveLickPath(indRec,:);
        fileName = listRecordingsPassiveLickFileName(indRec,:);
        mazeSess = mazeSessionPassiveLick(indRec);
    end
    
    fileNameInfo = [fileName '_Info.mat'];
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        disp('_Info.mat file does not exist.');
        return;
    end
    load(fullPath,'beh'); 
   
	fileNameThetaPhase = [fileName '_ThetaPhaseL_Run' num2str(onlyRun) '.mat'];
    fullPath = [path fileNameThetaPhase];
    if(exist(fullPath) == 0)
        disp('_ThetaPhaseL file does not exist.');
        return;
    end
    load(fullPath,'spikeThetaPhaseStruct','spikeThetaPhaseStructSess');
    if(length(beh.mazeSessAll) > 1)
        spikeThetaPhase = spikeThetaPhaseStructSess{mazeSess};
    else
        spikeThetaPhase = spikeThetaPhaseStruct;
    end
        
    [figNew,pos] = CreateFig();
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
%     colorArr = [163 207 98;...
%                 234 131 114]/255;
    
    ind = find(spikeThetaPhase.posPhase/180*pi >= 2*pi & spikeThetaPhase.posPhase/180*pi <= 4*pi +10/180*pi);  
    polarplot(spikeThetaPhase.posPhase(ind)/180*pi-2*pi,...
        spikeThetaPhase.histPhasePerNeuron{indNeu}(ind)/max(spikeThetaPhase.histPhasePerNeuron{indNeu}(ind)),...
        'Color',[17 114 186]/255,'LineWidth',2);
    title(ti);
    hold on;
    %%%%arrow head %%%%
    arrowhead_length    = phaseLen/10; % arrow head length relative to resultant_length
    num_arrowlines = 100;
    arrowhead_angle = deg2rad(30); % degrees
    %%%%arrow tip coordinates %%%%
    t1 = repmat(phaseMean,1,num_arrowlines);
    r1 = repmat(phaseLen,1,num_arrowlines);
    %%%%arrow base coordinates %%%%
    b = arrowhead_length.*tan(linspace(0,arrowhead_angle,num_arrowlines/2));
    theta = atan(b./(phaseLen-arrowhead_length));
    pre_t2 = [theta, -theta];
    r2 = (phaseLen-arrowhead_length)./cos(pre_t2);
    t2 = t1(1)+pre_t2;
    %%%%plot %%%%
    polarplot([t1; t2],[r1; r2],'r')
    polarplot([phaseMean; phaseMean],[0; phaseLen],'Color',[1 0 0],'LineWidth',2)
    
    fileName1 = [pathAnal fn '-task' num2str(task) '-rec' num2str(indRec) '-neu' num2str(indNeu)];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotPolarBurstNeu(task,indRec,indNeu,onlyRun,ti,fn,pathAnal)
    
    RecordingList;
    if(task == 1)
        path = listRecordingsNoCuePath(indRec,:);
        fileName = listRecordingsNoCueFileName(indRec,:);
        mazeSess = mazeSessionNoCue(indRec);
    elseif(task == 2)
        path = listRecordingsActiveLickPath(indRec,:);
        fileName = listRecordingsActiveLickFileName(indRec,:);
        mazeSess = mazeSessionActiveLick(indRec);
    else
        path = listRecordingsPassiveLickPath(indRec,:);
        fileName = listRecordingsPassiveLickFileName(indRec,:);
        mazeSess = mazeSessionPassiveLick(indRec);
    end
    
    fileNameInfo = [fileName '_Info.mat'];
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        disp('_Info.mat file does not exist.');
        return;
    end
    load(fullPath,'beh'); 
   
	fileNameBurst = [fileName '_burstAll_THL_Run' num2str(onlyRun) ...
                     '.mat'];
    fullPath = [path fileNameBurst];
    if(exist(fullPath) == 0)
        disp('_bustAll file does not exist.');
        return;
    end
    load(fullPath,'burstIsiPerNeuron','burstIsiPerNeuronSess');
    if(length(beh.mazeSessAll) > 1)
        burstIsi = burstIsiPerNeuronSess{mazeSess};
    else
        burstIsi = burstIsiPerNeuron;
    end
        
    [figNew,pos] = CreateFig();
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
%     colorArr = [163 207 98;...
%                 234 131 114]/255;
    
    stepPhase = 5;
    thetaSpikesMulti = getMultiCycles(burstIsi.phaseBurst{indNeu}', 3);
    [histPhasePerNeuron,posPhase] =...
            hist(thetaSpikesMulti,[stepPhase/2:stepPhase:1080-stepPhase/2]);
    ind = find(posPhase/180*pi >= 2*pi & posPhase/180*pi <= 4*pi +10/180*pi);  
    polarplot(posPhase(ind)/180*pi-2*pi,...
        histPhasePerNeuron(ind)/max(histPhasePerNeuron(ind)),...
        'Color',[17 114 186]/255,'LineWidth',2);
    title(ti);
    
    fileName1 = [pathAnal fn '-task' num2str(task) '-rec' num2str(indRec) '-neu' num2str(indNeu)];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function thetaSpikesMulti = getMultiCycles(thetaPhase, numCycles)
% construct spike trains which repete the neuronal activity for
% multiple phase cycles
% thetaTime:        time of spikes
% thetaPhase:       phases of spikes

    indexThetaPos = find(thetaPhase >= 0);
    indexThetaNeg = find(thetaPhase < 0);
    thetaSpikesDeg = thetaPhase*360/(2*pi);
    thetaPhaseTmp = [];
    for i = 1:numCycles-1
        thetaPhaseTmp = [thetaPhaseTmp; thetaSpikesDeg+i*360];
    end
    thetaSpikesMulti = [thetaSpikesDeg(indexThetaPos); thetaPhaseTmp; thetaSpikesDeg(indexThetaNeg)+numCycles*360];
        
end