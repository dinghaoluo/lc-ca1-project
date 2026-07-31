path = 'Z:\Raphael_tests\mice_expdata\ANM009\A009-20190128\A009-20190128-01\'; 
fileName = 'A009-20190128-01_DataStructure_mazeSection1_TrialType1.mat';
% load([path fileName],'trials');
trialNo = 89:91;
% trialNo = 50:52;
plotSpeedLickDistNoCue(path,fileName,trialNo)

path = 'Z:\Raphael_tests\mice_expdata\ANM011\A011-20190218\A011-20190218-01\'; 
fileName = 'A011-20190218-01_DataStructure_mazeSection1_TrialType1.mat';
% load([path fileName],'trials');
trialNo = 137:139;
plotSpeedLickDist(path,fileName,trialNo)

path = 'Z:\Raphael_tests\mice_expdata\ANM013\A013-20190428\A013-20190428-01\'; 
fileName = 'A013-20190428-01_DataStructure_mazeSection1_TrialType1.mat';
% load([path fileName],'trials');
trialNo = 124:126;
plotSpeedLickDist(path,fileName,trialNo)

path = 'Z:\Raphael_tests\mice_expdata\ANM016\A016-20190603\A016-20190603-01\'; 
fileName = 'A016-20190603-01_DataStructure_mazeSection1_TrialType1.mat';
% load([path fileName],'trials');
trialNo = 85:87;
plotSpeedLickTime(path,fileName,trialNo)

function plotSpeedLickTime(path,fileName,trialNo)
    load([path fileName],'trials');
    
    GlobalConst;

    speed = [];
    time = [];
    licktime = [];
    pumpTime  = [];
    totT = 0;
    timeStart = [];
    timeEnd = [];
    for i = trialNo  %Raphi A011-02190218-01
        ind = trials{i}.speedAll >= 0;
        speedTmp = smooth(trials{i}.speedAll(ind),125);
        speed = [speed speedTmp'];        
        timeStart = [timeStart totT];
        timeEnd = [timeEnd totT+1/2];
        timeCur = (trials{i}.lfpIndStart:trials{i}.lfpIndEnd)/1250;
        time = [time timeCur(ind)];  
        licktime = [licktime (trials{i}.lickLfpInd'+trials{i}.lfpIndStart-1)/1250];
        pumpTime = [pumpTime (trials{i}.pumpLfpInd(1)+trials{i}.lfpIndStart-1)/1250];
        totT = totT + trials{i}.Nsamples/1250;
    end

    speed = smooth(speed/10,125);
    startTime = time(1);
    time = time - startTime + 1;
    licktime = licktime - startTime + 1;
    pumpTime = pumpTime - startTime + 1;
    figure
    h = plot(time,speed,'k');
    set(h,'LineWidth',1.2);
    hold on
    for j = 1:length(licktime)
        h = plot([licktime(j) licktime(j)],[60 64],'m');
        set(h,'LineWidth',0.1);
    end
    for j = 1:length(trialNo)
        h = plot([pumpTime(j) pumpTime(j)], [0 70]);
        set(h,'LineWidth',1,'Color',[0.6 0.9 0.6]);
        h = plot([timeStart(j) timeStart(j)], [0 70]);
        set(h,'LineWidth',1,'Color',[0.6 0.6 0.6]);
        h = plot([timeEnd(j) timeEnd(j)], [0 70]);
        set(h,'LineWidth',1,'Color',[0.6 0.6 0.6]);
    end
    set(gca,'yLim',[0 70])
    xlabel('Time(s)')
    ylabel('Speed (cm/s)')
    ind = findstr(fileName,'_');
    print ('-painters', '-dpdf', ['SpeedAndLickT_' fileName(1:ind(1)-1)], '-r600')
end

function plotSpeedLickDist(path,fileName,trialNo)
    load([path fileName],'trials');
    
    GlobalConst;

    speed = [];
    xMM = [];
    lickDist = [];
    pumpDist = [];
    totXMM = 0;
    xMMStart = [];
    xMMEnd = [];
    for i = trialNo  %Raphi A011-02190218-01
        ind = trials{i}.speedAll < 0;
        trials{i}.speedAll(ind) = 0;
        speed = [speed trials{i}.speedAll];
        
        xMMTmp = xMMContinuous(trials{i}.xMMAll);
        xMM = [xMM xMMTmp'+totXMM];   
        xMMStart = [xMMStart totXMM];
        xMMEnd = [xMMEnd xMMTmp(floor(sampleFq/2))+totXMM];
        lickDist = [lickDist xMMTmp(trials{i}.lickLfpInd)'+totXMM];
        pumpDist = [pumpDist xMMTmp(trials{i}.pumpLfpInd(1))+totXMM];
        totXMM = totXMM + xMMTmp(end);
    end

    speed = smooth(speed/10,37);
    xMM = xMM/10;
    xMMStart = xMMStart/10;
    xMMEnd = xMMEnd/10;
    lickDist = lickDist/10;
    pumpDist = pumpDist/10;
    figure
    h = plot(xMM,speed,'k');
    set(h,'LineWidth',1.2);
    hold on
    for j = 1:length(lickDist)
        h = plot([lickDist(j) lickDist(j)],[60 62],'m');
        set(h,'LineWidth',0.1);
    end
    for j = 1:length(trialNo)
        h = plot([pumpDist(j) pumpDist(j)], [0 70]);
        set(h,'LineWidth',1,'Color',[0.6 0.9 0.6]);
        h = plot([xMMStart(j) xMMStart(j)], [0 70]);
        set(h,'LineWidth',1,'Color',[0.6 0.6 0.6]);
        h = plot([xMMEnd(j) xMMEnd(j)], [0 70]);
        set(h,'LineWidth',1,'Color',[0.6 0.6 0.6]);
    end
    set(gca,'yLim',[0 70])
    xlabel('Dist (cm)')
    ylabel('Speed (cm/s)')
    ind = findstr(fileName,'_');
    print ('-painters', '-dpdf', ['SpeedAndLick_' fileName(1:ind(1)-1)], '-r600')
end

function plotSpeedLickDistNoCue(path,fileName,trialNo)
    load([path fileName],'trials');
    
    GlobalConst;

    speed = [];
    xMM = [];
    lickDist = [];
    pumpDist = [];
    totXMM = 0;
    xMMStart = [];
    xMMEnd = [];
    for i = trialNo  %Raphi A011-02190218-01
        ind = trials{i}.speedAll < 0;
        trials{i}.speedAll(ind) = 0;
        speed = [speed trials{i}.speedAll'];
        
        xMMTmp = xMMContinuous(trials{i}.xMMAll);
        xMM = [xMM xMMTmp'+totXMM];   
        xMMStart = [xMMStart totXMM];
        lickDist = [lickDist xMMTmp(trials{i}.lickLfpInd)'+totXMM];
        pumpDist = [pumpDist xMMTmp(trials{i}.pumpLfpInd(1))+totXMM];
        totXMM = totXMM + xMMTmp(end);
    end

    speed = smooth(speed/10,37);
    xMM = xMM/10;
    xMMStart = xMMStart/10;
    xMMEnd = xMMEnd/10;
    lickDist = lickDist/10;
    pumpDist = pumpDist/10;
    figure
    h = plot(xMM,speed,'k');
    set(h,'LineWidth',1.2);
    hold on
    for j = 1:length(lickDist)
        h = plot([lickDist(j) lickDist(j)],[60 62],'m');
        set(h,'LineWidth',0.1);
    end
    for j = 1:length(trialNo)
        h = plot([pumpDist(j) pumpDist(j)], [0 70]);
        set(h,'LineWidth',1,'Color',[0.6 0.9 0.6]);
        h = plot([xMMStart(j) xMMStart(j)], [0 70]);
        set(h,'LineWidth',1,'Color',[0.6 0.6 0.6]);
    end
    set(gca,'yLim',[0 70])
    xlabel('Dist (cm)')
    ylabel('Speed (cm/s)')
    ind = findstr(fileName,'_');
    print ('-painters', '-dpdf', ['SpeedAndLick_' fileName(1:ind(1)-1)], '-r600')
end

function xMM1 = xMMContinuous(xMM)
    indResetXMM = find(diff(xMM) < 0);
    indResetXMM = [indResetXMM+1; length(xMM)+1];
    for n = 1:length(indResetXMM)-1
        xMM(indResetXMM(n):indResetXMM(n+1)-1) = ...
            xMM(indResetXMM(n):indResetXMM(n+1)-1) + xMM(indResetXMM(n)-1);
    end
    xMM1 = xMM-xMM(1);  
end