path = 'Z:\Raphael_tests\mice_expdata\ANM009\A009-20190128\A009-20190128-01\'; 
fileName = 'A009-20190128-01_DataStructure_mazeSection1_TrialType1.mat';
% load([path fileName],'trials');
% trialNo = 88:91;
trialNo = 50:52;
plotSpeedLickDistNoCue(path,fileName,trialNo)

function plotSpeedLickDistNoCue(path,fileName,trialNo)
    ind = strfind(fileName, '_');
    fileNameGen = [fileName(1:ind(1)) 'BehavElectrDataLFP.mat'];
    load([path fileNameGen],'Track','Laps');
    
    fileNameBeh = [fileName(1:ind(1)-1) 'BTDT.mat'];
    load([path fileNameBeh]);
    
    GlobalConst;

    speed = [];
    xMM = [];
    lickDist = [];
    pumpDist = [];
    totXMM = 0;
    xMMStart = [];
    xMMEnd = [];
    for i = trialNo  %Raphi A011-02190218-01
        indStart = Laps.startLfpInd(i);
        indEnd = Laps.startLfpInd(i+1)-1;
        speedAll = Track.speed_MMsecAll(indStart:indEnd);
        ind = speedAll < 0;
        speedAll(ind) = 0;
        speed = [speed speedAll'];
        
        indLick = behEventsTdt.lick(:,2) >= indStart & ...
            behEventsTdt.lick(:,2) <= indEnd;
        lickLFPInd = behEventsTdt.lick(indLick,2)-indStart+1;
        
        xMMTmp = xMMContinuous(Track.xMMAll(indStart:indEnd));
        xMM = [xMM xMMTmp'+totXMM];   
        xMMStart = [xMMStart totXMM];
        lickDist = [lickDist xMMTmp(lickLFPInd)'+totXMM];
        pumpDist = [pumpDist xMMTmp(Laps.pumpLfpInd{i}(1)-indStart)+totXMM];
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
    cd Z:\Yingxue\DataAnalysisRaphi
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