function plotGaChSpeedEg(path, filename, trialNo)
%% plot GaCh and corresponding speed signal
%% E.G.: plotGaChSpeedEg('Z:\Xiaoliang\mice-expdata\ANMZ583\A583-20211115\A583-20211115-04\','A583-20211115-04',54:57)

    GlobalConst2P;
    
    fullName = [path filename '_Behav2PDataLFP.mat'];
    load(fullName,'Track','Laps');

    startInd = Laps.startLfpInd(trialNo(1));
    endInd = Laps.endLfpInd(trialNo(end));
    
    figure
    speed = Track.speed_MMsecAll;
    idx = speed >= 0;
    speed(~idx) = 0;
    plot((1:endInd-startInd+1)/sampleFq,speed(startInd:endInd)/10,'k-');
    hold on;
    for i = 1:length(trialNo)
        for j = 1:length(Laps.lickLfpInd{trialNo(i)})
            t = (Laps.lickLfpInd{trialNo(i)}(j)-startInd+1)/sampleFq;
            plot([t t],[90 95],'m-');
        end
    end
    set(gca,'XLim',[0 (endInd-startInd+1)/sampleFq],'YLim',[0 100],'FontSize',12);
    xlabel('Time (s)')
    ylabel('Speed (cm/s)')
    
    fullpath = [path filename '_SpeedLick_Tr' num2str(trialNo(1)) '_' num2str(trialNo(end))];
    print('-painters', '-dpdf', fullpath, '-r600')
    savefig([fullpath '.fig']);
    
    figure
    plot((1:endInd-startInd+1)/sampleFq,Track.F(startInd:endInd)-min(Track.F(startInd:endInd)),'k-');
    set(gca,'XLim',[0 (endInd-startInd+1)/sampleFq],'YLim',[0 4],'FontSize',12);
    xlabel('Time (s)')
    ylabel('F')
    
    fullpath = [path filename '_F_Tr' num2str(trialNo(1)) '_' num2str(trialNo(end))];
    print('-painters', '-dpdf', fullpath, '-r600')
    savefig([fullpath '.fig']);
    
    figure
    F = Track.F(startInd:endInd);
    sp = speed >= 20;
    [continuousRun,stopRun] = numOfConsecutiveOnes(sp);
    idxLongStop = find(stopRun > 10*sampleFq);
    if(~isempty(idxLongStop))
        indFirstRun = find(sp == 1,1);
        
        if(indFirstRun > 1)
            indStart = sum(continuousRun(1:idxLongStop(1)-1)) + ...
                sum(stopRun(1:idxLongStop(1)-1))+1;
        else
            indStart = sum(continuousRun(1:idxLongStop(1))) + ...
                sum(stopRun(1:idxLongStop(1)-1))+1;
        end
        F0 = mean(Track.F(indStart+1:indStart+stopRun(idxLongStop(1))-1));
        dFF = (F-F0)/F0; 
    else
        return;
    end
    plot((1:endInd-startInd+1)/sampleFq,dFF,'k-');
    set(gca,'XLim',[0 (endInd-startInd+1)/sampleFq],'YLim',[0 0.2],'FontSize',12);
    xlabel('Time (s)')
    ylabel('dF/F')
    
    fullpath = [path filename '_dFF_Tr' num2str(trialNo(1)) '_' num2str(trialNo(end))];
    print('-painters', '-dpdf', fullpath, '-r600')
    savefig([fullpath '.fig']);
    
end

function [data,data1] = numOfConsecutiveOnes(arr)
    data = [];
    data1 = [];
    s = sprintf('%d', arr);
    %Reading the consequences of 1's from the string by using 0's as delimiters
    t1=textscan(s,'%s','delimiter','0','multipleDelimsAsOne',1);
    % Converting cell array of cell into a single cell array
    d = t1{:};
    % Computing the length of each run by going through the array and assigning
    % it into 
    for k = 1:length(d)
          data(k) = length(d{k});
    end
    
    t2=textscan(s,'%s','delimiter','1','multipleDelimsAsOne',1);
    % Converting cell array of cell into a single cell array
    f = t2{:};
    % Computing the length of each run by going through the array and assigning
    % it into 
    for k = 1:length(f)
          data1(k) = length(f{k});
    end
end
