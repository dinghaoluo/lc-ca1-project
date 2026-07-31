function plotEEGandSpeed(path,filename)
% plot speed, EEG and spikes from all the neurons

    load([path filename '.mat']);
    
    fileNameInfo = [filename '_Info.mat'];    
    fullPath = [path fileNameInfo];
    load(fullPath);
    indGoodLap = find(beh.pauseWithinTrial < 0.5 & ...
                    beh.pauseBefTrial > 0.2);
    indGoodLap = intersect(indGoodLap, 93:137);
    indGoodLap1 = intersect(indGoodLap, beh.indGoodLap);
    goodLapSess = beh.mazeSess(indGoodLap1);
    indSessBorder = find(diff(goodLapSess) == 1);   
    
    lfpSampleRate = 1250;
    % theta filtering parameter
    thetaFiltParam.FreqRange = [4 10]; 
        %originally [4 25], changed to [4 16] on 2016.07.14
    thetaFiltParam.FilterOrd = 4;
    thetaFiltParam.Ripple = 20;
    [b a] = cheby2(thetaFiltParam.FilterOrd, thetaFiltParam.Ripple, ...
                thetaFiltParam.FreqRange/(lfpSampleRate/2),'bandpass');
              
%     int = find(autoCorr.isInterneuron == 1);
%     pyr = find(autoCorr.isPyrneuron == 1);
%     indNeurons = [20 35 37 83 88 89 101 103 105]; % A004-20181028
%     indNeurons = [17 20 48 49 56 76 96 108];
%     indNeurons = [21 22 30 31 33 42 44 45 46 52 66 68 73 76 79 91 92 103 105 111 129 134 143 163];
    indNeurons = 85;
    
    for i = 1:35 %indGoodLap1    
        % plot spikes
        numNeurons = 0;
        figure(20)
        subplot(3,1,1)
        while(ishold)
            hold off;
        end
        
        % plot eeg
        Eegf = filtfilt(b,a,trials{i}.eeg);
        Eegf = Eegf - mean(Eegf);
        plot((1:trials{i}.Nsamples)/lfpSampleRate,Eegf);
        hold on;
        
        % plot speed
        plot((1:trials{i}.Nsamples)/lfpSampleRate, trials{i}.speed/1000-0.25,'r');
        set(gca,'XLim',[0 trials{i}.Nsamples/lfpSampleRate]);
        xlabel('time (s)')
        title(['Trial No.' num2str(i)]);
        
        for j = indNeurons
            h = plot(trials{i}.spikes{j}/lfpSampleRate,...
                    0.3+0.03*numNeurons*ones(1,length(trials{i}.spikes{j})),'k.');
            set(h,'MarkerSize',8);
            numNeurons = numNeurons + 1;
        end
        
        subplot(3,1,2)
        plot((1:trials{i}.Nsamples)/lfpSampleRate, trials{i}.xMM,'b');
        set(gca,'XLim',[0 trials{i}.Nsamples/lfpSampleRate]);
        
        subplot(3,1,3)
        while(ishold)
            hold off;
        end
        for j = indNeurons
            indSpikes = trials{i}.spikesThetaLin{j} < 0;
            spikesTheta = trials{i}.spikesThetaLin{j};
            spikesTheta(indSpikes) = spikesTheta(indSpikes) + 2*pi;
            h = plot(trials{i}.spikes{j}/lfpSampleRate,...
                        spikesTheta,'.');
            if(j == indNeurons(1))
                hold on;
            end
            set(gca,'XLim',[0 trials{i}.Nsamples/lfpSampleRate]);
        end
            
        pause;
    end
end
