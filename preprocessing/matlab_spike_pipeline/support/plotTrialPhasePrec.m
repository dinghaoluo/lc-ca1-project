function plotTrialPhasePrec(path, fileName, neuNo)
% plot phase precession for a single cell

    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end

    fullPath = [path fileName '.mat']; 
    if(exist(fullPath) == 0)
        disp('The file does not exist');
        return;
    end
    load(fullPath);
    
    GlobalConst;
    
    count = 0;
    for i = 1:length(trials)
        count = count + 1;

        if(mod(count-1,16) == 0)
            [figNew,pos] = CreateFig();
            set(0,'Units','pixels') 
            figTitle = 'Phase';
            set(figure(figNew),'OuterPosition',...
                [pos(1) pos(2) pos(3)*2 pos(4)*2.2],'Name',figTitle)
        end
        subplot(4,4,mod(count-1,16)+1)  
        ind = trials{i}.spikesThetaLin{neuNo} < 0;
        spikeThetaLin = trials{i}.spikesThetaLin{neuNo};
        spikeThetaLin(ind) = spikeThetaLin(ind) + 2*pi;
        plot(trials{i}.spikes{neuNo}*timeStep,spikeThetaLin,'.',...
            'MarkerSize',8);
        set(gca,'Xlim',[0 trials{i}.Nsamples*timeStep],...
            'Ylim',[0,2*pi]);
        hold on;
        if(mod(count-1,4) == 0)
            ylabel('Phase');
        end
        if(mod(count-1,16) > 11)
            xlabel('Dist (mm)');
        end
        yyaxis right

        plot((1:trials{i}.Nsamples)*timeStep,...
                trials{i}.speed,'r');
        set(gca,'Xlim',[0 trials{i}.Nsamples*timeStep]);
        
        title(['Neu' num2str(neuNo) 'Tr' num2str(i)]);        
    end
    
%         ax = axes('Position',get(gca,'Position'),...
%                            'XAxisLocation','top',...
%                            'YAxisLocation','right',...
%                            'XTickLabel',[],...
%                            'Xlim',[0 trials{i}.Nsamples*timeStep],...
%                            'YLim',[0 500],...
%                            'Color','none',...
%                            'XColor','k','YColor','k',...
%                            'FontSize',14);
%          axes(ax);