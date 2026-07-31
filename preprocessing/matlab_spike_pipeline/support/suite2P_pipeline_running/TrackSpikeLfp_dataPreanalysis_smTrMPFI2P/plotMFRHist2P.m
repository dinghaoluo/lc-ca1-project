function plotMFRHist(mFRArr)
% plot the mean firing rate histogram
%
% mFR:          mean firing rate array

    [figNew,pos] = CreateFig();
    set(figure(figNew),'OuterPosition',pos,'Name','Mean firing rate'); 

    [neuronCount,frPos] = hist(mFRArr,20);
    hist(mFRArr,20);
    
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','k','EdgeColor','k');
    peakNeuronCount = max(neuronCount);
    medianFiringRate = median(mFRArr);
    hold on
    h = plot([medianFiringRate,medianFiringRate],[peakNeuronCount*1.15,peakNeuronCount*1.05]);
    set(h,'Color',[0 0 0],'LineWidth',3.0);
    h = plot(medianFiringRate,peakNeuronCount*1.05,'v');
    set(h,'Color',[0 0 0],'MarkerFaceColor',[0 0 0]);

    set(gca,'FontSize',14.0,'Box','on','YLim',[0,ceil(peakNeuronCount*1.2)]);
    xlabel('Mean firing rate (Hz)')
    ylabel('N of neurons')