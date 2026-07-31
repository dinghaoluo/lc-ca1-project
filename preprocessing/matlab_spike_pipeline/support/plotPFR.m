function plotPFR(numNeurons,pFRArr,meanInstFRArr)
% plot peak firing rate of all the neurons
% numNeurons:       number of neurons
% pFRArr:           peak firing rate 
% meanInstFRArr:    mean instantaneous firing rate

    [figNew,pos] = CreateFig();
    set(figure(figNew),'OuterPosition',pos,'Name','Peak firing rate'); 

    fidsub = gca;
    h = plot(1:numNeurons,pFRArr,'ko');
    set(h,'LineWidth',2.0);
    hold on
    h = plot(1:numNeurons,meanInstFRArr,'r*');
    set(h,'LineWidth',2.0);
    xLim = [0.5 numNeurons+0.5];
    set(fidsub,'XLim',xLim, 'FontSize',14.0,'Box','on')
    ylabel('Firing rate (Hz)')
    xlabel('Neurons')