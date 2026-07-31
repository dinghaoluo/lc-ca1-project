function plotMFR(numNeurons,mFRArr,stdMFRArr)
% plot mean firing rate
%
% numNeurons:       number of neurons
% mFRArr:           mean firing rate 
% stdMFRArr:        std of mean firing rate
% mFRArrErr:        mean firing rate of the error trials
    
    [figNew,pos] = CreateFig();
    set(figure(figNew),'OuterPosition',pos,'Name','Mean firing rate'); 

    % plot the mean and std of the firing rate of each neuron
    h = errorbar(1:numNeurons,mFRArr,stdMFRArr,'ko');
    set(h,'LineWidth',2.0);

    xLim = [0.5 numNeurons+0.5];
    set(gca,'XLim',xLim, 'FontSize',14.0,'Box','on')
    xlabel('Neurons')
    ylabel('Mean firing rate (Hz)')