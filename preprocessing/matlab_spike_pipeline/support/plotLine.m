function plotLine(x,y,figTitle,xLabel,yLabel)    
% create a figuer and plot a line 
%
% Inputs:
% x and y:          2 vectors with equal length
% figTitle:         figure title
% xLabel:           x label of the figure
% yLabel:           y label of the figure

    [figNew,pos] = CreateFig();
    set(figure(figNew),'OuterPosition',pos,'Name',figTitle); 

    h = plot(x,y,'ko');
    set(h,'LineWidth',2.0);
    [rho, pVal] = corr(x',y','type','Spearman'); % spearman correlation between x and y

    set(gca,'FontSize',14.0,'Box','on');
    xlabel(xLabel)
    ylabel(yLabel)
    title(['r = ' num2str(rho) ' p = ' num2str(pVal)]);