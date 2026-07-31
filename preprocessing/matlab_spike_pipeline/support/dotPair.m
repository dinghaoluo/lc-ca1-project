function dotPair(x,xlab,ylab,tit)
% each column of x is data from one subsession
    figure;
    h = plot(x(:,1),x(:,2),'.');
    set(h,'MarkerSize',9,'Color','r');
    hold on;
    allX = x(:);
    h = plot([min(allX),max(allX)],[min(allX),max(allX)],'b-');
    set(h,'LineWidth',0.5);
    set(gca,'XLim',[min(allX) max(allX)],'YLim',[min(allX) max(allX)]);
	xlabel(xlab);
    ylabel(ylab);
    if(~isempty(tit))
        title(tit);
    end
    set(gca,'fontSize',12)
    
end
