function barPlot(ind,x,y,xlab,ylab,tit)
    figure;
    h = bar(ind,x);
    set(h,'faceColor',[0.7,0.7,0.7],'lineWidth',0.5);
    hold on;
    h = errorbar(ind,x,y,'.');
    set(h,'Color',[0 0 0],'lineWidth',1);
    xlabel(xlab);
    ylabel(ylab);
    if(~isempty(tit))
        title(tit);
    end
    set(gca,'fontSize',12)
    
end
