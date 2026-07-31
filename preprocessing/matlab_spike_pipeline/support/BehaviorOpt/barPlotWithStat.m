function barPlotWithStat(ind,x,y,data,xlab,ylab,tit,rankXY,rankXZ,rankYZ)
    figure;
    h = bar(ind,x);
    set(h,'faceColor',[0.7,0.7,0.7],'lineWidth',0.5);
    hold on;
    h = errorbar(ind,x,y,'.');
    set(h,'Color',[0 0 0],'lineWidth',1.5);
    h = plot([1.1 2.1],data(:,1:2),'k.-');
    set(h,'MarkerSize',6,'LineWidth',1,'Color',[0.5 0.5 0.5]);
    if length(ind) > 2
        h = plot([2.1 3.1],data(:,2:3),'k.-');
        set(h,'MarkerSize',6,'LineWidth',1,'Color',[0.5 0.5 0.5]);
    end
    text(1,1.1*max(max(data(:,1:2))),num2str(rankXY,'p = %f'));
    if length(ind) > 2
        text(2.2,1.1*max(max(data(:,2:3))),num2str(rankYZ,'p = %f'));
        text(1.5,1.25*max(data(:)),num2str(rankXZ,'p = %f'));
    end
    xlabel(xlab);
    ylabel(ylab);
    if(~isempty(tit))
        title(tit);
    end
    if(max(data(:)) > 0)
        set(gca,'fontSize',12,'YLim',[min(data(:)) 1.2*max(data(:))])
    elseif((max(data(:)) < 0))
        set(gca,'fontSize',12,'YLim',[1.5*min(data(:)) max(data(:))])
    else
        set(gca,'fontSize',12,'YLim',[min(data(:)) 1])
    end
    
end
