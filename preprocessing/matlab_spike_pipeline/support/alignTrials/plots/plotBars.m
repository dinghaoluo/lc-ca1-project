function plotBars(data1,data2,mean,std,x1,y1,t,pathAnal,fn)
    fig = figure;
    set(0,'Units','pixels') 
    set(figure(fig),'OuterPosition',...
        [500 500 210 280])
    fig.Renderer = 'Painters';
    
    h = bar([1,2],mean,0.5);
    set(h,'EdgeColor',[0.3 0.3 0.3],'FaceColor',[187 189 192]/255);
    hold
    
    h = errorbar([1,2],mean,std);
    set(h,'Marker','.','MarkerSize',0.1,'Color',[0 0 0],'LineStyle','none')
    
    h = plot(1+0.15*rand(1,length(data1)),data1,'o');
    set(h,'MarkerSize',3,'Color',[167 169 171]/255);
    
    h = plot(2+0.15*rand(1,length(data2)),data2,'o');
    set(h,'MarkerSize',3,'Color',[27 117 187]/255);
    set(gca,'XLim',[0.5 2.5],'YLim',[0 max([data1 data2])+0.01]);
    
    xlabel(x1);
    ylabel(y1);
    title(t);
    
    fileName1 = [pathAnal fn];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end
