function plotBarOnly(mean,std,x1,y1,t,pathAnal,fn)
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
    
    xlabel(x1);
    ylabel(y1);
    title(t);
    
    fileName1 = [pathAnal fn];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end
