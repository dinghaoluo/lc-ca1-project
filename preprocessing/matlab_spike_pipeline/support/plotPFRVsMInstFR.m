function plotPFRVsMInstFR(pFRArr,meanInstFRArr)
% plot peak firing rate vs mean instantaneous firing rate
% pFRArr:           peak firing rate 
% meanInstFRArr:    mean instantaneous firing rate

    [figNew,pos] = CreateFig();
    set(figure(figNew),'OuterPosition',pos,'Name','Peak vs mean inst FR'); 
    fid = gca;
    h = plot(meanInstFRArr,pFRArr,'*');
    set(h,'LineWidth',2.0,'Color',[0.502 0.502 0.502]);
    set(fid, 'FontSize',14.0,'Box','on')
    xlabel('Mean FR (Hz)')
    ylabel('Peak FR (Hz)')