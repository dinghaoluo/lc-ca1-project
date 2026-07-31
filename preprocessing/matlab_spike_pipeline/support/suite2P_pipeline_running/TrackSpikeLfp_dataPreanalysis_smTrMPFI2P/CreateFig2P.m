function [figNew,pos] = CreateFig2P()
% Create a new figure and display it at the defined location on the screen
% return
% fidNew:       figure handle
% pos:          figure position

    set(0, 'DefaultFigureRendererMode', 'manual')
    set(0,'DefaultFigureRenderer','zbuffer')
    figNew = figure;
    fid = gca;
    scnsize = get(0,'ScreenSize');

    % Obtain the entire window size and position
    outerpos = get(figNew,'OuterPosition');
    figurecount = length(findobj('Type','figure')); 
    horiShiftStep = scnsize(3)/16;
    verShiftStep = scnsize(4)/8;
    pos = [mod((figurecount-1),12)*horiShiftStep,...
            1600-504-floor((figurecount-1)/8)*verShiftStep,...
            outerpos(3),...
            outerpos(4)];
    set(gcf,'Doublebuffer','on')