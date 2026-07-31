function [sp,ac] = MazeSpeedAccel(whldata,sampleFq)

    % function [speed, accel] = MazeSpeedAccel(whldata,smoothlen)
    % Calculates running speed and acceleration from X position data. 

    %calculate speed for values that aren't -1 or distorded by filtering
    %---------------------------------
    timeLfp = (1:length(whldata))/sampleFq;
    timeLfp = timeLfp';
    speeddata = diff(whldata)./diff(timeLfp); % mm/s
    acdata = [diff(speeddata(1:2)); diff(speeddata)]./diff(timeLfp);
    sp = [speeddata(1); speeddata];
    ac = [acdata(1); acdata];
end