function [thetaTimeMulti,thetaSpikesMulti] = getMultiCycles(thetaTime, thetaPhase, numCycles)
% construct spike trains which repete the neuronal activity for
% multiple phase cycles
% thetaTime:        time of spikes
% thetaPhase:       phases of spikes

indexThetaPos = find(thetaPhase >= 0);
indexThetaNeg = find(thetaPhase < 0);
thetaSpikesDeg = thetaPhase*360/(2*pi);
thetaTimeTmp = [];
thetaPhaseTmp = [];
for i = 1:numCycles-1
    thetaTimeTmp = [thetaTimeTmp; thetaTime];
    thetaPhaseTmp = [thetaPhaseTmp; thetaSpikesDeg+i*360];
end
thetaTimeMulti = [thetaTime(indexThetaPos); thetaTimeTmp; thetaTime(indexThetaNeg)];
thetaSpikesMulti = [thetaSpikesDeg(indexThetaPos); thetaPhaseTmp; thetaSpikesDeg(indexThetaNeg)+numCycles*360];