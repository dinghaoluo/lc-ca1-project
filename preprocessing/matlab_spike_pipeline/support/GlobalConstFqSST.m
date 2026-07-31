sampleFq = 1250;

methodTheta = 1;
minFR = 0.15;
maxFR = 7;

spaceBin = 20;
intervalT = 10;
intervalTPopCorr = 20;

minFRInt = 3;

nSampBef = 3*sampleFq; % used in align to running onset

minNumGoodTr = 15;

anmNoInact = [49, 62, 63, 64, 67 75 77 81 82 86 87 92 95 96];
anmNoAct = [28 29 30]; % for pulse method 2
pulseMethod{1} = [2 3 4]; % inactivation
pulseMethod{2} = [2 3 4]; % activation

% significance probability
p = [99.9,99,95];
numShuffle = 1000;