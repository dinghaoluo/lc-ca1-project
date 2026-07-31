% initialize constants
sampleFq = 500; %1220.703125;
sampleFqOri = 30.01; %24414.0625;
timeStep = 1/sampleFq; 
timeStepOri = 1/sampleFqOri;

timeBin = 0.03; % s % smooth factor in spike convolution in time
spaceBin = 20; % mm % smooth factor in spike convolution in distance

minSpeed = 20; % the min running speed, mm/s
minSpeed1 = 100; % mm/s
spaceMergeBin = 1; % mm
meanSpeed = 100; % mean running speed, mm/s

minFR = 0.01;

trialLenT = 20; % sec
smoothSpan = 100;

startTrNo = 5; % do not consider the first 5 trials due to the imaging drift at the beginning of the session

minNumTrials = 10; % the minimum number of trials in which the neuron fires

% align to different behavioral parameters
nSampBef = 3*sampleFq; % used in align to running onset
intervalTSpInfo = 10; % sec
corrIntervalT = 20; % sec
corrIntervalTMin = -10;
corrIntervalD = 1800; % mm
tc = 500; % ms
cost = 2/(tc/1000*sampleFq);
trackLen = 2000; % track length used in convultion distance
