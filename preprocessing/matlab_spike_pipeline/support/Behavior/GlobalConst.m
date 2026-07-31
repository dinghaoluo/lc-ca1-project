% initialize constants
sampleFq = 1250; %1220.703125;
sampleFqOri = 20000; %24414.0625;
timeStep = 1/sampleFq; 
timeStepOri = 1/sampleFqOri;

nSampBef = 3*sampleFq; % used in align to running onset
nSampBefRew = 5*sampleFq;

minSpeed = 20; % the min running speed, mm/s
spaceMergeBin = 1; % mm
meanSpeed = 100; % mean running speed, mm/s

trialLenT = 20; % sec
smoothSpan = 100;

stopTimeRew2RunThrL = 1;
stopTimeRew2RunThrH = 3;

% firing rate criteria
minPeakFR = 3; % the minimum peak firing rate (Hz)
minSDTimes = 4.5; % the peak firing rate is larger than minSDTimes*SD+Mean 
minPeak2MeanInstRatio = 8; % the ratio between the peak and mean instantaneous firing rate (used as a measure of the existence of phase lock)
minPeak2MeanInstRatioFW = 3.5; % the ratio between the peak and mean instantaneous firing rate (used as a measure of the existence of a field)
                             % 4: when timeBin = 0.2s, and 6: when timeBin = 0.1s
minNumTrials = 10; % the minimum number of trials in which the neuron fires

% excitatory neuron firing rate
minFR = 0.15; % the minimum mean firing rate (Hz)
maxFR = 7; % the minimum mean firing rate (Hz)
maxFRPCA = 1.5; % the max mean firing rate of the neurons included in the PCA calculation
maxFREP = 6; % the max mean firing rate of episode cells

% cue location
cue = [82,132; 7,57; 57,82; 132,7];
switchMap = [2 1 4 3];

timeAftStim = 0.5;
timeBefStim = 0.5;

% theta phase
minSamplePerCycle = 0.05/timeStep; % the minimum length of a theta cycle
maxSamplePercycle = 0.2/timeStep; % the maximum length of a theta cycle
thetaPhaseJump = -4.5; % the min sudden change in theta phase which indicates the start of a new cycle

% field width
fieldBound = 0.1; % 10% of the peak defines the boundary of the field

% oscillation frequency related constants
stdOscFreq = 0.01*sampleFq;
dpssOrder = 4;
% theta band
upperCrossCorrTheta = 11; % the upper bound theta frequency for cross correlation between unit and LFP
lowerCrossCorrTheta = 5.5; % the lower bound theta frequency for cross correlation between unit and LFP
% gamma band
upperCrossCorrGamma = 100; % the upper bound theta frequency for cross correlation between unit and LFP
lowerCrossCorrGamma = 40; % the lower bound theta frequency for cross correlation between unit and LFP
% calculating the change of spectrum over time
lenSptmWin = 0.5; %(sec) the length of window 
numSampleSptmWin = floor(lenSptmWin*sampleFq);
numSampleOverlap = ceil(numSampleSptmWin/2); % the overlap between windows during spectrum calculation
nfftSample = 10*sampleFq; % nfft of the spectrum

% burst isi
burstIsi = 0.012*sampleFq; % the max isi between a pair of spikes which can be considered as a burst
burstIsi1st = 0.006*sampleFq; % the max isi between the first two spikes in a burst
burstSpeedThre = 100; % running speed limit 

% debug
dispDebug = 1; % display debug details

% mahala distance and refrectory period
mahalDistThre = 5; % 15  the threshold of the mahalDist of neuron clusters (last modification 04/22/2013)
refracViolPercentThre = 1; % 0.35 the threshold of the percentage of refractory period violation (last modification 04/22/2013)
centerMaxThre = -150; % threshold for the max of the waveshape center peak (last modification 04/22/2013)

% length of the trial used in calculation
fLenTrial = 1;     % 1: using the longest trial length 
                   % 0: using the shortest trial length
                   
% maximum speed in the wheel and in the maze
maxSpeedM = 3000;  % mm/sec
maxSpeedWhl = 3000; % clicks/sec

backward = 0; % for the return segment in the maze, align to the end of run = 1, align to the beginning of run = 0

threNumNeuronsWh = 5; %9 the min number of episode/place cells within a session before this session is considered a good session
threNumNeuronsMz = 5; %6 the min number of episode/place cells within a session before this session is considered a good session
commPrePostNumNeu = -1; % the min number of neurons which have fields in both pre and post conditions (10.13.2013)

MMperPix = 2.485; % mm per pixel in the maze (2014.10.16: this value is not correct. This value should be 4.2333 on  the floor, 3.9329 on the floor in the maze with wheel, 4.2333 on the top of the maze, See MMperPixTest folder)
occupThre = 3; % each spatial bin occupies how many samples

recAnm = [{1:3}, {4:9}, {10:13}]; % recording numbers for each animal
recAnmSaline = [{}, {1:4}, {5:7}]; % recording numbers for each animal (saline)

numCycleSeqThre = 100; %100; % the minimum number of cycles. if there aren't enough number of sequence cycles, then the recording is neglected
percCycleSeqThre = 0.005; % the minimum percent of cycles. if there aren't enough number of sequence cycles, then the recording is neglected
ratioPercSeqControlMusc = 7; % the percentage of sequences under muscimol should not be substantially larger than the control (pre/post), and it should also be true the other way