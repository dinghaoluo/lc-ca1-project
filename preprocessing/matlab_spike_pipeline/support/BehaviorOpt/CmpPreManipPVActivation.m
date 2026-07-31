cd Z:\Yingxue\DataAnalysisDongyan
load allRecData.mat

i = 20;

%% mean speed (seems to increase)
a = [recDataRunPre{i}.meanSpeedMean;recDataRunManip{i}.meanSpeedMean];
figure
plot(1:size(a,2),a)
figure
histogram(a(1,:))
figure
histogram(a(2,:))

%% num run segments
a = [recDataRunPre{i}.numRunMean;recDataRunManip{i}.numRunMean];
figure
plot(1:size(a,2),a)
figure
histogram(a(1,:),1:0.3:7)
figure
histogram(a(2,:),1:0.3:7)

%% acceleration
a = [recDataRunPre{i}.meanAccMean;recDataRunManip{i}.meanAccMean];
figure
plot(1:size(a,2),a)
figure;
histogram(a(1,:))
figure
histogram(a(2,:))

%% total stop time (increase?)
a = [recDataRunPre{i}.totStopLenTMean;recDataRunManip{i}.totStopLenTMean];
figure
plot(1:size(a,2),a)
figure;
histogram(a(1,:),0:0.1:1)
figure
histogram(a(2,:),0:0.1:1)

%% start cue to run time
a = [recDataRunPre{i}.startCueToRunMean;recDataRunManip{i}.startCueToRunMean];
figure
plot(1:size(a,2),a)
figure;
histogram(a(1,:),-1:0.2:2)
figure
histogram(a(2,:),-1:0.2:2)

%% number of lick before reward (decrease?)
a = [recDataRunPre{i}.numLicksBefRewMean;recDataRunManip{i}.numLicksBefRewMean];
figure
plot(1:size(a,2),a)
figure;
histogram(a(1,:),0:30:500)
figure
histogram(a(2,:),0:30:500)

%% number of lick reward (decrease)
a = [recDataRunPre{i}.numLicksRewMean;recDataRunManip{i}.numLicksRewMean];
figure
plot(1:size(a,2),a)
figure;
histogram(a(1,:),0:30:500)
figure
histogram(a(2,:),0:30:500)

%% med first 5 lick distance (increase)
a = [recDataRunPre{i}.med1stFiveLickDistMean;recDataRunManip{i}.med1stFiveLickDistMean];
figure
plot(1:size(a,2),a)
figure;
histogram(a(1,:),500:50:1800)
figure
histogram(a(2,:),500:50:1800)

%% med first 5 lick distance before reward (increase)
a = [recDataRunPre{i}.med1stFiveLickDistBefRewMean;recDataRunManip{i}.med1stFiveLickDistBefRewMean];
figure
plot(1:size(a,2),a)
figure;
histogram(a(1,:),500:50:1800)
figure
histogram(a(2,:),500:50:1800)

%% perc non-stopping
a = [recDataRunPre{i}.percNonStop;recDataRunManip{i}.percNonStop];
figure
plot(1:size(a,2),a)
figure;
histogram(a(1,:))
figure
histogram(a(2,:))

%% speed similarity (increase?)
a = [recDataRunPre{i}.speedSimMean;recDataRunManip{i}.speedSimMean];
figure
plot(1:size(a,2),a)
figure;
histogram(a(1,:))
figure
histogram(a(2,:))

%% lick similarity 
a = [recDataRunPre{i}.lickSimMean;recDataRunManip{i}.lickSimMean];
figure
plot(1:size(a,2),a)
figure;
histogram(a(1,:),0:0.02:0.5)
figure
histogram(a(2,:),0:0.02:0.5)

%% speed euclidean
a = [recDataRunPre{i}.speedEucMean;recDataRunManip{i}.speedEucMean];
figure
plot(1:size(a,2),a)
figure;
histogram(a(1,:))
figure
histogram(a(2,:))

%% lick euclidean 
a = [recDataRunPre{i}.lickEucMean;recDataRunManip{i}.lickEucMean];
figure
plot(1:size(a,2),a)
figure;
histogram(a(1,:))
figure
histogram(a(2,:))

%% lick 30 to 100 cm
figure;
histogram(recDataRunPre{i}.meanRun30to100,0:0.2:5)
figure
histogram(recDataRunManip{i}.meanRun30to100,0:0.2:5)

%% lick 100 to 150 cm
figure;
histogram(recDataRunPre{i}.meanRun100to150,0:0.2:10)
figure
histogram(recDataRunManip{i}.meanRun100to150,0:0.2:10)

%% lick 150 to 180 cm (less?)
figure;
histogram(recDataRunPre{i}.meanRun150to180,0:0.2:30)
figure
histogram(recDataRunManip{i}.meanRun150to180,0:0.2:30)

%% speed 0 to 100 cm (less?)
figure;
histogram(recDataRunPre{i}.meanSpeedOverDistRun0to100,0:5:150)
figure
histogram(recDataRunManip{i}.meanSpeedOverDistRun0to100,0:5:150)

%% speed after 100 cm (increase?)
figure;
histogram(recDataRunPre{i}.meanSpeedOverDistRunAfter100,0:5:150)
figure
histogram(recDataRunManip{i}.meanSpeedOverDistRunAfter100,0:5:150)