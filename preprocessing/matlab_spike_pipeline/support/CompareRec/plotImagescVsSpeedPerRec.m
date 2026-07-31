cd Z:\Yingxue\Draft\PV\Interneuron
load autoCorrIntAllRec.mat
ind = find(autoCorrIntAll.idxC2 == 1);
[autoCorrIntAll.task(ind);autoCorrIntAll.indRec(ind);autoCorrIntAll.indNeu(ind)]

% cd Z:\Raphael_tests\mice_expdata\ANM012\A012-20190224\A012-20190224-01\
% load('A012-20190224-01_DataStructure_mazeSection1_TrialType1_Info','autoCorr');
% load('A012-20190224-01_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat','trialsRun');
% % load ('A012-20190224-01_DataStructure_mazeSection1_TrialType1_alignedSpikesPerNPerT_msess1_Run0.mat',...
% %     'trialsRunSpikes');
% load('A012-20190224-01_DataStructure_mazeSection1_TrialType1_convSpikesAligned_msess1_BefRun0');
% load('A012-20190224-01_DataStructure_mazeSection1_TrialType1_PeakFR_msess1_RunOnset0','pFRNonStimGoodStruct');
% load('A012-20190224-01_DataStructure_mazeSection1_TrialType1_behPar_msess1')

cd Z:\Raphael_tests\mice_expdata\ANM037\A037-20201221\A037-20201221-01\
load('A037-20201221-01_DataStructure_mazeSection1_TrialType1_Info','autoCorr');
load('A037-20201221-01_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat','trialsRun');
load('A037-20201221-01_DataStructure_mazeSection1_TrialType1_convSpikesAligned_msess1_BefRun0');
load('A037-20201221-01_DataStructure_mazeSection1_TrialType1_PeakFR_msess1_RunOnset0','pFRNonStimGoodStruct');
load('A037-20201221-01_DataStructure_mazeSection1_TrialType1_FieldSpCorrAligned_Run1_Run1','fieldSpCorrSessNonStimGood');
load('A037-20201221-01_DataStructure_mazeSection1_TrialType1_behPar_msess1')

cd Z:\Raphael_tests\mice_expdata\ANM046\A046-20210423\A046-20210423-01
load('A046-20210423-01_DataStructure_mazeSection1_TrialType1_Info','autoCorr');
load('A046-20210423-01_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat','trialsRun');
load('A046-20210423-01_DataStructure_mazeSection1_TrialType1_convSpikesAligned_msess1_BefRun0');
load('A046-20210423-01_DataStructure_mazeSection1_TrialType1_PeakFR_msess1_RunOnset0','pFRNonStimGoodStruct');
load('A046-20210423-01_DataStructure_mazeSection1_TrialType1_FieldSpCorrAligned_Run1_Run1','fieldSpCorrSessNonStimGood');
load('A046-20210423-01_DataStructure_mazeSection1_TrialType1_behPar_msess1')

paramC.trialLenT = 20; %sec
paramC.timeBin = 0.0025; %sec
std = 0.15/paramC.timeBin;
paramC.gaussFilt = gaussFilter(12*std, std);
lenGaussKernel = length(paramC.gaussFilt);
normFactor = sum(paramC.gaussFilt);
paramC.gaussFilt = paramC.gaussFilt./normFactor;
avgFRProfile = conv2(pFRNonStimGoodStruct.avgFRProfile,paramC.gaussFilt,'same');
 
cosSim = zeros(1,length(filteredSpikeArrayRunOnSet));
indTime = timeStepRun/1250<7;
for i = 1:length(filteredSpikeArrayRunOnSet)
    tmp = pdist2(filteredSpikeArrayRunOnSet{i}(:,indTime), avgFRProfile(i,indTime),'spearman');
    cosSim(i) = mean(tmp(~isnan(tmp))); 
end

medTrLen = prctile(behPar.numSamplesRun,75)/1250;
indTime = find(timeStepRun/1250 >= medTrLen,1);
isPeakNeuArr = zeros(1,length(filteredSpikeArrayRunOnSet));
for i = 1:length(filteredSpikeArrayRunOnSet)
    isPeakNeuArr(i) = neuPeakDetection(filteredSpikeArrayRunOnSet{i}(:,1:indTime),paramC,i);
end

%% plot individual neurons
figure
for i = 1:length(filteredSpikeArrayRunOnSet)
    yyaxis left
    imagesc(timeStepRun/1250,1:size(filteredSpikeArrayRunOnSet{1},1),filteredSpikeArrayRunOnSet{i});
    title(['Neu ' num2str(i)]);
    xlabel('Time (s)')
    ylabel('Neuron no. ')
    
    yyaxis right
%     hold off;
%     plot(timeStepRun/1250,pFRNonStimGoodStruct.avgFRProfile(i,:),'g');
%     hold on;
    plot(timeStepRun/1250,avgFRProfile(i,:));
    ylabel('FR (Hz)')
    set(gca,'XLim',[-3 10]);
    
    pause;
end

%% PV interneurons
neu = 32;... %32    42    58 --- A012-20190224-01
neu = 69; % 69 --- A037-20201221-01
neu = 47; % 47   105   124 --- A046-20210423-01
figure(1);
for i = 2:length(trialsRun.speed_MMsecBef)
    hold off;
    speed = [trialsRun.speed_MMsecBef{i}; ...
        trialsRun.speed_MMsec{i}];
    speed(speed<0) = 0;
    plot((1:length(speed))/1250-3,speed/max(speed));
    hold on;
    
%     timeSp = [trialsRunSpikes.TimeBef{neu,i}; trialsRunSpikes.Time{neu,i}];
%     timeStep = -3:0.02:length(trialsRun.speed_MMsec{i})/1250;
%     histX = hist(timeSp/1250,timeStep); %91 119
%     plot(timeStep,histX/max(histX));
%     set(gca,'XLim',[0 6])

    plot(timeStepRun/1250,filteredSpikeArrayRunOnSet{neu}(i,:)./...
        max(filteredSpikeArrayRunOnSet{neu}(i,:)));
    pause;
end

% A012-20190224-01
% Neu 32 Tr 65 88 89
% Neu 42 Tr 16? 46? 53(sync start run)? 66 73 89
% Neu 58 Tr 15? 19? 30? 34? 53? 78? 79? 90?

% A037-20201221-01

%% Pyramidal neurons
figure(2)
figure(3)
figure(4)
for i = 2:length(trialsRun.speed_MMsecBef)
    figure(2)
    hold off
    speed = [trialsRun.speed_MMsecBef{i}; ...
        trialsRun.speed_MMsec{i}];
    speed(speed<0) = 0;
    plot((1:length(speed))/1250-3,speed);
    maxSpeed = max(speed);
    lick = [trialsRun.lickLfpInd{i}-trialsRun.startLfpInd(i)+1];
    hold on
    plot(lick/1250,maxSpeed*ones(1,length(lick)),'ro');
    plot(timeStepRun/1250,filteredSpikeArrayRunOnSet{neu}(i,:)*5);
    set(gca,'XLim',[-3 10])
    title(['Tr ' num2str(i) ' isStim = ' num2str(behPar.stimOn(i))])
    
    figure(3)
    indZero = find(timeStepRun >= 0,1);
%     indPyr = find(autoCorr.isPyrneuron == 1 & ...
%         pFRNonStimGoodStruct.meanInstFR > 0.1 & ...
%         pFRNonStimGoodStruct.meanInstFR < 5 & ...
%         pFRNonStimGoodStruct.p2MInstRatio > 7);
    indPyr = [2 6 8 12 14 16 17 25 27 32 37 38 43 45 54 ...
    57 59 61 62 64 69 70 73 81 85 90 92 94 99 101 103 108 109 115 128 133]; %A046-20210423-01
%     indPyr = [2 3 8 12 13 14 19 23 30 31 42 48 67 68 71 73 81 89 102]; %A037-20201221-01
    indPyr =  find(autoCorr.isPyrneuron == 1 & ...
        pFRNonStimGoodStruct.meanInstFR > 0.1 & ...
        pFRNonStimGoodStruct.meanInstFR < 7 & ...
        isPeakNeuArr == 1);
    pyrProfAll = zeros(length(indPyr),length(timeStepRun));
    m = 1;
    for n = indPyr
        pyrProfAll(m,:) = filteredSpikeArrayRunOnSet{n}(i,:);
        m = m+1;
    end
    indPeakPyr = pFRNonStimGoodStruct.peakFRInd(indPyr);
    [~,indPeakSorted] = sort(indPeakPyr); 
    pyrProfAll = pyrProfAll(indPeakSorted,:);
    [peakPyrBef] = max(pyrProfAll');
    peakPyrBef(peakPyrBef == 0) = 1;
    pyrProfAll = bsxfun(@rdivide,pyrProfAll,peakPyrBef(:));
    imagesc(timeStepRun/1250,1:length(indPyr),pyrProfAll);
    set(gca,'XLim',[-3 10])
    
    figure(4) % neurons with fields
    indZero = find(timeStepRun >= 0,1);
    indPyr = fieldSpCorrSessNonStimGood.indNeuron;
    pyrProfAll = zeros(length(indPyr),length(timeStepRun));
    m = 1;
    for n = indPyr
        pyrProfAll(m,:) = filteredSpikeArrayRunOnSet{n}(i,:);
        m = m+1;
    end
    indPeakPyr = pFRNonStimGoodStruct.peakFRInd(indPyr);
    [~,indPeakSorted] = sort(indPeakPyr); 
    pyrProfAll = pyrProfAll(indPeakSorted,:);
    [peakPyrBef] = max(pyrProfAll');
    peakPyrBef(peakPyrBef == 0) = 1;
    pyrProfAll = bsxfun(@rdivide,pyrProfAll,peakPyrBef(:));
    imagesc(timeStepRun/1250,1:length(indPyr),pyrProfAll);
    set(gca,'XLim',[-3 10])
    
%     figure(3)
%     indZero = find(timeStepRun >= 0,1);
%     indPyr = find(autoCorr.isPyrneuron == 1);
%     pyrProf = zeros(length(indPyr),length(timeStepRun)-indZero+1);
%     pyrProfAll = zeros(length(indPyr),length(timeStepRun));
%     m = 1;
%     for n = indPyr
%         pyrProf(m,:) = filteredSpikeArrayRunOnSet{n}(i,indZero:end);
%         pyrProfAll(m,:) = filteredSpikeArrayRunOnSet{n}(i,:);
%         m = m+1;
%     end
%     [peakPyr,indPeakPyr] = max(pyrProf');
%     [~,indPeakSorted] = sort(indPeakPyr); 
%     pyrProfAll = pyrProfAll(indPeakSorted,:);
%     [peakPyrBef] = max(pyrProfAll');
%     peakPyrBef(peakPyrBef == 0) = 1;
%     pyrProfAll = bsxfun(@rdivide,pyrProfAll,peakPyrBef(:));
%     imagesc(timeStepRun/1250,1:length(indPyr),pyrProfAll);
%     set(gca,'XLim',[-3 10])
    pause;
end