function orderNeuSingleTrial(path,fileName,onlyRun)
% order the neurons based on their activity within a single trial

    timeBin = 0.2; % s
    maxTrialLenLimit = 4.5; %s
    
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameFR = [fileName '_FR_Run' num2str(onlyRun) '.mat'];
    fileNamePeakFR = [fileName(1:end-1) '_RunOnSet_PeakFR' num2str(timeBin*1000) ...
                      'ms.mat']; 
    fileNamePeakFR = [fileName '_PeakFR' num2str(timeBin*1000) ...
                      'ms.mat']; 
    fileNameInfo = [fileName '_Info.mat'];
    fileNameConv = [fileName(1:end-1) '_RunOnSet_convSpikesTime' num2str(timeBin*1000) ...
                    'ms.mat'];
    fileNameOrig = [fileName '.mat'];
    fileNameRunOnsetOrig = [fileName(1:end-1) '_RunOnSet.mat'];
    load([path fileNameFR],'mFRStruct');
    load([path fileNamePeakFR],'pFRStruct','pFRStructGoodTr');
    load([path fileNameInfo],'autoCorr');
    load([path fileNameOrig],'cluList');
    load([path fileNameConv],'filteredSpikeArrayNormAmp');
    load([path fileNameRunOnsetOrig]);
    
    GlobalConst; 
    
    numTrials = length(trials);
    trialLenArr = [];
    for i = 1:numTrials
        trialLenArr = [trialLenArr trials{i}.Nsamples];
    end
    indGoodTrial = find(trialLenArr <= maxTrialLenLimit*sampleFq);
    
    minFR = 0.3;
    maxFR = 4;
    minPeakTo20ms = 3;
    indNeurons = mFRStruct.mFR > minFR & mFRStruct.mFR < maxFR;
    indNeuronsPyr1 = cluList.isIntern == 0 & autoCorr.peakTo20ms > 2 ...
                    & autoCorr.peakTo20ms < 5 & autoCorr.peakTime < 15;
    indNeuronsPyr2 = cluList.isIntern == 0 ...
                    & autoCorr.peakTo20ms >= 5;
    indNeurons1 = find(indNeurons == 1 & indNeuronsPyr1 == 1);
    indNeurons2 = find(indNeurons == 1 & indNeuronsPyr2 == 1);
    indNeurons = find(indNeurons == 1 & (indNeuronsPyr1 == 1 | ...
                    indNeuronsPyr2 == 1));
    
    [~,indNeuronsSorted1] = sort(pFRStructGoodTr.peakFRInd(indNeurons1));
    [~,indNeuronsSorted2] = sort(pFRStructGoodTr.peakFRInd(indNeurons2));
    [~,indNeuronsSorted] = sort(pFRStructGoodTr.peakFRInd(indNeurons));
    
    % within each trial, order neurons according to their peak firing rate 
    indOrder1 = zeros(numTrials,length(indNeurons1));
    for i = 1:numTrials
        peakFRIndPerTrial = pFRStruct.peakFRIndPerTrial(indNeurons1,i);
        indZeros = find(peakFRIndPerTrial == -1);
        peakFRIndPerTrial(indZeros) = Inf;
        [~,indOrder1(i,:)] = sort(peakFRIndPerTrial);
    end
    indOrder2 = zeros(numTrials,length(indNeurons2));
    for i = 1:numTrials
        peakFRIndPerTrial = pFRStruct.peakFRIndPerTrial(indNeurons2,i);
        indZeros = find(peakFRIndPerTrial == -1);
        peakFRIndPerTrial(indZeros) = Inf;
        [~,indOrder2(i,:)] = sort(peakFRIndPerTrial);
    end
    
    % calculate rank correlation on each pair
    peakFRIndPerTrial = pFRStruct.peakFRIndPerTrial(indNeurons,:); 
    rankCorrPeakPerTr = zeros(numTrials,numTrials);
    pRankCorrPeakPerTr = zeros(numTrials,numTrials);
    sigPairs = [];
    for i = 1:numTrials
        for j = 1:numTrials
            [rankCorrPeakPerTr(i,j),pRankCorrPeakPerTr(i,j)]...
                = corr(peakFRIndPerTrial(:,i),...
                       peakFRIndPerTrial(:,j),'Type','spearman');
            if(rankCorrPeakPerTr(i,j) > 0.4 & ...
                pRankCorrPeakPerTr(i,j) < 0.05 & i~=j)
                sigPairs = [sigPairs; i,j];
            end
        end
    end
    
    for i = 1:length(indNeurons)
        ind = find(sigPairs(:,1) == i)
        
    end
    
    count = 0;
    for i = 1:length(indGoodTrial)
        count = count + 1;

        if(mod(count-1,16) == 0)
            [figNew,pos] = CreateFig();
            set(0,'Units','pixels') 
            figTitle = 'FR per trial';
%             figTitle = 'Spikes vs Dist';
            set(figure(figNew),'OuterPosition',...
                [pos(1) pos(2) pos(3)*2 pos(4)*2.2],'Name',figTitle)
        end
        subplot(4,4,mod(count-1,16)+1)  
        imagesc(filteredSpikeArrayNormAmp{indGoodTrial(i)}...
            (indNeurons(indNeuronsSorted),:));
        if(mod(count-1,4) == 0)
            ylabel('Neuron');
        end
        if(mod(count-1,16) > 11)
            xlabel('Time (s)');
        end
        figTitle = ['Trial ' num2str(i)];
        title(figTitle);
    end
    
% function LCS(seqA,seqB,m,n) 
% 
%     L = zeros(m+1,n+1);
% 
%     % Following steps build L[m+1][n+1] in bottom up fashion. Note 
%     % that L[i][j] contains length of LCS of X[0..i-1] and Y[0..j-1]
%     % int L[m+1][n+1]; 
%     for i = 1:m+1
%         for j = 1:n+1
%             if(i == 1 || j == 1)
%                 L(i,j) = 0;
%             elseif(seqA(i-1) == seqB(j-1))
%                 L(i,j) = L(i-1,j-1)+1;
%             else
%                 L(i,j) = max(L(i-1,j),L(i,j-1));
%             end
%         end
%     end
% 
%     index = 
% 
% /* Following steps build L[m+1][n+1] in bottom up fashion. Note 
%   that L[i][j] contains length of LCS of X[0..i-1] and Y[0..j-1] */
% for (int i=0; i<=m; i++) 
% { 
%  for (int j=0; j<=n; j++) 
%  { 
%    if (i == 0 || j == 0) 
%      L[i][j] = 0; 
%    else if (X[i-1] == Y[j-1]) 
%      L[i][j] = L[i-1][j-1] + 1; 
%    else
%      L[i][j] = max(L[i-1][j], L[i][j-1]); 
%  } 
% } 
% 
% // Following code is used to print LCS 
% int index = L[m][n]; 
% 
% // Create a character array to store the lcs string 
% char lcs[index+1]; 
% lcs[index] = ''; // Set the terminating character 
% 
% // Start from the right-most-bottom-most corner and 
% // one by one store characters in lcs[] 
% int i = m, j = n; 
% while (i > 0 && j > 0) 
% { 
%   // If current character in X[] and Y are same, then 
%   // current character is part of LCS 
%   if (X[i-1] == Y[j-1]) 
%   { 
%       lcs[index-1] = X[i-1]; // Put current character in result 
%       i--; j--; index--;     // reduce values of i, j and index 
%   } 
% 
%   // If not same, then find the larger of two and 
%   // go in the direction of larger value 
%   else if (L[i-1][j] > L[i][j-1]) 
%      i--; 
%   else
%      j--; 
% } 
% 
% // Print the lcs 
% cout << "LCS of " << X << " and " << Y << " is " << lcs; 
% } 
% end
