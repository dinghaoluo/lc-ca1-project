function CCGPeakTrough(path,fileName)
% finding out whether there is significant cofiring or mutual inhibition
% through CCG. Here we used the method reported in Diba.JNeurosci.2014, and
% Stark.jNeurosciMethods.2009

    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '_');
    if(~isempty(indexFileName))
        fileNameBehGenCCG = [fileName(1:indexFileName(1)) ...
            'BehavElectrDataLFP_CCG.mat'];
        fileNameBehGen = [fileName(1:indexFileName(1)) ...
            'BehavElectrDataLFP.mat'];
    else
        fileNameBehGenCCG = [];
        
    end 
    fullPath = [path fileNameBehGenCCG];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath); 
    
    fullPath = [path fileNameBehGen];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath,'Spike'); 
    
    GlobalConst;
    
    numNeurons = size(ccgVal,2);
    lenT = size(ccgVal,1);
    
    peakTroughThetaWin = 50;
    binSize = sampleFqOri*(ccgT(2)-ccgT(1))/1000; 
        % number of samples in a bin
    peakTroughThetaInd = find(ccgT >= -peakTroughThetaWin & ...
                ccgT <= peakTroughThetaWin);
    lenInd = length(peakTroughThetaInd);
    
    wind = [ 1 1 1 0.42 1 1 1; ...
            -3 -2 -1 0 1 2 3];
    wind(1,:) = wind(1,:)/sum(wind(1,:));
    lenWind = size(wind,2);
    ccgPredict = zeros(lenInd,numNeurons,numNeurons);
    pCCG = zeros(lenInd,numNeurons,numNeurons);
    ccgPredictSig99 = zeros(lenInd,numNeurons,numNeurons);
    ccgPredictSig95 = zeros(lenInd,numNeurons,numNeurons);
    ccgPredictSig05 = zeros(lenInd,numNeurons,numNeurons);
    ccgPredictSig01 = zeros(lenInd,numNeurons,numNeurons);
    ccgPredictSigInd99 = zeros(lenInd,numNeurons,numNeurons);
    ccgPredictSigInd95 = zeros(lenInd,numNeurons,numNeurons);
    ccgPredictSigInd05 = zeros(lenInd,numNeurons,numNeurons);
    ccgPredictSigInd01 = zeros(lenInd,numNeurons,numNeurons);
    sigPairs99 = [];
    sigPairs95 = [];
    sigPairs05 = [];
    sigPairs01 = [];
    sigPairs99Interval = [];
    sigPairs95Interval = [];
    sigPairs05Interval = [];
    sigPairs01Interval = [];
    nPairs99 = 1;
    nPairs95 = 1;
    nPairs05 = 1;
    nPairs01 = 1;
    for i = 1:numNeurons
        numSpikes = sum(Spike.totclu == i);  
        for j = i+1:numNeurons
            % remove normalization
            ccgTmp = ccgVal(:,i,j) ...
                    / sampleFqOri * binSize * numSpikes;
            if(mean(ccgTmp(peakTroughThetaInd)) > 2.5)
                % calculate the ccg prediction
                ccgPredictTmp = conv(ccgTmp,wind(1,:));
                ccgPredictTmp = ccgPredictTmp(round(lenWind/2):...
                                        end-round(lenWind/2)+1);
                % remove the normalization
                ccgPredict(:,i,j) = ccgPredictTmp(peakTroughThetaInd);
                % propability of occurance based on poisson distribution
                pCCG(:,i,j) = 1-poisscdf(ccgTmp(peakTroughThetaInd),...
                        ccgPredict(:,i,j));
                    
                % 99% upper bound confidence level of the ccg prediction
                ccgPredictSig99(:,i,j) = ...
                    poissinv(0.99*ones(lenInd,1),ccgPredict(:,i,j));
                % find indices that are significantly different from 
                % prediction
                ccgPredictSigInd99(:,i,j) = ccgTmp(peakTroughThetaInd)...
                    > ccgPredictSig99(:,i,j);
                % rule out that the significant difference occurs purely
                % by chance
                [sigPair,sigPairsInterval] = testSignificance(...
                    ccgPredictSigInd99(:,i,j)',i,j);
                if(~isempty(sigPair))
                    sigPairs99 = [sigPairs99;sigPair];
                    sigPairs99Interval{nPairs99} = sigPairsInterval;
                    nPairs99 = nPairs99+1;
                end
                
                % 95% upper bound confidence level of the ccg prediction
                ccgPredictSig95(:,i,j) = ...
                    poissinv(0.95*ones(lenInd,1),ccgPredict(:,i,j));
                ccgPredictSigInd95(:,i,j) = ccgTmp(peakTroughThetaInd)...
                    > ccgPredictSig95(:,i,j);
                [sigPair,sigPairsInterval] = testSignificance(...
                    ccgPredictSigInd95(:,i,j)',i,j);
                if(~isempty(sigPair))
                    sigPairs95 = [sigPairs95;sigPair];
                    sigPairs95Interval{nPairs95} = sigPairsInterval;
                    nPairs95 = nPairs95+1;
                end
                
                % 5% lower bound confidence level of the ccg prediction
                ccgPredictSig05(:,i,j) = ...
                    poissinv(0.05*ones(lenInd,1),ccgPredict(:,i,j));
                ccgPredictSigInd05(:,i,j) = ccgTmp(peakTroughThetaInd)...
                    < ccgPredictSig05(:,i,j);
                [sigPair,sigPairsInterval] = testSignificance(...
                    ccgPredictSigInd05(:,i,j)',i,j);
                if(~isempty(sigPair))
                    sigPairs05 = [sigPairs05;sigPair];
                    sigPairs05Interval{nPairs05} = sigPairsInterval;
                    nPairs05 = nPairs05+1;
                end
                
                % 1% lower bound confidence level of the ccg prediction
                ccgPredictSig01(:,i,j) = ...
                    poissinv(0.01*ones(lenInd,1),ccgPredict(:,i,j)); 
                ccgPredictSigInd01(:,i,j) = ccgTmp(peakTroughThetaInd)...
                    < ccgPredictSig01(:,i,j);
                [sigPair,sigPairsInterval] = testSignificance(...
                    ccgPredictSigInd01(:,i,j)',i,j);
                if(~isempty(sigPair))
                    sigPairs01 = [sigPairs01;sigPair];
                    sigPairs01Interval{nPairs01} = sigPairsInterval;
                    nPairs01 = nPairs01+1;
                end          
            end            
        end
    end
    
    peakTrough.ccgPredict = ccgPredict;
    peakTrough.pCCG = pCCG;
    peakTrough.ccgPredictSig99 = ccgPredictSig99;
    peakTrough.ccgPredictSig95 = ccgPredictSig95;
    peakTrough.ccgPredictSig05 = ccgPredictSig05;
    peakTrough.ccgPredictSig01 = ccgPredictSig01;
    peakTrough.ccgPredictSigInd99 = ccgPredictSigInd99;
    peakTrough.ccgPredictSigInd95 = ccgPredictSigInd95;
    peakTrough.ccgPredictSigInd05 = ccgPredictSigInd05;
    peakTrough.ccgPredictSigInd01 = ccgPredictSigInd01;
    peakTrough.sigPairs99 = sigPairs99;
    peakTrough.sigPairs95 = sigPairs95;
    peakTrough.sigPairs05 = sigPairs05;
    peakTrough.sigPairs01 = sigPairs01;
    peakTrough.sigPairs99Interval = sigPairs99Interval;
    peakTrough.sigPairs95Interval = sigPairs95Interval;
    peakTrough.sigPairs05Interval = sigPairs05Interval;
    peakTrough.sigPairs01Interval = sigPairs01Interval;
    
    indexFileName = findstr(fileName, '.mat');
    if(isempty(indexFileName))
        fileNamePeakTrough = [fileName '_CCGPeakTrough' ...
            num2str(peakTroughThetaWin) 'ms.mat'];
    else
        fileNamePeakTrough = [fileName(1:indexFileName(end)-1) '_CCGPeakTrough' ...
            num2str(peakTroughThetaWin) 'ms.mat'];
    end 
    fullPath = [path fileName];
    save(fullPath, 'peakTrough');
end

function [sigPair,sigInterval] = testSignificance(sigInd,i,j)
    numRepeats = countOnes(sigInd);
    indOne = find(numRepeats(2,:) == 1);
    if(length(indOne) > 0 && max(numRepeats(1,indOne)) >= 2) 
        sigPairs = [i,j];
        repeatIndex = [0 cumsum(numRepeats(1,:))];
        sigInterval = [];
        for n = 1:length(indOne)
            if(numRepeats(1,indOne(n)) >= 2)
                siginterval = [...
                    siginterval;...
                    repeatIndex(indOne(n))+1 ...
                    repeatIndex(indOne(n)+1)];
            end
        end
    end
end

function y = countOnes(x)
    i = find(diff(x)); 
    n = [i numel(x)] - [0 i];
    c = arrayfun(@(X) X-1:-1:0, n , 'un',0);
    y = zeros(2,length(c));
    cumCount = 0;
    for i = 1:length(c)
        cumCount = cumCount + length(c{i});
        y(1,i) = length(c{i});
        y(2,i) = x(cumCount);
    end
%     y = cat(2,c{:});
end