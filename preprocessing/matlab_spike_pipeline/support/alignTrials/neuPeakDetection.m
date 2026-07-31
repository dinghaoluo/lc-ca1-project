function isPeakNeu = neuPeakDetection(filteredSpikeArray,paramC,neuNo,p)
    isPeakNeu = zeros(length(p),1);
    [avgShuf,sigShuf] = neuActivityShuffle(filteredSpikeArray,p);
    avgProfile = mean(filteredSpikeArray);
    avgProfile = conv(avgProfile,paramC.gaussFilt,'same');
    
    for i = 1:length(p)
        diffAvgShuf = avgProfile - sigShuf(i,:);
        indDiffAvgShuf = diffAvgShuf > 0;
        [continuousHigh,continuousLow] = numOfConsecutiveOnes(indDiffAvgShuf);
        if(sum(continuousHigh > 0.5/paramC.timeBin) > 0)
            isPeakNeu(i) = 1;
        end
    end
        
    figure(1);
    hold off;
    plot(avgProfile);
    hold on;
    plot(avgShuf);
    for i = 1:length(p)
        plot(sigShuf(i,:));
    end
    title(['isPeakNeu = ' num2str(isPeakNeu') ' neuNo = ' num2str(neuNo)]);
%     pause;
    
end

function [data,data1] = numOfConsecutiveOnes(arr)
    data = [];
    data1 = [];
    s = sprintf('%d', arr);
    %Reading the consequences of 1's from the string by using 0's as delimiters
    t1=textscan(s,'%s','delimiter','0','multipleDelimsAsOne',1);
    % Converting cell array of cell into a single cell array
    d = t1{:};
    % Computing the length of each run by going through the array and assigning
    % it into 
    for k = 1:length(d)
          data(k) = length(d{k});
    end
    
    t2=textscan(s,'%s','delimiter','1','multipleDelimsAsOne',1);
    % Converting cell array of cell into a single cell array
    f = t2{:};
    % Computing the length of each run by going through the array and assigning
    % it into 
    for k = 1:length(f)
          data1(k) = length(f{k});
    end
end