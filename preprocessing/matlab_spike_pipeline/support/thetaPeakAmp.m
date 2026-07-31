function [thetaPeakAmpMeanArr, thetaPeakAmpStdArr, ...
    thetaPeakAmpMean, thetaPeakAmpStd] = ...
    thetaPeakAmp(thetaPeakAmpArr,indSpeed,indLap)
% this function is to calculate the variance of theta peak amplitude
% thetaPeakAmpStdArr:      theta peak amplitude from all the trials
% indLap:                  indices of laps

numEegTrace = 0;
if(~isempty(thetaPeakAmpArr{indLap(1)}))
    if(iscell(thetaPeakAmpArr{indLap(1)}))
        numEegTrace = size(thetaPeakAmpArr{indLap(1)},1);
    else
        numEegTrace = 1;
    end
end

thetaPeakAmpStdArr = zeros(length(indLap),numEegTrace);
thetaPeakAmpMeanArr = zeros(length(indLap),numEegTrace);
for j = 1:numEegTrace
    thetaPeakArr = [];    
    for i = 1:length(indLap)   
        if(~isempty(thetaPeakAmpArr{indLap(i)})) 
            if iscell(thetaPeakAmpArr{indLap(i)})
                [~,indTheta] = intersect(thetaPeakAmpArr{indLap(i)}{j}(:,1),...
                                        indSpeed{i});
                thetaPeakAmpMeanArr(i,j) = ...
                    mean(thetaPeakAmpArr{indLap(i)}{j}(indTheta,2));        
                    % first column: time; second column: amplitude
                thetaPeakAmpStdArr(i,j) = ...
                    std(thetaPeakAmpArr{indLap(i)}{j}(indTheta,2));
                thetaPeakArr = ...
                    [thetaPeakArr; thetaPeakAmpArr{indLap(i)}{j}(indTheta,2)];
            else
                [~,indTheta] = intersect(thetaPeakAmpArr{indLap(i)}(:,1),...
                                        indSpeed{i});
                thetaPeakAmpMeanArr(i,j) = ...
                    mean(thetaPeakAmpArr{indLap(i)}(indTheta,2));
                thetaPeakAmpStdArr(i,j) = ...
                    std(thetaPeakAmpArr{indLap(i)}(indTheta,2));
                thetaPeakArr = ...
                    [thetaPeakArr; thetaPeakAmpArr{indLap(i)}(indTheta,2)];
            end
        end
    end
    thetaPeakAmpMean(j) = mean(thetaPeakArr);
    thetaPeakAmpStd(j) = std(thetaPeakArr);
end
