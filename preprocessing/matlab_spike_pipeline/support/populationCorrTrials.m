function populationCorrTrials(path,fileName,spaceBin,onlyRun)
% calculate the correlation of population activity between trials

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
     elseif nargin == 2
        spaceBin = 2; % cm
        onlyRun = 1;
    elseif nargin == 3
        onlyRun = 1;
    elseif nargin > 4
        disp('Too many input arguments.');
        return;
    end
    
    fullPath = [path fileName '.mat'];
    if(exist(fullPath,'file') == 0)
        disp('recording file does not exist');
    else
        load(fullPath,'cluList');
    end
    
    fileNameConv = [fileName '_convSpikesDist' num2str(spaceBin) ...
                    'mm_Run' num2str(onlyRun) '.mat'];          
    fullPath = [path fileNameConv];
    if(exist(fullPath,'file') == 0)
        if(fileState == 0)
            disp('File does not exist.');
        else
            disp(['The firing profile file does not exist. Try to run the',...
                    'function again with fileState = 0.']);
        end
        return;
    end
    load(fullPath,'filteredSpikeArrayNormTNormAmp');
    
    fileNameInfo = [fileName '_Info.mat'];     
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath,'autoCorr');
    
    ind = autoCorr.isPyrneuron;
%     ind = abs(cluList.centerMax) > 500;
    numTrials = length(filteredSpikeArrayNormTNormAmp);
    activityTrialCorr = zeros(numTrials,numTrials);
    for i = 1:numTrials
        for j = i:numTrials
            activityTrialCorr(i,j) = corr2(filteredSpikeArrayNormTNormAmp{i}(ind,:),...
                        filteredSpikeArrayNormTNormAmp{j}(ind,:));
            activityTrialCorr(j,i) = activityTrialCorr(i,j);
        end
    end
    
    figure;
    imagesc(activityTrialCorr);
    xlabel('Trial No.');
    ylabel('Trial No.');
    
    activityTrialCorr = ...
        bsxfun(@minus,activityTrialCorr,mean(activityTrialCorr));
       
    [coeff,score,latent,tsquared,explained,mu] = ...
        pca(activityTrialCorr);
    
    proj1 = activityTrialCorr*coeff(:,1);
    [~,indProj1] = sort(proj1);
    figure;
    imagesc(activityTrialCorr(indProj1,:))
    
    save(fileNameFull,'activityTrialCorr','-v7.3','-append')
end
