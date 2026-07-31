function calACG_mpfi(FileNameBase,sampleRate)
% Calculate ACG for each cluster
    
    if exist([FileNameBase '_BehavElectrDataLFP.mat'], 'file') == 2
        fprintf('Check whether %s file already contains Spike and Clu.\n', ...
            [FileNameBase '_BehavElectrDataLFP.mat']);
        load([FileNameBase '_BehavElectrDataLFP.mat'], 'Processing');
        if(sum(Processing == 2) == 0)
            fprintf(['\nError: you need to generate Spike and Clu structures', ...
                ' before calculating autocorrelogram.']);
            return;
        end 
    end
    
    % check whether CCG has already been calculated and saved
    if(sum(Processing == 3) == 0)
        fprintf(...
            'Load Spike from %s file.\n', ...
            [FileNameBase '_BehavElectrDataLFP.mat']);
        load([FileNameBase '_BehavElectrDataLFP.mat'], ...
                'Spike');
        % get ACG of each cluster
        HalfBins = 5000;     % total number of bins on each side
        samplePerBin = 1 * (sampleRate/1000);      
            % number of samples within each bin
        [ccgVal, ccgT] = CCG_wang(Spike.res20kHz, Spike.totclu, samplePerBin, ...
            HalfBins, sampleRate);
        ccgT = ccgT'; 

        Processing = [Processing 3]; % processing stage three, calculating CCG
        fprintf('CCG saved into the structure file: %s....\n',...
                [FileNameBase '_BehavElectrDataLFP_CCG.mat']);
        save([FileNameBase '_BehavElectrDataLFP_CCG.mat'], ...
                'ccgVal','ccgT','-v7.3');
        save([FileNameBase '_BehavElectrDataLFP.mat'], ...
                'Processing', '-append'); 
    end
    
    if(sum(Processing == 4) == 0)
        fprintf(...
            'Load Spike from %s file.\n', ...
            [FileNameBase '_BehavElectrDataLFP_CCG.mat']);
        load([FileNameBase '_BehavElectrDataLFP_CCG.mat'], ...
                'ccgVal','ccgT');
        for IDclu = 1:size(ccgVal,2)
            CCGExt.peaktoMean(IDclu) = ...
                max(ccgVal(ccgT > -30 & ccgT < 30,IDclu,IDclu))...
                ./mean(ccgVal(ccgT < -500 | ccgT > 500,IDclu,IDclu));
            idxTmp = find(ccgT == -50);
            [~,idxPeak] = max(ccgVal(ccgT > -50 & ccgT < 50,IDclu,IDclu));
            CCGExt.peakTime(IDclu) = ccgT(idxPeak + idxTmp);                                 
        end   
    
        Processing = [Processing 4]; % processing stage four, calculating CCG Peak2mean
        fprintf('Update CCG saved into the structure file: %s....\n',...
                [FileNameBase '_BehavElectrDataLFP.mat']);
        save([FileNameBase '_BehavElectrDataLFP_CCG.mat'], ...
                'CCGExt', '-append');
        save([FileNameBase '_BehavElectrDataLFP.mat'], ...
                'Processing', '-append');
    end
end
