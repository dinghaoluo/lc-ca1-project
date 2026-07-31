function ThetaPower_Hilbert(path,fileName)
%% estimate theta phase and power using hilbert transform

    ind = findstr(fileName, '_');
    fullPath = [path fileName(1:ind(1)) 'BehavElectrDataLFP.mat'];
    if(exist(fullPath) == 0)
        disp('The _BehavElectrDataLFP file does not exist');
        return;
    end    
    load(fullPath,'Track');
    eeg = Track.eeg;
   
    GlobalConst;
    thetaFiltParam.FreqRange = [6 10]; 
        %originally [4 25], changed to [4 16] on 2016.07.14
    thetaFiltParam.FilterOrd = 3;
    thetaFiltParam.Ripple = 20;
    [thetaPower.ThetaPhase, thetaPower.ThetaAmp, thetaPower.TotPhase, ...
        ThetaFreqTmp, thetaPower.Eegf] = ...
        GetThetaPhase_hilbert(eeg, sampleFq, thetaFiltParam);
    thetaPower.ThetaFreq = [ThetaFreqTmp(1);ThetaFreqTmp];
    
    fullPath = [path fileName '_thetaPower.mat'];
    save(fullPath,'thetaPower');
end
