cd Z:\Raphael_tests\mice_expdata\ANM016\A016-20190603\A016-20190603-01\

load A016-20190603-01_DataStructure_mazeSection1_TrialType1_alignRun_msess3.mat

load A016-20190603-01_DataStructure_mazeSection1_TrialType1_behPar.mat

indGoodTr = find(behPar.indTrBadBeh == 0 & behPar.indTrNoArt == 1);
for i = indGoodTr
    subplot(4,1,1)
    hold off
    plot([trialsRun.thetaAmpBef{i}',trialsRun.thetaAmp{i}']);
    hold on
    plot([trialsRun.thetaAmpBef{i}']);
    
    subplot(4,1,2)
    hold off
    plot([trialsRun.ThetaFreqBef{i}',trialsRun.ThetaFreq{i}']);
    hold on
    plot([trialsRun.ThetaFreqBef{i}']);
    
    subplot(4,1,3)
    hold off
    plot([trialsRun.eegBef{i}',trialsRun.eeg{i}']);
    hold on
    plot([trialsRun.eegBef{i}']);
    
    subplot(4,1,4)
    hold off
    plot([trialsRun.thetaPhHilbBef{i}',trialsRun.thetaPhHilb{i}']);
    hold on
    plot([trialsRun.thetaPhHilbBef{i}']);
        
    pause;
end
