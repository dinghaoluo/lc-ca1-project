function FRProfileMeanStat = accumMeanPyrIntInitPeakStatC(FRProfileMean)
      
    if(isempty(FRProfileMean.meanAvgFRProfileBaseline))
        FRProfileMeanStat = [];
        return;
    end
    FRProfileMeanStat.pRS0to1VsBL = ranksum(FRProfileMean.meanAvgFRProfileBaseline,...
                FRProfileMean.meanAvgFRProfile0to1);    
    FRProfileMeanStat.pRSBefRunVsBL = ranksum(FRProfileMean.meanAvgFRProfileBaseline,...
                FRProfileMean.meanAvgFRProfileBefRun);    
    FRProfileMeanStat.pRS3to5VsBL = ranksum(FRProfileMean.meanAvgFRProfileBaseline,...
                FRProfileMean.meanAvgFRProfile3to5);
    FRProfileMeanStat.pRSBefRunVs0to1 = ranksum(FRProfileMean.meanAvgFRProfile0to1,...
                FRProfileMean.meanAvgFRProfileBefRun);
    FRProfileMeanStat.pRS3to5Vs0to1 = ranksum(FRProfileMean.meanAvgFRProfile0to1,...
                FRProfileMean.meanAvgFRProfile3to5);
    FRProfileMeanStat.pRS3to5VsBefRun = ranksum(FRProfileMean.meanAvgFRProfileBefRun,...
                FRProfileMean.meanAvgFRProfile3to5);

    [~,FRProfileMeanStat.pTTPercChange0to1VsBL] = ttest(FRProfileMean.percChange0to1VsBL);
    [~,FRProfileMeanStat.pTTPercChangeBefRunVsBL] = ttest(FRProfileMean.percChangeBefRunVsBL);
    [~,FRProfileMeanStat.pTTPercChangeBefRunVs0to1] = ttest(FRProfileMean.percChangeBefRunVs0to1);
    [~,FRProfileMeanStat.pTTPercChange0to1Vs3to5] = ttest(FRProfileMean.percChange0to1Vs3to5);
    [~,FRProfileMeanStat.pTTPercChangeBefRunVs3to5] = ttest(FRProfileMean.percChangeBefRunVs3to5);
    
end