function plotExamples()
    
    cd Z:\Raphael_tests\mice_expdata\ANM016\A016-20190603\A016-20190603-01
    
    %% plot example sequence aligned to different behavioral parameters (neurons with field)
    plotSequenceRunOnT('./','A016-20190603-01_DataStructure_mazeSection1_TrialType1',1,3,20,3);
    %% plot example sequence aligned to different behavioral parameters (neurons corrT > threshold)
    plotSequenceRunOnT('./','A016-20190603-01_DataStructure_mazeSection1_TrialType1',1,3,20,1);
    %% plot example sequence in distance
    plotSequenceDist('./','A016-20190603-01_DataStructure_mazeSection1_TrialType1',1,3,20,20,3);
    %% plot example raster as a function of distance
    load('A016-20190603-01_DataStructure_mazeSection1_TrialType1_PeakFRAligned_msess3_Run1.mat',...
        'pFRNonStimGoodStruct','pFRNonStimBadStruct');
    trialNo1 = pFRNonStimGoodStruct.indLapList;
    trialNo1 = trialNo1(trialNo1 ~= 1);
    plotSpikeRasterEg('./','A016-20190603-01_DataStructure_mazeSection1_TrialType1',1,trialNo1,68);
    
    %% plot example raster of PV interneuron
    plotSpikeRaster_AlignedEg('./','A016-20190603-01_DataStructure_mazeSection1_TrialType1',0,3,[],49)
    
    %% plot example raster of PV interneuron
    cd Z:\Raphael_tests\mice_expdata\ANM016\A016-20190531\A016-20190531-01
    plotSpikeRaster_AlignedEg('./','A016-20190531-01_DataStructure_mazeSection1_TrialType1',0,5,[],19)
    
    %% plot example raster of MS neurons
    cd Z:\Raphael_tests\mice_expdata\ANM035\A035-20201125\A035-20201125-03
    plotSpikeRaster_AlignedEg('./','A035-20201125-03_DataStructure_mazeSection1_TrialType1',0,1,[1:50]',1:13)
    
    cd Z:\Raphael_tests\mice_expdata\ANM035\A035-20201127\A035-20201127-02
    plotSpikeRaster_AlignedEg('./','A035-20201127-02_DataStructure_mazeSection1_TrialType1',0,1,[1:50]',1:13)
    plotSpikeRaster_AlignedEg('./','A035-20201127-02_DataStructure_mazeSection1_TrialType1',0,1,[20:70]',12)
    
    %% plot example running profiles over 3 trials
    PlotSpeedAndLickEg
    
    %% plot spike rasters for individual neurons aligned to different behavioral parameters
    load('A016-20190603-01_DataStructure_mazeSection1_TrialType1_FieldSpCorrAligned_Run3_Run1.mat');
    plotSpikeRasterRunOnT_GoodVsBadBeh('./','A016-20190603-01_DataStructure_mazeSection1_TrialType1',1,3,fieldSpCorrSessNonStimGood.indNeuron)

    cd Z:\Raphael_tests\mice_expdata\ANM026\A026-20200323\A026-20200323-01
    %% plot example sequence aligned to different behavioral parameters
    plotSequenceRunOnT('./','A026-20200323-01_DataStructure_mazeSection1_TrialType1',1,1,20,3);
    %% plot example sequence in distance
    plotSequenceDist('./','A026-20200323-01_DataStructure_mazeSection1_TrialType1',1,1,20,20,3);
    %% plot example sequence over multiple trials
    plotSequenceDistOverTrials('./','A026-20200323-01_DataStructure_mazeSection1_TrialType1',1,1,20,20,3,2:4);
    plotSequenceDistOverTrials('./','A026-20200323-01_DataStructure_mazeSection1_TrialType1',1,1,20,20,3,33:35); %selected
    plotSequenceDistOverTrials('./','A026-20200323-01_DataStructure_mazeSection1_TrialType1',1,1,20,20,3,36:38);
    plotSequenceDistOverTrials('./','A026-20200323-01_DataStructure_mazeSection1_TrialType1',1,1,20,20,3,44:46);
    plotSequenceDistOverTrials('./','A026-20200323-01_DataStructure_mazeSection1_TrialType1',1,1,20,20,3,45:47);
    
    %% plot example raster of PV interneuron
    cd Z:\Raphael_tests\mice_expdata\ANM023\A023-20191219\A023-20191219-01
    plotSpikeRaster_AlignedEg('./','A023-20191219-01_DataStructure_mazeSection1_TrialType1',0,1,[],65)
    
    %% plot example raster as a function of distance
    cd Z:\Raphael_tests\mice_expdata\ANM007\A007-20190116\A007-20190116-01
    load('A007-20190116-01_DataStructure_mazeSection1_TrialType1_PeakFRAligned_msess1_Run1.mat',...
        'pFRNonStimGoodStruct','pFRNonStimBadStruct');
    trialNo1 = pFRNonStimGoodStruct.indLapList;
    trialNo1 = trialNo1(trialNo1 ~= 1);
%     plotSpikeRasterEg('./','A007-20190116-01_DataStructure_mazeSection1_TrialType1',1,[1:49],36)
    plotSpikeRasterEg('./','A007-20190116-01_DataStructure_mazeSection1_TrialType1',1,trialNo1,36)