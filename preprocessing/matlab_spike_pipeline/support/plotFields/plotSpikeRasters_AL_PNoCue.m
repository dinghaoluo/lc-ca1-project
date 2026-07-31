function plotSpikeRasters_AL_PNoCue()
    
    %% plot an example neuron in passive no cue task 
    pathPNoCue = 'Z:\Raphael_tests\mice_expdata\ANM007\A007-20190116\A007-20190116-01\';
    fileNamePNoCue = 'A007-20190116-01_DataStructure_mazeSection1_TrialType1';
    plotSpikeRasterEg(pathPNoCue,fileNamePNoCue,1,[1:49],36);
    
    %% plot an example neuron in active licking task
    pathAL = 'Z:\Raphael_tests\mice_expdata\ANM016\A016-20190603\A016-20190603-01\';
    fileNameAL = 'A016-20190603-01_DataStructure_mazeSection1_TrialType1';
    mazeSess = 3;
    fileNameInfo = [fileNameAL '_Info.mat'];
        
    fullPath = [pathAL fileNameInfo];
    if(exist(fullPath,'file') == 0)
        BasicInfo(pathAL,fileNameAL);
    end
    load(fullPath);
    indTr = find(beh.mazeSess == mazeSess);
    
    fullPath = [pathAL fileNameAL '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    indTr = indTr(behPar.indTrBadBeh == 0);
    
%     plotSpikeRasterEg(pathAL,fileNameAL,1,indTr(1:49),47);
    
    plotSpikeRasterEg(pathAL,fileNameAL,1,indTr(end-48:end),68);
end
