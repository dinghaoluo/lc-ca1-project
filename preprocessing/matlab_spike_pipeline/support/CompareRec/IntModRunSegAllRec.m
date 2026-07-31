function IntModRunSegAllRec(onlyRun)

    methodTheta = 1;
    minFRInt = 2;
    
    RecordingList;
    pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\interneuron\';
    
    load([pathAnal 'autoCorrIntAllRec.mat']);
    
    pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\interneuronRunSeg\';
    
    load([pathAnal 'modRunSegIntAllRec.mat']);
    
%     %% interneurons in no cue passive task
%     disp('No cue')
%     modSegIntNoCue = accumIntNeurons1(listRecordingsNoCuePath,...
%         listRecordingsNoCueFileName,mazeSessionNoCue,minFRInt,1,onlyRun);
%     
%     disp('Active licking')
%     modSegIntAL = accumIntNeurons1(listRecordingsActiveLickPath,...
%         listRecordingsActiveLickFileName,mazeSessionActiveLick,minFRInt,2,onlyRun);
%     
%     disp('Passive licking')
%     modSegIntPL = accumIntNeurons1(listRecordingsPassiveLickPath,...
%         listRecordingsPassiveLickFileName,mazeSessionPassiveLick,minFRInt,3,onlyRun);
%     
%     save([pathAnal 'modRunSegIntAllRec.mat'],'modSegIntNoCue','modSegIntAL','modSegIntPL','-append'); 
    
    mod.task = [autoCorrIntNoCue.task autoCorrIntAL.task autoCorrIntPL.task];
    mod.indRec = [autoCorrIntNoCue.indRec autoCorrIntAL.indRec autoCorrIntPL.indRec];
    mod.nNeuWithField = [modIntNoCue.nNeuWithField modIntAL.nNeuWithField modIntPL.nNeuWithField];
    
    mod.thetaModHistBefRun = [modSegIntNoCue.thetaModHistBefRun modSegIntAL.thetaModHistBefRun modSegIntPL.thetaModHistBefRun];
    mod.phaseMeanDireBefRun = [modSegIntNoCue.phaseMeanDireBefRun modSegIntAL.phaseMeanDireBefRun modSegIntPL.phaseMeanDireBefRun];
    mod.maxPhaseArrBefRun = [modSegIntNoCue.maxPhaseFilBefRun modSegIntAL.maxPhaseFilBefRun modSegIntPL.maxPhaseFilBefRun];
    mod.minPhaseArrBefRun = [modSegIntNoCue.minPhaseFilBefRun modSegIntAL.minPhaseFilBefRun modSegIntPL.minPhaseFilBefRun];
    mod.phaseMeanResultantLenBefRun = [modSegIntNoCue.phaseMeanResultantLenBefRun modSegIntAL.phaseMeanResultantLenBefRun modSegIntPL.phaseMeanResultantLenBefRun];
    
    mod.thetaModHist0to1 = [modSegIntNoCue.thetaModHist0to1 modSegIntAL.thetaModHist0to1 modSegIntPL.thetaModHist0to1];
    mod.thetaModHistH0to1 = [modSegIntNoCue.thetaModHistH0to1 modSegIntAL.thetaModHistH0to1 modSegIntPL.thetaModHistH0to1];
    mod.phaseMeanDire0to1 = [modSegIntNoCue.phaseMeanDire0to1 modSegIntAL.phaseMeanDire0to1 modSegIntPL.phaseMeanDire0to1];
    mod.phaseMeanDireH0to1 = [modSegIntNoCue.phaseMeanDireH0to1 modSegIntAL.phaseMeanDireH0to1 modSegIntPL.phaseMeanDireH0to1];
    mod.maxPhaseArr0to1 = [modSegIntNoCue.maxPhaseFil0to1 modSegIntAL.maxPhaseFil0to1 modSegIntPL.maxPhaseFil0to1];
    mod.maxPhaseArrH0to1 = [modSegIntNoCue.maxPhaseFilH0to1 modSegIntAL.maxPhaseFilH0to1 modSegIntPL.maxPhaseFilH0to1];
    mod.minPhaseArr0to1 = [modSegIntNoCue.minPhaseFil0to1 modSegIntAL.minPhaseFil0to1 modSegIntPL.minPhaseFil0to1];
    mod.minPhaseArrH0to1 = [modSegIntNoCue.minPhaseFilH0to1 modSegIntAL.minPhaseFilH0to1 modSegIntPL.minPhaseFilH0to1];
    mod.phaseMeanResultantLen0to1 = [modSegIntNoCue.phaseMeanResultantLen0to1 modSegIntAL.phaseMeanResultantLen0to1 modSegIntPL.phaseMeanResultantLen0to1];
    
    mod.thetaModHist3to5 = [modSegIntNoCue.thetaModHist3to5 modSegIntAL.thetaModHist3to5 modSegIntPL.thetaModHist3to5];
    mod.thetaModHistH3to5 = [modSegIntNoCue.thetaModHistH3to5 modSegIntAL.thetaModHistH3to5 modSegIntPL.thetaModHistH3to5];
    mod.phaseMeanDire3to5 = [modSegIntNoCue.phaseMeanDire3to5 modSegIntAL.phaseMeanDire3to5 modSegIntPL.phaseMeanDire3to5];
    mod.phaseMeanDireH3to5 = [modSegIntNoCue.phaseMeanDireH3to5 modSegIntAL.phaseMeanDireH3to5 modSegIntPL.phaseMeanDireH3to5];
    mod.maxPhaseArr3to5 = [modSegIntNoCue.maxPhaseFil3to5 modSegIntAL.maxPhaseFil3to5 modSegIntPL.maxPhaseFil3to5];
    mod.maxPhaseArrH3to5 = [modSegIntNoCue.maxPhaseFilH3to5 modSegIntAL.maxPhaseFilH3to5 modSegIntPL.maxPhaseFilH3to5];
    mod.minPhaseArr3to5 = [modSegIntNoCue.minPhaseFil3to5 modSegIntAL.minPhaseFil3to5 modSegIntPL.minPhaseFil3to5];
    mod.minPhaseArrH3to5 = [modSegIntNoCue.minPhaseFilH3to5 modSegIntAL.minPhaseFilH3to5 modSegIntPL.minPhaseFilH3to5];
    mod.phaseMeanResultantLen3to5 = [modSegIntNoCue.phaseMeanResultantLen3to5 modSegIntAL.phaseMeanResultantLen3to5 modSegIntPL.phaseMeanResultantLen3to5];
    
    mod.thetaModHist3to4 = [modSegIntNoCue.thetaModHist3to4 modSegIntAL.thetaModHist3to4 modSegIntPL.thetaModHist3to4];
    mod.thetaModHistH3to4 = [modSegIntNoCue.thetaModHistH3to4 modSegIntAL.thetaModHistH3to4 modSegIntPL.thetaModHistH3to4];
    mod.phaseMeanDire3to4 = [modSegIntNoCue.phaseMeanDire3to4 modSegIntAL.phaseMeanDire3to4 modSegIntPL.phaseMeanDire3to4];
    mod.phaseMeanDireH3to4 = [modSegIntNoCue.phaseMeanDireH3to4 modSegIntAL.phaseMeanDireH3to4 modSegIntPL.phaseMeanDireH3to4];
    mod.maxPhaseArr3to4 = [modSegIntNoCue.maxPhaseFil3to4 modSegIntAL.maxPhaseFil3to4 modSegIntPL.maxPhaseFil3to4];
    mod.maxPhaseArrH3to4 = [modSegIntNoCue.maxPhaseFilH3to4 modSegIntAL.maxPhaseFilH3to4 modSegIntPL.maxPhaseFilH3to4];
    mod.minPhaseArr3to4 = [modSegIntNoCue.minPhaseFil3to4 modSegIntAL.minPhaseFil3to4 modSegIntPL.minPhaseFil3to4];
    mod.minPhaseArrH3to4 = [modSegIntNoCue.minPhaseFilH3to4 modSegIntAL.minPhaseFilH3to4 modSegIntPL.minPhaseFilH3to4];
    mod.phaseMeanResultantLen3to4 = [modSegIntNoCue.phaseMeanResultantLen3to4 modSegIntAL.phaseMeanResultantLen3to4 modSegIntPL.phaseMeanResultantLen3to4];
    
%     idxC = [autoCorrIntAll.idxC(autoCorrIntAll.task == autoCorrIntNoCue.task(1))' ...
%         autoCorrIntAll.idxC(autoCorrIntAll.task == autoCorrIntAL.task(1))' ...
%         autoCorrIntAll.idxC(autoCorrIntAll.task == autoCorrIntPL.task(1))'];
    mod.idxC = [autoCorrIntAll.idxC1(autoCorrIntAll.task == autoCorrIntNoCue.task(1))' ...
        autoCorrIntAll.idxC1(autoCorrIntAll.task == autoCorrIntAL.task(1))' ...
        autoCorrIntAll.idxC1(autoCorrIntAll.task == autoCorrIntPL.task(1))'];
%     idxC = [autoCorrIntAll.idxC2(autoCorrIntAll.task == autoCorrIntNoCue.task(1))' ...
%         autoCorrIntAll.idxC2(autoCorrIntAll.task == autoCorrIntAL.task(1))' ...
%         autoCorrIntAll.idxC2(autoCorrIntAll.task == autoCorrIntPL.task(1))'];

    % for each cluster, compare neurons in the recordings with fields vs.
    % recordings without field
    modSegIntStatsField = modIntStatsFieldPerC(mod);
    
    % compare neurons between the two clusters
    modSegIntStatsSeg = modIntStatsSegments(mod);
    
    save([pathAnal 'modRunSegIntAllRec.mat'],'modSegIntStatsField','modSegIntStatsSeg',...
        '-append'); 
        
    %% compare two clusters
    colorSel = 1;
    cluster = 3;
    idxC = mod.idxC;
    indCurCField = idxC == cluster & mod.nNeuWithField >= 2;
    indCurCNoField = idxC == cluster & mod.nNeuWithField < 1;
    plotModRunSeg(mod,modSegIntStatsField,modSegIntStatsSeg,colorSel,...
        cluster,indCurCField,indCurCNoField,pathAnal);
    
    %% 
    colorSel = 6;
    cluster = 5;
    indCurCField = idxC == cluster & mod.nNeuWithField >= 2;
    indCurCNoField = idxC == cluster & mod.nNeuWithField < 1;
    plotModRunSeg(mod,modSegIntStatsField,modSegIntStatsSeg,colorSel,...
        cluster,indCurCField,indCurCNoField,pathAnal);
    
    %% 
    colorSel = 0;
    cluster = 2;
    indCurCField = idxC == cluster & mod.nNeuWithField >= 2;
    indCurCNoField = idxC == cluster & mod.nNeuWithField < 1;
    plotModRunSeg(mod,modSegIntStatsField,modSegIntStatsSeg,colorSel,...
        cluster,indCurCField,indCurCNoField,pathAnal);
end

function plotModRunSeg(mod,modSegIntStatsField,modSegIntStatsSeg,colorSel,...
        cluster,indCurCField,indCurCNoField,pathAnal)
    idxC = mod.idxC;
    % before run, field vs no field 
    plotBoxPlot(mod.phaseMeanDireBefRun(indCurCField),...
        mod.phaseMeanDireBefRun(indCurCNoField),'Mean theta phase (bef run)',...
        ['Int_ThetaMeanC' num2str(cluster) 'BefRunField-NoFieldBox'],...
        pathAnal,[0 2*pi],modSegIntStatsField.pWWPhaseMeanDireBefRun(cluster),colorSel);
    
    plotBoxPlot(mod.phaseMeanResultantLenBefRun(indCurCField),...
        mod.phaseMeanResultantLenBefRun(indCurCNoField),'Theta phase resultant length (bef run)',...
        ['Int_PhaseMeanResultantLen' num2str(cluster) 'BefRunField-NoFieldBox'],...
        pathAnal,[-0.1 0.8],modSegIntStatsField.pRSPhaseMeanResultantLenBefRun(cluster),colorSel);
    
    plotBoxPlot(mod.thetaModHistBefRun(indCurCField),...
        mod.thetaModHistBefRun(indCurCNoField),'Theta modulation (bef run)',...
        ['Int_thetaModHist' num2str(cluster) 'BefRunField-NoFieldBox'],...
        pathAnal,[-0.1 1.2],modSegIntStatsField.pRSThetaModHistBefRun(cluster),colorSel);
    
    % 0-1s, field vs no field 
    plotBoxPlot(mod.phaseMeanDire0to1(indCurCField),...
        mod.phaseMeanDire0to1(indCurCNoField),'Mean theta phase (0to1s)',...
        ['Int_ThetaMeanC' num2str(cluster) '0to1sField-NoFieldBox'],...
        pathAnal,[0 2*pi],modSegIntStatsField.pWWPhaseMeanDire0to1(cluster),colorSel);
    
    plotBoxPlot(mod.phaseMeanResultantLen0to1(indCurCField),...
        mod.phaseMeanResultantLen0to1(indCurCNoField),'Theta phase resultant length (0to1s)',...
        ['Int_PhaseMeanResultantLen' num2str(cluster) '0to1sField-NoFieldBox'],...
        pathAnal,[-0.1 0.6],modSegIntStatsField.pRSPhaseMeanResultantLen0to1(cluster),colorSel);
    
    plotBoxPlot(mod.thetaModHist0to1(indCurCField),...
        mod.thetaModHist0to1(indCurCNoField),'Theta modulation (0to1s)',...
        ['Int_thetaModHist' num2str(cluster) '0to1sField-NoFieldBox'],...
        pathAnal,[-0.1 1.2],modSegIntStatsField.pRSThetaModHist0to1(cluster),colorSel);
    
    % 3-4s, field vs no field 
    plotBoxPlot(mod.phaseMeanDire3to4(indCurCField),...
        mod.phaseMeanDire3to4(indCurCNoField),'Mean theta phase (3to4s)',...
        ['Int_ThetaMeanC' num2str(cluster) '3to4sField-NoFieldBox'],...
        pathAnal,[0 2*pi],modSegIntStatsField.pWWPhaseMeanDire3to4(cluster),colorSel);
    
    plotBoxPlot(mod.phaseMeanResultantLen3to4(indCurCField),...
        mod.phaseMeanResultantLen3to4(indCurCNoField),'Theta phase resultant length (3to4s)',...
        ['Int_PhaseMeanResultantLen' num2str(cluster) '3to4sField-NoFieldBox'],...
        pathAnal,[-0.1 0.6],modSegIntStatsField.pRSPhaseMeanResultantLen3to4(cluster),colorSel);
    
    plotBoxPlot(mod.thetaModHist3to4(indCurCField),...
        mod.thetaModHist3to4(indCurCNoField),'Theta modulation (3to4s)',...
        ['Int_thetaModHist' num2str(cluster) '3to4sField-NoFieldBox'],...
        pathAnal,[-0.1 1.2],modSegIntStatsField.pRSThetaModHist3to4(cluster),colorSel);
   
    % 0-1s vs 3-4s 
    plotBoxPlot(mod.phaseMeanDire0to1(idxC == cluster),...
        mod.phaseMeanDire3to4(idxC == cluster),'Mean theta phase (0to1 vs 3to4)',...
        ['Int_ThetaMeanC' num2str(cluster) '0to1sVs3to4sBox'],...
        pathAnal,[0 2*pi],modSegIntStatsSeg.pWWPhaseMeanDireC0to1vs3to4(cluster),colorSel);
    
    plotBoxPlot(mod.phaseMeanResultantLen0to1(idxC == cluster),...
        mod.phaseMeanResultantLen3to4(idxC == cluster),'Theta phase resultant length (0to1 vs 3to4)',...
        ['Int_PhaseMeanResultantLen' num2str(cluster) '0to1sVs3to4sBox'],...
        pathAnal,[-0.1 0.6],modSegIntStatsSeg.pRSPhaseMeanResultantLenC0to1vs3to4(cluster),colorSel);
    
    plotBoxPlot(mod.thetaModHist0to1(idxC == cluster),...
        mod.thetaModHist3to4(idxC == cluster),'Theta modulation (0to1 vs 3to4)',...
        ['Int_thetaModHist' num2str(cluster) '0to1sVs3to4sBox'],...
        pathAnal,[-0.1 1.2],modSegIntStatsSeg.pRSThetaModHistC0to1vs3to4(cluster),colorSel);
    
    % 0-1s vs 3-4s no field
    plotBoxPlot(mod.phaseMeanDire0to1(indCurCNoField),...
        mod.phaseMeanDire3to4(indCurCNoField),'Mean theta phase (0to1 vs 3to4 no field)',...
        ['Int_ThetaMeanC' num2str(cluster) '0to1sVs3to4sNoFieldBox'],...
        pathAnal,[0 2*pi],modSegIntStatsSeg.pWWPhaseMeanDireCNoField0to1vs3to4(cluster),colorSel);
    
    plotBoxPlot(mod.phaseMeanResultantLen0to1(indCurCNoField),...
        mod.phaseMeanResultantLen3to4(indCurCNoField),'Theta phase resultant length (0to1 vs 3to4 no field)',...
        ['Int_PhaseMeanResultantLen' num2str(cluster) '0to1sVs3to4sNoFieldBox'],...
        pathAnal,[-0.1 0.8],modSegIntStatsSeg.pRSPhaseMeanResultantLenCNoField0to1vs3to4(cluster),colorSel);
    
    plotBoxPlot(mod.thetaModHist0to1(indCurCNoField),...
        mod.thetaModHist3to4(indCurCNoField),'Theta modulation (0to1 vs 3to4 no field)',...
        ['Int_thetaModHist' num2str(cluster) '0to1sVs3to4sNoFieldBox'],...
        pathAnal,[-0.1 1.2],modSegIntStatsSeg.pRSThetaModHistCNoField0to1vs3to4(cluster),colorSel);
   
    % 0-1s vs 3-4s field    
    plotBoxPlot(mod.phaseMeanDire0to1(indCurCField),...
        mod.phaseMeanDire3to4(indCurCField),'Mean theta phase (0to1 vs 3to4 field)',...
        ['Int_ThetaMeanC' num2str(cluster) '0to1sVs3to4sFieldBox'],...
        pathAnal,[0 2*pi],modSegIntStatsSeg.pWWPhaseMeanDireCField0to1vs3to4(cluster),colorSel);
    
    plotBoxPlot(mod.phaseMeanResultantLen0to1(indCurCField),...
        mod.phaseMeanResultantLen3to4(indCurCField),'Theta phase resultant length (0to1 vs 3to4 field)',...
        ['Int_PhaseMeanResultantLen' num2str(cluster) '0to1sVs3to4sFieldBox'],...
        pathAnal,[-0.1 0.8],modSegIntStatsSeg.pRSPhaseMeanResultantLenCField0to1vs3to4(cluster),colorSel);
    
    plotBoxPlot(mod.thetaModHist0to1(indCurCField),...
        mod.thetaModHist3to4(indCurCField),'Theta modulation (0to1 vs 3to4 field)',...
        ['Int_thetaModHist' num2str(cluster) '0to1sVs3to4sFieldBox'],...
        pathAnal,[-0.1 1.2],modSegIntStatsSeg.pRSThetaModHistCField0to1vs3to4(cluster),colorSel);
    
    % 0-1s vs 3-5s field
    plotBoxPlot(mod.phaseMeanDire0to1(indCurCField),...
        mod.phaseMeanDire3to5(indCurCField),'Mean theta phase (0to1 vs 3to5 field)',...
        ['Int_ThetaMeanC' num2str(cluster) '0to1sVs3to5sFieldBox'],...
        pathAnal,[0 2*pi],modSegIntStatsSeg.pWWPhaseMeanDireCField0to1vs3to5(cluster),colorSel);
    
    plotBoxPlot(mod.phaseMeanResultantLen0to1(indCurCField),...
        mod.phaseMeanResultantLen3to5(indCurCField),'Theta phase resultant length (0to1 vs 3to5 field)',...
        ['Int_PhaseMeanResultantLen' num2str(cluster) '0to1sVs3to5sFieldBox'],...
        pathAnal,[-0.1 0.8],modSegIntStatsSeg.pRSPhaseMeanResultantLenCField0to1vs3to5(cluster),colorSel);
    
    plotBoxPlot(mod.thetaModHist0to1(indCurCField),...
        mod.thetaModHist3to5(indCurCField),'Theta modulation (0to1 vs 3to5 field)',...
        ['Int_thetaModHist' num2str(cluster) '0to1sVs3to5sFieldBox'],...
        pathAnal,[-0.1 1.2],modSegIntStatsSeg.pRSThetaModHistCField0to1vs3to5(cluster),colorSel);
    
    % bef run vs 3-4s field
    plotBoxPlot(mod.phaseMeanDireBefRun(indCurCField),...
        mod.phaseMeanDire3to4(indCurCField),'Mean theta phase (bef run vs 3to4 field)',...
        ['Int_ThetaMeanC' num2str(cluster) 'BefRunVs3to4sFieldBox'],...
        pathAnal,[0 2*pi],modSegIntStatsSeg.pWWPhaseMeanDireCFieldBefRunvs3to4(cluster),colorSel);
    
    plotBoxPlot(mod.phaseMeanResultantLenBefRun(indCurCField),...
        mod.phaseMeanResultantLen3to4(indCurCField),'Theta phase resultant length (bef run vs 3to4 field)',...
        ['Int_PhaseMeanResultantLen' num2str(cluster) 'BefRunVs3to4sFieldBox'],...
        pathAnal,[-0.1 0.8],modSegIntStatsSeg.pRSPhaseMeanResultantLenCFieldBefRunvs3to4(cluster),colorSel);
    
    plotBoxPlot(mod.thetaModHistBefRun(indCurCField),...
        mod.thetaModHist3to4(indCurCField),'Theta modulation (bef run vs 3to4 field)',...
        ['Int_thetaModHist' num2str(cluster) 'BefRunVs3to4sFieldBox'],...
        pathAnal,[-0.1 1.2],modSegIntStatsSeg.pRSThetaModHistCFieldBefRunvs3to4(cluster),colorSel);
    
    % bef run vs 3-4s no field
    plotBoxPlot(mod.phaseMeanDireBefRun(indCurCNoField),...
        mod.phaseMeanDire3to4(indCurCNoField),'Mean theta phase (bef run vs 3to4 no field)',...
        ['Int_ThetaMeanC' num2str(cluster) 'BefRunVs3to4sNoFieldBox'],...
        pathAnal,[0 2*pi],modSegIntStatsSeg.pWWPhaseMeanDireCNoFieldBefRunvs3to4(cluster),colorSel);
    
    plotBoxPlot(mod.phaseMeanResultantLenBefRun(indCurCNoField),...
        mod.phaseMeanResultantLen3to4(indCurCNoField),'Theta phase resultant length (bef run vs 3to4 no field)',...
        ['Int_PhaseMeanResultantLen' num2str(cluster) 'BefRunVs3to4sNoFieldBox'],...
        pathAnal,[-0.1 0.8],modSegIntStatsSeg.pRSPhaseMeanResultantLenCNoFieldBefRunvs3to4(cluster),colorSel);
    
    plotBoxPlot(mod.thetaModHistBefRun(indCurCNoField),...
        mod.thetaModHist3to4(indCurCNoField),'Theta modulation (bef run vs 3to4 no field)',...
        ['Int_thetaModHist' num2str(cluster) 'BefRunVs3to4sNoFieldBox'],...
        pathAnal,[-0.1 1.2],modSegIntStatsSeg.pRSThetaModHistCNoFieldBefRunvs3to4(cluster),colorSel);
end

function modInt = accumIntNeurons1(paths,filenames,mazeSess,minFRInt,task,onlyRun)

    %% Ints in no cue passive task
    numRec = size(paths,1);
    modInt = struct('task',[],... % no cue - 1, AL - 2, PL - 3
                              'indRec',[],... % recording index
                              'indNeu',[],... % neuron indices
                              ...  
                              'minPhaseFilBefRun',[],... % the phase which fires the least number of spikes 
                              'maxPhaseFilBefRun',[],... % the phase which fires the largest number of spikes 
                              'thetaModHistBefRun',[],... % theta modulation calculated based on theta phase histogram
                              'phaseMeanDireBefRun',[],... % the mean phase direction
                              'phaseMeanResultantLenBefRun',[],... % the mean resultant length of the mean phase direction
                              ...
                              'minPhaseFilH0to1',[],... % the phase which fires the least number of spikes (hilbert)
                              'maxPhaseFilH0to1',[],... % the phase which fires the largest number of spikes (hilbert)
                              'thetaModHistH0to1',[],... % theta modulation calculated based on theta phase histogram (hilbert)
                              'phaseMeanDireH0to1',[],... % the mean phase direction (hilbert)
                              'phaseMeanResultantLenH0to1',[],... % the mean resultant length of the mean phase direction (hilbert)
                              ...
                              'minPhaseFil0to1',[],... % the phase which fires the least number of spikes 
                              'maxPhaseFil0to1',[],... % the phase which fires the largest number of spikes 
                              'thetaModHist0to1',[],... % theta modulation calculated based on theta phase histogram
                              'phaseMeanDire0to1',[],... % the mean phase direction
                              'phaseMeanResultantLen0to1',[],... % the mean resultant length of the mean phase direction
                              ...
                              'minPhaseFilH3to5',[],... % the phase which fires the least number of spikes (hilbert)
                              'maxPhaseFilH3to5',[],... % the phase which fires the largest number of spikes (hilbert)
                              'thetaModHistH3to5',[],... % theta modulation calculated based on theta phase histogram (hilbert)
                              'phaseMeanDireH3to5',[],... % the mean phase direction (hilbert)
                              'phaseMeanResultantLenH3to5',[],... % the mean resultant length of the mean phase direction (hilbert)
                              ...
                              'minPhaseFil3to5',[],... % the phase which fires the least number of spikes 
                              'maxPhaseFil3to5',[],... % the phase which fires the largest number of spikes 
                              'thetaModHist3to5',[],... % theta modulation calculated based on theta phase histogram
                              'phaseMeanDire3to5',[],... % the mean phase direction
                              'phaseMeanResultantLen3to5',[],... % the mean resultant length of the mean phase direction
                              ...
                              'minPhaseFilH3to4',[],... % the phase which fires the least number of spikes (hilbert)
                              'maxPhaseFilH3to4',[],... % the phase which fires the largest number of spikes (hilbert)
                              'thetaModHistH3to4',[],... % theta modulation calculated based on theta phase histogram (hilbert)
                              'phaseMeanDireH3to4',[],... % the mean phase direction (hilbert)
                              'phaseMeanResultantLenH3to4',[],... % the mean resultant length of the mean phase direction (hilbert)
                              ...
                              'minPhaseFil3to4',[],... % the phase which fires the least number of spikes 
                              'maxPhaseFil3to4',[],... % the phase which fires the largest number of spikes 
                              'thetaModHist3to4',[],... % theta modulation calculated based on theta phase histogram
                              'phaseMeanDire3to4',[],... % the mean phase direction
                              'phaseMeanResultantLen3to4',[]); % the mean resultant length of the mean phase direction
                          
    for i = 1:numRec
        disp(filenames(i,:));
        if(i == 7)
            a = 1;
        end
        
        fullPath = [paths(i,:) filenames(i,:) '.mat'];
        if(exist(fullPath) == 0)
            disp('File does not exist.');
            return;
        end
        load(fullPath,'cluList'); 
        
        fileNameInfo = [filenames(i,:) '_Info.mat'];
        fullPath = [paths(i,:) fileNameInfo];
        if(exist(fullPath) == 0)
            disp('_Info.mat file does not exist.');
            return;
        end
        load(fullPath,'autoCorr'); 
        
        fileNameThetaPhaseSeg = [filenames(i,:) '_ThetaPhaseLAlignedSeg_msess'...
            num2str(mazeSess(i)) '_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameThetaPhaseSeg];
        if(exist(fullPath) == 0)
            disp('_ThetaPhaseLAlignedSeg file does not exist.');
            return;
        end
        load(fullPath,'spikeThetaPhaseRunNoStimGood0to1','spikeThetaPhaseRunNoStimGood3to5',...
            'spikeThetaPhaseRunNoStimGood3to4');
        spikeThetaPhaseRunNoStimGood0to1L = spikeThetaPhaseRunNoStimGood0to1;
        spikeThetaPhaseRunNoStimGood3to5L = spikeThetaPhaseRunNoStimGood3to5;
        spikeThetaPhaseRunNoStimGood3to4L = spikeThetaPhaseRunNoStimGood3to4;
        
        fileNameThetaPhaseSeg = [filenames(i,:) '_ThetaPhaseHAlignedSeg_msess'...
            num2str(mazeSess(i)) '_Run' num2str(onlyRun) '.mat'];
        fullPath = [paths(i,:) fileNameThetaPhaseSeg];
        if(exist(fullPath) == 0)
            disp('_ThetaPhaseHAlignedSeg file does not exist.');
            return;
        end
        load(fullPath,'spikeThetaPhaseRunNoStimGood0to1','spikeThetaPhaseRunNoStimGood3to5',...
            'spikeThetaPhaseRunNoStimGood3to4');
        
        fileNameThetaPhaseSegBefRun = [filenames(i,:) '_ThetaPhaseLSeg_msess' ...
            num2str(mazeSess(i)) '_RunOnset0.mat'];
        fullPath = [paths(i,:) fileNameThetaPhaseSegBefRun];
        if(exist(fullPath) == 0)
            disp('_ThetaPhaseLSeg_RunOnset file does not exist.');
            return;
        end
        load(fullPath,'spikeThetaPhaseRunNoStimGoodbefRun');
        
        indNeu = cluList.firingRate > minFRInt &...
                    autoCorr.isInterneuron == 1;
        modInt.task = [modInt.task task*ones(1,sum(indNeu))];
        modInt.indRec = [modInt.indRec i*ones(1,sum(indNeu))];
        modInt.indNeu = [modInt.indNeu find(indNeu == 1)]; 
        
        modInt.minPhaseFil0to1 = [modInt.minPhaseFil0to1 spikeThetaPhaseRunNoStimGood0to1L.minPhaseFilArr(indNeu)];
        modInt.maxPhaseFil0to1 = [modInt.maxPhaseFil0to1 spikeThetaPhaseRunNoStimGood0to1L.maxPhaseFilArr(indNeu)];
        modInt.phaseMeanDire0to1 = [modInt.phaseMeanDire0to1 spikeThetaPhaseRunNoStimGood0to1L.meanDire(indNeu)];
        modInt.phaseMeanResultantLen0to1 = [modInt.phaseMeanResultantLen0to1 spikeThetaPhaseRunNoStimGood0to1L.meanResultantLen(indNeu)];
        modInt.thetaModHist0to1 = [modInt.thetaModHist0to1 spikeThetaPhaseRunNoStimGood0to1L.thetaMod(indNeu)]; 
        
        modInt.minPhaseFilH0to1 = [modInt.minPhaseFilH0to1 spikeThetaPhaseRunNoStimGood0to1.minPhaseFilArr(indNeu)];
        modInt.maxPhaseFilH0to1 = [modInt.maxPhaseFilH0to1 spikeThetaPhaseRunNoStimGood0to1.maxPhaseFilArr(indNeu)];
        modInt.phaseMeanDireH0to1 = [modInt.phaseMeanDireH0to1 spikeThetaPhaseRunNoStimGood0to1.meanDire(indNeu)];
        modInt.phaseMeanResultantLenH0to1 = [modInt.phaseMeanResultantLenH0to1 spikeThetaPhaseRunNoStimGood0to1.meanResultantLen(indNeu)];
        modInt.thetaModHistH0to1 = [modInt.thetaModHistH0to1 spikeThetaPhaseRunNoStimGood0to1.thetaMod(indNeu)]; 
        
        modInt.minPhaseFil3to5 = [modInt.minPhaseFil3to5 spikeThetaPhaseRunNoStimGood3to5L.minPhaseFilArr(indNeu)];
        modInt.maxPhaseFil3to5 = [modInt.maxPhaseFil3to5 spikeThetaPhaseRunNoStimGood3to5L.maxPhaseFilArr(indNeu)];
        modInt.phaseMeanDire3to5 = [modInt.phaseMeanDire3to5 spikeThetaPhaseRunNoStimGood3to5L.meanDire(indNeu)];
        modInt.phaseMeanResultantLen3to5 = [modInt.phaseMeanResultantLen3to5 spikeThetaPhaseRunNoStimGood3to5L.meanResultantLen(indNeu)];
        modInt.thetaModHist3to5 = [modInt.thetaModHist3to5 spikeThetaPhaseRunNoStimGood3to5L.thetaMod(indNeu)]; 
        
        modInt.minPhaseFilH3to5 = [modInt.minPhaseFilH3to5 spikeThetaPhaseRunNoStimGood3to5.minPhaseFilArr(indNeu)];
        modInt.maxPhaseFilH3to5 = [modInt.maxPhaseFilH3to5 spikeThetaPhaseRunNoStimGood3to5.maxPhaseFilArr(indNeu)];
        modInt.phaseMeanDireH3to5 = [modInt.phaseMeanDireH3to5 spikeThetaPhaseRunNoStimGood3to5.meanDire(indNeu)];
        modInt.phaseMeanResultantLenH3to5 = [modInt.phaseMeanResultantLenH3to5 spikeThetaPhaseRunNoStimGood3to5.meanResultantLen(indNeu)];
        modInt.thetaModHistH3to5 = [modInt.thetaModHistH3to5 spikeThetaPhaseRunNoStimGood3to5.thetaMod(indNeu)]; 
        
        modInt.minPhaseFil3to4 = [modInt.minPhaseFil3to4 spikeThetaPhaseRunNoStimGood3to4L.minPhaseFilArr(indNeu)];
        modInt.maxPhaseFil3to4 = [modInt.maxPhaseFil3to4 spikeThetaPhaseRunNoStimGood3to4L.maxPhaseFilArr(indNeu)];
        modInt.phaseMeanDire3to4 = [modInt.phaseMeanDire3to4 spikeThetaPhaseRunNoStimGood3to4L.meanDire(indNeu)];
        modInt.phaseMeanResultantLen3to4 = [modInt.phaseMeanResultantLen3to4 spikeThetaPhaseRunNoStimGood3to4L.meanResultantLen(indNeu)];
        modInt.thetaModHist3to4 = [modInt.thetaModHist3to4 spikeThetaPhaseRunNoStimGood3to4L.thetaMod(indNeu)]; 
        
        modInt.minPhaseFilH3to4 = [modInt.minPhaseFilH3to4 spikeThetaPhaseRunNoStimGood3to4.minPhaseFilArr(indNeu)];
        modInt.maxPhaseFilH3to4 = [modInt.maxPhaseFilH3to4 spikeThetaPhaseRunNoStimGood3to4.maxPhaseFilArr(indNeu)];
        modInt.phaseMeanDireH3to4 = [modInt.phaseMeanDireH3to4 spikeThetaPhaseRunNoStimGood3to4.meanDire(indNeu)];
        modInt.phaseMeanResultantLenH3to4 = [modInt.phaseMeanResultantLenH3to4 spikeThetaPhaseRunNoStimGood3to4.meanResultantLen(indNeu)];
        modInt.thetaModHistH3to4 = [modInt.thetaModHistH3to4 spikeThetaPhaseRunNoStimGood3to4.thetaMod(indNeu)]; 
        
        modInt.minPhaseFilBefRun = [modInt.minPhaseFilBefRun spikeThetaPhaseRunNoStimGoodbefRun.minPhaseFilArr(indNeu)];
        modInt.maxPhaseFilBefRun = [modInt.maxPhaseFilBefRun spikeThetaPhaseRunNoStimGoodbefRun.maxPhaseFilArr(indNeu)];
        modInt.phaseMeanDireBefRun = [modInt.phaseMeanDireBefRun spikeThetaPhaseRunNoStimGoodbefRun.meanDire(indNeu)];
        modInt.phaseMeanResultantLenBefRun = [modInt.phaseMeanResultantLenBefRun spikeThetaPhaseRunNoStimGoodbefRun.meanResultantLen(indNeu)];
        modInt.thetaModHistBefRun = [modInt.thetaModHistBefRun spikeThetaPhaseRunNoStimGoodbefRun.thetaMod(indNeu)]; 
    end
    
    indDire = find(modInt.phaseMeanDireBefRun < 0);
    modInt.phaseMeanDireBefRun(indDire) = modInt.phaseMeanDireBefRun(indDire) + 2*pi;
    
    indDire = find(modInt.phaseMeanDire0to1 < 0);
    modInt.phaseMeanDire0to1(indDire) = modInt.phaseMeanDire0to1(indDire) + 2*pi;
    
    indDire = find(modInt.phaseMeanDire3to5 < 0);
    modInt.phaseMeanDire3to5(indDire) = modInt.phaseMeanDire3to5(indDire) + 2*pi;
    
    indDire = find(modInt.phaseMeanDire3to4 < 0);
    modInt.phaseMeanDire3to4(indDire) = modInt.phaseMeanDire3to4(indDire) + 2*pi;
    
    indDire = find(modInt.phaseMeanDireH0to1 < 0);
    modInt.phaseMeanDireH0to1(indDire) = modInt.phaseMeanDireH0to1(indDire) + 2*pi;
    
    indDire = find(modInt.phaseMeanDireH3to5 < 0);
    modInt.phaseMeanDireH3to5(indDire) = modInt.phaseMeanDireH3to5(indDire) + 2*pi;
    
    indDire = find(modInt.phaseMeanDireH3to4 < 0);
    modInt.phaseMeanDireH3to4(indDire) = modInt.phaseMeanDireH3to4(indDire) + 2*pi;
end

function modIntStatsField = modIntStatsFieldPerC(mod)
    modIntStatsField = [];
    idxC = mod.idxC;
    for i = 1:max(idxC)
        indCurCField = idxC == i & mod.nNeuWithField > 1;
        indCurCNoField = idxC == i & mod.nNeuWithField < 1;
        modIntStatsField.indCurCField{i} = indCurCField;
        modIntStatsField.indCurCNoField{i} = indCurCNoField;
        
        %% before run
        modIntStatsField.meanPhaseMeanDireBefRun(i,:) = [circ_mean(mod.phaseMeanDireBefRun(indCurCField)'),circ_mean(mod.phaseMeanDireBefRun(indCurCNoField)')];
        modIntStatsField.meanMaxPhaseBefRun(i,:) = [circ_mean(mod.maxPhaseArrBefRun(indCurCField)'/180*pi),circ_mean(mod.maxPhaseArrBefRun(indCurCNoField)'/180*pi)];
        modIntStatsField.meanminPhaseBefRun(i,:) = [circ_mean(mod.minPhaseArrBefRun(indCurCField)'/180*pi),circ_mean(mod.minPhaseArrBefRun(indCurCNoField)'/180*pi)];
        modIntStatsField.meanThetaModHistBefRun(i,:) = [mean(mod.thetaModHistBefRun(indCurCField)),mean(mod.thetaModHistBefRun(indCurCNoField))];
        modIntStatsField.meanPhaseMeanResultantLenBefRun(i,:) = [mean(mod.phaseMeanResultantLenBefRun(indCurCField)),mean(mod.phaseMeanResultantLenBefRun(indCurCNoField))];
        
        modIntStatsField.pKPhaseMeanDireBefRun(i) = circ_ktest(mod.phaseMeanDireBefRun(indCurCField)',mod.phaseMeanDireBefRun(indCurCNoField)');
        modIntStatsField.pKMaxPhaseBefRun(i) = circ_ktest(mod.maxPhaseArrBefRun(indCurCField)'/180*pi,mod.maxPhaseArrBefRun(indCurCNoField)'/180*pi);
        modIntStatsField.pKMinPhaseBefRun(i) = circ_ktest(mod.minPhaseArrBefRun(indCurCField)'/180*pi,mod.minPhaseArrBefRun(indCurCNoField)'/180*pi);
        
        modIntStatsField.pWWPhaseMeanDireBefRun(i) = circ_wwtest(mod.phaseMeanDireBefRun(indCurCField)',mod.phaseMeanDireBefRun(indCurCNoField)');
        modIntStatsField.pWWMaxPhaseBefRun(i) = circ_wwtest(mod.maxPhaseArrBefRun(indCurCField)'/180*pi,mod.maxPhaseArrBefRun(indCurCNoField)'/180*pi);
        modIntStatsField.pWWMinPhaseBefRun(i) = circ_wwtest(mod.minPhaseArrBefRun(indCurCField)'/180*pi,mod.minPhaseArrBefRun(indCurCNoField)'/180*pi);
        
        modIntStatsField.pRSThetaModHistBefRun(i) = ranksum(mod.thetaModHistBefRun(indCurCField),mod.thetaModHistBefRun(indCurCNoField));
        
        modIntStatsField.pRSPhaseMeanResultantLenBefRun(i) = ranksum(mod.phaseMeanResultantLenBefRun(indCurCField),mod.phaseMeanResultantLenBefRun(indCurCNoField));
        
        %% 0-1s
        modIntStatsField.meanPhaseMeanDire0to1(i,:) = [circ_mean(mod.phaseMeanDire0to1(indCurCField)'),circ_mean(mod.phaseMeanDire0to1(indCurCNoField)')];
        modIntStatsField.meanPhaseMeanDireH0to1(i,:) = [circ_mean(mod.phaseMeanDireH0to1(indCurCField)'),circ_mean(mod.phaseMeanDireH0to1(indCurCNoField)')];
        modIntStatsField.meanMaxPhase0to1(i,:) = [circ_mean(mod.maxPhaseArr0to1(indCurCField)'/180*pi),circ_mean(mod.maxPhaseArr0to1(indCurCNoField)'/180*pi)];
        modIntStatsField.meanMaxPhaseH0to1(i,:) = [circ_mean(mod.maxPhaseArrH0to1(indCurCField)'/180*pi),circ_mean(mod.maxPhaseArrH0to1(indCurCNoField)'/180*pi)];
        modIntStatsField.meanminPhase0to1(i,:) = [circ_mean(mod.minPhaseArr0to1(indCurCField)'/180*pi),circ_mean(mod.minPhaseArr0to1(indCurCNoField)'/180*pi)];
        modIntStatsField.meanminPhaseH0to1(i,:) = [circ_mean(mod.minPhaseArrH0to1(indCurCField)'/180*pi),circ_mean(mod.minPhaseArrH0to1(indCurCNoField)'/180*pi)];
        modIntStatsField.meanThetaModHist0to1(i,:) = [mean(mod.thetaModHist0to1(indCurCField)),mean(mod.thetaModHist0to1(indCurCNoField))];
        modIntStatsField.meanThetaModHistH0to1(i,:) = [mean(mod.thetaModHistH0to1(indCurCField)),mean(mod.thetaModHistH0to1(indCurCNoField))];       
        modIntStatsField.meanPhaseMeanResultantLen0to1(i,:) = [mean(mod.phaseMeanResultantLen0to1(indCurCField)),mean(mod.phaseMeanResultantLen0to1(indCurCNoField))];
        
        modIntStatsField.pKPhaseMeanDire0to1(i) = circ_ktest(mod.phaseMeanDire0to1(indCurCField)',mod.phaseMeanDire0to1(indCurCNoField)');
        modIntStatsField.pKPhaseMeanDireH0to1(i) = circ_ktest(mod.phaseMeanDireH0to1(indCurCField)',mod.phaseMeanDireH0to1(indCurCNoField)');
        modIntStatsField.pKMaxPhase0to1(i) = circ_ktest(mod.maxPhaseArr0to1(indCurCField)'/180*pi,mod.maxPhaseArr0to1(indCurCNoField)'/180*pi);
        modIntStatsField.pKMaxPhaseH0to1(i) = circ_ktest(mod.maxPhaseArrH0to1(indCurCField)'/180*pi,mod.maxPhaseArrH0to1(indCurCNoField)'/180*pi);
        modIntStatsField.pKMinPhase0to1(i) = circ_ktest(mod.minPhaseArr0to1(indCurCField)'/180*pi,mod.minPhaseArr0to1(indCurCNoField)'/180*pi);
        modIntStatsField.pKMinPhaseH0to1(i) = circ_ktest(mod.minPhaseArrH0to1(indCurCField)'/180*pi,mod.minPhaseArrH0to1(indCurCNoField)'/180*pi);
        
        modIntStatsField.pWWPhaseMeanDire0to1(i) = circ_wwtest(mod.phaseMeanDire0to1(indCurCField)',mod.phaseMeanDire0to1(indCurCNoField)');
        modIntStatsField.pWWPhaseMeanDireH0to1(i) = circ_wwtest(mod.phaseMeanDireH0to1(indCurCField)',mod.phaseMeanDireH0to1(indCurCNoField)');
        modIntStatsField.pWWMaxPhase0to1(i) = circ_wwtest(mod.maxPhaseArr0to1(indCurCField)'/180*pi,mod.maxPhaseArr0to1(indCurCNoField)'/180*pi);
        modIntStatsField.pWWMaxPhaseH0to1(i) = circ_wwtest(mod.maxPhaseArrH0to1(indCurCField)'/180*pi,mod.maxPhaseArrH0to1(indCurCNoField)'/180*pi);
        modIntStatsField.pWWMinPhase0to1(i) = circ_wwtest(mod.minPhaseArr0to1(indCurCField)'/180*pi,mod.minPhaseArr0to1(indCurCNoField)'/180*pi);
        modIntStatsField.pWWMinPhaseH0to1(i) = circ_wwtest(mod.minPhaseArrH0to1(indCurCField)'/180*pi,mod.minPhaseArrH0to1(indCurCNoField)'/180*pi);
        
        modIntStatsField.pRSThetaModHist0to1(i) = ranksum(mod.thetaModHist0to1(indCurCField),mod.thetaModHist0to1(indCurCNoField));
        modIntStatsField.pRSThetaModHistH0to1(i) = ranksum(mod.thetaModHistH0to1(indCurCField),mod.thetaModHistH0to1(indCurCNoField));
        
        modIntStatsField.pRSPhaseMeanResultantLen0to1(i) = ranksum(mod.phaseMeanResultantLen0to1(indCurCField),mod.phaseMeanResultantLen0to1(indCurCNoField));
        
        %% 3-5s
        modIntStatsField.meanPhaseMeanDire3to5(i,:) = [circ_mean(mod.phaseMeanDire3to5(indCurCField)'),circ_mean(mod.phaseMeanDire3to5(indCurCNoField)')];
        modIntStatsField.meanPhaseMeanDireH3to5(i,:) = [circ_mean(mod.phaseMeanDireH3to5(indCurCField)'),circ_mean(mod.phaseMeanDireH3to5(indCurCNoField)')];
        modIntStatsField.meanMaxPhase3to5(i,:) = [circ_mean(mod.maxPhaseArr3to5(indCurCField)'/180*pi),circ_mean(mod.maxPhaseArr3to5(indCurCNoField)'/180*pi)];
        modIntStatsField.meanMaxPhaseH3to5(i,:) = [circ_mean(mod.maxPhaseArrH3to5(indCurCField)'/180*pi),circ_mean(mod.maxPhaseArrH3to5(indCurCNoField)'/180*pi)];
        modIntStatsField.meanminPhase3to5(i,:) = [circ_mean(mod.minPhaseArr3to5(indCurCField)'/180*pi),circ_mean(mod.minPhaseArr3to5(indCurCNoField)'/180*pi)];
        modIntStatsField.meanminPhaseH3to5(i,:) = [circ_mean(mod.minPhaseArrH3to5(indCurCField)'/180*pi),circ_mean(mod.minPhaseArrH3to5(indCurCNoField)'/180*pi)];
        modIntStatsField.meanThetaModHist3to5(i,:) = [mean(mod.thetaModHist3to5(indCurCField)),mean(mod.thetaModHist3to5(indCurCNoField))];
        modIntStatsField.meanThetaModHistH3to5(i,:) = [mean(mod.thetaModHistH3to5(indCurCField)),mean(mod.thetaModHistH3to5(indCurCNoField))];       
        modIntStatsField.meanPhaseMeanResultantLen3to5(i,:) = [mean(mod.phaseMeanResultantLen3to5(indCurCField)),mean(mod.phaseMeanResultantLen3to5(indCurCNoField))];
        
        modIntStatsField.pKPhaseMeanDire3to5(i) = circ_ktest(mod.phaseMeanDire3to5(indCurCField)',mod.phaseMeanDire3to5(indCurCNoField)');
        modIntStatsField.pKPhaseMeanDireH3to5(i) = circ_ktest(mod.phaseMeanDireH3to5(indCurCField)',mod.phaseMeanDireH3to5(indCurCNoField)');
        modIntStatsField.pKMaxPhase3to5(i) = circ_ktest(mod.maxPhaseArr3to5(indCurCField)'/180*pi,mod.maxPhaseArr3to5(indCurCNoField)'/180*pi);
        modIntStatsField.pKMaxPhaseH3to5(i) = circ_ktest(mod.maxPhaseArrH3to5(indCurCField)'/180*pi,mod.maxPhaseArrH3to5(indCurCNoField)'/180*pi);
        modIntStatsField.pKMinPhase3to5(i) = circ_ktest(mod.minPhaseArr3to5(indCurCField)'/180*pi,mod.minPhaseArr3to5(indCurCNoField)'/180*pi);
        modIntStatsField.pKMinPhaseH3to5(i) = circ_ktest(mod.minPhaseArrH3to5(indCurCField)'/180*pi,mod.minPhaseArrH3to5(indCurCNoField)'/180*pi);
        
        modIntStatsField.pWWPhaseMeanDire3to5(i) = circ_wwtest(mod.phaseMeanDire3to5(indCurCField)',mod.phaseMeanDire3to5(indCurCNoField)');
        modIntStatsField.pWWPhaseMeanDireH3to5(i) = circ_wwtest(mod.phaseMeanDireH3to5(indCurCField)',mod.phaseMeanDireH3to5(indCurCNoField)');
        modIntStatsField.pWWMaxPhase3to5(i) = circ_wwtest(mod.maxPhaseArr3to5(indCurCField)'/180*pi,mod.maxPhaseArr3to5(indCurCNoField)'/180*pi);
        modIntStatsField.pWWMaxPhaseH3to5(i) = circ_wwtest(mod.maxPhaseArrH3to5(indCurCField)'/180*pi,mod.maxPhaseArrH3to5(indCurCNoField)'/180*pi);
        modIntStatsField.pWWMinPhase3to5(i) = circ_wwtest(mod.minPhaseArr3to5(indCurCField)'/180*pi,mod.minPhaseArr3to5(indCurCNoField)'/180*pi);
        modIntStatsField.pWWMinPhaseH3to5(i) = circ_wwtest(mod.minPhaseArrH3to5(indCurCField)'/180*pi,mod.minPhaseArrH3to5(indCurCNoField)'/180*pi);
        
        modIntStatsField.pRSThetaModHist3to5(i) = ranksum(mod.thetaModHist3to5(indCurCField),mod.thetaModHist3to5(indCurCNoField));
        modIntStatsField.pRSThetaModHistH3to5(i) = ranksum(mod.thetaModHistH3to5(indCurCField),mod.thetaModHistH3to5(indCurCNoField));
        
        modIntStatsField.pRSPhaseMeanResultantLen3to5(i) = ranksum(mod.phaseMeanResultantLen3to5(indCurCField),mod.phaseMeanResultantLen3to5(indCurCNoField));
        
        %% 3-4s
        modIntStatsField.meanPhaseMeanDire3to4(i,:) = [circ_mean(mod.phaseMeanDire3to4(indCurCField)'),circ_mean(mod.phaseMeanDire3to4(indCurCNoField)')];
        modIntStatsField.meanPhaseMeanDireH3to4(i,:) = [circ_mean(mod.phaseMeanDireH3to4(indCurCField)'),circ_mean(mod.phaseMeanDireH3to4(indCurCNoField)')];
        modIntStatsField.meanMaxPhase3to4(i,:) = [circ_mean(mod.maxPhaseArr3to4(indCurCField)'/180*pi),circ_mean(mod.maxPhaseArr3to4(indCurCNoField)'/180*pi)];
        modIntStatsField.meanMaxPhaseH3to4(i,:) = [circ_mean(mod.maxPhaseArrH3to4(indCurCField)'/180*pi),circ_mean(mod.maxPhaseArrH3to4(indCurCNoField)'/180*pi)];
        modIntStatsField.meanminPhase3to4(i,:) = [circ_mean(mod.minPhaseArr3to4(indCurCField)'/180*pi),circ_mean(mod.minPhaseArr3to4(indCurCNoField)'/180*pi)];
        modIntStatsField.meanminPhaseH3to4(i,:) = [circ_mean(mod.minPhaseArrH3to4(indCurCField)'/180*pi),circ_mean(mod.minPhaseArrH3to4(indCurCNoField)'/180*pi)];
        modIntStatsField.meanThetaModHist3to4(i,:) = [mean(mod.thetaModHist3to4(indCurCField)),mean(mod.thetaModHist3to4(indCurCNoField))];
        modIntStatsField.meanThetaModHistH3to4(i,:) = [mean(mod.thetaModHistH3to4(indCurCField)),mean(mod.thetaModHistH3to4(indCurCNoField))];       
        modIntStatsField.meanPhaseMeanResultantLen3to4(i,:) = [mean(mod.phaseMeanResultantLen3to4(indCurCField)),mean(mod.phaseMeanResultantLen3to4(indCurCNoField))];
        
        modIntStatsField.pKPhaseMeanDire3to4(i) = circ_ktest(mod.phaseMeanDire3to4(indCurCField)',mod.phaseMeanDire3to4(indCurCNoField)');
        modIntStatsField.pKPhaseMeanDireH3to4(i) = circ_ktest(mod.phaseMeanDireH3to4(indCurCField)',mod.phaseMeanDireH3to4(indCurCNoField)');
        modIntStatsField.pKMaxPhase3to4(i) = circ_ktest(mod.maxPhaseArr3to4(indCurCField)'/180*pi,mod.maxPhaseArr3to4(indCurCNoField)'/180*pi);
        modIntStatsField.pKMaxPhaseH3to4(i) = circ_ktest(mod.maxPhaseArrH3to4(indCurCField)'/180*pi,mod.maxPhaseArrH3to4(indCurCNoField)'/180*pi);
        modIntStatsField.pKMinPhase3to4(i) = circ_ktest(mod.minPhaseArr3to4(indCurCField)'/180*pi,mod.minPhaseArr3to4(indCurCNoField)'/180*pi);
        modIntStatsField.pKMinPhaseH3to4(i) = circ_ktest(mod.minPhaseArrH3to4(indCurCField)'/180*pi,mod.minPhaseArrH3to4(indCurCNoField)'/180*pi);
        
        modIntStatsField.pWWPhaseMeanDire3to4(i) = circ_wwtest(mod.phaseMeanDire3to4(indCurCField)',mod.phaseMeanDire3to4(indCurCNoField)');
        modIntStatsField.pWWPhaseMeanDireH3to4(i) = circ_wwtest(mod.phaseMeanDireH3to4(indCurCField)',mod.phaseMeanDireH3to4(indCurCNoField)');
        modIntStatsField.pWWMaxPhase3to4(i) = circ_wwtest(mod.maxPhaseArr3to4(indCurCField)'/180*pi,mod.maxPhaseArr3to4(indCurCNoField)'/180*pi);
        modIntStatsField.pWWMaxPhaseH3to4(i) = circ_wwtest(mod.maxPhaseArrH3to4(indCurCField)'/180*pi,mod.maxPhaseArrH3to4(indCurCNoField)'/180*pi);
        modIntStatsField.pWWMinPhase3to4(i) = circ_wwtest(mod.minPhaseArr3to4(indCurCField)'/180*pi,mod.minPhaseArr3to4(indCurCNoField)'/180*pi);
        modIntStatsField.pWWMinPhaseH3to4(i) = circ_wwtest(mod.minPhaseArrH3to4(indCurCField)'/180*pi,mod.minPhaseArrH3to4(indCurCNoField)'/180*pi);
        
        modIntStatsField.pRSThetaModHist3to4(i) = ranksum(mod.thetaModHist3to4(indCurCField),mod.thetaModHist3to4(indCurCNoField));
        modIntStatsField.pRSThetaModHistH3to4(i) = ranksum(mod.thetaModHistH3to4(indCurCField),mod.thetaModHistH3to4(indCurCNoField));
        
        modIntStatsField.pRSPhaseMeanResultantLen3to4(i) = ranksum(mod.phaseMeanResultantLen3to4(indCurCField),mod.phaseMeanResultantLen3to4(indCurCNoField));
        
    end
end

function modIntStatsField = modIntStatsSegments(mod)
    idxC = mod.idxC;
    for i = 1:max(idxC)
        %% befRun s vs 3-4 s
        modIntStatsField.pWWPhaseMeanDireCBefRunvs3to4(i) = circ_wwtest(mod.phaseMeanDireBefRun(idxC == i)',mod.phaseMeanDire3to4(idxC == i)');
        modIntStatsField.pKPhaseMeanDireCBefRunvs3to4(i) = circ_ktest(mod.phaseMeanDireBefRun(idxC == i)',mod.phaseMeanDire3to4(idxC == i)');
        modIntStatsField.pWWMaxPhaseCBefRunvs3to4(i) = circ_wwtest(mod.maxPhaseArrBefRun(idxC == i)'/180*pi,mod.maxPhaseArr3to4(idxC == i)'/180*pi);
        modIntStatsField.pKMaxPhaseCBefRunvs3to4(i) = circ_ktest(mod.maxPhaseArrBefRun(idxC == i)'/180*pi,mod.maxPhaseArr3to4(idxC == i)'/180*pi);
        modIntStatsField.pWWMinPhaseCBefRunvs3to4(i) = circ_wwtest(mod.minPhaseArrBefRun(idxC == i)'/180*pi,mod.minPhaseArr3to4(idxC == i)'/180*pi);
        modIntStatsField.pKMinPhaseCBefRunvs3to4(i) = circ_ktest(mod.minPhaseArrBefRun(idxC == i)'/180*pi,mod.minPhaseArr3to4(idxC == i)'/180*pi);

        modIntStatsField.pRSThetaModHistCBefRunvs3to4(i) = ranksum(mod.thetaModHistBefRun(idxC == i),mod.thetaModHist3to4(idxC == i));
        
        modIntStatsField.pRSPhaseMeanResultantLenCBefRunvs3to4(i) = ranksum(mod.phaseMeanResultantLenBefRun(idxC == i),mod.phaseMeanResultantLen3to4(idxC == i));

        %% 0-1 s vs 3-4 s
        modIntStatsField.pWWPhaseMeanDireC0to1vs3to4(i) = circ_wwtest(mod.phaseMeanDire0to1(idxC == i)',mod.phaseMeanDire3to4(idxC == i)');
        modIntStatsField.pKPhaseMeanDireC0to1vs3to4(i) = circ_ktest(mod.phaseMeanDire0to1(idxC == i)',mod.phaseMeanDire3to4(idxC == i)');
        modIntStatsField.pWWPhaseMeanDireHC0to1vs3to4(i) = circ_wwtest(mod.phaseMeanDireH0to1(idxC == i)',mod.phaseMeanDireH3to4(idxC == i)');
        modIntStatsField.pKPhaseMeanDireHC0to1vs3to4(i) = circ_ktest(mod.phaseMeanDireH0to1(idxC == i)',mod.phaseMeanDireH3to4(idxC == i)');
        modIntStatsField.pWWMaxPhaseC0to1vs3to4(i) = circ_wwtest(mod.maxPhaseArr0to1(idxC == i)'/180*pi,mod.maxPhaseArr3to4(idxC == i)'/180*pi);
        modIntStatsField.pKMaxPhaseC0to1vs3to4(i) = circ_ktest(mod.maxPhaseArr0to1(idxC == i)'/180*pi,mod.maxPhaseArr3to4(idxC == i)'/180*pi);
        modIntStatsField.pWWMinPhaseC0to1vs3to4(i) = circ_wwtest(mod.minPhaseArr0to1(idxC == i)'/180*pi,mod.minPhaseArr3to4(idxC == i)'/180*pi);
        modIntStatsField.pKMinPhaseC0to1vs3to4(i) = circ_ktest(mod.minPhaseArr0to1(idxC == i)'/180*pi,mod.minPhaseArr3to4(idxC == i)'/180*pi);

        modIntStatsField.pRSThetaModHistC0to1vs3to4(i) = ranksum(mod.thetaModHist0to1(idxC == i),mod.thetaModHist3to4(idxC == i));
        modIntStatsField.pRSThetaModHistHC0to1vs3to4(i) = ranksum(mod.thetaModHistH0to1(idxC == i),mod.thetaModHistH3to4(idxC == i));

        modIntStatsField.pRSPhaseMeanResultantLenC0to1vs3to4(i) = ranksum(mod.phaseMeanResultantLen0to1(idxC == i),mod.phaseMeanResultantLen3to4(idxC == i));

        %% 0-1 s vs 3-5 s
        modIntStatsField.pWWPhaseMeanDireC0to1vs3to5(i) = circ_wwtest(mod.phaseMeanDire0to1(idxC == i)',mod.phaseMeanDire3to5(idxC == i)');
        modIntStatsField.pKPhaseMeanDireC0to1vs3to5(i) = circ_ktest(mod.phaseMeanDire0to1(idxC == i)',mod.phaseMeanDire3to5(idxC == i)');
        modIntStatsField.pWWPhaseMeanDireHC0to1vs3to5(i) = circ_wwtest(mod.phaseMeanDireH0to1(idxC == i)',mod.phaseMeanDireH3to5(idxC == i)');
        modIntStatsField.pKPhaseMeanDireHC0to1vs3to5(i) = circ_ktest(mod.phaseMeanDireH0to1(idxC == i)',mod.phaseMeanDireH3to5(idxC == i)');
        modIntStatsField.pWWMaxPhaseC0to1vs3to5(i) = circ_wwtest(mod.maxPhaseArr0to1(idxC == i)'/180*pi,mod.maxPhaseArr3to5(idxC == i)'/180*pi);
        modIntStatsField.pKMaxPhaseC0to1vs3to5(i) = circ_ktest(mod.maxPhaseArr0to1(idxC == i)'/180*pi,mod.maxPhaseArr3to5(idxC == i)'/180*pi);
        modIntStatsField.pWWMinPhaseC0to1vs3to5(i) = circ_wwtest(mod.minPhaseArr0to1(idxC == i)'/180*pi,mod.minPhaseArr3to5(idxC == i)'/180*pi);
        modIntStatsField.pKMinPhaseC0to1vs3to5(i) = circ_ktest(mod.minPhaseArr0to1(idxC == i)'/180*pi,mod.minPhaseArr3to5(idxC == i)'/180*pi);

        modIntStatsField.pRSThetaModHistC0to1vs3to5(i) = ranksum(mod.thetaModHist0to1(idxC == i),mod.thetaModHist3to5(idxC == i));
        modIntStatsField.pRSThetaModHistHC0to1vs3to5(i) = ranksum(mod.thetaModHistH0to1(idxC == i),mod.thetaModHistH3to5(idxC == i));

        modIntStatsField.pRSPhaseMeanResultantLenC0to1vs3to5(i) = ranksum(mod.phaseMeanResultantLen0to1(idxC == i),mod.phaseMeanResultantLen3to5(idxC == i));

        %% befRun s vs 3-4 s field recordings
        idxCField = mod.nNeuWithField > 1 & idxC == i;
        modIntStatsField.pWWPhaseMeanDireCFieldBefRunvs3to4(i) = circ_wwtest(mod.phaseMeanDireBefRun(idxCField)',mod.phaseMeanDire3to4(idxCField)');
        modIntStatsField.pKPhaseMeanDireCFieldBefRunvs3to4(i) = circ_ktest(mod.phaseMeanDireBefRun(idxCField)',mod.phaseMeanDire3to4(idxCField)');
        modIntStatsField.pWWMaxPhaseCFieldBefRunvs3to4(i) = circ_wwtest(mod.maxPhaseArrBefRun(idxCField)'/180*pi,mod.maxPhaseArr3to4(idxCField)'/180*pi);
        modIntStatsField.pKMaxPhaseCFieldBefRunvs3to4(i) = circ_ktest(mod.maxPhaseArrBefRun(idxCField)'/180*pi,mod.maxPhaseArr3to4(idxCField)'/180*pi);
        modIntStatsField.pWWMinPhaseCFieldBefRunvs3to4(i) = circ_wwtest(mod.minPhaseArrBefRun(idxCField)'/180*pi,mod.minPhaseArr3to4(idxCField)'/180*pi);
        modIntStatsField.pKMinPhaseCFieldBefRunvs3to4(i) = circ_ktest(mod.minPhaseArrBefRun(idxCField)'/180*pi,mod.minPhaseArr3to4(idxCField)'/180*pi);

        modIntStatsField.pRSThetaModHistCFieldBefRunvs3to4(i) = ranksum(mod.thetaModHistBefRun(idxCField),mod.thetaModHist3to4(idxCField));
        
        modIntStatsField.pRSPhaseMeanResultantLenCFieldBefRunvs3to4(i) = ranksum(mod.phaseMeanResultantLenBefRun(idxCField),mod.phaseMeanResultantLen3to4(idxCField));

        %% 0-1 s vs 3-4 s field recordings
        modIntStatsField.pWWPhaseMeanDireCField0to1vs3to4(i) = circ_wwtest(mod.phaseMeanDire0to1(idxCField)',mod.phaseMeanDire3to4(idxCField)');
        modIntStatsField.pKPhaseMeanDireCField0to1vs3to4(i) = circ_ktest(mod.phaseMeanDire0to1(idxCField)',mod.phaseMeanDire3to4(idxCField)');
        modIntStatsField.pWWPhaseMeanDireHCField0to1vs3to4(i) = circ_wwtest(mod.phaseMeanDireH0to1(idxCField)',mod.phaseMeanDireH3to4(idxCField)');
        modIntStatsField.pKPhaseMeanDireHCField0to1vs3to4(i) = circ_ktest(mod.phaseMeanDireH0to1(idxCField)',mod.phaseMeanDireH3to4(idxCField)');
        modIntStatsField.pWWMaxPhaseCField0to1vs3to4(i) = circ_wwtest(mod.maxPhaseArr0to1(idxCField)'/180*pi,mod.maxPhaseArr3to4(idxCField)'/180*pi);
        modIntStatsField.pKMaxPhaseCField0to1vs3to4(i) = circ_ktest(mod.maxPhaseArr0to1(idxCField)'/180*pi,mod.maxPhaseArr3to4(idxCField)'/180*pi);
        modIntStatsField.pWWMinPhaseCField0to1vs3to4(i) = circ_wwtest(mod.minPhaseArr0to1(idxCField)'/180*pi,mod.minPhaseArr3to4(idxCField)'/180*pi);
        modIntStatsField.pKMinPhaseCField0to1vs3to4(i) = circ_ktest(mod.minPhaseArr0to1(idxCField)'/180*pi,mod.minPhaseArr3to4(idxCField)'/180*pi);

        modIntStatsField.pRSThetaModHistCField0to1vs3to4(i) = ranksum(mod.thetaModHist0to1(idxCField),mod.thetaModHist3to4(idxCField));
        modIntStatsField.pRSThetaModHistHCField0to1vs3to4(i) = ranksum(mod.thetaModHistH0to1(idxCField),mod.thetaModHistH3to4(idxCField));

        modIntStatsField.pRSPhaseMeanResultantLenCField0to1vs3to4(i) = ranksum(mod.phaseMeanResultantLen0to1(idxCField),mod.phaseMeanResultantLen3to4(idxCField));

        %% 0-1 s vs 3-5 s field recordings
        modIntStatsField.pWWPhaseMeanDireCField0to1vs3to5(i) = circ_wwtest(mod.phaseMeanDire0to1(idxCField)',mod.phaseMeanDire3to5(idxCField)');
        modIntStatsField.pKPhaseMeanDireCField0to1vs3to5(i) = circ_ktest(mod.phaseMeanDire0to1(idxCField)',mod.phaseMeanDire3to5(idxCField)');
        modIntStatsField.pWWPhaseMeanDireHCField0to1vs3to5(i) = circ_wwtest(mod.phaseMeanDireH0to1(idxCField)',mod.phaseMeanDireH3to5(idxCField)');
        modIntStatsField.pKPhaseMeanDireHCField0to1vs3to5(i) = circ_ktest(mod.phaseMeanDireH0to1(idxCField)',mod.phaseMeanDireH3to5(idxCField)');
        modIntStatsField.pWWMaxPhaseCField0to1vs3to5(i) = circ_wwtest(mod.maxPhaseArr0to1(idxCField)'/180*pi,mod.maxPhaseArr3to5(idxCField)'/180*pi);
        modIntStatsField.pKMaxPhaseCField0to1vs3to5(i) = circ_ktest(mod.maxPhaseArr0to1(idxCField)'/180*pi,mod.maxPhaseArr3to5(idxCField)'/180*pi);
        modIntStatsField.pWWMinPhaseCField0to1vs3to5(i) = circ_wwtest(mod.minPhaseArr0to1(idxCField)'/180*pi,mod.minPhaseArr3to5(idxCField)'/180*pi);
        modIntStatsField.pKMinPhaseCField0to1vs3to5(i) = circ_ktest(mod.minPhaseArr0to1(idxCField)'/180*pi,mod.minPhaseArr3to5(idxCField)'/180*pi);

        modIntStatsField.pRSThetaModHistCField0to1vs3to5(i) = ranksum(mod.thetaModHist0to1(idxCField),mod.thetaModHist3to5(idxCField));
        modIntStatsField.pRSThetaModHistHCField0to1vs3to5(i) = ranksum(mod.thetaModHistH0to1(idxCField),mod.thetaModHistH3to5(idxCField));

        modIntStatsField.pRSPhaseMeanResultantLenCField0to1vs3to5(i) = ranksum(mod.phaseMeanResultantLen0to1(idxCField),mod.phaseMeanResultantLen3to5(idxCField));

        %% 0-1 s vs 3-4 s no field recordings
        idxCNoField = mod.nNeuWithField < 1 & idxC == i;
        modIntStatsField.pWWPhaseMeanDireCNoField0to1vs3to4(i) = circ_wwtest(mod.phaseMeanDire0to1(idxCNoField)',mod.phaseMeanDire3to4(idxCNoField)');
        modIntStatsField.pKPhaseMeanDireCNoField0to1vs3to4(i) = circ_ktest(mod.phaseMeanDire0to1(idxCNoField)',mod.phaseMeanDire3to4(idxCNoField)');
        modIntStatsField.pWWPhaseMeanDireHCNoField0to1vs3to4(i) = circ_wwtest(mod.phaseMeanDireH0to1(idxCNoField)',mod.phaseMeanDireH3to4(idxCNoField)');
        modIntStatsField.pKPhaseMeanDireHCNoField0to1vs3to4(i) = circ_ktest(mod.phaseMeanDireH0to1(idxCNoField)',mod.phaseMeanDireH3to4(idxCNoField)');
        modIntStatsField.pWWMaxPhaseCNoField0to1vs3to4(i) = circ_wwtest(mod.maxPhaseArr0to1(idxCNoField)'/180*pi,mod.maxPhaseArr3to4(idxCNoField)'/180*pi);
        modIntStatsField.pKMaxPhaseCNoField0to1vs3to4(i) = circ_ktest(mod.maxPhaseArr0to1(idxCNoField)'/180*pi,mod.maxPhaseArr3to4(idxCNoField)'/180*pi);
        modIntStatsField.pWWMinPhaseCNoField0to1vs3to4(i) = circ_wwtest(mod.minPhaseArr0to1(idxCNoField)'/180*pi,mod.minPhaseArr3to4(idxCNoField)'/180*pi);
        modIntStatsField.pKMinPhaseCNoField0to1vs3to4(i) = circ_ktest(mod.minPhaseArr0to1(idxCNoField)'/180*pi,mod.minPhaseArr3to4(idxCNoField)'/180*pi);

        modIntStatsField.pRSThetaModHistCNoField0to1vs3to4(i) = ranksum(mod.thetaModHist0to1(idxCNoField),mod.thetaModHist3to4(idxCNoField));
        modIntStatsField.pRSThetaModHistHCNoField0to1vs3to4(i) = ranksum(mod.thetaModHistH0to1(idxCNoField),mod.thetaModHistH3to4(idxCNoField));

        modIntStatsField.pRSPhaseMeanResultantLenCNoField0to1vs3to4(i) = ranksum(mod.phaseMeanResultantLen0to1(idxCNoField),mod.phaseMeanResultantLen3to4(idxCNoField));

        %% 0-1 s vs 3-5 s no field recordings
        modIntStatsField.pWWPhaseMeanDireCNoField0to1vs3to5(i) = circ_wwtest(mod.phaseMeanDire0to1(idxCNoField)',mod.phaseMeanDire3to5(idxCNoField)');
        modIntStatsField.pKPhaseMeanDireCNoField0to1vs3to5(i) = circ_ktest(mod.phaseMeanDire0to1(idxCNoField)',mod.phaseMeanDire3to5(idxCNoField)');
        modIntStatsField.pWWPhaseMeanDireHCNoField0to1vs3to5(i) = circ_wwtest(mod.phaseMeanDireH0to1(idxCNoField)',mod.phaseMeanDireH3to5(idxCNoField)');
        modIntStatsField.pKPhaseMeanDireHCNoField0to1vs3to5(i) = circ_ktest(mod.phaseMeanDireH0to1(idxCNoField)',mod.phaseMeanDireH3to5(idxCNoField)');
        modIntStatsField.pWWMaxPhaseCNoField0to1vs3to5(i) = circ_wwtest(mod.maxPhaseArr0to1(idxCNoField)'/180*pi,mod.maxPhaseArr3to5(idxCNoField)'/180*pi);
        modIntStatsField.pKMaxPhaseCNoField0to1vs3to5(i) = circ_ktest(mod.maxPhaseArr0to1(idxCNoField)'/180*pi,mod.maxPhaseArr3to5(idxCNoField)'/180*pi);
        modIntStatsField.pWWMinPhaseCNoField0to1vs3to5(i) = circ_wwtest(mod.minPhaseArr0to1(idxCNoField)'/180*pi,mod.minPhaseArr3to5(idxCNoField)'/180*pi);
        modIntStatsField.pKMinPhaseCNoField0to1vs3to5(i) = circ_ktest(mod.minPhaseArr0to1(idxCNoField)'/180*pi,mod.minPhaseArr3to5(idxCNoField)'/180*pi);

        modIntStatsField.pRSThetaModHistCNoField0to1vs3to5(i) = ranksum(mod.thetaModHist0to1(idxCNoField),mod.thetaModHist3to5(idxCNoField));
        modIntStatsField.pRSThetaModHistHCNoField0to1vs3to5(i) = ranksum(mod.thetaModHistH0to1(idxCNoField),mod.thetaModHistH3to5(idxCNoField));

        modIntStatsField.pRSPhaseMeanResultantLenCNoField0to1vs3to5(i) = ranksum(mod.phaseMeanResultantLen0to1(idxCNoField),mod.phaseMeanResultantLen3to5(idxCNoField));

        %% befRun s vs 3-4 s no field recordings
        modIntStatsField.pWWPhaseMeanDireCNoFieldBefRunvs3to4(i) = circ_wwtest(mod.phaseMeanDireBefRun(idxCNoField)',mod.phaseMeanDire3to4(idxCNoField)');
        modIntStatsField.pKPhaseMeanDireCNoFieldBefRunvs3to4(i) = circ_ktest(mod.phaseMeanDireBefRun(idxCNoField)',mod.phaseMeanDire3to4(idxCNoField)');
        modIntStatsField.pWWMaxPhaseCNoFieldBefRunvs3to4(i) = circ_wwtest(mod.maxPhaseArrBefRun(idxCNoField)'/180*pi,mod.maxPhaseArr3to4(idxCNoField)'/180*pi);
        modIntStatsField.pKMaxPhaseCNoFieldBefRunvs3to4(i) = circ_ktest(mod.maxPhaseArrBefRun(idxCNoField)'/180*pi,mod.maxPhaseArr3to4(idxCNoField)'/180*pi);
        modIntStatsField.pWWMinPhaseCNoFieldBefRunvs3to4(i) = circ_wwtest(mod.minPhaseArrBefRun(idxCNoField)'/180*pi,mod.minPhaseArr3to4(idxCNoField)'/180*pi);
        modIntStatsField.pKMinPhaseCNoFieldBefRunvs3to4(i) = circ_ktest(mod.minPhaseArrBefRun(idxCNoField)'/180*pi,mod.minPhaseArr3to4(idxCNoField)'/180*pi);

        modIntStatsField.pRSThetaModHistCNoFieldBefRunvs3to4(i) = ranksum(mod.thetaModHistBefRun(idxCNoField),mod.thetaModHist3to4(idxCNoField));
        
        modIntStatsField.pRSPhaseMeanResultantLenCNoFieldBefRunvs3to4(i) = ranksum(mod.phaseMeanResultantLenBefRun(idxCNoField),mod.phaseMeanResultantLen3to4(idxCNoField));

    end
end

function plotBoxPlot(x1,x2,yl,fn,pathAnal,ylimit,p,colorSel)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 200 400])
    if(colorSel == 6)
        colorArr = [...
                127 127 229;...
                125 125 168]/255;
    elseif(colorSel == 1)            
        colorArr = [...
                230 127 128;...
                153 85 86]/255;
    else        
        colorArr = [163 207 98;... 
            63 79 37]/255;
    end
    x = [x1';x2'];
    g = [repmat({'C1'},length(x1),1);...
        repmat({'C2'},length(x2),1)];
    boxplot(x,g,'Notch','on','Widths',0.3,'Symbol','');
    h = findobj(gca,'Tag','Box');
    for j = 1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colorArr(j,:),'FaceAlpha',0.5);
    end
    if(~isempty(ylimit))
        set(gca,'YLim',ylimit);
    end
    ylabel(yl);
    title(['p = ' num2str(p)]);
    
    fileName1 = [pathAnal fn];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotBarPlot(x1,x2,yl,fn,pathAnal,ylimit,p,colorSel)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 200 400])
    if(colorSel == 0)
        colorArr = [163 207 98;...
                234 131 114]/255;
    elseif(colorSel == 1)            
        colorArr = [234 131 114;...
                116 53 61]/255;
    else        
        colorArr = [163 207 98;... 
            63 79 37]/255;
    end
    meanx1 = mean(x1);
    meanx2 = mean(x2);
    sem1 = std(x1)/sqrt(size(x1,1));
    sem2 = std(x2)/sqrt(size(x2,1));
    hold on;
    h = bar(1,meanx1,0.6,'FaceColor',colorArr(2,:));
    set(h,'FaceAlpha',0.5);
    h = bar(2,meanx2,0.6,'FaceColor',colorArr(1,:));
    set(h,'FaceAlpha',0.5);
    h = errorbar(1,meanx1,sem1,'LineWidth',1,'Marker','.','Color','k');    
    h = errorbar(2,meanx2,sem2,'LineWidth',1,'Marker','.','Color','k');
    if(~isempty(ylimit))
        set(gca,'YLim',ylimit);
    end
    ylabel(yl);
    title(['p = ' num2str(p)]);
    
    fileName1 = [pathAnal fn];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotPolarPlot(x1,x2,ti,fn,pathAnal,p)
    [figNew,pos] = CreateFig();
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    colorArr = [163 207 98;...
                234 131 114]/255;
    
    polarhistogram(x1,0:pi/18:2*pi,'Normalization','probability','DisplayStyle','bar',....
        'FaceAlpha',0.5,'FaceColor',colorArr(2,:),'EdgeColor',[0.5 0.5 0.5]);
    hold on;
    polarhistogram(x2,0:pi/18:2*pi,'Normalization','probability','DisplayStyle','bar',....
        'FaceAlpha',0.5,'FaceColor',colorArr(1,:),'EdgeColor',[0.5 0.5 0.5]);
    title([ti ' p = ' num2str(p)])    
    
    fileName1 = [pathAnal fn];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotDistri(x,y,histBins,xl,leg,pathAnal,fileN,xlimit)
    
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    colorArr = [163 207 98;...
                234 131 114]/255;
    hold on;
    hx = hist(x,histBins);
    hy = hist(y,histBins);
    h = plot(histBins,hx/sum(hx),'-');
    set(h,'LineWidth',2,'Color',colorArr(1,:));
    h = plot(histBins,hy/sum(hy),'-');
    set(h,'LineWidth',2,'Color',colorArr(2,:));
    if(~isempty(xlimit))
        set(gca,'xLim',xlimit);
    end
    xlabel(xl);
    ylabel('Probability');
    legend(leg);
    
    fileName1 = [pathAnal fileN];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotConditions(x1,y1,x2,y2,x3,y3,xl,yl,ti)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    colorArr = [0.5 0.5 0.9;...
                0.9 0.5 0.5;...
                0.2 0.8 0.5];
    h = plot(x1,y1,'o');
    set(h,'MarkerSize',6,'Color',colorArr(1,:));
    hold on;
    h = plot(x2,y2,'o');
    set(h,'MarkerSize',6,'Color',colorArr(2,:));
    h = plot(x3,y3,'o');
    set(h,'MarkerSize',6,'Color',colorArr(3,:));
    maxX = max([x1,x2,x3]);
    maxY = max([y1,y2,y3]);
    minX = min([x1,x2,x3]);
    minY = min([y1,y2,y3]);
    set(gca,'XLim',[minX maxX],'YLim',[minY maxY]);
    xlabel(xl)
    ylabel(yl)
    title(ti)
    legend('NoCue','AL','PL')
end

function plotClusters(x,y,idx,xl,yl,ti)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])

    % only for plotting two clusters with fields
    colorArr = [...
                234 131 114;...
                163 207 98]/255;
            
%     colorArr = [0.5 0.5 0.9;...
%         0.9 0.5 0.5;...
%         0.3 0.3 0.7;...
%         0.7 0.3 0.3;...
%         0.5 0.9 0.5;...
%         0.2 0.5 0.8;...
%         0.2 0.8 0.5;...
%         0.8 0.5 0.2;...
%         0.3 0.7 0.3];
    hold on;
    for i = 1:max(idx)
        indTmp = idx == i;
        h = plot(x(indTmp),y(indTmp),'.');
        set(h,'MarkerSize',8,'Color',colorArr(mod(i,2)+1,:));
    end
%     h = plot(xt,yt,'k+');
%     set(h,'MarkerSize',7);
    maxX = max(x);
    maxY = max(y);
    minX = min(x);
    minY = min(y);
    set(gca,'XLim',[minX maxX],'YLim',[minY maxY]);
    xlabel(xl)
    ylabel(yl)
    title(ti);
end

function plotClustersFieldsPhase(pathAnal,x,y,idx,indField,xl,yl,ti,fileN)
    plotClusters(x(indField),y(indField),...
        idx(indField),xl,yl,[ti ' - field']);
    plot([0 2*pi],[0 2*pi],'k-')
    
    fileName1 = [pathAnal fileN 'Field'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    if(sum(~indField) > 0)
        plotClusters(x(~indField),y(~indField),...
            idx(~indField),xl,yl,[ti ' - no field']);
        plot([0 2*pi],[0 2*pi],'k-')

        fileName1 = [pathAnal fileN 'NoField'];
        saveas(gcf,fileName1);
        print('-painters', '-dpdf', fileName1, '-r600')
    end
end

function plotClustersFieldsPhaseLabelF(pathAnal,x,y,idx,xfF,yfF,xfNoF,yfNoF,indField,xl,yl,ti,fileN)
    plotClusters(x(indField),y(indField),...
        idx(indField),xl,yl,[ti ' - field']);
%     plot([0 2*pi],[0 2*pi],'k-')
    h = plot(xfF,yfF,'k+');
    set(h,'MarkerSize',9);
    
    fileName1 = [pathAnal fileN 'Field'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    if(sum(~indField) > 0)
        plotClusters(x(~indField),y(~indField),...
            idx(~indField),xl,yl,[ti ' - no field']);
    %     plot([0 2*pi],[0 2*pi],'k-')
        h = plot(xfNoF,yfNoF,'k+');
        set(h,'MarkerSize',9);

        fileName1 = [pathAnal fileN 'NoField'];
        saveas(gcf,fileName1);
        print('-painters', '-dpdf', fileName1, '-r600')
    end
end

function plotClustersFields(pathAnal,x,y,idx,indField,xl,yl,ti,fileN)
    plotClusters(x(indField),y(indField),...
        idx(indField),xl,yl,[ti ' - field']);
    
    fileName1 = [pathAnal fileN 'Field'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotClusters(x(~indField),y(~indField),...
        idx(~indField),xl,yl,[ti ' - no field']);
        
    fileName1 = [pathAnal fileN 'NoField'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end

function plotClustersFieldsLabelF(pathAnal,x,y,idx,xfF,yfF,xfNoF,yfNoF,indField,xl,yl,ti,fileN)
    plotClusters(x(indField),y(indField),...
        idx(indField),xl,yl,[ti ' - field']);
    h = plot(xfF,yfF,'k+');
    set(h,'MarkerSize',9);
    
    fileName1 = [pathAnal fileN 'Field'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotClusters(x(~indField),y(~indField),...
        idx(~indField),xl,yl,[ti ' - no field']);
    h = plot(xfNoF,yfNoF,'k+');
    set(h,'MarkerSize',9);
    
    fileName1 = [pathAnal fileN 'NoField'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end
