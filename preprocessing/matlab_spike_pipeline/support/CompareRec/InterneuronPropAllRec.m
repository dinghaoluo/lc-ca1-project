function InterneuronPropAllRec()

    RecordingList;
    pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\interneuron\';
    
    minFRInt = 3;
    methodTheta = 1;
    onlyRun = 1;
    
    if(exist([pathAnal 'autoCorrIntAllRec.mat']))
        load([pathAnal 'autoCorrIntAllRec.mat']);
    end

%     %% interneurons in no cue passive task
%     disp('No cue')
%     autoCorrIntNoCue = accumInterneurons(listRecordingsNoCuePath,listRecordingsNoCueFileName,mazeSessionNoCue,minFRInt,1,methodTheta,onlyRun);
%     
%     %% interneurons in active licking task
%     disp('Active licking')
%     autoCorrIntAL = accumInterneurons(listRecordingsActiveLickPath,listRecordingsActiveLickFileName,mazeSessionActiveLick,minFRInt,2,methodTheta,onlyRun);
%     
%     %% interneurons in passive licking task with start cues
%     disp('Passive licking')
%     autoCorrIntPL = accumInterneurons(listRecordingsPassiveLickPath,listRecordingsPassiveLickFileName,mazeSessionPassiveLick,minFRInt,3,methodTheta,onlyRun);
%     
%     autoCorrIntAll.ccgVal = [autoCorrIntNoCue.ccgVal;autoCorrIntAL.ccgVal;autoCorrIntPL.ccgVal];
%     autoCorrIntAll.peakTo40ms = [autoCorrIntNoCue.peakTo40ms,autoCorrIntAL.peakTo40ms,autoCorrIntPL.peakTo40ms];
%     autoCorrIntAll.mean = [autoCorrIntNoCue.mean,autoCorrIntAL.mean,autoCorrIntPL.mean];
%     autoCorrIntAll.refract = [autoCorrIntNoCue.refract,autoCorrIntAL.refract,autoCorrIntPL.refract];
%     autoCorrIntAll.burstInd = [autoCorrIntNoCue.burstInd,autoCorrIntAL.burstInd,autoCorrIntPL.burstInd];
%     autoCorrIntAll.peakTime = [autoCorrIntNoCue.peakTime,autoCorrIntAL.peakTime,autoCorrIntPL.peakTime];
%     
%     autoCorrIntAll.phaseMeanDire = [autoCorrIntNoCue.phaseMeanDire,autoCorrIntAL.phaseMeanDire,autoCorrIntPL.phaseMeanDire]; 
%     autoCorrIntAll.thetaModHist = [autoCorrIntNoCue.thetaModHist,autoCorrIntAL.thetaModHist,autoCorrIntPL.thetaModHist]; 
%     autoCorrIntAll.maxPhaseArr = [autoCorrIntNoCue.maxPhaseArr,autoCorrIntAL.maxPhaseArr,autoCorrIntPL.maxPhaseArr]; 
%     autoCorrIntAll.minPhaseArr = [autoCorrIntNoCue.minPhaseArr,autoCorrIntAL.minPhaseArr,autoCorrIntPL.minPhaseArr]; 
%     autoCorrIntAll.maxPhaseOArr = [autoCorrIntNoCue.maxPhaseOArr,autoCorrIntAL.maxPhaseOArr,autoCorrIntPL.maxPhaseOArr]; 
%     autoCorrIntAll.minPhaseOArr = [autoCorrIntNoCue.minPhaseOArr,autoCorrIntAL.minPhaseOArr,autoCorrIntPL.minPhaseOArr]; 
%     autoCorrIntAll.histPhaseFil = [autoCorrIntNoCue.histPhaseFil;autoCorrIntAL.histPhaseFil;autoCorrIntPL.histPhaseFil];
%     autoCorrIntAll.histPhase = [autoCorrIntNoCue.histPhase;autoCorrIntAL.histPhase;autoCorrIntPL.histPhase];
%     
%     autoCorrIntAll.isSpikeHighAmp = [autoCorrIntNoCue.isSpikeHighAmp,autoCorrIntAL.isSpikeHighAmp,autoCorrIntPL.isSpikeHighAmp];
%     autoCorrIntAll.isSpikeHighAmp200 = [autoCorrIntNoCue.isSpikeHighAmp200,autoCorrIntAL.isSpikeHighAmp200,autoCorrIntPL.isSpikeHighAmp200];
%     autoCorrIntAll.relDepthNeuHDef = [autoCorrIntNoCue.relDepthNeuHDef,autoCorrIntAL.relDepthNeuHDef,autoCorrIntPL.relDepthNeuHDef];      
%     
%     autoCorrIntAll.task = [autoCorrIntNoCue.task,autoCorrIntAL.task,autoCorrIntPL.task];
%     autoCorrIntAll.indRec = [autoCorrIntNoCue.indRec,autoCorrIntAL.indRec,autoCorrIntPL.indRec];
%     autoCorrIntAll.indNeu = [autoCorrIntNoCue.indNeu,autoCorrIntAL.indNeu,autoCorrIntPL.indNeu];
%     
%     RecordingListTagging;
%     autoCorrIntTag = accumInterneuronsTag(listRecordingsTaggingPath,...
%         listRecordingsTaggingFileName,listRecordingsTaggingMazeSess,...
%         minFRInt,listRecordingsTaggingTask,...
%         listRecordingsTaggingIndRec,listRecordingsTaggingCellType,...
%         methodTheta,onlyRun);
    
    [autoCorrIntAll.idxC,autoCorrIntAll.klust] = kmeansInterneurons(autoCorrIntAll);
    
    [autoCorrIntAll.idxC1,autoCorrIntAll.klust1,...
        autoCorrIntAll.scorekm11,autoCorrIntAll.explainedkm11,...
        autoCorrIntAll.scorekm12,autoCorrIntAll.explainedkm12...
        ] = kmeansInterneurons1(autoCorrIntAll);
    
    [autoCorrIntAll.idxC2,autoCorrIntAll.klust2,...
        autoCorrIntAll.scorekm21,autoCorrIntAll.explainedkm21,...
        autoCorrIntAll.scorekm22,autoCorrIntAll.explainedkm22...
        ] = kmeansInterneurons2(autoCorrIntAll);

%     [autoCorrIntAll.idxC3,autoCorrIntAll.klust3,...
%         autoCorrIntAll.TSNE1,autoCorrIntAll.TSNE2...
%         ] = kmeansInterneurons3(autoCorrIntAll);
    
    autoCorrIntAll.relDepthNeuHDefC = [];
    autoCorrIntAll.relDepthNeuHDefMean = zeros(1,length(unique(autoCorrIntAll.idxC1)));
    for i = 1:max(autoCorrIntAll.idxC1)
        idxCurC = autoCorrIntAll.idxC1 == i;
        autoCorrIntAll.relDepthNeuHDefC{i} = autoCorrIntAll.relDepthNeuHDef(idxCurC);
        autoCorrIntAll.relDepthNeuHDefMean(i) = mean(autoCorrIntAll.relDepthNeuHDefC{i});
    end
        
    autoCorrIntTag.idxC = zeros(1,length(autoCorrIntTag.task));
    autoCorrIntTag.idxC1 = zeros(1,length(autoCorrIntTag.task));
    autoCorrIntTag.idxC2 = zeros(1,length(autoCorrIntTag.task));
%     autoCorrIntTag.idxC3 = zeros(1,length(autoCorrIntTag.task));
    autoCorrIntTag.scorekm11 = ...
        zeros(length(autoCorrIntTag.task),size(autoCorrIntAll.scorekm11,2));
    autoCorrIntTag.scorekm12 = ...
        zeros(length(autoCorrIntTag.task),size(autoCorrIntAll.scorekm12,2));
%     autoCorrIntTag.TSNE1 = ...
%         zeros(length(autoCorrIntTag.task),size(autoCorrIntAll.TSNE1,2));
%     autoCorrIntTag.TSNE2 = ...
%         zeros(length(autoCorrIntTag.task),size(autoCorrIntAll.TSNE2,2));
    for i = 1:length(autoCorrIntTag.task)
        indNeu = autoCorrIntAll.task == autoCorrIntTag.task(i) &...
            autoCorrIntAll.indRec == autoCorrIntTag.indRec(i) &...
            autoCorrIntAll.indNeu == autoCorrIntTag.indNeu(i);
        autoCorrIntTag.idxC(i) = autoCorrIntAll.idxC(indNeu);
        autoCorrIntTag.idxC1(i) = autoCorrIntAll.idxC1(indNeu);
        autoCorrIntTag.idxC2(i) = autoCorrIntAll.idxC2(indNeu);
%         autoCorrIntTag.idxC3(i) = autoCorrIntAll.idxC3(indNeu);
        autoCorrIntTag.scorekm11(i,:) = autoCorrIntAll.scorekm11(indNeu,:);
        autoCorrIntTag.scorekm12(i,:) = autoCorrIntAll.scorekm12(indNeu,:);
        autoCorrIntTag.scorekm21(i,:) = autoCorrIntAll.scorekm21(indNeu,:);
        autoCorrIntTag.scorekm22(i,:) = autoCorrIntAll.scorekm22(indNeu,:);
%         autoCorrIntTag.TSNE1(i,:) = autoCorrIntAll.TSNE1(indNeu,:);
%         autoCorrIntTag.TSNE2(i,:) = autoCorrIntAll.TSNE2(indNeu,:);
    end
    
    if(exist([pathAnal 'autoCorrIntAllRec.mat']))
        save([pathAnal 'autoCorrIntAllRec.mat'],'autoCorrIntNoCue','autoCorrIntAL',...
            'autoCorrIntPL','autoCorrIntTag','autoCorrIntAll','-append');
    else
        save([pathAnal 'autoCorrIntAllRec.mat'],'autoCorrIntNoCue','autoCorrIntAL',...
            'autoCorrIntPL','autoCorrIntTag','autoCorrIntAll');
    end
    
    plotCompCPCA(autoCorrIntAll.scorekm11(:,1),autoCorrIntAll.scorekm12(:,1),...
        autoCorrIntAll.scorekm12(:,2),'PC1','PC2','PC3',autoCorrIntAll.idxC1,...
        'PCA space');
    fileName1 = [pathAnal 'Interneuron_PCAClustersAng12'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
        
    plotCompCPCA(autoCorrIntAll.scorekm21(:,1),autoCorrIntAll.scorekm22(:,1),...
        autoCorrIntAll.scorekm22(:,2),'PC1','PC2','PC3',autoCorrIntAll.idxC2,...
        'PCA space');
    fileName1 = [pathAnal 'Interneuron_PCAClustersAng21'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    idxTagPV = autoCorrIntTag.cellType == 1;
    plotCompCPCATagC5(autoCorrIntAll.scorekm11(:,1),autoCorrIntAll.scorekm12(:,1),...
        autoCorrIntAll.scorekm12(:,2),...
        autoCorrIntTag.scorekm11(idxTagPV,1),autoCorrIntTag.scorekm12(idxTagPV,1),...
        autoCorrIntTag.scorekm12(idxTagPV,2),'PC1','PC2','PC3',autoCorrIntAll.idxC1,...
        'PCA space PV');
    fileName1 = [pathAnal 'Interneuron_PCAClustersScoreKm11TagPV'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    idxTagSST = autoCorrIntTag.cellType == 2;
    plotCompCPCATagC5(autoCorrIntAll.scorekm11(:,1),autoCorrIntAll.scorekm12(:,1),...
        autoCorrIntAll.scorekm12(:,2),...
        autoCorrIntTag.scorekm11(idxTagSST,1),autoCorrIntTag.scorekm12(idxTagSST,1),...
        autoCorrIntTag.scorekm12(idxTagSST,2),'PC1','PC2','PC3',autoCorrIntAll.idxC1,...
        'PCA space SST');
    fileName1 = [pathAnal 'Interneuron_PCAClustersScoreKm11TagSST'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompCPCATagC5(autoCorrIntAll.scorekm21(:,1),autoCorrIntAll.scorekm22(:,1),...
        autoCorrIntAll.scorekm22(:,2),...
        autoCorrIntTag.scorekm21(idxTagPV,1),autoCorrIntTag.scorekm22(idxTagPV,1),...
        autoCorrIntTag.scorekm22(idxTagPV,2),'PC1','PC2','PC3',autoCorrIntAll.idxC2,...
        'PCA space PV');
    fileName1 = [pathAnal 'Interneuron_PCAClustersScoreKm21TagPV'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompCPCATagC5(autoCorrIntAll.scorekm21(:,1),autoCorrIntAll.scorekm22(:,1),...
        autoCorrIntAll.scorekm22(:,2),...
        autoCorrIntTag.scorekm21(idxTagSST,1),autoCorrIntTag.scorekm22(idxTagSST,1),...
        autoCorrIntTag.scorekm22(idxTagSST,2),'PC1','PC2','PC3',autoCorrIntAll.idxC2,...
        'PCA space SST');
    fileName1 = [pathAnal 'Interneuron_PCAClustersScoreKm21TagSST'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotKmeanResult(autoCorrIntAll,autoCorrIntTag,idxTagPV,'All Interneurons','PV',pathAnal);
    plotKmeanResult(autoCorrIntAll,autoCorrIntTag,idxTagSST,'All Interneurons','SST',pathAnal);
    
%     plotCompCNoTagC5(autoCorrIntAll.peakTime,autoCorrIntAll.phaseMeanDire/pi*180,...
%         'peak time (ms)','phase mean dire',autoCorrIntAll.idxC1,'All Interneurons',[],[]);
%     fileName1 = [pathAnal 'Interneuron_PeakTVsPhaseMeanDireC5NoTag'];
%     saveas(gcf,fileName1);
%     print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompCNoTagC5(autoCorrIntAll.peakTime,autoCorrIntAll.phaseMeanDire/pi*180,...
        'peak time (ms)','phase mean dire',autoCorrIntAll.idxC2,'All Interneurons',[],[]);
    fileName1 = [pathAnal 'Interneuron_PeakTVsPhaseMeanDireNoTag'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
   
    cluPV = [3 4];
    plotCompCTagC5(autoCorrIntAll.peakTime,autoCorrIntAll.phaseMeanDire/pi*180,...
        autoCorrIntTag.peakTime(idxTagPV),autoCorrIntTag.phaseMeanDire(idxTagPV)/pi*180,...
        'peak time (ms)','phase mean dire',autoCorrIntAll.idxC2,'All Interneurons',[],[],cluPV);
    fileName1 = [pathAnal 'Interneuron_PeakTVsPhaseMeanDirePVTag'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    cluSST = [2 3];
    plotCompCTagC5(autoCorrIntAll.peakTime,autoCorrIntAll.phaseMeanDire/pi*180,...
        autoCorrIntTag.peakTime(idxTagSST),autoCorrIntTag.phaseMeanDire(idxTagSST)/pi*180,...
        'peak time (ms)','phase mean dire',autoCorrIntAll.idxC2,'All Interneurons',[],[],cluSST);
    fileName1 = [pathAnal 'Interneuron_PeakTVsPhaseMeanDireSSTTag'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    cluSST1 = [2];
    idxTagSST1 = autoCorrIntTag.cellType == 2 & autoCorrIntTag.idxC2 == cluSST1;
    plotCompCTagC5(autoCorrIntAll.peakTime(autoCorrIntAll.idxC2 == cluSST1),...
        autoCorrIntAll.phaseMeanDire(autoCorrIntAll.idxC2 == cluSST1)/pi*180,...
        autoCorrIntTag.peakTime(idxTagSST1),autoCorrIntTag.phaseMeanDire(idxTagSST1)/pi*180,...
        'peak time (ms)','phase mean dire',autoCorrIntAll.idxC2(autoCorrIntAll.idxC2 == cluSST1),...
        'SST O_LM neurons',[],[],cluSST1);
    fileName1 = [pathAnal 'Interneuron_PeakTVsPhaseMeanDireSSTTagClu2'];
    saveas(gcf,fileName1);
    
    cluPV1 = [4];
    idxTagPV1 = autoCorrIntTag.cellType == 1 & autoCorrIntTag.idxC2 == cluPV1;
    plotCompCTagC5(autoCorrIntAll.peakTime(autoCorrIntAll.idxC2 == cluPV1),...
        autoCorrIntAll.phaseMeanDire(autoCorrIntAll.idxC2 == cluPV1)/pi*180,...
        autoCorrIntTag.peakTime(idxTagPV1),autoCorrIntTag.phaseMeanDire(idxTagPV1)/pi*180,...
        'peak time (ms)','phase mean dire',autoCorrIntAll.idxC2(autoCorrIntAll.idxC2 == cluPV1),...
        'PV basket neurons',[],[],cluPV1);
    fileName1 = [pathAnal 'Interneuron_PeakTVsPhaseMeanDirePVTagClu4'];
    saveas(gcf,fileName1);
    
    plotCompCNoTagC5(autoCorrIntAll.burstInd,autoCorrIntAll.phaseMeanDire/pi*180,...
        'burst index','phase mean dire',autoCorrIntAll.idxC2,'All Interneurons',[],[]);
    fileName1 = [pathAnal 'Interneuron_BurstIndVsPhaseMeanDireC5NoTag'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompCTagC5(autoCorrIntAll.burstInd,autoCorrIntAll.phaseMeanDire/pi*180,...
        autoCorrIntTag.burstInd(idxTagPV),autoCorrIntTag.phaseMeanDire(idxTagPV)/pi*180,...
        'burst index','phase mean dire',autoCorrIntAll.idxC2,'All Interneurons',[],[],cluPV);
    fileName1 = [pathAnal 'Interneuron_burstIndVsPhaseMeanDireC5PVTag'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompCTagC5(autoCorrIntAll.burstInd,autoCorrIntAll.phaseMeanDire/pi*180,...
        autoCorrIntTag.burstInd(idxTagSST),autoCorrIntTag.phaseMeanDire(idxTagSST)/pi*180,...
        'burst index','phase mean dire',autoCorrIntAll.idxC2,'All Interneurons',[],[],cluSST);
    fileName1 = [pathAnal 'Interneuron_burstIndVsPhaseMeanDireC5SSTTag'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompCTagC5BG(autoCorrIntAll.burstInd,autoCorrIntAll.phaseMeanDire/pi*180,...
        autoCorrIntTag.burstInd(idxTagPV),autoCorrIntTag.phaseMeanDire(idxTagPV)/pi*180,...
        'burst index','phase mean dire',autoCorrIntAll.idxC2,'All Interneurons',[],[],cluPV(2));
    fileName1 = [pathAnal 'Interneuron_burstIndVsPhaseMeanDireC5PVTagBG'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompCTagC5BG(autoCorrIntAll.burstInd,autoCorrIntAll.phaseMeanDire/pi*180,...
        autoCorrIntTag.burstInd(idxTagSST),autoCorrIntTag.phaseMeanDire(idxTagSST)/pi*180,...
        'burst index','phase mean dire',autoCorrIntAll.idxC2,'All Interneurons',[],[],cluSST(1));
    fileName1 = [pathAnal 'Interneuron_burstIndVsPhaseMeanDireC5SSTTagBG'];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    %% plot figures
%     %% peak time vs. burst index
%     plotComp(autoCorrIntNoCue.peakTime,autoCorrIntNoCue.burstInd,...
%         'PL-NoCue peak time (ms)','PL-NoCue burst index');
%     
%     plotComp(autoCorrIntAL.peakTime,autoCorrIntAL.burstInd,...
%         'AL peak time (ms)','AL burst index');
%     
%     plotComp(autoCorrIntPL.peakTime,autoCorrIntPL.burstInd,...
%         'PL peak time (ms)','PL burst index');
%     
%     %% peak time vs. peak to 40 ms mean
%     plotComp(autoCorrIntNoCue.peakTime,autoCorrIntNoCue.peakTo40ms,...
%         'PL-NoCue peak time (ms)','PL-NoCue peak to 40-50 ms');
%     
%     plotComp(autoCorrIntAL.peakTime,autoCorrIntAL.peakTo40ms,...
%         'AL peak time (ms)','AL peak to 40-50 ms');
%     
%     plotComp(autoCorrIntPL.peakTime,autoCorrIntPL.peakTo40ms,...
%         'PL peak time (ms)','PL peak to 40-50 ms');
%     
%     %% peak time vs. peak to 40 ms mean
%     plotComp(autoCorrIntNoCue.peakTime,autoCorrIntNoCue.peakToMean,...
%         'PL-NoCue peak time (ms)','PL-NoCue peak to mean');
%     
%     plotComp(autoCorrIntAL.peakTime,autoCorrIntAL.peakToMean,...
%         'AL peak time (ms)','AL peak to mean');
%     
%     plotComp(autoCorrIntPL.peakTime,autoCorrIntPL.peakToMean,...
%         'PL peak time (ms)','PL peak to mean');
%     
%     %% refractory vs. peak to 40 ms mea
%     plotComp(autoCorrIntNoCue.refract,autoCorrIntNoCue.burstInd,...
%         'PL-NoCue refractory (ms)','PL-NoCue burst index');
%     
%     plotComp(autoCorrIntAL.refract,autoCorrIntAL.burstInd,...
%         'AL refractory (ms)','AL burst index');
%     
%     plotComp(autoCorrIntPL.refract,autoCorrIntPL.burstInd,...
%         'PL refractory (ms)','PL burst index');
    
    %% refractory vs. peak time
    plotComp(autoCorrIntAL.refract,autoCorrIntAL.peakTime,...
        'AL refractory (ms)','AL peak time (ms)');
end

function [idxInt,cInt] = gmInterneurons(autoCorr)
%     X = [autoCorr.peakTo40ms' autoCorr.refract' autoCorr.burstInd' autoCorr.peakTime' autoCorr.phaseMeanDire'];
    X = [autoCorr.peakTo40ms' autoCorr.refract' autoCorr.burstInd' autoCorr.peakTime' autoCorr.phaseMeanDire' autoCorr.thetaModHist']; % autoCorr.maxPhaseArr'];
    klist = [2:10];
    len = 500;
    options = statset('MaxIter',500);
    AIC = zeros(1,length(klist));
    GMModels = zeros(1,length(klist));
    for i = klist
        GMModels{i} = fitgmdist(X,i,'Options',options,'CovarianceType','diagonal');
        AIC(i) = GMModels{i}.AIC;
    end
    [minAIC,numC] = min(AIC);
    [~,optimalK] = max(binCount);
    optimalK = klist(optimalK);
    [idxInt,cInt,sumD] = kmeans(X,optimalK,'MaxIter',5000,'Display','final','Replicates',100);
    
end

function [idxInt,cInt] = kmeansInterneurons(autoCorr)
%     X = [autoCorr.peakTo40ms' autoCorr.refract' autoCorr.burstInd' autoCorr.peakTime' autoCorr.phaseMeanDire'];
    X = [normData1(autoCorr.peakTo40ms') ...
        normData1(autoCorr.refract') ...
        normData1(autoCorr.burstInd') ...
        1.5*normData1(autoCorr.peakTime') ... 1.5
        1.5*normData1(autoCorr.phaseMeanDire') ...
        normData1(autoCorr.thetaModHist') ...
        normData1(autoCorr.maxPhaseArr') ... 2
        normData1(autoCorr.minPhaseArr') ...
        normData1(autoCorr.relDepthNeuHDef')*2];
    klist = [6:10];
    len = 10;
    evaC = [];
    kn = zeros(len,1);
    for i = 1:len
        evaC{i} = evalclusters(X,'kmeans','CalinskiHarabasz','kList',klist);
        kn(i) = evaC{i}.OptimalK;
    end
    binCount = histc(kn,klist);
    [~,optimalK] = max(binCount);
    optimalK = klist(optimalK);
    [idxInt,cInt,sumD] = kmeans(X,optimalK,'MaxIter',5000,'Display','final','Replicates',100);
    
end

function [idxInt,cInt,score1,explained1,score2,explained2] = kmeansInterneurons1(autoCorr)
%     X = [autoCorr.peakTo40ms' autoCorr.refract' autoCorr.burstInd' autoCorr.peakTime' autoCorr.phaseMeanDire'];
    std1 = 1.5;
    paramT.gaussFilt = gaussFilter(12*std1, std1);
    lenGaussKernel = length(paramT.gaussFilt);
    normFactor = sum(paramT.gaussFilt);
    paramT.gaussFilt = paramT.gaussFilt./normFactor;
    
    [r,c] = size(autoCorr.ccgVal);
    X = zeros(r,c);
    for i = 1:r
        autocor = conv(autoCorr.ccgVal(i,:),paramT.gaussFilt);
        if(mod(lenGaussKernel,2) == 0)
            autocor = ...
                autocor(floor(lenGaussKernel/2)+1:...
                    (end-floor(lenGaussKernel/2))+1); 
        else
            autocor = ...
                autocor(floor(lenGaussKernel/2)+1:...
                    (end-floor(lenGaussKernel/2))); 
        end
        X(i,:) = normData(autocor);        
    end
    
    [coeff,score1,latent,ts,explained1] = pca(X);
    numComp = find(cumsum(explained1) > 99,1); %98
    scoreNorm = score1(:,1:numComp);
    scoreNorm = scoreNorm/max(abs(scoreNorm(:,1))); % normalization based on PC1
%     for i = 1:numComp
%         scoreNorm(:,i) = normData(scoreNorm(:,i));
%     end
    [r,c] = size(autoCorr.histPhaseFil);
    X1 = zeros(r,c);
    for i = 1:r
        X1(i,:) = normData(autoCorr.histPhaseFil(i,:));
    end
    [~,score2,latent2,ts2,explained2] = pca(X1);
    numComp2 = find(cumsum(explained2) > 99,1); %98
    scoreNorm2 = score2(:,1:numComp2);
    scoreNorm2 = scoreNorm2/max(abs(scoreNorm2(:,1))); % normalization based on PC1
%     for i = 1:numComp2
%         scoreNorm2(:,i) = normData(scoreNorm2(:,i));
%     end
%     XAll = [scoreNorm(:,1)*1.3 scoreNorm(:,2:end) scoreNorm2]; %1.5

    depth = autoCorr.relDepthNeuHDef/max(abs(autoCorr.relDepthNeuHDef));
    
    XAll = [scoreNorm scoreNorm2 depth'*2];
    klist = [6:10];
    len = 10;
    evaC = [];
    kn = zeros(len,1);
    for i = 1:len
        evaC{i} = evalclusters(XAll,'kmeans','CalinskiHarabasz','kList',klist);
        kn(i) = evaC{i}.OptimalK;
    end
    binCount = histc(kn,klist);
    [~,optimalK] = max(binCount);
    optimalK = klist(optimalK);
    [idxInt,cInt,sumD] = kmeans(XAll,optimalK,'MaxIter',5000,'Display','final','Replicates',100);
    
end

function [idxInt,cInt,score,explained,score2,explained2] = kmeansInterneurons2(autoCorr)
    
    [r,c] = size(autoCorr.ccgVal);
    X = zeros(r,c);
    for i = 1:r
        X(i,:) = normData(autoCorr.ccgVal(i,:));
    end
    [coeff,score,latent,ts,explained] = pca(X);
    numComp = find(cumsum(explained) > 95,1);
    if(numComp > 30)
        numComp = 30;
    end
    scoreNorm = score(:,1:numComp);
    scoreNorm = scoreNorm/max(abs(scoreNorm(:,1))); % normalization based on PC1
    
    [r,c] = size(autoCorr.histPhase);
    X1 = zeros(r,c);
    for i = 1:r
        X1(i,:) = normData(autoCorr.histPhase(i,:));
    end
    [coeff2,score2,latent2,ts2,explained2] = pca(X1);
    numComp2 = find(cumsum(explained2) > 95,1);
    if(numComp2 > 30)
        numComp2 = 30;
    end
    scoreNorm2 = score2(:,1:numComp2);
    scoreNorm2 = scoreNorm2/max(abs(scoreNorm2(:,1))); % normalization based on PC1
  
    depth = autoCorr.relDepthNeuHDef/max(abs(autoCorr.relDepthNeuHDef));
%     peakTime = autoCorr.peakTime/max(abs(autoCorr.peakTime));

    XAll = [scoreNorm scoreNorm2 depth'*2]; %2
    klist = [6:10];
    len = 10;
    evaC = [];
    kn = zeros(len,1);
    for i = 1:len
        evaC{i} = evalclusters(XAll,'kmeans','CalinskiHarabasz','kList',klist);
        kn(i) = evaC{i}.OptimalK;
    end
    binCount = histc(kn,klist);
    [~,optimalK] = max(binCount);
    optimalK = klist(optimalK);
    [idxInt,cInt,sumD] = kmeans(XAll,optimalK,'MaxIter',5000,'Display','final','Replicates',100);
    
end

% function [idxInt,cInt,Y,Y1] = kmeansInterneurons3(autoCorr)
%     X0 = [autoCorr.peakTo40ms'/max(abs(autoCorr.peakTo40ms))...
%         autoCorr.refract'/max(abs(autoCorr.refract))...
%         autoCorr.burstInd'/max(abs(autoCorr.burstInd))...
%         autoCorr.peakTime'/max(abs(autoCorr.peakTime))...
%         autoCorr.phaseMeanDire'/max(abs(autoCorr.phaseMeanDire))...
%         autoCorr.relDepthNeuHDef'/max(abs(autoCorr.relDepthNeuHDef))];
%     
%     [r,c] = size(autoCorr.ccgVal);
%     X = zeros(r,c);
%     for i = 1:r
%         X(i,:) = normData(autoCorr.ccgVal(i,:));
%     end
%     Y = tsne(X,'Algorithm','exact','Distance','cosine','NumDimensions',3);
%     for i = 1:size(Y,2)
%         Y(:,i) = Y(:,i)/max(abs(Y(:,i))); % normalization 
%     end
% %     Y = tsne(autoCorr.ccgVal,'Algorithm','exact','Distance','cosine','NumDimensions',3);
% %     plot3(Y(:,1),Y(:,2),Y(:,3),'.')
% %     gscatter(Y(:,1),Y(:,2),ones(r,1));
% 
%     [r,c] = size(autoCorr.histPhase);
%     X1 = zeros(r,c);
%     for i = 1:r
%         X1(i,:) = normData(autoCorr.histPhase(i,:));
%     end
%     Y1 = tsne(X1,'Algorithm','exact','Distance','cosine','NumDimensions',3);
%     for i = 1:size(Y1,2)
%         Y1(:,i) = Y1(:,i)/max(abs(Y1(:,i))); % normalization 
%     end
% %     plot3(Y1(:,1),Y1(:,2),Y1(:,3),'.')
% %     gscatter(Y1(:,1),Y1(:,2),ones(r,1));
%     
%     YAll = [Y Y1 X0(:,end)*2];
%     klist = [6:10];
%     len = 10;
%     evaC = [];
%     kn = zeros(len,1);
%     for i = 1:len
%         evaC{i} = evalclusters(YAll,'kmeans','CalinskiHarabasz','kList',klist);
%         kn(i) = evaC{i}.OptimalK;
%     end
%     binCount = histc(kn,klist);
%     [~,optimalK] = max(binCount);
%     optimalK = klist(optimalK);
%     [idxInt,cInt,sumD] = kmeans(YAll,optimalK,'MaxIter',5000,'Display','final','Replicates',100);
%     
% end

% function [idxInt,cInt,score,explained,score2,explained2] = kmeansInterneurons2(autoCorr)
% 
%     std1 = 1.5;
%     paramT.gaussFilt = gaussFilter(12*std1, std1);
%     lenGaussKernel = length(paramT.gaussFilt);
%     normFactor = sum(paramT.gaussFilt);
%     paramT.gaussFilt = paramT.gaussFilt./normFactor;
%     
%     [r,c] = size(autoCorr.ccgVal);
%     X = zeros(r,c);
%     for i = 1:r
%         autocor = conv(autoCorr.ccgVal(i,:),paramT.gaussFilt);
%         if(mod(lenGaussKernel,2) == 0)
%             autocor = ...
%                 autocor(floor(lenGaussKernel/2)+1:...
%                     (end-floor(lenGaussKernel/2))+1); 
%         else
%             autocor = ...
%                 autocor(floor(lenGaussKernel/2)+1:...
%                     (end-floor(lenGaussKernel/2))); 
%         end     
%         X(i,:) = normData(autocor);     
%     end
%     
%     [coeff,score,latent,ts,explained] = pca(X);
%     numComp = find(cumsum(explained) > 99,1);
%     scoreNorm = score(:,1:numComp);
%     scoreNorm = scoreNorm/std(scoreNorm(:,1)); % normalization based on PC1
% %     for i = 1:numComp
% %         scoreNorm(:,i) = normData1(scoreNorm(:,i));
% %     end
% 
%     [r,c] = size(autoCorr.histPhaseFil);
%     X1 = zeros(r,c);
%     for i = 1:r
%         X1(i,:) = normData(autoCorr.histPhaseFil(i,:));
%     end
%     [coeff2,score2,latent2,ts2,explained2] = pca(X1);
%     numComp2 = find(cumsum(explained2) > 99,1);
%     scoreNorm2 = score2(:,1:numComp2);
%     scoreNorm2 = scoreNorm2/std(scoreNorm2(:,1)); % normalization based on PC1
% %     for i = 1:numComp2
% %         scoreNorm2(:,i) = normData1(scoreNorm2(:,i));
% %     end
%     XAll = [scoreNorm scoreNorm2]; %2
%     klist = [6:10];
%     len = 100;
%     evaC = [];
%     kn = zeros(len,1);
%     for i = 1:len
%         evaC{i} = evalclusters(XAll,'kmeans','CalinskiHarabasz','kList',klist);
%         kn(i) = evaC{i}.OptimalK;
%     end
%     binCount = histc(kn,klist);
%     [~,optimalK] = max(binCount);
%     optimalK = klist(optimalK);
%     [idxInt,cInt,sumD] = kmeans(XAll,optimalK,'MaxIter',5000,'Display','final','Replicates',100);
%     
% end

function result = normData(data)
    result = (data-min(data))/(max(data)-min(data));
end

function result = normData1(data)
    result = data/max(data);
end

function plotKmeanResult(autoCorr,autoCorrTag,idxTag,titleX,Tagtype,pathAnal)
%     plotCompC(autoCorr.refract,autoCorr.burstInd,...
%         autoCorrTag.refract(idxTag),autoCorrTag.burstInd(idxTag),...
%         'refractory (ms)','burst index',autoCorr.idxC,titleX);  
%     fileName1 = [pathAnal 'Interneuron_RefrVsBurst' Tagtype];
%     saveas(gcf,fileName1);
%     print('-painters', '-dpdf', fileName1, '-r600')
%     
%     plotCompC(autoCorr.peakTime,autoCorr.peakTo40ms,...
%         autoCorrTag.peakTime(idxTag),autoCorrTag.peakTo40ms(idxTag),...
%         'peak time (ms)','peak to 40-50 ms',autoCorr.idxC,titleX);
%     fileName1 = [pathAnal 'Interneuron_PeakTVsPeakTo40ms' Tagtype];
%     saveas(gcf,fileName1);
%     print('-painters', '-dpdf', fileName1, '-r600')
%     
%     plotCompC(autoCorr.peakTime,autoCorr.phaseMeanDire,...
%         autoCorrTag.peakTime(idxTag),autoCorrTag.phaseMeanDire(idxTag),...
%         'peak time (ms)','phase mean dire',autoCorr.idxC,titleX);
%     fileName1 = [pathAnal 'Interneuron_PeakTVsPhaseMeanDire' Tagtype];
%     saveas(gcf,fileName1);
%     print('-painters', '-dpdf', fileName1, '-r600')
%     
%     plotCompC(autoCorr.peakTime,autoCorr.minPhaseArr,...
%         autoCorrTag.peakTime(idxTag),autoCorrTag.minPhaseArr(idxTag),...
%         'peak time (ms)','min phase theta hist',autoCorr.idxC,titleX);
%     fileName1 = [pathAnal 'Interneuron_PeakTVsMinPhaseArr' Tagtype];
%     saveas(gcf,fileName1);
%     print('-painters', '-dpdf', fileName1, '-r600')
%     
%     plotCompC(autoCorr.peakTime,autoCorr.maxPhaseArr,...
%         autoCorrTag.peakTime(idxTag),autoCorrTag.maxPhaseArr(idxTag),...
%         'peak time (ms)','max phase theta hist',autoCorr.idxC,titleX);
%     fileName1 = [pathAnal 'Interneuron_PeakTVsMaxPhaseArr' Tagtype];
%     saveas(gcf,fileName1);
%     print('-painters', '-dpdf', fileName1, '-r600')
%     
%     isHighAmp = autoCorr.isSpikeHighAmp == 1;
%     isHighAmpT = autoCorrTag.isSpikeHighAmp == 1 & idxTag;
%     plotCompC(autoCorr.peakTime(isHighAmp),autoCorr.relDepthNeuHDef(isHighAmp),...
%         autoCorrTag.peakTime(isHighAmpT),autoCorrTag.relDepthNeuHDef(isHighAmpT),...
%         'peak time (ms)','depth',autoCorr.idxC(isHighAmp),titleX); 
%     fileName1 = [pathAnal 'Interneuron_PeakTVsDepth' Tagtype];
%     saveas(gcf,fileName1);
%     print('-painters', '-dpdf', fileName1, '-r600')
%     
%     plotCompC(autoCorr.peakTime,autoCorr.relDepthNeuHDef,...
%         autoCorrTag.peakTime(idxTag),autoCorrTag.relDepthNeuHDef(idxTag),...
%         'peak time (ms)','depth',autoCorr.idxC,titleX); 
%     fileName1 = [pathAnal 'Interneuron_PeakTVsDepthAllInt' Tagtype];
%     saveas(gcf,fileName1);
%     print('-painters', '-dpdf', fileName1, '-r600')
    
    
    %% klustering method 1
    plotCompC(autoCorr.refract,autoCorr.burstInd,...
        autoCorrTag.refract(idxTag),autoCorrTag.burstInd(idxTag),...
        'refractory (ms)','burst index',autoCorr.idxC1,titleX);  
    fileName1 = [pathAnal 'Interneuron_RefrVsBurst1' Tagtype];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompC(autoCorr.peakTime,autoCorr.peakTo40ms,...
        autoCorrTag.peakTime(idxTag),autoCorrTag.peakTo40ms(idxTag),...
        'peak time (ms)','peak to 40-50 ms',autoCorr.idxC1,titleX);
    fileName1 = [pathAnal 'Interneuron_PeakTVsPeakTo40ms1' Tagtype];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompC(autoCorr.peakTime,autoCorr.phaseMeanDire,...
        autoCorrTag.peakTime(idxTag),autoCorrTag.phaseMeanDire(idxTag),...
        'Peak time (ms)','Mean theta phase',autoCorr.idxC1,titleX,...
        [0 max(autoCorr.peakTime)],[0 2*pi]);
    fileName1 = [pathAnal 'Interneuron_PeakTVsPhaseMeanDire1' Tagtype];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompC(autoCorr.peakTime,autoCorr.minPhaseArr,...
        autoCorrTag.peakTime(idxTag),autoCorrTag.minPhaseArr(idxTag),...
        'peak time (ms)','min phase theta hist',autoCorr.idxC1,titleX);
    fileName1 = [pathAnal 'Interneuron_PeakTVsMinPhaseArr1' Tagtype];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompC(autoCorr.peakTime,autoCorr.maxPhaseArr,...
        autoCorrTag.peakTime(idxTag),autoCorrTag.maxPhaseArr(idxTag),...
        'peak time (ms)','max phase theta hist',autoCorr.idxC1,titleX);
    fileName1 = [pathAnal 'Interneuron_PeakTVsMaxPhaseArr1' Tagtype];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    isHighAmp = autoCorr.isSpikeHighAmp == 1;
    isHighAmpT = autoCorrTag.isSpikeHighAmp == 1 & idxTag;
    plotCompC(autoCorr.peakTime(isHighAmp),autoCorr.relDepthNeuHDef(isHighAmp),...
        autoCorrTag.peakTime(isHighAmpT),autoCorrTag.relDepthNeuHDef(isHighAmpT),...
        'peak time (ms)','depth',autoCorr.idxC1(isHighAmp),titleX); 
    fileName1 = [pathAnal 'Interneuron_PeakTVsDepth1' Tagtype];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompC(autoCorr.peakTime,autoCorr.relDepthNeuHDef,...
        autoCorrTag.peakTime(idxTag),autoCorrTag.relDepthNeuHDef(idxTag),...
        'peak time (ms)','depth',autoCorr.idxC1,titleX); 
    fileName1 = [pathAnal 'Interneuron_PeakTVsDepthAllInt1' Tagtype];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    %% klustering method 2
    plotCompC(autoCorr.refract,autoCorr.burstInd,...
        autoCorrTag.refract(idxTag),autoCorrTag.burstInd(idxTag),...
        'refractory (ms)','burst index',autoCorr.idxC2,titleX);  
    fileName1 = [pathAnal 'Interneuron_RefrVsBurst2' Tagtype];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompC(autoCorr.peakTime,autoCorr.peakTo40ms,...
        autoCorrTag.peakTime(idxTag),autoCorrTag.peakTo40ms(idxTag),...
        'peak time (ms)','peak to 40-50 ms',autoCorr.idxC2,titleX);
    fileName1 = [pathAnal 'Interneuron_PeakTVsPeakTo40ms2' Tagtype];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompC(autoCorr.peakTime,autoCorr.phaseMeanDire,...
        autoCorrTag.peakTime(idxTag),autoCorrTag.phaseMeanDire(idxTag),...
        'peak time (ms)','phase mean dire',autoCorr.idxC2,titleX);
    fileName1 = [pathAnal 'Interneuron_PeakTVsPhaseMeanDire2' Tagtype];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompC(autoCorr.peakTime,autoCorr.minPhaseArr,...
        autoCorrTag.peakTime(idxTag),autoCorrTag.minPhaseArr(idxTag),...
        'peak time (ms)','min phase theta hist',autoCorr.idxC2,titleX);
    fileName1 = [pathAnal 'Interneuron_PeakTVsMinPhaseArr2' Tagtype];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompC(autoCorr.peakTime,autoCorr.maxPhaseArr,...
        autoCorrTag.peakTime(idxTag),autoCorrTag.maxPhaseArr(idxTag),...
        'peak time (ms)','max phase theta hist',autoCorr.idxC2,titleX);
    fileName1 = [pathAnal 'Interneuron_PeakTVsMaxPhaseArr2' Tagtype];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    isHighAmp = autoCorr.isSpikeHighAmp == 1;
    isHighAmpT = autoCorrTag.isSpikeHighAmp == 1 & idxTag;
    plotCompC(autoCorr.peakTime(isHighAmp),autoCorr.relDepthNeuHDef(isHighAmp),...
        autoCorrTag.peakTime(isHighAmpT),autoCorrTag.relDepthNeuHDef(isHighAmpT),...
        'peak time (ms)','depth',autoCorr.idxC2(isHighAmp),titleX); 
    fileName1 = [pathAnal 'Interneuron_PeakTVsDepth2' Tagtype];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
    plotCompC(autoCorr.peakTime,autoCorr.relDepthNeuHDef,...
        autoCorrTag.peakTime(idxTag),autoCorrTag.relDepthNeuHDef(idxTag),...
        'peak time (ms)','depth',autoCorr.idxC2,titleX); 
    fileName1 = [pathAnal 'Interneuron_PeakTVsDepthAllInt2' Tagtype];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
    
end

function plotComp(x,y,xl,yl)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    h = plot(x,y,'o');
    set(h,'MarkerSize',6,'Color',[0.5 0.5 0.9]);
    maxX = max(x);
    maxY = max(y);
    minX = min(x);
    minY = min(y);
    set(gca,'XLim',[minX maxX],'YLim',[minY maxY]);
    xlabel(xl)
    ylabel(yl)
end

function plotCompC(x,y,xt,yt,xl,yl,idx,ti,xlimit,ylimit)
    if(nargin == 8)
        xlimit = [];
        ylimit = [];
    end
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    colorArr = [51/255 102/255 255/255;...      
                218/255 179/255 255/255;...                
                255/255 100/255 100/255;...   
                0.9 0.9 0.3;...
                0.3 0.9 0.3;...
                59/255 113/255 86/255;...           
                0.2 0.8 0.5;...
                0.8 0.5 0.2;...
                0.3 0.7 0.3];
    hold on;
    for i = 1:max(idx)
        indTmp = idx == i;
        disp(['Cluster' num2str(i) ' has ' num2str(sum(indTmp)) ' components']);
        h = plot(x(indTmp),y(indTmp),'.');
        set(h,'MarkerSize',11,'Color',colorArr(mod(i,max(idx))+1,:));
    end
    h = plot(xt,yt,'k+');
    set(h,'MarkerSize',9);
    if(isempty(xlimit))
        maxX = max(x);
        minX = min(x);        
    else
        maxX = xlimit(2);
        minX = xlimit(1);
    end
    if(isempty(ylimit))
        maxY = max(y);
        minY = min(y);        
    else
        maxY = ylimit(2);
        minY = ylimit(1);
    end
    set(gca,'XLim',[minX maxX],'YLim',[minY maxY]);
    xlabel(xl)
    ylabel(yl)
    title(ti);
end

function plotCompCNoTagC5(x,y,xl,yl,idx,ti,xlimit,ylimit)
    if(nargin == 8)
        xlimit = [];
        ylimit = [];
    end
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    colorArr = [51/255 102/255 255/255;...      
                218/255 179/255 255/255;...                
                255/255 100/255 100/255;...   
                0.9 0.9 0.3;...
                0.3 0.9 0.3;...
                59/255 113/255 86/255;...           
                0.2 0.8 0.5;...
                0.8 0.5 0.2;...
                0.3 0.7 0.3];
    hold on;
    for i = 1:max(idx)
        indTmp = idx == i;
        disp(['Cluster' num2str(i) ' has ' num2str(sum(indTmp)) ' components']);
        h = plot(x(indTmp),y(indTmp),'.');
        set(h,'MarkerSize',11,'Color',colorArr(mod(i,max(idx))+1,:));
    end
%     indTmp = idx == 4;
%     h = plot(x(indTmp),y(indTmp),'.');
%     set(h,'MarkerSize',11,'Color',colorArr(5,:));
    if(isempty(xlimit))
        maxX = max(x);
        minX = min(x);        
    else
        maxX = xlimit(2);
        minX = xlimit(1);
    end
    if(isempty(ylimit))
        maxY = max(y);
        minY = min(y);        
    else
        maxY = ylimit(2);
        minY = ylimit(1);
    end
    set(gca,'XLim',[minX maxX],'YLim',[minY maxY]);
    xlabel(xl)
    ylabel(yl)
    title(ti);
end

function plotCompCTagC5(x,y,xt,yt,xl,yl,idx,ti,xlimit,ylimit,clu)
    if(nargin == 8)
        xlimit = [];
        ylimit = [];
    end
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    colorArr = [51/255 102/255 255/255;...      
                218/255 179/255 255/255;...                
                255/255 100/255 100/255;...   
                0.9 0.9 0.3;...
                0.3 0.9 0.3;...
                59/255 113/255 86/255;...           
                0.2 0.8 0.5;...
                0.8 0.5 0.2;...
                0.3 0.7 0.3];
    hold on;
    for i = 1:max(idx)
        indTmp = idx == i;
        disp(['Cluster' num2str(i) ' has ' num2str(sum(indTmp)) ' components']);
        h = plot(x(indTmp),y(indTmp),'.');
        set(h,'MarkerSize',11,'Color',colorArr(mod(i,max(idx))+1,:));
    end
    
    indTmp = idx == clu(1);
    h = plot(x(indTmp),y(indTmp),'.');
    set(h,'MarkerSize',11,'Color',colorArr(mod(clu(1),max(idx))+1,:));
    if(length(clu) > 1)
        indTmp = idx == clu(2);
        h = plot(x(indTmp),y(indTmp),'.');
        set(h,'MarkerSize',11,'Color',colorArr(mod(clu(2),max(idx))+1,:));
    end
    h = plot(xt,yt,'k+');
    set(h,'MarkerSize',9);
    if(isempty(xlimit))
        maxX = max(x);
        minX = min(x);        
    else
        maxX = xlimit(2);
        minX = xlimit(1);
    end
    if(isempty(ylimit))
        maxY = max(y);
        minY = min(y);        
    else
        maxY = ylimit(2);
        minY = ylimit(1);
    end
    set(gca,'XLim',[minX maxX],'YLim',[minY maxY]);
    xlabel(xl)
    ylabel(yl)
    title(ti);
end

function plotCompCTagC5BG(x,y,xt,yt,xl,yl,idx,ti,xlimit,ylimit,clu)
    if(nargin == 8)
        xlimit = [];
        ylimit = [];
    end
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    colorArr = [51/255 102/255 255/255;...      
                218/255 179/255 255/255;...                
                255/255 100/255 100/255;...   
                0.9 0.9 0.3;...
                0.3 0.9 0.3;...
                59/255 113/255 86/255;...           
                0.2 0.8 0.5;...
                0.8 0.5 0.2;...
                0.3 0.7 0.3];
    hold on;
    for i = 1:max(idx)
        indTmp = idx == i;
        disp(['Cluster' num2str(i) ' has ' num2str(sum(indTmp)) ' components']);
        h = plot(x(indTmp),y(indTmp),'.');
        set(h,'MarkerSize',11,'Color',[0.7 0.7 0.7]);
    end
    indTmp = idx == clu;
    h = plot(x(indTmp),y(indTmp),'.');
    set(h,'MarkerSize',15,'Color',[0 0 1]);
    h = plot(xt,yt,'+');
    set(h,'MarkerSize',9,'Color',[230 127 128]/255);
    if(isempty(xlimit))
        maxX = max(x);
        minX = min(x);        
    else
        maxX = xlimit(2);
        minX = xlimit(1);
    end
    if(isempty(ylimit))
        maxY = max(y);
        minY = min(y);        
    else
        maxY = ylimit(2);
        minY = ylimit(1);
    end
    set(gca,'XLim',[minX maxX],'YLim',[minY maxY]);
    xlabel(xl)
    ylabel(yl)
    title(ti);
end

function plotCompCPCA(x,y,z,xl,yl,zl,idx,ti,xlimit,ylimit,zlimit)
    if(nargin == 8)
        xlimit = [];
        ylimit = [];
        zlimit = [];
    end
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    colorArr = [51/255 102/255 255/255;...      
                218/255 179/255 255/255;...                
                255/255 100/255 100/255;...   
                0.9 0.9 0.3;...
                0.3 0.9 0.3;...
                59/255 113/255 86/255;...           
                0.2 0.8 0.5;...
                0.8 0.5 0.2;...
                0.3 0.7 0.3];
    
    for i = 1:max(idx)
        indTmp = idx == i;
        disp(['Cluster' num2str(i) ' has ' num2str(sum(indTmp)) ' components']);
        h = plot3(x(indTmp),y(indTmp),z(indTmp),'.');
        set(h,'MarkerSize',11,'Color',colorArr(mod(i,max(idx))+1,:));
        if(i == 1)
            hold on;
        end
    end
    if(isempty(xlimit))
        maxX = max(x);
        minX = min(x);        
    else
        maxX = xlimit(2);
        minX = xlimit(1);
    end
    if(isempty(ylimit))
        maxY = max(y);
        minY = min(y);        
    else
        maxY = ylimit(2);
        minY = ylimit(1);
    end
    if(isempty(zlimit))
        maxZ = max(z);
        minZ = min(z);        
    else
        maxZ = zlimit(2);
        minZ = zlimit(1);
    end
    set(gca,'XLim',[minX maxX],'YLim',[minY maxY]);
    xlabel(xl)
    ylabel(yl)
    zlabel(zl)
    title(ti);
end

function plotCompCPCATagC5(x,y,z,xt,yt,zt,xl,yl,zl,idx,ti,xlimit,ylimit,zlimit)
    if(nargin == 11)
        xlimit = [];
        ylimit = [];
        zlimit = [];
    end
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    colorArr = [51/255 102/255 255/255;...      
                218/255 179/255 255/255;...                
                255/255 100/255 100/255;...   
                0.9 0.9 0.3;...
                0.3 0.9 0.3;...
                59/255 113/255 86/255;...           
                0.2 0.8 0.5;...
                0.8 0.5 0.2;...
                0.3 0.7 0.3];
    
    for i = 1:max(idx)
        indTmp = idx == i;
        disp(['Cluster' num2str(i) ' has ' num2str(sum(indTmp)) ' components']);
        h = plot3(x(indTmp),y(indTmp),z(indTmp),'.');
        set(h,'MarkerSize',11,'Color',colorArr(mod(i,max(idx))+1,:));
        if(i == 1)
            hold on;
        end
    end
%     indTmp = idx == 5;
%     h = plot3(x(indTmp),y(indTmp),z(indTmp),'.');
%     set(h,'MarkerSize',13,'Color',[0 0 1]);
    h = plot3(xt,yt,zt,'k+');
    set(h,'MarkerSize',11);
    
    if(isempty(xlimit))
        maxX = max(x);
        minX = min(x);        
    else
        maxX = xlimit(2);
        minX = xlimit(1);
    end
    if(isempty(ylimit))
        maxY = max(y);
        minY = min(y);        
    else
        maxY = ylimit(2);
        minY = ylimit(1);
    end
    if(isempty(zlimit))
        maxZ = max(z);
        minZ = min(z);        
    else
        maxZ = zlimit(2);
        minZ = zlimit(1);
    end
    set(gca,'XLim',[minX maxX],'YLim',[minY maxY]);
    xlabel(xl)
    ylabel(yl)
    zlabel(zl)
    title(ti);
end

function plotCompCPCATag(x,y,z,xt,yt,zt,xl,yl,zl,idx,ti,xlimit,ylimit,zlimit)
    if(nargin == 8)
        xlimit = [];
        ylimit = [];
        zlimit = [];
    end
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 400])
    colorArr = [51/255 102/255 255/255;...      
                218/255 179/255 255/255;...                
                255/255 100/255 100/255;...   
                0.9 0.9 0.3;...
                0.3 0.9 0.3;...
                59/255 113/255 86/255;...           
                0.2 0.8 0.5;...
                0.8 0.5 0.2;...
                0.3 0.7 0.3];
    
    for i = 1:max(idx)
        indTmp = idx == i;
        disp(['Cluster' num2str(i) ' has ' num2str(sum(indTmp)) ' components']);
        h = plot3(x(indTmp),y(indTmp),z(indTmp),'.');
        set(h,'MarkerSize',11,'Color',colorArr(mod(i,max(idx))+1,:));
        if(i == 1)
            hold on;
        end
    end
    h = plot3(xt,yt,zt,'k+');
    set(h,'MarkerSize',9);
    if(isempty(xlimit))
        maxX = max(x);
        minX = min(x);        
    else
        maxX = xlimit(2);
        minX = xlimit(1);
    end
    if(isempty(ylimit))
        maxY = max(y);
        minY = min(y);        
    else
        maxY = ylimit(2);
        minY = ylimit(1);
    end
    if(isempty(zlimit))
        maxZ = max(z);
        minZ = min(z);        
    else
        maxZ = zlimit(2);
        minZ = zlimit(1);
    end
    set(gca,'XLim',[minX maxX],'YLim',[minY maxY]);
    xlabel(xl)
    ylabel(yl)
    zlabel(zl)
    title(ti);
end
