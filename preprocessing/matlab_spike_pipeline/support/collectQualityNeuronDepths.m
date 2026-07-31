function collectQualityNeuronDepths(path, fileName)
%     Select pyramidal neurons that meet criteria based on firing rate, amplitude,
%     and refractory period contamination    
%     Save the layercenterHdef, relDepthHdef, and absolute index of these Neurons into a file
%
% example call: 
%              collectQualityNeuronDepths('./','A022-20191107-01_DataStructure_mazeSection1_TrialType1')
% 

    %---- Define Criteria for High Quality Neurons ---%
    
    onlyRun = 1;
    minFR = 0.2;
    maxFR = 7;
    minRelAmp = 4;
    maxRefractViol = 1;
    maxCenterMax = -300;
    
    %---- Load Relevant Files and Variables ---%
    
    fullPathcluList = [path fileName '.mat'];
    fullPathFR = [path fileName '_FR_Run' num2str(onlyRun) '.mat'];
    fullPathDepth = [path fileName '_Depth' '.mat'];
    fullPathField = []
    fullPathInfo = [path fileName '_Info' '.mat'];
    
    if(exist(fullPathcluList) == 0)
        disp('cannot find TrialType1 file');
        return;
    end
    
    if(exist(fullPathFR) == 0)
        disp('cannot find FR_Run file');
        return;
    end
    
    if(exist(fullPathDepth) == 0)
        disp('cannot find Depth file');
        return;
    end
        
    load(fullPathcluList, 'cluList');
    load(fullPathFR, 'mFRStruct');
    load(fullPathDepth,'depthNeu');
    load(fullPathInfo,'autoCorr','fieldStructSess');
    
    indexFileName = strfind(fileName,'_');
    fileNameCCGT = [fileName(1:indexFileName(1)-1)...
                        '_BehavElectrDataLFP_CCG.mat'];
    fullPath = [path fileNameCCGT];
    if(exist(fullPath) == 0)
        disp('_BehavElectrDataLFP_CCG.mat file does not exist.');
        return;
    end
    load(fullPath,'ccgT','ccgVal'); 
    
    fileNameCCGT = [fileName(1:indexFileName(1)-1)...
                        '_BehavElectrDataLFP.mat'];
    fullPath = [path fileNameCCGT];
    if(exist(fullPath) == 0)
        disp('_BehavElectrDataLFP.mat file does not exist.');
        return;
    end
    load(fullPath,'Clu'); 
    
    GlobalConst;
    
    mFR = mFRStruct.mFR;
    relAmp = cluList.spatLocalRelAmpl;
    refViol = cluList.refracViolPercent;
    isInt = cluList.isIntern;
    centerMax = cluList.centerMax;
    leftMax = cluList.leftMax;
    isPyr = autoCorr.isPyrneuron;
    
    %--- Determine index of clusters that meet criteria ---%
    
    indGoodFR = mFR > minFR & mFR < maxFR;
    indGoodAmp = relAmp > minRelAmp;
    indGoodCenterMax = autoCorr.isSpikeHighAmp == 1;
    indGoodRef = refViol < maxRefractViol;
    
    indFinal = find((indGoodFR & indGoodRef & indGoodCenterMax & isPyr) == 1); 
    
    %--- Extract relevant information about depth and cluster id from
    %--- selected clusters
    
    layerCenterHDef = depthNeu.layerCenterHDef;
    relDepthNeuHDef = depthNeu.relDepthNeuHDef;
    shank = cluList.shank;
    localClu = cluList.localClu;
    
    QualityNeurons.clusterID = indFinal;
    QualityNeurons.shankID = shank(indFinal);
    QualityNeurons.localClu = localClu(indFinal);
    QualityNeurons.relDepthNeuHDef = relDepthNeuHDef(indFinal);
    QualityNeurons.layerCenterHDef = layerCenterHDef;
    QualityNeurons.peakTo20ms = autoCorr.peakTo20ms(indFinal);
    withField = zeros(1,length(isPyr));
    withField(fieldStructSess.indNeuron) = 1;
    QualityNeurons.withField = withField(indFinal);
    
    save([path fileName '_qualityNeuronDepths' '.mat'], 'QualityNeurons');
     
end
