function [indCluGood,indCluBad] = HandSelectedCluFromBehKMean(fileName)
% hand selected best and worst cluster based on the kmean classification of behavior

    if(~isempty(strfind(fileName,'A011-20190218')))
        indCluGood = 2; 
        indCluBad = 0;
    elseif(~isempty(strfind(fileName,'A011-20190219')))
        indCluGood = 4; 
        indCluBad = 0;
    elseif(~isempty(strfind(fileName,'A012-20190224')))
        indCluGood = 5; 
        indCluBad = 0;
    elseif(~isempty(strfind(fileName,'A012-20190221')))
        indCluGood = 4; 
        indCluBad = 3; %(? all the clusters seem to have good behavior trials)
    end
end
