function accumRecDataRunSpeedLick()
    RecordingListNT;
    
    recDataLabel{1} = 'Muscimol active licking';
    [recDataRunPre1{1},recDataRunManip1{1},recDataRunPost1{1}] = ...
                accumRecDataSpeedLick(activeLickMuscPath,recSessionsALMusc,...
                                ALMuscSessions,ALMuscMazeSession);
end

function [recDataPre,recDataManip,recDataPost] = ...
                accumRecDataSpeedLick(path,recSess,...
                                MSess,MazeSess)
    nRec = size(path,1);
    
    tmp.k
    
    for i = 1:nRec
        ind = findstr(path(i,:),'\');
        recName = path(i,ind(end-1)+1:ind(end)-1);
        disp(recName);
        fullPath = [path(i,:) recName '_lickDist.mat'];
        if(~exist(fullPath,'file'))
            disp([fullPath ' does not exist.']);
            return;
        end
        load(fullPath,'lickOverDist');
        
        indPre = find(recSess{i} < manipSess{i}(1));
        if(isempty(indPre))
            indPre = NaN;
        end
        indPost = find(recSess{i} > manipSess{i}(end));
        if(isempty(indPost))
            indPost = NaN;
        end
        [~,indManip] = intersect(recSess{i}, manipSess{i});
        nSessManip = length(manipSess{i});
        recDataRunPre = concatData(recDataRunPre,sessDataRun,sessDataRecTime,...
                    indPre*ones(1,nSessManip));
    end
end