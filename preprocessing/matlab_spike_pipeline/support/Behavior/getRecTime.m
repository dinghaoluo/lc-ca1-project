function getRecTime()

    RecordingListPerAnm

    for i = 13:14
        numRec = size(recFile{i},1);
        for j = 1:numRec
            disp(recFile{i}(j,:));
            fInfo = [];
            for k = 1:length(recSess{i}{j})
                fullPath = [recPath{i} recFile{i}(j,:) '-0' num2str(recSess{i}{j}(k)) 'T.txt'];
                fileInfo = dir(fullPath);
                addDate = 1;
                while(isempty(fileInfo))
                    if(addDate > 2)
                        disp('There is no recording file found in the next two days.')
                        break;
                    end
                    ind = strfind(recFile{i}(j,:),'-');
                    dateRec = recFile{i}(j,ind(1)+1:end);
                    dateRec = [dateRec(1:4) '-' dateRec(5:6) '-' dateRec(7:8)];
                    actDate = datetime(dateRec) + days(addDate);     
                    str = datestr(actDate,30);
                    fullPath = [recPath{i} recFile{i}(j,1:ind) str(1:8) '-01T.txt'];
                    fileInfo = dir(fullPath);
                    addDate = addDate + 1;
                end
                fInfo.time{k} = fileInfo.date;
                fInfo.recFile{k} = [recFile{i}(j,:) '-0' num2str(recSess{i}{j}(k))];
                fInfo.sessNo(k) = recSess{i}{j}(k);
            end
            savePath = [recSavePath{i} recFile{i}(j,:) '\' recFile{i}(j,:)  '_recTime.mat'];
            save(savePath,'fInfo');
        end
    end
end
