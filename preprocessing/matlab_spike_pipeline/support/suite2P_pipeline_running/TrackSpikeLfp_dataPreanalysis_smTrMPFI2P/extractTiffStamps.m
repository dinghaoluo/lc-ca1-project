function extractTiffStamps(fileName)

    info = imfinfo([fileName '.tif']);
    timeStamp = zeros(1,length(info));
    for i = 1:length(info)
        strImage = info(i).ImageDescription;
        indStrStart = strfind(info(i).ImageDescription,'frameTimestamps_sec');
        indStrEnd = strfind(info(i).ImageDescription,'acqTriggerTimestamps_sec');
        timeStamp(i) = str2num(strImage(indStrStart+22:indStrEnd-1));
    end
