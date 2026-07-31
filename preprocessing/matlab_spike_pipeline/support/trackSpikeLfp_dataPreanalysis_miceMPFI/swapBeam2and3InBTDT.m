function swapBeam2and3InBTDT(baseFileName)
% In some cases, the beam break 2 and 3 appear in the wrong order.
% Therefore, the order need to be corrected before running SortTrials_2armMaze_v2.m 

    fileName = [baseFileName 'BTDT.mat'];
    if exist(fileName, 'file') ~= 2
        disp('BTDT file does not exist.')
        return;
    else
        load(fileName);
        ind3 = behEventsTdt.beam(:,3) == 3;
        ind2 = behEventsTdt.beam(:,3) == 2;
        behEventsTdt.beam(ind3,3) = 2;
        behEventsTdt.beam(ind2,3) = 3;
        save(fileName,'behEventsTdt');
    end
end