load('A561-20210422-02B.mat');

msec = behEvents.TwoPSyncMsec;
arr = [arr(1)];
for i = 2:length(behEvents.TwoPSyncMsec)
   
    if(msec(i) - msec(i-1) < 50)
       arr(end+1) = i;
    end
    
end
length(arr)
% arr
disp('done')