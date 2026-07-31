function calCorrectedFluo (FileNameBase, data_2p, isInt)
% calculate dF/F from the fluorescence data

    % extract cluster information 
    Clu = extractClu(data_2p);
    
    % calculate corrected fluorescence
    Fc = data_2p.F(Clu.localClu,:) - data_2p.ops.neucoeff*data_2p.Fneu(Clu.localClu,:); % corrected fluorescence
    F0 = zeros(size(Fc,1),size(Fc,2));
    
    % generate gaussian filter
    sig_baseline = 10; % sec
    std = round(sig_baseline*data_2p.ops.fs);
    h = gaussFilter2P(6*std,std);
    lenGaussKernel = length(h);
    
    % calculate the baseline using minmax filter
    win = round(data_2p.ops.win_baseline*data_2p.ops.fs);
    for i = 1:length(Clu.localClu)
        filtFcTmp = gauss_filter(Fc(i,:),h,lenGaussKernel);
        filtFcTmp = minimum_filter1d(filtFcTmp,win);
        filtFcTmp = maximum_filter1d(filtFcTmp,win);
        F0(i,:) = filtFcTmp;  % baseline fluorescence F0
        Fc(i,:) = Fc(i,:) - filtFcTmp; % F-F0
    end
    
    % calculate dFF
    dFF = Fc./F0;  
    
    save([FileNameBase '_corrFluo.mat'],'Clu', 'F0', 'Fc', 'dFF');
end

function output = gauss_filter(input,h,lenh)
    % get a row vector
    if(size(input,1) ~= 1)
        input = input';
    end
    
    % reflect the input vector at the beginning and end
    halfWin = ceil(lenh/2);
    inReflectStart = input(halfWin:-1:1);
    inReflectEnd = input(end:-1:end-halfWin+1);
    inTmp = [inReflectStart input inReflectEnd];
    
    output = conv(inTmp,h);
    
    output = output(halfWin+floor(lenh/2)+1:end-halfWin-lenh+floor(lenh/2)+1); 
                    % cut the convolution result to be the same length as the original data
end

function output = minimum_filter1d(input,win)
    % get a row vector
    if(size(input,1) ~= 1)
        input = input';
    end
    
    % reflect the input vector at the beginning and end
    halfWin = ceil(win/2);
    inReflectStart = input(halfWin:-1:1);
    inReflectEnd = input(end:-1:end-halfWin+1);
    inTmp = [inReflectStart input inReflectEnd];
    
    outTmp = movmin(inTmp,win);
    output = outTmp(halfWin+1:end-halfWin);
end

function output = maximum_filter1d(input,win)
    % get a row vector
    if(size(input,1) ~= 1)
        input = input';
    end
    
    % reflect the input vector at the beginning and end
    halfWin = ceil(win/2);
    inReflectStart = input(halfWin:-1:1);
    inReflectEnd = input(end:-1:end-halfWin+1);
    inTmp = [inReflectStart input inReflectEnd];
    
    outTmp = movmax(inTmp,win);
    output = outTmp(halfWin+1:end-halfWin);
end
