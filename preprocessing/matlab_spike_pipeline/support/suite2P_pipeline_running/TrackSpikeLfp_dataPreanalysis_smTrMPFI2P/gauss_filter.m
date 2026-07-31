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