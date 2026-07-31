function h = gaussFilter(length, std)
% generate Gaussian kernel
% length:       kernel length
% std:          the std of Gaussian
    if(mod(length,2) == 0)
        siz = length/2;
    else
        siz   = (length-1)/2;
    end
    x = -siz:siz;
    arg   = -(x.*x)/(2*std*std);

    h     = exp(arg);
    h(h<eps*max(h(:))) = 0;

    sumh = sum(h(:)); % normalize energy to 1
    if sumh ~= 0,
        h  = h/sumh;
    end;
    
%     sumh = sum(h(:).*h(:)); % normalize energy to 1
%     if sumh ~= 0,
%         h  = h/sqrt(sumh);
%     end;