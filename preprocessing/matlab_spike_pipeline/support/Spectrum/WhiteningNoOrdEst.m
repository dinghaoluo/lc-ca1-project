function [whiteV] = WhiteningNoOrdEst(v, p, arType)
% whiten the signal represented by vector v. This function does not
% estimate the order of AR model
%
% v:            vector to be whitened
% p:            the order of AR model
% arType:       1: using ARFIT for AR estimation
%               2: using BURG method for AR estimation
%
% whiteV:       whitened vector

%%%%%%%%% check arguments
if nargin<2
    disp('At least two arguments are needed for this function.');
    return;
elseif nargin == 2
    arType = 1; % by default, use ARFIT to estimate the parameters
elseif nargin > 3
    disp('Too many input arguments.');
    return;
end    

if(arType == 1) % using arfit
    [w,Atmp] = arfit(v,p,p);
    A = [1 -Atmp];
else % using burg method
    [w,A] = arburg(v,p);
end

whiteV = filter(A,1,v);
