function [A] = rcumsum(A,D)
%RCUMSUM cumulative sum of elements, restarted after every zero.
% For vectors, RCUMSUM(X) is a vector containing the cumulative sum of
% the elements of X, with the summation restarting after a zero is 
% encountered. For matrices, RCUMSUM(X) is a matrix the same size
% as X containing the cumulative restarted sums over each column.  For N-D
% arrays, RCUMSUM(X) operates along the first non-singleton dimension.
%  
% RCUMSUM(X,DIM) works along the dimension DIM.
%  
%     Example: If X = [0 1 2 
%                      3 0 4 
%                      1 1 3]
%  
%     then rcumsum(X,1) is [0 1 2  and rcumsum(X,2) is [0 1  3
%                           3 0 6                       3 0 4
%                           4 1 9]                      1 2 5]
%  
% This function has a MEX file helper function which will run only if
% compiled.  If the MEX function is not compiled, then a MATLAB
% implementation will be used.  The MEX function is RCUMSUMC.CPP
%
%    Class support: double, single and (logical - will convert to double).
%
%
%    See also cumsum, cumprod, sum, prod.
%
% Author:  Matt Fig
% Date:  9/2/2010
S = size(A);
if nargin<2
    D = find(S>1,1,'first'); % Default is the same as CUMSUM.
end
if D>ndims(A) || ~isnumeric(D) || ~isscalar(D) || isempty(A) || D<1
    return  % Just return A as it is passed.
end
CLA = strcmp(class(A),{'single';'double';'logical'}); % Get input class.
if ~any(CLA)
    error('Input A must be a single, double or logical array.')
end
try  % See if the MEX helper function has been compiled.
    if CLA(3)
        A = double(A);  % Turn a logical into a double for the MEX.
    end
    A = rcumsumc(A,double(D));  % Call helper function.
catch
%   This is the available pure MATLAB fallback when the helper MEX function is
%   unavailable.
%   Note there is the possibility of overflow here.  The MEX function avoids
%   this much better than this code.  If for some reason you cannot use the
%   MEX function and are experiencing overflow, uncomment the code below
%   and comment out this first section of code.
    SD = S(D);
    L = 1:length(S);
    S(D) = [];
    L(D) = [];
    A = [permute(A,[D L]);zeros([1,S],class(A))];  % Pad with zeros.
    A = A(:);
    I = ~A;  % Find the zeros.
    H = cumsum(A);  % Here is where the overflow could occur.
    H = H(I);  % We want the cumsummed numbers in the zero locations.
    A(I) = -[H(1); diff(H)]; % Subtract them out.
    A = cumsum(A);  % And redo it.
    A(SD+1:SD+1:numel(A)) = [];  % Remove the zero padding.
    A = ipermute(reshape(A,[SD S]),[D L]);  % Get back into shape.
    
%   The alternative, (more-so) overflow avoiding code.  To use, comment
%   the above code and uncomment the following.   
%     L = 1:length(S);
%     L(D) = [];
%     A = permute(A,[D L]);
%     N = numel(A);
%     H = 2:N;
%     H = H(logical(mod(H-1,S(D))));
% 
%     for ii = H  % O.K., this will be slow.
%         if A(ii)
%             A(ii) = A(ii-1) + A(ii);
%         end
%     end
%     A = ipermute(A,[D L]);
end