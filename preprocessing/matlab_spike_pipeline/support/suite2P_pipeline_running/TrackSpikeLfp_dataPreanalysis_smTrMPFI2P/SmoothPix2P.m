function smoothed = SmoothPix(data,smooth)

% smoothed = SmoothPix(data,smooth)
%
% Smooth using a Gaussian kernel
%  smooth:  horizontal and vertical standard deviations [Sh Sv]
%           for Gaussian kernel (0 = no smoothing)
%

vector = min(size(data)) == 1;

if nargin < 2,
  error('Incorrect number of parameters (type ''help Bin'' for details).');
end

%
% Required parameters:
%  data:    the data to smooth
if ~isa(smooth,'numeric'),
  error('Incorrect value for property ''smooth'' (type ''help Smooth'' for details).');
end
if (vector && length(smooth) ~= 1) || length(smooth) > 2 || sum(smooth < 0) > 0,
  error('Incorrect value for property ''smooth'' (type ''help Smooth'' for details).');
end

% If Sh = Sv = 0, no smoothing required
if sum(smooth) == 0,
	smoothed = data;
	return
end

% For 2D data, providing only one value S for the std is interpreted as Sh = Sv = S
if ~vector && length(smooth) == 1,
	smooth = [smooth smooth];
end

% Vectors must be vertical
if vector,
	if size(data,2) ~= 1,
		data = data';
	end
	smooth = [0 smooth];
end

[vSize,hSize] = size(data);
maxSize = max([hSize vSize]);
r = (-maxSize:maxSize)/maxSize;
if smooth(1) ~= 0,
  hSmoothFactor = smooth(1)/hSize;           % determines size of e(det) -> the smaller the sharper finall template 
  hSmoother = exp(-r.^2/hSmoothFactor^2/2);  % min= very large negative # -> exp turnes that to 0, max=peak = 0 -> exp trunes to 1 => creates a sharply peaked function = template
  hSmoother = hSmoother/sum(hSmoother);
else
  hSmoother = zeros(size(r));
  hSmoother((end-1)/2+1) = 1;
end
if smooth(2) ~= 0,
  vSmoothFactor = smooth(2)/vSize;      
  vSmoother = exp(-r.^2/vSmoothFactor^2/2);
  vSmoother = vSmoother/sum(vSmoother);
else
  vSmoother = zeros(size(r));
  vSmoother((end-1)/2+1) = 1;
end

if(smooth(1) ~= 0)
    smoothed = conv2(vSmoother,hSmoother,data,'same');
else
    smoothed = conv(data,vSmoother,'same');
end
