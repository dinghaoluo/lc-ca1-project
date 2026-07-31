% if nargout>0, clear Eeg; end;% [ThetaPhase ThetaAmp TotPhase ThetaFreq] = GetThetaPhase_hilbert(Eeg, lfpSampleRate, thetaFiltParam)
%
% takes a 1-channel Eeg file (assumes 1250 Hz) and produces
% instantaneous theta phase and amplitude.  TotPhase is
% unwrapped phase.
%
% Theta is filtered with filtfilt and a cheby2 filter with parameters
% FreqRange, FilterOrd, Ripple defautls [4 10], 4, 20.
% good values for gamma are [40 100], 8, 20.
%
% if no args are provided, will plot some diagnostics
%   
% thetaFiltParam.FreqRange = [4 15];
% thetaFiltParam.FilterOrd = 4;
% thetaFiltParam.Ripple = 20;

function [ThetaPhase, ThetaAmp, TotPhase, ThetaFreq, Eegf] = ...
        GetThetaPhase_hilbert(Eeg, lfpSampleRate, thetaFiltParam)    

if min(size(Eeg)>1)
	error('Eeg should be 1 channel only!');
elseif size(Eeg,1)==1
    Eeg = Eeg(:);
end

[b a] = cheby2(thetaFiltParam.FilterOrd, thetaFiltParam.Ripple, ...
                thetaFiltParam.FreqRange/(lfpSampleRate/2),'bandpass');
Eegf = filtfilt(b,a,Eeg);
% remove constant term to avoid bias
Eegf = Eegf - mean(Eegf);

Hilb = hilbert(Eegf);
% if nargout>0, clear Eegf; end;
ThetaPhase = angle(Hilb);

if nargout>=2
	ThetaAmp = abs(Hilb);
end
if nargout>=3
	TotPhase = unwrap(ThetaPhase);
end

if nargout>=4
     ThetaFreq = diff(TotPhase)./(1/lfpSampleRate)./(2*pi);     
     % what is the change of phase per time bin with respect to one cycle -> freqency
end

if nargout==0
    subplot(4,1,1);
    [h w s] = freqz(b, a, 2048, 1250);
    plot(w,abs(h));
    grid on
    title('frequncy response of filter');
    
    subplot(4,1,2)
    xr = (1:length(Eeg))*1000/1250;
    plot(xr, [Eeg, Eegf]);
    title('eeg');
    legend('raw', 'filtered');
    
    subplot(4,1,3);
    plot(xr, [Eegf, ThetaPhase*std(Eegf)]);
    clear ThetaPhase
    title('extracted phase');
    legend('filtered wave', 'phase');
    
    subplot(4,1,4);
    plot(xr(2:end), ThetaFreq);
    clear ThetaFreq
    title('extracted frequency');
    legend('frequency');    
end
