function scores = my_simple_gcc(s1,s2,wf,sw)
% s1 and s2 are the two signals to compute cross correlations for. wf is
% the weighting function, wf = 1 (can use wf=[]) for regular cross 
% correlation. sw is the clipping size for the search window.

if isempty(wf)
    wf = @(x) 1;
end

% make sure the signals are columns
s1 = s1(:);
s2 = s2(:);

% calculate the DFT of the two signals
S1 = fft(s1,2*length(s1)-1);
S2 = fft(s2,2*length(s2)-1);

% compute the cross-correlation with the given weight function
tmp = S1.*conj(S2);
% tmp = ifft(tmp.*wf(tmp));
tmp = fftshift(ifft(tmp.*wf(tmp)),1);

if nargin == 4
    scores = tmp(round(end/2)-sw:round(end/2)+sw);
else
    scores = tmp;
end