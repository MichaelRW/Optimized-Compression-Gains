function [ data ] = NALR_old(data,FS,FREQS,H_IMP)
% This Code has been verified on July 20th, 2007.
%# PROBLEM ON LINE 21 and with INPUT PARAMETERS!!!
if (nargin == 1)
    disp('Need Data Sampling Rate');
    return;
elseif (nargin == 2)
    disp('Using arbitraty Threshold Shift');
    FREQS = [250 500 750 1000 1500 2000 3000 4000 6000 8000];
    H_IMP = [0 4 6 8 12 16 18 24 28 40];
end

data_length = length(data);
window = 256;
shift = 128;
fft_length = 2^(ceil(log2(window))+1);
zero_padding = fft_length - window;
f = -0.5:1/fft_length:0.5; f = f*FS; f(end) = [];

% H_IMP = H_IMP/2;
k = [-17 -8 -3 1 1 -1 -2 -2 -2 -2];
H3FA = (H_IMP(2) + H_IMP(4) + H_IMP(6))/3;
X = 0.15*H3FA;
IG_dB = X + 0.31*H_IMP + k;
IG_SC = 10.^(IG_dB/20);

% scale = fftshift(spline(FREQS,IG_SC,abs(f)));
% scale = fftshift(interp1([0 FREQS],[0 IG_SC],abs(f),'cubic'));
IG_SC_F = fftshift(interp1([0 FREQS],[0 IG_SC],abs(f),'linear'));
% Append zero at end to prevent inaccuracies when estimating gain near 0 Hz

data = stratify(data, window, 'shift', shift);
for i = 1:size(data,1);
    x = [data(i,:).*triang(window)' zeros(1,zero_padding)];  
    
    X = fft(x).*IG_SC_F;
    x = ifft(X);
    data(i,:) = x(1:window);
end
data = overlapandadd(data,window,shift);
data = data(1:data_length);