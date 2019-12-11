function [ data ] = DSL_old(data,FS,FREQS,H_IMP)
% Verified Friday July 20th, 2007.
% Implements half of the DSL alogrithm

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

%%%%%%%%%%DSL%%%%%%%%%%

[x y] = meshgrid(1:23);

x = x(:,1:10);
F = [250 500 750 1000 1500 2000 3000 4000 6000 8000];
F = F(x);

y = y(:,1:10);
H = [0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110];
H = H(y);

Z = [[0 3 5 7 9 12 14 17 20 22 25 29 32 36 39 43 47 51 55 59 62 66 68]'...
[2 4 6 8 11 13 15 18 20 23 26 29 32 35 38 42 45 48 52 55 59 62 66]'...
[3 5 7 10 12 14 17 19 22 25 28 31 34 37 40 43 47 50 54 57 61 64 68]'...
[3 5 8 10 13 15 18 21 24 27 30 33 36 40 43 46 50 53 57 60 64 68 71]'...
[5 8 10 13 15 18 20 23 26 29 32 35 38 42 45 48 52 55 59 62 66 70 73]'...
[12 15 17 19 22 24 27 30 33 36 39 42 46 49 52 56 59 63 66 70 73 77 80]'...
[16 18 20 23 25 28 30 33 36 39 42 45 48 52 55 59 62 66 69 73 76 80 83]'...
[14 17 19 21 24 27 29 32 35 38 41 45 48 51 55 58 62 65 69 73 76 80 84]'...
[8 11 14 17 20 23 26 29 32 36 39 43 46 50 54 58 61 65 69 73 76 81 85]'...
[4 5 9 11 14 16 19 22 25 28 31 34 37 41 44 47 51 54 58 61 65 69 72]'];

GAIN_dB = interp2(F,H,Z,FREQS,H_IMP);
GAIN_SC = 10.^(GAIN_dB/20);
GAIN_SC_F = fftshift(interp1([0 FREQS],[0 GAIN_SC],abs(f),'linear'));

data = stratify(data, window, 'shift', shift);
for i = 1:size(data,1);
    x = [data(i,:).*triang(window)' zeros(1,zero_padding)];  
    
    X = fft(x).*GAIN_SC_F;
    x = ifft(X);
    data(i,:) = x(1:window);
end
data = overlapandadd(data,window,shift);
data = data(1:data_length);