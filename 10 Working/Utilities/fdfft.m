function [f,output] = fdfft( input, FS );
% [f,output] = fdfft( input, FS );

input = reshape(input,1,length(input));     % Make Row.

len = length(input);
padding = 2^ceil(log(2*len)/log(2)) - len;

input = [input zeros(1,padding)];

% output = fftshift(abs(fft(input)/len));
output = abs(fftshift(fft(input)/len));

f = [-0.5:1/length(output):0.5]*FS;
f(end) = [];

% figure;
% plot(f,output)