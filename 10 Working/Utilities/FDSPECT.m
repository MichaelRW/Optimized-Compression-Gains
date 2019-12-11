function [ X,f,t ] = FDSPECT( data, FS, window, shift )

data_len = length(data);

fft_length = 2^(ceil(log2(window))+1);
zero_padding = fft_length - window;

f = -0.5:1/fft_length:0.5; f = f*FS; f(end) = [];
f_index = find(f>=0 & f<=8000);
f = f(f_index);

data = stratify(data, window, 'shift', shift);

for i = 1:size(data,1);
    x = [data(i,:).*FD_HAMMING(window) zeros(1,zero_padding)];
    x_fft = 20*log10(abs(fft(x)));
    X(:,i) = x_fft(f_index);
end

t = [0:1/size(X,2):1]*data_len/FS;