function [ t, x ] = sine( time, amplitude, frequency, sampling_rate );
% [ t, x ] = sine( time, freqeuncy, sampling_rate );

t = 0:1/sampling_rate:time;
t(end) = [];

x = zeros(size(t));
for i = 1:length(frequency)
    x = x + amplitude(i)*cos(2*pi*frequency(i)*t);
end