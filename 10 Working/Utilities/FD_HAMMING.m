function [ x ] = FD_HAMMING( N )

for i = 0:N-1
    x(i+1) = 0.53836 - 0.46164*cos(2*pi*i/(N-1));
end