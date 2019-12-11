function [ x ] = FD_HANNING( N )
i = 0:N;
x(i+1) = 0.5 - 0.5*cos(2*pi*i/(N));
x(end) = [];