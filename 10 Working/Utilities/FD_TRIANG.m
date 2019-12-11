function [ x ] = FD_TRIANG( N )

for i = 0:N-1
    x(i+1) = (2/N)*( (N/2) - abs(i - (N-1)/2) );
end