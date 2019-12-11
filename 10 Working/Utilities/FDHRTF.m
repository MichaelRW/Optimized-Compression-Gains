function [ y ] = FDHRTF( x)

load HRTF16K;
y = conv(ir16k,x);
y(end-100+2:end) = [];