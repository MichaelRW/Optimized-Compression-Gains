function [ spl ] = get_spl( data )
% function to read out sound pressure level from data vector using
% root-mean-square (rms)
spl = 20*log10(sqrt(sum(data.^2)/length(data))/20e-6);