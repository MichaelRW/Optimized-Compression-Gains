function [ data ] = set_spl( data, in_spl )
% function for rms normalization of input data input data and
% setting level by to desired SPL (in_spl).

rms = sqrt(sum(data.^2)/length(data));

spl = 20*log10(rms/20e-6);

dB_gain = in_spl - spl;

k = 10.^(dB_gain/20);

data = k.*data;
