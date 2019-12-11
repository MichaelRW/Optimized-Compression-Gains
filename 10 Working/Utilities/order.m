function [ out ] = order( in )
% Estimates the scientific order of the number.

out = sprintf('%1.0e', in);
out = str2num(out(strfind(out,'e')+1:end));