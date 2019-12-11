function [ spl ] = get_spl( data )
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here

spl = 20*log10(sqrt(sum(data.^2)/length(data))/20e-6);