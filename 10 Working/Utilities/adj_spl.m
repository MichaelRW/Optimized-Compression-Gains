function [ data ] = adj_spl( data, adj_spl )
% function to read out spl from data and apply input gain (adj_spl)
spl  = get_spl(data);
data = set_spl(data,spl + adj_spl);