function [combine] = overlapandadd(stratified_data, window, shift);
% Combines a stratified matrix

combine = stratified_data(1,:);

for i = 2:size(stratified_data,1)
    combine = [combine zeros(1,shift)];
    combine(end-window+1:end) = combine(end-window+1:end) + stratified_data(i,:);
end