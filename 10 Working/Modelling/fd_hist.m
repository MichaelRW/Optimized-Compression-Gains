function [ bins, psth_hist ] = fd_hist( psth_struct, varargin )
% FD_HIST computes the count distribution of rates in the neurogram.
% [ bins, histogram ] = fd_hist( Neurogram_Struct, varargin )
TTTT=nan;
if ~strcmp(psth_struct.type, 'AVG')
    bins = [];
    psth_hist = [];
    return
end

fibre = size(psth_struct.psth,1);
psth_hist = [];

for i = 1:fibre
    [a, bins] = hist(psth_struct.psth(i,:), 0:5:400);
    psth_hist(i,:) = smooth_Octave(a);
    %psth_hist_mat(i,:) = smooth(a);
  
end
 ttt=nan;
if (nargout == 0) | strncmpi(varargin,'y',1)
%     figure;
%     imagesc(bins, psth_struct.psth_freq, psth_hist); axis xy
%     set(gca, 'YScale', 'log');
%     set(gca, 'YTick', [250 500 750 1000 1500 2000 3000 4000 6000 8000]);
%     set(gca, 'YTickLabel', [250 500 750 1000 1500 2000 3000 4000 6000 8000]);
%     ylabel('Center Frequency (Hz)')
%     xlabel('Rate Bins (spikes/sec/fiber)')
%     title('Histogram Response')
end