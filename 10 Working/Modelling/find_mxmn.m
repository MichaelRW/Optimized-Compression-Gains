function [ out, b ] = find_mxmn( key, values, look_up, mxmn )
% example
% [error.psth_opti error.SPL_uniq] = find_mxmn( error.SPL, error.psth, error.ADJ, 'min' );

% comments adjusted to example

% [C,ia,ic] = unique() also returns index vectors ia and ic using any of the previous syntaxes.
% n therefore should be those indices where the elements of b occur in key.  

% find unique values in the original SPL of all phonemes and the indices 
[b, m, n] = unique(key);

out = [];
% go through 1 till maximum of unique indices % number of unique elements
for i = 1:max(n)
    % posv are the indices in original SPL, where values of the unique are the same 
    posv = find(key == b(i));
    if strcmp(mxmn, 'max')
        [val pos2] = max(values(:,posv), [], 2);
    elseif strcmp(mxmn, 'min')
        % find minimum of psth difference for this SPL 
        [val pos2] = min(values(:,posv), [], 2);
    end
    
    % find the adjustment of both SPLs where the psth error is
    % minimal/maximal 
    out = [out look_up(posv(pos2))'];
end
tttt=1