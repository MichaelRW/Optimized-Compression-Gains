function [ out, b ] = find_mxmn( key, values, look_up, mxmn )
% function to find ADJ, respective djustment step, at which the element of
% look_up, the smallest/largest (defined through mxmn) element of values 
% occurs
% 
% Input Arguments:
%                   key: orginal SPL values
%                   values: values, which should be minimized (usually 
%                           error metrics)
%                   look_up: parameters of the error metrics, which should
%                            be optimized (usually the ADJor adjustment 
%                            steps)
% Output Arguments: 
%                   error.psth_opti: element of look_up, at which values is
%                                    smallest (usually ADJ, adjustment step,
%                                    at which psth-error is smallest)
%                   error.SPL_uniq: SPL vector without redundant elements
%
% Example: 
% [error.psth_opti error.SPL_uniq] = find_mxmn( error.SPL, error.psth, error.ADJ, 'min' );

% comments adjusted to example

% [C,ia,ic] = unique() also returns index vectors ia and ic using any of the previous syntaxes.
% n therefore should be those indices where the elements of b occur in key.  

% find unique values in the original SPL of all phonemes and the indices 
[b, m, n] = unique(key);
ttt=nan;
out = [];
% go through all SPL values that are unique
for i = 1:max(n)
    posv = find(key == b(i));
    if strcmp(mxmn, 'max')
        % find maximum of values ( usually psth difference) for this SPL 
        [val pos2] = max(values(:,posv), [], 2);
    elseif strcmp(mxmn, 'min')
        % find minimum of values ( usually psth difference) for this SPL 
        [val pos2] = min(values(:,posv), [], 2);
    end
    
    % find the adjustment of both SPLs where the psth error is
    % minimal/maximal 
    out = [out look_up(posv(pos2))'];
end
tttt=1