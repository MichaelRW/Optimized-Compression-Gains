function [ error ] = psth_err_mean( psth_struct1, psth_struct2, window, varargin )
% function to compute error metrics between two input psths
%
% used functions:
% find_mxmn()
%
% Input Arguments:
%      psth_struct1
%      psth_struct2
%      window
%
% Output Arguments:
%      error
%      struct with fields:
%                         hist: difference between rate histograms of the
%                               two psths
%                         psth: mean absolute difference betwee the two
%                               psths
%                         SPL:  unchanged input SPL
%                         SPL_uniq: input SPL after applying unique()
%                         ADJ:      unchanged input SPL
%                         hist_opti: ADJ or adjustment step, at which
%                                    difference between rate histograms
%                                    is minimal
%                         psth_opti: ADJ or adjustment step, at which
%                                    difference between psths is minimal
%                         type:     'FINE' or 'AVG', fine timing
%                                   information or average discharge rate
%                         meth:     default: 'mean', no other method
%                                   implemented
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Fine Timing
if ( strcmp(psth_struct1.type, 'FINE') && strcmp(psth_struct2.type, 'FINE') )
    
    %     if (nargin == 2), additional one because of window
    
    if (nargin == 3)
        disp('Not enough input variables')
    elseif strcmp(varargin{1}, 'make') && strncmpi(psth_struct1.calc_details, 'detailed', 6)
        error = struct( 'pwrr', [], 'psth', [], 'boxp', [], 'SPL', [], ...
                        'SPL_uniq', [], 'ADJ', [], 'pwrr_opti', [], ...
                        'psth_opti', [], 'boxp_opti', []);
        error.type = psth_struct1.type;
        error.meth = 'mean';
        return
    elseif strcmp(varargin{1}, 'make')
        error = struct('psth', [], 'SPL', [], 'SPL_uniq', [], 'ADJ', [], ...
                       'psth_opti', []);
        error.type = psth_struct1.type;
        error.meth = 'mean';
        return
    else
        error = varargin{1};
    end
    
    %=====================================================================%
    %                                PSTH                                 %
    %=====================================================================%
    
    % error.psth is difference between both psths, possibly changes with 
    % every adjustment step, because psth_impaired (usually psth_struct2) 
    % changes
    
    % original without window
    % error.psth = [error.psth mean(abs(psth_struct1.psth(:) - psth_struct2.psth(:)))];
    
    % multiply with window, dim num(freq) x length(dat)
    p1w = psth_struct1.psth;%.*window;
    
    %     MAX_Y_VALUE = 500;
    %        figure; ...
    %            subplot(2, 1, 1); plot( p1w(1, :) ); xlim( [1 length(psth_struct1.psth(1, :))] ); grid on; title( 'Normal PSTH' );
    %            subplot(2, 1, 2); plot( psth_struct1.psth(1, :)); hold on; plot(window(1, :).*MAX_Y_VALUE, 'r' ); axis( [1 length(psth_struct1.psth(1, :)) -50 1e3] ); grid on;
    %
    p2w = psth_struct2.psth;%.*window;
    %        figure; ...
    %            subplot(2, 1, 1); plot( p2w(1, :) ); xlim( [1 length(psth_struct2.psth(1, :))] ); grid on; title( sprintf('Impaired PSTH at %d dB Gain', varargin{2} ) );
    %            subplot(2, 1, 2); plot( psth_struct2.psth(1, :) ); hold on; plot(window(1, :).*MAX_Y_VALUE, 'r' ); axis( [1 length(psth_struct2.psth(1, :)) -50 1e3] ); grid on;
    %
    
    
    % mean absolute error of the two psth
    error.psth = [error.psth mean(abs(p1w(:) - p2w(:)))];
    
    
    % find_mxmn find the ADJ (third input param), at which the psth error
    % (second input param) is smallest 'min' or largest 'max' (fourth input
    % param)
    [ error.psth_opti, error.SPL_uniq ] = find_mxmn( error.SPL, error.psth, error.ADJ, 'min' );
    
    %=====================================================================%
    
    % adjustment error is the difference in SPL of the normal and changed
    % data
    
    % SPL error is the input SPL of the normal data, changes everytime when
    % function psth_err_mean is called. New element is spl of the whole
    % stimulus at the newly added phoneme (stimulus phonemes are appended
    % after another), combination might have different overall SPLs.
    error.SPL  = [error.SPL psth_struct1.data_orig.SPL];
% -------------------------------------------------------------------------
% mean rate neurogram     
elseif strcmp(psth_struct1.type, 'AVG') && strcmp(psth_struct2.type, 'AVG')
    if (nargin == 2)
        disp('Not enough input variables')
    elseif strcmp(varargin{1}, 'make') && strncmpi(psth_struct1.calc_details, 'detailed', 6)
        error = struct('hist', [], 'psth', [], 'SPL', [], 'SPL_uniq', [],...
                        'ADJ', [], 'hist_opti', [], 'psth_opti', []);
        %error.approx_ants = psth_struct1.data_orig.approx_formants;
        error.type = psth_struct1.type;
        error.meth = 'mean';
        return
    elseif strcmp(varargin{1}, 'make')
        error = struct('psth', [], 'SPL', [], 'SPL_uniq', [], 'ADJ', [],...
                       'psth_opti', []);
        error.type = psth_struct1.type;
        error.meth = 'mean';
        return
    else
        error = varargin{1};
    end
    
    %=====================================================================%
    %                           PR, PH, and BOX                           %
    %=====================================================================%
    
    % for the AVG type, additionally to the psth error, a histogram error
    % is computed. In PSTHmay, a histogram over average discharge is
    % computed, which is compared between the histograms here and then
    % optimized (optimal ADJ found)
    if strncmpi(psth_struct1.calc_details, 'detailed', 6)
        
        error.hist = [error.hist mean(abs(psth_struct1.hist(:) - psth_struct2.hist(:)))];
        [ error.hist_opti, error.SPL_uniq ] = find_mxmn( error.SPL, error.hist, error.ADJ, 'min' );
        
    end
    
    
    %=====================================================================%
    %                                PSTH                                 %
    %=====================================================================%
    
    % windowing
    p1w = psth_struct1.psth;%.*window;
    
    p2w = psth_struct2.psth;%.*window;
    %         figure; plot(psth_struct1.psth(1, :)); hold on; plot(window(1, :), 'r'); grid on; title( 'Normal PSTH' );
    
    %         figure; plot(psth_struct2.psth(1, :)); hold on; plot(window(1, :), 'r'); grid on; title( 'Impaired PSTH' );
    
    % computing the mean absolute error
    error.psth = [error.psth mean(abs(p1w(:) - p2w(:)))];
    
    % finding the ADJ, at which psth error is minimal
    [ error.psth_opti, error.SPL_uniq ] = find_mxmn( error.SPL, error.psth, error.ADJ, 'min' );
    
    
    
    %=====================================================================%
    % appending newly added SPL
    error.SPL  = [error.SPL psth_struct1.data_orig.SPL];
    
else
    disp('Unknown Structure Type');
end


