
function [ error ] = psth_err_mean( psth_struct1, psth_struct2, window, varargin )


% Fine Timing
if ( strcmp(psth_struct1.type, 'FINE') && strcmp(psth_struct2.type, 'FINE') )
    
    %     if (nargin == 2), additional one because of window
    
    if (nargin == 3)
        disp('Not enough input variables')
    elseif strcmp(varargin{1}, 'make') && strncmpi(psth_struct1.calc_details, 'detailed', 6)
        error = struct( 'pwrr', [], 'psth', [], 'boxp', [], 'SPL', [], 'SPL_uniq', [], 'ADJ', [], 'pwrr_opti', [], 'psth_opti', [], 'boxp_opti', []);
        error.type = psth_struct1.type;
        error.meth = 'mean';
        return
    elseif strcmp(varargin{1}, 'make')
        error = struct('psth', [], 'SPL', [], 'SPL_uniq', [], 'ADJ', [], 'psth_opti', []);
        error.type = psth_struct1.type;
        error.meth = 'mean';
        return
    else
        error = varargin{1};
    end
    
    %=====================================================================%
    %                                PSTH                                 %
    %=====================================================================%
    
    % error.psth is difference between both psths, changes everytime
    % function psth_error_mean is called, because psth_impaired
    % (psth_struct2) changes
    
    % original without window
    % error.psth = [error.psth mean(abs(psth_struct1.psth(:) - psth_struct2.psth(:)))];
    
    % multiply with window, dim num(freq) x length(dat)
    
    MAX_Y_VALUE = 500;
    %
    p1w = psth_struct1.psth.*window;
        figure; ...
            subplot(2, 1, 1); plot( p1w(1, :) ); xlim( [1 length(psth_struct1.psth(1, :))] ); grid on; title( 'Normal PSTH' );
            subplot(2, 1, 2); plot( psth_struct1.psth(1, :)); hold on; plot(window(1, :).*MAX_Y_VALUE, 'r' ); axis( [1 length(psth_struct1.psth(1, :)) 0 1e3] ); grid on;
        
    p2w = psth_struct2.psth.*window;
        figure; ...
            subplot(2, 1, 1); plot( p2w(1, :) ); xlim( [1 length(psth_struct2.psth(1, :))] ); grid on; title( sprintf('Impaired PSTH - %d dB Gain', varargin{2} ) );
            subplot(2, 1, 2); plot( psth_struct2.psth(1, :) ); hold on; plot(window(1, :).*MAX_Y_VALUE, 'r' ); axis( [1 length(psth_struct2.psth(1, :)) 0 1e3] ); grid on;
    
    % mean absolute error of the two psth,
    error.psth = [error.psth mean(abs(p1w(:) - p2w(:)))];
    
    [error.psth_opti error.SPL_uniq] = find_mxmn( error.SPL, error.psth, error.ADJ, 'min' );
    
    %=====================================================================%
    
    % adjustment error is the difference in SPL of the normal and changed
    % data
    
    % SPL error is the SPL of the normal data, changes everytime when
    % function psth_err_mean is called. new element is spl of the whole
    % stimulus (stimulus phonemes are appended after another), combination
    % might have different overall SPLs.
    error.SPL  = [error.SPL psth_struct1.data_orig.SPL];
    
elseif strcmp(psth_struct1.type, 'AVG') && strcmp(psth_struct2.type, 'AVG')
    if (nargin == 2)
        disp('Not enough input variables')
    elseif strcmp(varargin{1}, 'make') && strncmpi(psth_struct1.calc_details, 'detailed', 6)
        error = struct('hist', [], 'psth', [], 'SPL', [], 'SPL_uniq', [], 'ADJ', [], 'hist_opti', [], 'psth_opti', []);
        %error.approx_ants = psth_struct1.data_orig.approx_formants;
        error.type = psth_struct1.type;
        error.meth = 'mean';
        return
    elseif strcmp(varargin{1}, 'make')
        error = struct('psth', [], 'SPL', [], 'SPL_uniq', [], 'ADJ', [], 'psth_opti', []);
        error.type = psth_struct1.type;
        error.meth = 'mean';
        return
    else
        error = varargin{1};
    end
    
    %=====================================================================%
    %                           PR, PH, and BOX                           %
    %=====================================================================%
    
    if strncmpi(psth_struct1.calc_details, 'detailed', 6)
        error.hist = [error.hist mean(abs(psth_struct1.hist(:) - psth_struct2.hist(:)))];
        [error.hist_opti error.SPL_uniq] = find_mxmn( error.SPL, error.hist, error.ADJ, 'min' );
    end
    
    
    %=====================================================================%
    %                                PSTH                                 %
    %=====================================================================%
    p1w = psth_struct1.psth.*window;
        figure; plot(psth_struct1.psth(1, :)); hold on; plot(window(1, :), 'r'); grid on; title( 'Normal PSTH' );
    
    p2w = psth_struct2.psth.*window;
        figure; plot(psth_struct2.psth(1, :)); hold on; plot(window(1, :), 'r'); grid on; title( 'Impaired PSTH' );
    
    error.psth = [error.psth mean(abs(p1w(:) - p2w(:)))];
    
    [error.psth_opti error.SPL_uniq] = find_mxmn( error.SPL, error.psth, error.ADJ, 'min' );
    %=====================================================================%
    
    error.SPL  = [error.SPL psth_struct1.data_orig.SPL];
    
else
    disp('Unknown Structure Type');
end


