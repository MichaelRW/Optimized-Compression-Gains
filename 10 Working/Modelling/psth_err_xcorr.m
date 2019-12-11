function [ error ] = psth_err_xcorr( psth_struct1, psth_struct2, varargin )

if strcmp(psth_struct1.type, 'FINE') && strcmp(psth_struct2.type, 'FINE')
    if (nargin == 2)
        disp('Not enough input variables')
    elseif strcmp(varargin{1}, 'make')
        error = struct('phsr', [], 'pwrr', [], 'psth', [], 'boxp', [], 'SPL', [], 'SPL_uniq', [], 'ADJ', [], 'phsr_opti', [], 'pwrr_opti', [], 'psth_opti', [], 'boxp_opti', []);
        error.approx_formants = psth_struct1.data_orig.approx_formants;
        error.type = psth_struct1.type;
        error.meth = 'xcorr';
        return
    else
        error = varargin{1};
    end

    % Make the ranges of phase the same %
    %[ phsr_freq1, phsr_freq2, phsr1, phsr2 ] = same_phsr_range( psth_struct1, psth_struct2 );
    
##    for i = 1:length(psth_struct1.phsr)
##        phsr(i) = xcorr(phsr1{i}, phsr2{i},0,'coeff');
##        pwrr(i) = xcorr(psth_struct1.pwrr(i,:), psth_struct2.pwrr(i,:),0,'coeff');
##    end
##    phsr(isnan(phsr)) = 0;
##    error.phsr = [error.phsr phsr'];
##    error.pwrr = [error.pwrr pwrr'];
    error.psth = [error.psth xcorr(psth_struct1.psth(:), psth_struct2.psth(:),0,'coeff')];
##   error.boxp = [error.boxp xcorr(psth_struct1.boxp(:), psth_struct2.boxp(:),0,'coeff')];
    error.ADJ  = [error.ADJ (psth_struct2.data_orig.SPL - psth_struct1.data_orig.SPL)];
    error.SPL  = [error.SPL psth_struct1.data_orig.SPL];
    
##    [error.phsr_opti error.SPL_uniq] = find_mxmn( error.SPL, error.phsr, error.ADJ, 'max' );
##    [error.pwrr_opti error.SPL_uniq] = find_mxmn( error.SPL, error.pwrr, error.ADJ, 'max' );
    [error.psth_opti error.SPL_uniq] = find_mxmn( error.SPL, error.psth, error.ADJ, 'max' );
##    [error.boxp_opti error.SPL_uniq] = find_mxmn( error.SPL, error.boxp, error.ADJ, 'max' );

elseif strcmp(psth_struct1.type, 'AVG') && strcmp(psth_struct2.type, 'AVG')
    if (nargin == 2)
        disp('Not enough input variables')
    elseif strcmp(varargin{1}, 'make')
        error = struct('hist', [], 'psth', [], 'SPL', [], 'SPL_uniq', [], 'ADJ', [], 'hist_opti', [], 'psth_opti', []);
        error.approx_formants = psth_struct1.data_orig.approx_formants;
        error.type = psth_struct1.type;
        error.meth = 'xcorr';
        return
    else
        error = varargin{1};
    end
    
    error.hist = [error.hist xcorr(psth_struct1.hist(:), psth_struct2.hist(:),0,'coeff')];
    error.psth = [error.psth xcorr(psth_struct1.psth(:), psth_struct2.psth(:),0,'coeff')];
    error.ADJ  = [error.ADJ (psth_struct2.data_orig.SPL - psth_struct1.data_orig.SPL)];
    error.SPL  = [error.SPL psth_struct1.data_orig.SPL];
    
    [error.hist_opti error.SPL_uniq] = find_mxmn( error.SPL, error.hist, error.ADJ, 'max' );
    [error.psth_opti error.SPL_uniq] = find_mxmn( error.SPL, error.psth, error.ADJ, 'max' );
else
    disp('Unknown Structure Type');
end