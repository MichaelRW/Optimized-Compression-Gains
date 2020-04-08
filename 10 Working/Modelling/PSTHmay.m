function [ psth_struct ] = PSTHmay( S, aud, binwidth, varargin )
% PSTHmay Computes the AN Neurogram.
%
% [ Neurogram_Structure ] =
% PSTHmay(data_structure, audi_structure, binwidth, plot)
%
% This function computes the Neurogram at 'binwidth' resolution using the
% spike timing information from the Zilany-Bruce Cat Auditory Model. The AN
% response is taken at 'CF_Points', different center frequency locations,
% spaced logarithmically along the basilar membrane. At each CF point, the
% function calculates the summed response of 50 nerve fibers organized as
% 60% high, 20% medium, and 20% low spontaneous rate fibers. The Audiogram
% structure specifies the hearing loss profile, as +dB loss, at specific
% frequencies. In addition to calculating the Neurogram, phase response,
% power ratio, box plots, and histogram analysis are calculated where
% possible. An optional potting parameter, if 'y' plots the Neurogram.
%
% Input Parameters:
% S:                input data structure
% aud:              audiogram structure
% binwidth:         specified binwidth (see binwidth below/output parameter)
%
% Output Parameters: psth_struct with fields:               
% type:             either 'FINE' or 'AVG', depending on the input binwidth
%                   'FINE' fine-timing neurogram, binwidth = 10e-6 [sec]
%                   'AVG' mean-rate-neurogram, binwidth = 100e-6 [sec]
% psth_time:        time vector of the psth
% psth:             output psth size
% psth_freq:        CFs of the used fibers
% psth_mnmx:        minimum and maximum value of computed psth size [1 x 2]
% binwidth:         binwidth used to compute the psth
%                   function, fine-time: 10e-6 [sec], mean-rate: 100e-6
%                   [sec]
% data_orig:        unamplified original data
% audiogram_struct: identical with input audigram struct 'aud'
% F:                frequency axis to phase and magnitude spectrum,
%                   provided by quickfft()
% PSTH:             magnitude spectrum of psth computed with quickfft() 
% PHASE:            phase spectrum computed from psth with quickfft()
% calc_details:     from input, options: 'simple' and 'detailed', if
%                   'detailed', then additional information is infered 
%                   from the psth (psth_norm for 'FINE' and hist for 'AVG') 
% hist_bins:        output of fd_hist(), bins of the rate histogram
%                   size [1 x number of histogram bins]
% hist:             output of fd_hist(), rate histogram computed from 
%                   average discharge rate neurogram,
%                   size [length(psth_freq)x number of histogram bins]
% hist_mnmx:        minimum and maximum value of average rate histogram 
%                   (computed from average discharge neurogram)
%                   [1 x 2]

% See Also 
%           make_data_struct()
%           audiograms()
%           fd_phsr()
%           fd_boxp()
%           fd_hist()
%           psth_plot()
%
%--------------------------------------------------------------------------
% Version 0: Faheem Dinath. June 7th 2008.
% Modifications and change to 2018 model version:
%   H.T.Heinermann April 2020, helen@heinermann.net
%--------------------------------------------------------------------------
%% Constants: Time-related
% sampling frequency
fs = 1e5;  % May change if FS is greater than fs %

% time step
dt = 1/fs;

% number of bins per second for psth, 10e-6 in case of FT and 100e-6 in 
% case of AVG neurogram
number_bins = double(single(1/binwidth));

% How many samples are in each bin.
number_samples_in_bin = round(binwidth*fs);

%% Frequency Specifics
% frequency for psth generation, also CF, size = 1x number of CF
psth_freq = aud.F;

%% Check for proper FS and Binwidth
if ( number_bins >= fs) && (S.FS < number_bins )
    type = 'FINE';
    fs = number_bins;
    dat = resample(S.data,fs,S.FS);
    T = length(dat)/fs;
elseif ( number_bins < fs ) && ( S.FS < fs )
    type = 'AVG';
    if ~iswhole(binwidth*fs)
        disp('Make sure that binwidth*100e3 is a whole number')
        return;
    end   
    dat = resample(S.data,fs,S.FS);
    T = length(dat)/fs;
    disp(['Returning Average Discharge Rate With Binwidth of ' num2str(binwidth) 's'])
elseif ( S.FS >= fs ) && ( number_bins < fs )
    type = 'AVG';
    fs = S.FS;
    if ~iswhole(binwidth*fs)
        disp('Make sure that binwidth*FS is a whole number')
        return;
    end
    dat = S.data;
    T = length(dat)/fs;
    disp(['Returning Average Discharge Rate With Binwidth of ' num2str(binwidth) 's'])
elseif ( S.FS >= number_bins ) && ( number_bins >=fs )
    type = 'FINE';
    fs = S.FS;
    dat = S.data;
    T = length(dat)/fs;
    disp('Returning Fine Timing Information')
else
    disp('Something is wrong with FS and binwidth.')
    return;
end
%% save copy of data input struct
if  isfield(S, 'data_orig') % If there's been amplification.
    data_orig = S.data_orig;
    S = rmfield(S,'data_orig');
    data_pres = S;
else                        % Otherwise there hasn't been amplification.
    data_orig = S;
end
%% Impairment Parameters
% vectors, size = 1 x number of CF
Cohc = aud.Cohc;  Cihc = aud.Cihc;

% number of repetitions for each fiber type (LSR,MSR,HSR) for healthy system
number_rep_healthy = round( 50*[0.2 0.2 0.6] );
%% Auditory nerve fiber population is generated or loaded
if exist('ANpopulation.mat','file')
    load('ANpopulation.mat','sponts', 'tabss', 'trels');
    disp('Loading existing population of AN fibers saved in ANpopulation.mat')
    if (size(sponts.LS,2)<number_rep_healthy(1))||(size(sponts.MS,2)<number_rep_healthy(2))||(size(sponts.HS,2)<number_rep_healthy(3))||(size(sponts.HS,1)<length(psth_freq)||~exist('tabss','var'))
        disp('Saved population of AN fibers in ANpopulation.mat is too small - generating a new population');
        [sponts,tabss,trels] = generateANpopulation(length(psth_freq),number_rep_healthy);
    end
else
    [sponts,tabss,trels] = generateANpopulation(length(psth_freq),number_rep_healthy);
    disp('Generating population of AN fibers, saved in ANpopulation.mat')
end
%% synaptopathy conditions:
% depending on the synaptopathy input string, the number of healthy LSR,
% MSR and HSR - auditory nerve fibers is changed.
if strcmp(S.SYNAPTOPATHY, 'healthy')
    nrep = round([1 1 1].*number_rep_healthy); % Healthy AN
elseif strcmp(S.SYNAPTOPATHY, '50all')
    nrep = round([0.5 0.5 0.5].*number_rep_healthy); % 50% fiber loss of all types
elseif strcmp(S.SYNAPTOPATHY, 'lossall_low')
    nrep = round([0 1 1].*number_rep_healthy); % Loss of all LS fibers
elseif strcmp(S.SYNAPTOPATHY, 'IHCproportional')
    nrep = round([mean(aud.Cihc) 1 mean(aud.Cihc)].*number_rep_healthy); % loss of LS and HS fibers proportional to IHC impairment
else
    warning('synaptopathy type unknown, healthy AN assumed');
end
%% duration of simulated psth
% define simulation duration to integer multiple of binwidth to avoid
% additional samples for later reshaping
simdur = ceil( T * 1.2 / binwidth) * binwidth;
%% Filter properties
% lines 52 to 56 from generate_neurogram_BEZ2018_octaveparallel()
windur_ft = 32;  % Size of window for fine-timing neurogram.
smw_ft = hamming(windur_ft);
%
windur_mr = 128;  % Size of window for mean-rate neurogram.
smw_mr = hamming(windur_mr);
%% Pre-allocation of LSR, MSR, HSR fiber psth, and psth of all fibers

% just a template for pre-allocation in the following lines
% size: [1x approx.simdur*fs]
vihc_temp = model_IHC_BEZ2018( dat, psth_freq(1), 1, dt, simdur, Cohc(1), Cihc(1), 2 );

% pre-allocate psth matrices, size: [number of fibers x length of IHC response]
psth_LSR_single_fibers = NaN( length(nrep(1)), length(vihc_temp) );
psth_MSR_single_fibers = NaN( length(nrep(2)), length(vihc_temp) );
psth_HSR_single_fibers = NaN( length(nrep(3)), length(vihc_temp) );

psth_all_fibers = NaN( length(psth_freq), length(vihc_temp) );

psth_ft = NaN( length(psth_freq), length(vihc_temp) );
psth_mr = NaN( length(psth_freq), length(vihc_temp)/number_samples_in_bin );

neurogram_ft = NaN( length(psth_freq), length(vihc_temp) );
neurogram_mr = NaN( length(psth_freq), length(vihc_temp)/number_samples_in_bin );
% ----- loop over all CFs--------------------------------------------------
for control_freq = 1:length(psth_freq)
    
    fprintf( 1, '\n\tCF: %d of %d - %d Hz', control_freq, numel(psth_freq), psth_freq(control_freq) );
    
    % spontaneous rates, absolute and relative refractory rates of all
    % fibers (LSR, MSR and HSR fibers) concatenated
    sponts_concat = [sponts.LS(control_freq,1:nrep(1)) sponts.MS(control_freq,1:nrep(2)) sponts.HS(control_freq,1:nrep(3))];
    tabss_concat = [tabss.LS(control_freq,1:nrep(1)) tabss.MS(control_freq,1:nrep(2)) tabss.HS(control_freq,1:nrep(3))];%sec
    trels_concat = [trels.LS(control_freq,1:nrep(1)) trels.MS(control_freq,1:nrep(2)) trels.HS(control_freq,1:nrep(3))];%sec
    
    % Generate inner hair cell response for each frequency.
    vihc_temp = model_IHC_BEZ2018( dat, psth_freq(control_freq), 1, dt, simdur, Cohc(control_freq), Cihc(control_freq), 2 );
    
    % Accumulate synapse responses for low spontaneous fibers.
    %     fprintf( 1, '\n\t\tComputing low spont. fiber responses.' );
    for ma = 1:nrep(1)
        spont = sponts_concat(ma);  tabs = tabss_concat(ma);  trel = trels_concat(ma);
        
        % Took length(dat)*ts because saving the variable caused an accuracy error and length for reshape did not work.
        [ psth_LSR_temp, ~, ~, ~, ~, ~ ] = model_Synapse_BEZ2018( vihc_temp, psth_freq(control_freq), 1, dt, 1, 0, spont, tabs, trel );
        psth_LSR_single_fibers(ma, :) = psth_LSR_temp;
    end
    % sum up psth of all LSR fibers for this control frequency
    psth_LSR = sum(psth_LSR_single_fibers,1);
    
    % Accumulate synapse responses for medium spontaneous fibers.
    %     fprintf( 1, '\n\t\tComputing medium spont. fiber responses.' );
    %
    for mb = 1:nrep(2)
        spont = sponts_concat(mb+nrep(1));  tabs = tabss_concat(mb+nrep(1));  trel = trels_concat(mb+nrep(1));
        
        [ psth_MSR_temp, ~, ~, ~, ~, ~ ] = model_Synapse_BEZ2018( vihc_temp, psth_freq(control_freq), 1, dt, 1, 0, spont, tabs, trel );
        psth_MSR_single_fibers(mb, :) = psth_MSR_temp;
        
        
    end
    % sum up psth of all MSR fibers for this control frequency
    psth_MSR = sum(psth_MSR_single_fibers,1);
    
    
    % Accumulate synapse responses for high spontaneous fibers.
    %     fprintf( 1, '\n\t\tComputing high spont. fiber responses.' );
    %
    for mc = 1:nrep(3)
        spont = sponts_concat(mc+nrep(1)+nrep(2));  tabs = tabss_concat(mc+nrep(1)+nrep(2));  trel = trels_concat(mc+nrep(1)+nrep(2));
        
        [ psth_HSR_temp, ~, ~, ~, ~, ~ ] = model_Synapse_BEZ2018( vihc_temp, psth_freq(control_freq), 1, dt, 1, 0, spont, tabs, trel );
        psth_HSR_single_fibers(mc, :) = psth_HSR_temp;
    end
    % sum up psth of all HSR fibers for this control frequency
    psth_HSR = sum(psth_HSR_single_fibers,1);
    
    % sum up psth of all fiber types and save for this control frequency
    psth_all_fibers(control_freq,:) = psth_LSR + psth_MSR + psth_HSR;  % Complete PSTH response.

    
    psth_ft(control_freq,:) =      psth_all_fibers(control_freq,:); %round(10e-6*fs);%sum( reshape(psth500k, round(10e-6*fs), length(psth500k) / round(10e-6*fs) ));
    psth_mr(control_freq,:) =      sum( reshape(psth_all_fibers(control_freq,:), number_samples_in_bin,            length(psth_all_fibers(control_freq,:)) / number_samples_in_bin ));
    
   % pre-allocate neurogram
    neurogram_ft(control_freq,:) = zeros(1, size(vihc_temp, 2) );
    neurogram_mr(control_freq,:) = zeros(1, size(vihc_temp, 2) / number_samples_in_bin);
    %figure; plot(psth_ft, 'x')
    
    % filter psth with hamming window defined in ll. 164-168 
    neurogram_ft(control_freq,:) = neurogram_ft(control_freq,:) + filter( smw_ft', 1, psth_ft(control_freq,:));
    neurogram_mr(control_freq,:) = neurogram_mr(control_freq,:) + filter( smw_mr', 1, psth_mr(control_freq,:) );
    
end


pr_ft = NaN(length(psth_freq), length(1:windur_ft/2:length(neurogram_ft)));
pr_mr = NaN(length(psth_freq), length(1:windur_mr/2:length(neurogram_mr)));

% resampling after filtering
for control_freq = 1:length(psth_freq)
    pr_ft(control_freq,:) = neurogram_ft(control_freq, 1:windur_ft/2:end );    
    pr_mr(control_freq,:) = neurogram_mr(control_freq, 1:windur_mr/2:end );
end

t_ft = 0:windur_ft/2/fs:( size( pr_ft, 2 ) - 1 ) * windur_ft / 2 / fs;
t_mr = 0:windur_mr/2*binwidth:( size( pr_mr, 2 ) - 1 ) * windur_mr / 2 *binwidth;
%% saving everything in the output structure
psth_struct.type = type;

if strcmp(type, 'FINE')
    psth=pr_ft;
fprintf(1, '\n\n');

    psth_struct.psth_time = t_ft;%[0:(size(psth,2)-1)]*binwidth;
    
elseif strcmp(type, 'AVG')
    psth=pr_mr;
    fprintf(1, '\n\n');

    psth_struct.psth_time = t_mr;%[0:(size(psth,2)-1)]*binwidth;
end

psth_struct.psth = psth(:,1:size(psth,2));

psth_struct.psth_freq = psth_freq;
psth_struct.psth_mnmx = [min(psth_struct.psth(:)) max(psth_struct.psth(:))];
psth_struct.binwidth = binwidth;
% psth_struct.data_struct = data;  % FIXME
psth_struct.data_orig = data_orig;

if exist( 'data_pres' )
    psth_struct.data_pres = data_pres;
end

psth_struct.audiogram_struct = aud;

%% PR, PH, BOX, HIST and Plots

if ( strncmpi(data_orig.calc_details, 'detailed', 6) )
    
    [psth_struct.F, psth_struct.PSTH, psth_struct.PHASE] = quickfft( psth_struct.psth, 1/binwidth );
    psth_struct.calc_details = 'detailed';
    
    if strcmp(psth_struct.type, 'FINE')
        
        psth_struct.psth_norm = 1/sum(nrep)/binwidth;   % Normalization Coeff.
        %disp(['Calculating the Box plot and the Power Ratio and Phase Responses @ ' num2str(psth_struct.data_struct.approx_formants) ' Hz']);
        %psth_struct.phsr_freq = {}; psth_struct.phsr = {}; psth_struct.pwrr = [];
        %for i = 1:length(psth_struct.data_struct.approx_formants)
        
        % calculate the phase of a partiuclar speech
        % frequency component in the region around a fiber with the same CF
        %[psth_struct.phsr_freq{i}, psth_struct.phsr{i}, psth_struct.pwrr(i,:)] = fd_phsr( psth_struct, psth_struct.data_struct.approx_formants(i) );
        %end
        % FD_BOXP computes the strength of AN fiber phase locking to individual
        % frequency components of a vowel.
        %[ psth_struct.boxp_freq, psth_struct.boxp] = fd_boxp( psth_struct );
        %psth_struct.boxp_mnmx = [min(psth_struct.boxp(:)) max(psth_struct.boxp(:))];
        
    elseif strcmp(psth_struct.type, 'AVG')
        % histogram over rate bins is computed here, but is only used as
        % additional information and not for neurogram comparison
        disp('Calculating the Histogram Response')
        [ psth_struct.hist_bins, psth_struct.hist ] = fd_hist( psth_struct );
        psth_struct.hist_mnmx = [min(psth_struct.hist(:)) max(psth_struct.hist(:))];
        
    end
    
else
    psth_struct.calc_details = 'simple';
end
%% Is a Whole Number?
function [ flag ] = iswhole( number )
flag = (number == round(number));


