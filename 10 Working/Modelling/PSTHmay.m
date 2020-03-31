
function [ psth_struct ] = PSTHmay( S, aud, binwidth, varargin )

%PSTHmay Computes the AN Neurogram.
%
% [ Neurogram_Strucutre ] =
% PSTHmay(data_structure, audi_structure, binwidth, plot) 
%
% This function computes the Neurogram at 'binwidth' resolution using the
% spike timing information from the Zilany-Bruce Cat Auditory Model. The AN
% response is taken at 'CF_Points' different center frequency locations,
% spaced logarithmically along the basilar membrane. At each CF point, the
% function calculates the summed response of 50 nerve fibers organized as
% 60%, high, 20% medium, and 20% low spontaneous rate fibers. The Audiogram
% structure specifies the hearing loss profile, as +dB loss, at specific
% frequencies. In addition to calculating the Neurogram, phase response,
% power ratio, box plots, and histogram analysis are calculated where
% possible. An optional potting parameter, if 'y' plots the Neurogram.
% 
% See Also make_data_struct, audiograms, fd_phsr, fd_boxp, fd_hist,
% psth_plot.

%disp('PSTHmay. Faheem Dinath. June 7th 2008.')



%% Constants

fs = 100e3;  % May change if FS is greater than fs %
ts = 1/fs;

b = double(single(1/binwidth));

nrep_healthy = round( 50*[0.6 0.2 0.2] );

psth_freq = aud.F;

Cohc = aud.Cohc;  Cihc = aud.Cihc;



%% Check for proper FS and Binwidth

if ( b >= fs) && (S.FS < b )
    type = 'FINE';
    fs = b;
    dat = resample(S.data,fs,S.FS);
    T = length(dat)/fs;
elseif ( b < fs ) && ( S.FS < fs )
    type = 'AVG';
    if ~iswhole(binwidth*fs)
        disp('Make sure that binwidth*100e3 is a whole number')
        return;
    end
    
    dat = resample(S.data,fs,S.FS);
    T = length(dat)/fs;
    disp(['Returning Average Discharge Rate With Binwidth of ' num2str(binwidth) 's'])
elseif ( S.FS >= fs ) && ( b < fs )
    type = 'AVG';
    fs = S.FS;
    if ~iswhole(binwidth*fs)
        disp('Make sure that binwidth*FS is a whole number')
        return;
    end
    dat = S.data;
    T = length(dat)/fs;
    disp(['Returning Average Discharge Rate With Binwidth of ' num2str(binwidth) 's'])
elseif ( S.FS >= b ) && ( b >=fs )
    type = 'FINE';
    fs = S.FS;
    dat = S.data;
    T = length(dat)/fs;
    disp('Returning Fine Timing Information')
else
    disp('Something is wrong with FS and binwidth.')
    return;
end

% data = make_data_struct( dat, fs, S.SPL, S.calc_details );  % FIXME?

if  isfield(S, 'data_orig') % If there's been amplification.
    data_orig = S.data_orig;
    S = rmfield(S,'data_orig');
    data_pres = S;
else                        % Otherwise there hasn't been amplification.
    data_orig = S;
end



%% Filter Design

if strcmp(type, 'AVG')
%    Ap = 0.1;       % Passband Ripple
%    Ast = 20;      % Amplitude below passband
%    Fp = 250*binwidth;    % Cut-off at 250
%    Fst = 350*binwidth;   % to 350 Hz
%    %disp('Low-Pass Filtering psth with cut-off @ 250 - 350 Hz')
elseif strcmp(type, 'FINE')
%    Ap = 0.1;
%    Ast = 20;
%    Fp = 6000*binwidth;
%    Fst = 8000*binwidth;
%    %disp('Low-Pass Filtering psth with cut-off @ 6000 - 8000 Hz')
end
% d = fdesign.lowpass(Fp,Fst,Ap,Ast);
% hd = design(d,'butter','matchexactly','passband');
% [h,t] = impz(hd);
%[B,A] = butter(6,Fp,'low');
% freqz(B,A)



%%

if exist('ANpopulation.mat','file')
    load('ANpopulation.mat');
    %disp('Loading existing population of AN fibers saved in ANpopulation.mat')
    if (size(sponts.LS,2)<nrep_healthy(1))||(size(sponts.MS,2)<nrep_healthy(2))||(size(sponts.HS,2)<nrep_healthy(3))||(size(sponts.HS,1)<length(psth_freq)||~exist('tabss','var'))
        disp('Saved population of AN fibers in ANpopulation.mat is too small - generating a new population');
        [sponts,tabss,trels] = generateANpopulation(length(psth_freq),nrep_healthy);
    end
else
    [sponts,tabss,trels] = generateANpopulation(length(psth_freq),nrep_healthy);
    disp('Generating population of AN fibers, saved in ANpopulation.mat')
end

if strcmp(S.SYNAPTOPATHY, 'healthy')
    nrep = round([1 1 1].*nrep_healthy); % Healthy AN
elseif strcmp(S.SYNAPTOPATHY, '50all')
    nrep = round([0.5 0.5 0.5].*nrep_healthy); % 50% fiber loss of all types
elseif strcmp(S.SYNAPTOPATHY, 'lossall_low')
    nrep = round([0 1 1].*nrep_healthy); % Loss of all LS fibers
elseif strcmp(S.SYNAPTOPATHY, 'IHCproportional')
    nrep = round([mean(aud.Cihc) 1 mean(aud.Cihc)].*nrep_healthy); % loss of LS and HS fibers proportional to IHC impairment
else
    warning('synaptopathy type unknown, healthy AN assumed');
end



%% Main Loop
    
binw = round(binwidth*fs);             % How many samples are in each bin.


    remainder = mod(length(dat),binw);      % Number of samples to cut off end.
    
dat = dat(1:end-remainder);             % Make length multiple of binwidth.

simdur = ceil( T * 1.2 / binwidth) * binwidth;
len_fact = 1.2;
    len = ceil(len_fact*length(dat)/binw);              % New length of binned data.

%simdur = len_fact*length(dat)/fs; 


		% lines 52 to 56
		    windur_ft = 32;  % Size of window for fine-timing neurogram.
        smw_ft = hamming(windur_ft);
     
        windur_mr = 128;  % Size of window for mean-rate neurogram.
        smw_mr = hamming(windur_mr);
        
        
##        figure; plot(smw_ft)
##        
##        figure; plot(smw_mr)
        
t=1
% Pre-allocation of a) LSR, b) MSR, c) HSR fibers
%psth500k_a_single_fibers = [];  psth500k_b_single_fibers = [];  psth500k_c_single_fibers = [];


vihc_temp = model_IHC_BEZ2018( dat, psth_freq(1), 1, ts, simdur, Cohc(1), Cihc(1), 2 );
    
psth500k_a_single_fibers = NaN( length(psth_freq), length(vihc_temp) );
psth500k_b_single_fibers = NaN( length(psth_freq), length(vihc_temp) );
psth500k_c_single_fibers = NaN( length(psth_freq), length(vihc_temp) ); 
 dat_len = length(dat)
%figure; hold on
for i = 1:length(psth_freq)
    
    fprintf( 1, '\n\tCF: %d of %d - %d Hz', i, numel(psth_freq), psth_freq(i) );
    
    sponts_concat = [sponts.LS(i,1:nrep(1)) sponts.MS(i,1:nrep(2)) sponts.HS(i,1:nrep(3))];
    tabss_concat = [tabss.LS(i,1:nrep(1)) tabss.MS(i,1:nrep(2)) tabss.HS(i,1:nrep(3))];
    trels_concat = [trels.LS(i,1:nrep(1)) trels.MS(i,1:nrep(2)) trels.HS(i,1:nrep(3))];
    
    % Generate inner hair cell response for each frequency.
    vihc_temp = model_IHC_BEZ2018( dat, psth_freq(i), 1, ts, simdur, Cohc(i), Cihc(i), 2 );
   
    
    % Accumulate synapse responses for low spontaneous fibers.
%     fprintf( 1, '\n\t\tComputing low spont. fiber responses.' );
    %
    for ma = 1:nrep(1)        
        spont = sponts_concat(ma);  tabs = tabss_concat(ma);  trel = trels_concat(ma);
        
        % Took length(dat)*ts because saving the variable caused an accuracy error and length for reshape did not work.
        [ psth500k_a_temp, ~, ~, ~, ~, ~ ] = model_Synapse_BEZ2018( vihc_temp, psth_freq(i), 1, ts, 1, 0, spont, tabs, trel );
            psth500k_a_single_fibers(ma, :) = psth500k_a_temp;
    end
    %
    psth500k_a = sum(psth500k_a_single_fibers);
    
    
    
    % Accumulate synapse responses for medium spontaneous fibers.
%     fprintf( 1, '\n\t\tComputing medium spont. fiber responses.' );
    %
    for mb = 1:nrep(2)        
        spont = sponts_concat(mb+nrep(1));  tabs = tabss_concat(mb+nrep(1));  trel = trels_concat(mb+nrep(1));        
        
        [ psth500k_b_temp, ~, ~, ~, ~, ~ ] = model_Synapse_BEZ2018( vihc_temp, psth_freq(i), 1, ts, 1, 0, spont, tabs, trel );
            psth500k_b_single_fibers(mb, :) = psth500k_b_temp;
    end
    %
    psth500k_b = sum(psth500k_b_single_fibers);
    
   
    % Accumulate synapse responses for high spontaneous fibers.
%     fprintf( 1, '\n\t\tComputing high spont. fiber responses.' );
    %
    for mc = 1:nrep(3)
        spont = sponts_concat(mc+nrep(1)+nrep(2));  tabs = tabss_concat(mc+nrep(1)+nrep(2));  trel = trels_concat(mc+nrep(1)+nrep(2));        
        
        [ psth500k_c_temp, ~, ~, ~, ~, ~ ] = model_Synapse_BEZ2018( vihc_temp, psth_freq(i), 1, ts, 1, 0, spont, tabs, trel );
            psth500k_c_single_fibers(mc, :) = psth500k_c_temp;
    end
    
    %figure; plot(psth500k_c_temp)
    
    %
    psth500k_c = sum(psth500k_c_single_fibers);
    
    %figure; plot(psth500k_c)
    
    
    
    len_psth = floor(length(psth500k_a)/10)*10;
    
    psth500k_a = psth500k_a(1: len_psth);
    psth500k_b = psth500k_b(1: len_psth);
    psth500k_c = psth500k_c(1: len_psth);
    
    
    
    
    
    %figure; plot(psth500k_c)
    
    
    psth500k = psth500k_a + psth500k_b + psth500k_c;  % Complete PSTH response.
        clear psth500k_a psth500k_b psth500k_c

    
    
        
    %pr = sum( reshape(psth500k, binw, len), 1 ) / sum(nrep) / binwidth; % psth in units of spikes/s/fiber
    %psth_sec = sum( psth500k, 1 ) / sum(nrep) * fs;
    
   %% l_sec = length(psth_sec)
	%figure; plot(psth_sec)
	
  
##   psth_len =   length(psth_ft)
## psthb = binw
## resh_fact = length(psth_ft) / binw
  
  
  
  ft_binw = round(10e-6*fs)
  mr_binw = binw
  
  psth_ft =      psth500k; %round(10e-6*fs);%sum( reshape(psth500k, round(10e-6*fs), length(psth500k) / round(10e-6*fs) ));
	psth_mr =      sum( reshape(psth500k, mr_binw,            length(psth500k) / mr_binw ));
  
  ps_mr_size = size(psth_mr)
	
  %figure; plot(psth_ft, 'x')
	
	
        neurogram_ft = psth_ft + filter( smw_ft, 1, psth_ft );
        neurogram_mr = psth_mr + filter( smw_mr, 1, psth_mr );
        
##   figure; plot(neurogram_ft, 'x')     
##     
##   figure; plot(neurogram_mr, 'o')     
        
        
        
        sz_ng_mr = size(neurogram_mr)
        
        %neurogram_Sout = neurogram_Sout + synout;  

	% ?
    %neurogram_ft = cell2mat( neurogram_ft' );    
    %neurogram_mr = cell2mat( neurogram_mr' );
    %neurogram_Sout = cell2mat( neurogram_Sout' );
    % ?
	
  
  % make sure neurogram is a column vector
  neurogram_ft = neurogram_ft';
  
  neurogram_mr = neurogram_mr';
  
  
    neurogram_ft = neurogram_ft(:, 1:windur_ft/2:end ); % this should be the overlap, I guess?????
    t_ft = 0:windur_ft/2/fs:( size( neurogram_ft, 2 ) - 1 ) * windur_ft / 2 / fs;        
    
    
    
    
    
    neurogram_mr = neurogram_mr(:, 1:windur_mr/2:end );
    t_mr = 0:windur_mr/2*(binw):( length(neurogram_mr) - 1 ) * windur_mr/2 * (binw);

size(neurogram_mr)
t_mr_disp = t_mr(end)   
    sz_ng_mr_after = size(neurogram_mr)

    
    %t_Sout = 0:1/Fs:( size(neurogram_Sout, 2) - 1 ) / Fs;    

    if strcmp(type, 'FINE')==1
      pr = neurogram_ft;
      %figure; plot(pr, 'o')
      
      % ... = t_ft
    elseif strcmp(type, 'AVG')==1
      t_mr = 0:windur_mr/2*(1/fs):( length(neurogram_mr) - 1 ) * windur_mr/2 * (1/fs);
      pr = neurogram_mr;
      len_pr = size(pr)
      figure; hold on; plot(pr, '-')
      xlabel('bin')
      ylabel('average spike count')
      hold off;
      
    end
    
    
    
    t=1
    
    %    pr = filtfilt( B, A, pr );
            psth(i, :) = single(pr);
            
% 	fprintf(1, '\n');
    
end

fprintf(1, '\n\n');


psth_struct.type = type;
psth_struct.psth = psth(:,1:length(psth));
psth_struct.psth_time = [0:(length(psth)-1)]*binwidth;

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
        
        disp('Calculating the Histogram Response')
        [ psth_struct.hist_bins, psth_struct.hist ] = fd_hist( psth_struct );
        psth_struct.hist_mnmx = [min(psth_struct.hist(:)) max(psth_struct.hist(:))];
        
    end

    
    if strncmpi(varargin{1},'y',1)
        %psth_plot( psth_struct ); % check again
    end
    
else
    psth_struct.calc_details = 'simple';
end



%% Is a Whole Number?

function [ flag ] = iswhole( number )
flag = (number == round(number));


