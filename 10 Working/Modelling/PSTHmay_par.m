function [ psth_struct ] = PSTHmay_par(S, aud, binwidth, varargin)
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
%%
%=========================================================================%
%                               Constants                                 %
%=========================================================================%

fs = 100e3  % May change if FS is greater than fs %
ts = 1/fs;
b = double(single(1/binwidth));
%spont = [50 5 .1];
%nrep = round(50*[0.6 0.2 0.2]);%Original!
nrep_healthy = round(50*[0.6 0.2 0.2]);

%nrep = 1;
%aud;

psth_freq = aud.F;

Cohc = aud.Cohc;
Cihc = aud.Cihc;


%%
%=========================================================================%
%                   Check for proper FS and Binwidth                      %
%=========================================================================%

if ( b >= fs) && (S.FS < b )
    type = 'FINE';
    fs = b;
    dat = resample(S.data,fs,S.FS);
    %disp(['Data resampled to fs : ' num2str(fs) 'Hz'])
    %disp('Returning Fine Timing Information')
elseif ( b < fs ) && ( S.FS < fs )
    type = 'AVG';
    if ~iswhole(binwidth*fs)
        disp('Make sure that binwidth*100e3 is a whole number')
        return;
    end
    
    dat = resample(S.data,fs,S.FS);
    %disp(['Data resampled to fs : ' num2str(fs) 'Hz'])
    disp(['Returning Average Discharge Rate With Binwidth of ' num2str(binwidth) 's'])
elseif ( S.FS >= fs ) && ( b < fs )
    type = 'AVG';
    fs = S.FS;
    if ~iswhole(binwidth*fs)
        disp('Make sure that binwidth*FS is a whole number')
        return;
    end
    dat = S.data;
    %disp(['Data presented at : ' num2str(fs) 'Hz'])
    disp(['Returning Average Discharge Rate With Binwidth of ' num2str(binwidth) 's'])
elseif ( S.FS >= b ) && ( b >=fs )
    type = 'FINE';
    fs = S.FS;
    dat = S.data;
    %disp(['Data presented at : ' num2str(fs) 'Hz'])
    disp('Returning Fine Timing Information')
else
    disp('Something is wrong with FS and binwidth.')
    return;
end

data = make_data_struct( dat, fs, S.SPL, S.calc_details );

%if strncmpi(S.calc_details, 'detailed', 6)
    %data.formants = S.formants;                % Keep the same
    %data.approx_formants = S.approx_formants;  % Keep the same
%end

if  isfield(S, 'data_orig') % If there's been amplification.
    data_orig = S.data_orig;
    S = rmfield(S,'data_orig');
    data_pres = S;
else                        % Otherwise there hasn't been amplification.
    data_orig = S;
end

%T = length(dat)*ts;

%%
%=========================================================================%
%                             Filter Design                               %
%=========================================================================%

if strcmp(type, 'AVG')
    Ap = 0.1;       % Passband Ripple
    Ast = 20;      % Amplitude below passband
    Fp = 250*binwidth;    % Cut-off at 250
    Fst = 350*binwidth;   % to 350 Hz
    %disp('Low-Pass Filtering psth with cut-off @ 250 - 350 Hz')
elseif strcmp(type, 'FINE')
    Ap = 0.1;
    Ast = 20;
    Fp = 6000*binwidth;
    Fst = 8000*binwidth;
    %disp('Low-Pass Filtering psth with cut-off @ 6000 - 8000 Hz')
end
% d = fdesign.lowpass(Fp,Fst,Ap,Ast);
% hd = design(d,'butter','matchexactly','passband');
% [h,t] = impz(hd);
[B,A] = butter(6,Fp,'low');
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
  warning('synaptopathy type unknown, healthy AN assumed')      
endif
%%
%=========================================================================%
%                               Main  Loop                                %
%=========================================================================%
    
%=++++++++++++++++++++++++++===% Clip Data %===++++++++++++++++++++++++++=%
binw = round(binwidth*fs);              % How many samples are in each bin.
remainder = mod(length(dat),binw);      % Number of samples to cut off end.
dat = dat(1:end-remainder);             % Make length multiple of binwidth.
len = 2*length(dat)/binw;                 % New length of binned data.

reptime = 2*length(dat)/fs; 
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%


% first approach: use new model with same input params, also for nrep spont

% after that go on ...

    %psth500k_a_temp = zeros(1,length(dat))
    %psth500k_a_single_fibers = zeros(nrep(1),length(dat));
    psth500k_a_single_fibers = [];
    %psth500k_b_temp = zeros(1,length(dat))
    %psth500k_b_single_fibers = zeros(nrep(2),length(dat));
    psth500k_b_single_fibers = [];
    %psth500k_c_temp = zeros(1,length(dat))
    %psth500k_c_single_fibers = zeros(nrep(3),length(dat));
    psth500k_c_single_fibers = [];

%tic 
disp(['length(psth_freq) = ', num2str(length(psth_freq))]);
disp(['size(dat) = ', num2str(size(dat))]);
psth=zeros(2,length(dat));

%try either this in parall:
tic;
% 				       (i,sponts,tabss,trels,nrep,dat,psth_freq,reptime, ts,Cohc,Cihc,binw,len,binwidth,B,A)


psth_c=pararrayfun(30,@(i)psth_freq_fun(i,sponts,tabss,trels,nrep,dat,psth_freq,reptime, ts,Cohc,Cihc,binw,len,binwidth,B,A),[1:length(psth_freq)], 'UniformOutput', false);
psth=cell2mat(psth_c');
toc
size(psth)

nn= nrep




%... or this in serial:
%tic;
%for i = 1:length(psth_freq)
%  psth(i,:)=psth_freq_fun(i,sponts,tabss,trels,nrep,dat,psth_freq,ts,Cohc,Cihc,binw,len,binwidth,B,A);
%end
%toc

%... instead of the original:
%for i = 1:length(psth_freq)
%for i = 1
%     
%    sponts_concat = [sponts.LS(i,1:nrep(1)) sponts.MS(i,1:nrep(2)) sponts.HS(i,1:nrep(3))];
%    tabss_concat = [tabss.LS(i,1:nrep(1)) tabss.MS(i,1:nrep(2)) tabss.HS(i,1:nrep(3))];
%    trels_concat = [trels.LS(i,1:nrep(1)) trels.MS(i,1:nrep(2)) trels.HS(i,1:nrep(3))];
%
%    %[a,a,a,a,a,a,a,a,psth500k_a] = zbcatmodel(dat,psth_freq(i),nrep(1),ts,length(dat)*ts,Cohc(i),Cihc(i),spont(1)); % High sr
%    %[a,a,a,a,a,a,a,a,psth500k_b] = zbcatmodel(dat,psth_freq(i),nrep(2),ts,length(dat)*ts,Cohc(i),Cihc(i),spont(2)); % Medium sr
%    %[a,a,a,a,a,a,a,a,psth500k_c] = zbcatmodel(dat,psth_freq(i),nrep(3),ts,length(dat)*ts,Cohc(i),Cihc(i),spont(3)); % Low sr
%                                  % zbcatmodel(pin,CF          ,nrep,  binwidth,stimtime,  cohc, cihc,   spont);
% 
%                                  
%    %vihc = model_IHC_BEZ2018(pin,CF,nrep,dt,2*T,cohc,cihc,species);
%    %[psth,~,~,~, ~,~] = model_Synapse_BEZ2018(vihc,CF,nrep,dt,noiseType,implnt,spont,tabs,trel);
%    %loop for single fibers
%%    vihc_temp = model_IHC_BEZ2018(dat,psth_freq(i),1,ts,length(dat)*ts,Cohc(i),Cihc(i),1);
%    
%%    figure;hold on;
%    disp(['nrep(1) = ', num2str(nrep(1))]);
%    for ma = 1:nrep(1)
%    spont = sponts_concat(ma);
%    tabs = tabss_concat(ma);
%    trel = trels_concat(ma);
%    
%   % took length(dat)*ts because saving the variable caused an accuracy
%   % error and length for reshape didn't work
%    vihc_temp = model_IHC_BEZ2018(dat,psth_freq(i),1,ts,length(dat)*ts,Cohc(i),Cihc(i),1);
%%    plot(vihc_temp);
%    [psth500k_a_temp,~,~,~,~,~] = model_Synapse_BEZ2018(vihc_temp,psth_freq(i),1,ts,1,0,spont,tabs,trel);
%    
%    psth500k_a_single_fibers(ma,:) = psth500k_a_temp;
%    
%    end
%%    hold off
%
%    psth500k_a = sum(psth500k_a_single_fibers);
%    
%
%%    figure;hold on;
%    disp(['nrep(2) = ', num2str(nrep(2))]);
%    for mb = 1:nrep(2)
%        
%    spont = sponts_concat(mb+nrep(1));
%    tabs = tabss_concat(mb+nrep(1));
%    trel = trels_concat(mb+nrep(1));
%       
%    vihc_temp = model_IHC_BEZ2018(dat,psth_freq(i),1,ts,length(dat)*ts,Cohc(i),Cihc(i),1);
%%    plot(vihc_temp);
%    [psth500k_b_temp,~,~,~, ~,~] = model_Synapse_BEZ2018(vihc_temp,psth_freq(i),1,ts,1,0,spont,tabs,trel);
%    
%    psth500k_b_single_fibers(mb,:) = psth500k_b_temp;
%    
%    end
%%    hold off
%
%    
%    psth500k_b = sum(psth500k_b_single_fibers);
%    
%   
%%    figure;hold on;
%    disp(['nrep(3) = ', num2str(nrep(3))]);
%    for mc = 1:nrep(3)
%    spont = sponts_concat(mc+nrep(1)+nrep(2));
%    tabs = tabss_concat(mc+nrep(1)+nrep(2));
%    trel = trels_concat(mc+nrep(1)+nrep(2));
%        
%        
%    vihc_temp = model_IHC_BEZ2018(dat,psth_freq(i),1,ts,length(dat)*ts,Cohc(i),Cihc(i),1);
%%    plot(vihc_temp);
%    [psth500k_c_temp,~,~,~, ~,~] = model_Synapse_BEZ2018(vihc_temp,psth_freq(i),1,ts,1,0,spont,tabs,trel);
%    
%    psth500k_c_single_fibers(mc,:) = psth500k_c_temp;
%    end
%%    hold off
% 
%    psth500k_c = sum(psth500k_c_single_fibers);
%    
%    
%    psth500k = psth500k_a + psth500k_b + psth500k_c;
%    clear psth500k_a psth500k_b psth500k_c
%
%    pr = sum(reshape(psth500k,binw,len),1)/sum(nrep)/binwidth; % psth in units of spikes/s/fiber
%
%    
%    %Psth = sum(reshape(psth,psthbins,length(psth)/psthbins)); %
%    %     pr = filter(hd, pr);
%    pr = filtfilt(B,A,pr);
%    psth(i,:) = single(pr);
%end
%toc
psth_struct.type = type;
psth_struct.psth = psth(:,1:len);
psth_struct.psth_time = [0:(len-1)]*binwidth;
psth_struct.psth_freq = psth_freq;
psth_struct.psth_mnmx = [min(psth_struct.psth(:)) max(psth_struct.psth(:))];
psth_struct.binwidth = binwidth;
psth_struct.data_struct = data;
psth_struct.data_orig = data_orig;
if exist('data_pres')
    psth_struct.data_pres = data_pres;
end
psth_struct.audiogram_struct = aud;


%%
%=========================================================================%
%                       PR, PH, BOX, HIST & Plots                         %
%=========================================================================%
if strncmpi(data_orig.calc_details, 'detailed', 6)
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
        psth_plot( psth_struct );
    end
else
    psth_struct.calc_details = 'simple';
end


%%
%=========================================================================%
%                          Is a Whole Number?                             %
%=========================================================================%
function [ flag ] = iswhole( number )
flag = (number == round(number));

