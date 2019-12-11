function [ psth, psth_freq, psth_time, clim ] = PSTH( segment,freq,H,FS,num,psth_compare_type);

% Constants                                                               %
%=========================================================================%
% hammer = FD_HAMMING(128)/128;
% Fs = 500e3;

hammer = hamming(128)/128;
Fs = 500e3;

spont = [50 5 .1];
%=========================================================================%

% Variables for determining the average discharge rate response           %
%=========================================================================%
psthbinwidth = 1/FS;
psthbins = round(psthbinwidth*Fs);  % number of psth500k bins per psth bin
remainder = rem(length(segment),psthbins);
%=========================================================================%


% Interpolates frequnecy/threshold curve, and determines hair cell coef's %
%=========================================================================%
psth_freq = logspace(log10(freq(1)),log10(freq(end)),num);           % Configure number of CF's Here!
H_interp = interp1(freq,H,psth_freq,'cubic');                        % Changing this to linear causes problems

[Cohc,Cihc]=fitaudiogram(psth_freq,H_interp);
%=========================================================================%


% Main loop                                                               %
%=========================================================================%
for i = 1:length(psth_freq)
%     nrep = ceil(50*[0.6*Cihc(i) 0.2 0.2]);
    nrep = round(50*[0.6 0.2 0.2]);
    
    if nrep(1) < 1
        nrep(1) = 1;
    end
    
    [a,a,a,a,a,a,a,a,psth500k_a] = zbcatmodel(segment,psth_freq(i),nrep(1),1/Fs,length(segment)/Fs,Cohc(i),Cihc(i),spont(1)); % High spont-rate fibers die.
    [a,a,a,a,a,a,a,a,psth500k_b] = zbcatmodel(segment,psth_freq(i),nrep(2),1/Fs,length(segment)/Fs,Cohc(i),Cihc(i),spont(2));   % medium spont-rate remain.
    [a,a,a,a,a,a,a,a,psth500k_c] = zbcatmodel(segment,psth_freq(i),nrep(3),1/Fs,length(segment)/Fs,Cohc(i),Cihc(i),spont(3));   % low spont-rate remain.
    
    psth500k = psth500k_a + psth500k_b + psth500k_c;
    clear a psth500k_a psth500k_b psth500k_c

    if strncmpi(psth_compare_type, 'Avg',3)
        psth500k = psth500k(1:end-remainder);
        pr = sum(reshape(psth500k,psthbins,length(psth500k)/psthbins))/sum(nrep); % pr of spike in each bin
        len = length(pr);
        psth(i,:) = conv(hammer, pr/psthbinwidth); % psth in units of spikes/s
        psth_time = [0:(length(pr)-1)]/FS;
    elseif strncmpi(psth_compare_type, 'Fine',3)
        psth(i,:) = conv(hammer, psth500k);
        len = length(psth500k);
        psth_time = [0:(length(psth500k)-1)]/500000;
    else
        disp('Improper Save Specification');
        return;
    end
    psth = single(psth);
end
psth = psth(:,1:len);
clim = maxmin(psth);
%=========================================================================%