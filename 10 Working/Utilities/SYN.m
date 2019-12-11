function [ synout, syn_freq, syn_time ] = SYN( segment,freq,H,FS,num,comp_type);

%   Hamming Window Length: N = (FS*Sec)/(1 - 2*acos(2/23)/(2*pi));

if ~(strncmpi(comp_type, 'Fine',3) | strncmpi(comp_type, 'Avg',3))
    disp('Improper Compare Specification');
    return;
end

% Constants                                                               %
%=========================================================================%
Fs = 500e3;
spont = [50 5 .1];
N = 1895;
hammer = FD_HAMMING(N)/N;
%=========================================================================%

% Interpolates frequnecy/threshold curve, and determines hair cell coef's %
%=========================================================================%
syn_freq = logspace(log10(freq(1)),log10(freq(end)),num);           % Configure number of CF's Here!
H_interp = interp1(freq,H,syn_freq,'cubic');                        % Changing this to linear causes problems

[Cohc,Cihc]=fitaudiogram(syn_freq,H_interp);
%=========================================================================%

% Main loop                                                               %
%=========================================================================%
for i = 1:length(syn_freq)

    nrep = round(50*[0.6 0.2 0.2]);
    
    if nrep(1) < 1
        nrep(1) = 1;
    end
    
    [a,a,a,a,a,a,a,synout_a,a] = zbcatmodel(segment,syn_freq(i),1,1/Fs,length(segment)/Fs,Cohc(i),Cihc(i),spont(1));  % High spont-rate fibers die.
    [a,a,a,a,a,a,a,synout_b,a] = zbcatmodel(segment,syn_freq(i),1,1/Fs,length(segment)/Fs,Cohc(i),Cihc(i),spont(2));  % medium spont-rate remain.
    [a,a,a,a,a,a,a,synout_c,a] = zbcatmodel(segment,syn_freq(i),1,1/Fs,length(segment)/Fs,Cohc(i),Cihc(i),spont(3));  % low spont-rate remain.
    
    synout_cf = nrep(1)*synout_a + nrep(2)*synout_b + nrep(3)*synout_c;
    clear a synout_a synout_b synout_c
    
    if strncmpi(comp_type, 'Avg',3)
        synout(i,:) = conv(hammer, synout_cf);
    else
        synout(i,:) = synout_cf;
    end
    
    synout = single(synout);
    
end

syn_time = [0:(size(synout,2)-1)]/FS;
%=========================================================================%