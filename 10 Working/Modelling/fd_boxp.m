function [ f, ft_f0 ] = fd_boxp(psth_struct, varargin);
% FD_BOXP computes the strength of AN fiber phase locking to individual
% frequency components of a vowel.
% 
% [ f, box_data] = fd_boxp(Neurogram_Strucutre) receives a neurogram data
% structure which it uses to compute the box plot. An 80ms Hamming window
% is isolates a portion of the neurogram in 20 ms from the stimulus onset.
% An offset is requied to avoid any adaptation irregularities and a
% windowed response is taken because it is assumed that the stimulus is
% perioid, ie, a vowel. The Fourier transform of the windowed neurogram is
% computed and normalized to express the transform components in units of
% spikes/s. The magnitude of the resulting transform components are the
% synchonized rates.
%  
% REF: Roger L. Miller et al. Effects of acoustic trauma on the
% representation of the vowel /e/ in cat auditory nerve fibers. J. Acoust
% Soc. Am. 1997.

%disp('fd_boxp. Faheem Dinath. May 29th 2008.')

% Make sure that we're looking at fine timing PSTH
if ~strcmp(psth_struct.type, 'FINE')
    f = [];
    ft_f0 = [];
    return
end

%=========================================================================%
if (nargout == 0) | strncmpi(varargin,'y',1)
    h1 = figure;
    h3 = fill([0.1 10 10 0.1],[0.1*2^(1/2) 10*2^(1/2) 10/2^(1/2) 0.1/2^(1/2)],0.9*ones(1,3));
    set(h3,'edgecolor','none')
    set(gca,'xscale','log','yscale','log')
    xlabel('BF (kHz)')
    ylabel('Frequency (kHz)')
    axis([0.1 10 0.095 4.5])
    set(gca,'xticklabel','0.1|1.0|10')
    set(gca,'yticklabel','0.1|1.0')
    hold on
    plot([0.5 0.5],[0.095 4.5],'k--')
    plot([1.7 1.7],[0.095 4.5],'k--')
    plot([2.5 2.5],[0.095 4.5],'k--')
    plot([0.1 10],[0.5 0.5],'k--')
    plot([0.1 10],[1.7 1.7],'k--')
    plot([0.1 10],[2.5 2.5],'k--')
    wysiwyg
end
%=========================================================================%

psth_freq = psth_struct.psth_freq;
psth_time = psth_struct.psth_time;
binwidth = psth_struct.binwidth;

for cfslp=1:length(psth_freq)
        
	cf=psth_freq(cfslp);
    
%     disp(['CF = ' num2str(cf) ' (' int2str(cfslp) '/' int2str(length(psth_freq)) ')'])
    
    if psth_time(end) < 0.1
        disp('Data run-time too short; must be greater than 100 ms');
        return
    end
    
    [dummyy, tonset] =min(abs(psth_time-20e-3));

    toffset = tonset+ 80e-3/binwidth;
    toffset = round(toffset) - 1;
     
%     Hamming window starting 20ms after STIMULUS onset
	p = psth_struct.psth(cfslp, tonset:toffset);
	w = hamming(length(p))';

    [f, MX, P] = quickfft( w.*p, 1/binwidth );
    ft(cfslp,:)=MX/sqrt(sum(w.^2)/length(w));
    
    f0 = psth_struct.data_struct.fundamental;
   
    [mn f_min] = min(abs(f - f0));
    [mn f_max] = min(abs(f - 4e3));
    f_ind = f_min:f_min-1:f_max;
    [ freq, val ] = max_value_around( psth_struct.data_struct, f(f_ind));
    
    for i = 1:length(freq)
        [v add] = min( abs( f - freq(i) ) );
        f_ind(i) = add;
    end

%     f_ind = find(f==100):find(f==100)-1:find(f==4e3);

    if (nargout == 0) | strncmpi(varargin,'y',1)
        figure(h1)
        h2 = loglog(cf/1e3,f(f_ind)/1e3,'ks','markerfacecolor','k');
        for lp=1:length(h2)
            if ft(cfslp,f_ind(lp))<15
                set(h2(lp),'marker','none')
            elseif ft(cfslp,f_ind(lp))<30
                set(h2(lp),'markersize',2)
            elseif ft(cfslp,f_ind(lp))<45
                set(h2(lp),'markersize',4)
            elseif ft(cfslp,f_ind(lp))<60
                set(h2(lp),'markersize',6)
            elseif ft(cfslp,f_ind(lp))<75
                set(h2(lp),'markersize',8)
            elseif ft(cfslp,f_ind(lp))<90
                set(h2(lp),'markersize',10)
            elseif ft(cfslp,f_ind(lp))<105
                set(h2(lp),'markersize',12)
            elseif ft(cfslp,f_ind(lp))<120
                set(h2(lp),'markersize',14)
            else
                set(h2(lp),'markersize',16)
            end
        end
    end
    
    ft_f0(cfslp,:) = ft(cfslp,f_ind);
end

ft_f0(ft_f0<=15) = 0;
ft_f0((ft_f0>15) & (ft_f0<=30)) = 1;
ft_f0((ft_f0>30) & (ft_f0<=45)) = 2;
ft_f0((ft_f0>45) & (ft_f0<=60)) = 3;
ft_f0((ft_f0>60) & (ft_f0<=75)) = 4;
ft_f0((ft_f0>75) & (ft_f0<=90)) = 5;
ft_f0((ft_f0>90) & (ft_f0<=105)) = 6;
ft_f0((ft_f0>105) & (ft_f0<=120)) = 7;
ft_f0(ft_f0>120) = 8;
f = f(f_ind);