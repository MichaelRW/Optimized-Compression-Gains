function pr = fd_pwrr(psth_struct, fx, varargin)
%[PR, H] = fd_pwrr(Neurogram_Struct, fx) computes the
%strength of AN fiber phase locking to individual frequency components of a
%vowel.
% 
%   fdpowerratio recieves a peritimulus time historgram, psth, its 
%   corresponding center frequency axis, psth_freq, and the binwidth used
%   to compute the PSTH from the spike timing information. The variable fx
%   is the frequency component of the vowel to which phase locking is
%   measured. fdpowerratio returns the handle of the figure plot, H, and
%   figure data, PR.
%
%   The power ratio is defined as the sum of power in the AN fft response
%   at the frequency fx and its harmonics, divided by the total power in
%   the response. Because phase locking in cats is not observed above 5kHz,
%   the summations are limited to frequency componenets below 5kHz.
%   
%  
%                     u
%                    --- 
%                    \     ( R^2(m*fx) )
%                    /
%                    ---
%                    m=1
%   PR(fx) =      ----------------------
%                     v
%                    --- 
%                    \     ( R^2(n*f0) )
%                    /
%                    ---
%                    n=1
% 
%   with u<4 and (u*fx) <= 5kHz
%   and v=50 and (v*f0) <= 5kHz
% 
%   NOTE: fdpowerration assumes that vowels have fundamental harmonics, f0,
%   equalto 100hz and that fx is a multiple of the fundamental harmonic.
% 
%   REF: Roger L. Miller et al. Effects of acoustic trauma on the
%   representation of the vowel /e/ in cat auditory nerve fibers.
%   J. Acoust Soc. Am. 1997.

%disp('fd_pwrr. Faheem Dinath. May 29th 2008.')

% Make sure that we're looking at fine timing PSTH
if ~strcmp(psth_struct.type, 'FINE')
    pr = [];
    return
end

f0 = psth_struct.data_struct.fundamental;
f0_har = floor(5000/f0);

fx_har = floor(5000/fx);
if fx_har > 3
    fx_har = 3;
end


for i = 1:f0_har
    [ freq_f0(i), val ] = max_value_around( psth_struct.data_struct, i*f0);
end

for i = 1:fx_har
    [ freq_fx(i), val ] = max_value_around( psth_struct.data_struct, i*fx);
end

% [PSTH_F, PSTH, PSTH_P] = quickfft( psth, 1/binwidth );
% PSTH = abs(PSTH).^2;    % Looking at the Power!

PSTH_F = psth_struct.F;
PSTH = abs(psth_struct.PSTH).^2;

for i = 1:size(PSTH,1)
    
    sum_n(i) = 0;
    for j = 1:fx_har
        [m1 m2] = min( abs( PSTH_F - freq_fx(j) ) );
        sum_n(i) = sum_n(i) + PSTH(i,m2);
    end
    
    sum_d(i) = 0;
    for j = 1:f0_har
        [m1 m2] = min( abs( PSTH_F - freq_f0(j) ) );
        sum_d(i) = sum_d(i) + PSTH(i,m2);
    end
    
end

pr = sum_n./sum_d;

if (nargout == 0) | strncmpi(varargin,'y',1)
%     figure;
%     plot(gca,psth_struct.psth_freq,pr);
%     drawnow;
%     xlabel('Centre Frequencies (Hz)')
%     title(['Power Ratio at ' num2str(fx) ' Hz'])
%     ylim([0 1]);
end

%=========================================================================%