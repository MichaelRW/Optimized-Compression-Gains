function [freq, phase, pr] = fd_phsr(psth_struct, CF, varargin)
% FD_PHSR computes the phase response.
% 
% [freq, phase, pr] = fd_phsr(Neurogram_Struct, CF, varargin)
% uses the neurogram to calculate the phase of a partiuclar speech
% frequency component in the region around a fiber with the same CF. The
% region is determined using the power ratio, as a continuois adjacent
% segment where the power ratio is greater than 0.1. This function
% calculates both the phase and power ratio repsonse for the specified CF.

% Make sure that we're looking at fine timing PSTH
if ~strcmp(psth_struct.type, 'FINE')
    freq = [];
    phase = [];
    return
end

pr = fd_pwrr(psth_struct, CF, varargin);

%disp('fd_phsr. Faheem Dinath. May 29th 2008.')

psth_freq = psth_struct.psth_freq;

% / Find Continuous Region Greater than 0.1 \
%=========================================================================%
[temp pos] = min(abs(psth_freq - CF));

if pr(pos) <=0.1
    freq = [];
    phase = [];
    return;
end

BW = pr>0.1;
[BW,num] = bwlabel(BW);
num = BW(pos);
BW = find(BW == num);     % Binary Vector indicating where region lie.

%=========================================================================%

for i = 1:length(BW)
    [temp pos] = min( abs( psth_struct.F - CF ) );
    phase(i) = psth_struct.PHASE(BW(i), pos);
end

% / Unwrap Phase \
%=========================================================================%
freq = psth_freq(BW);
[temp pos] = min(abs(freq-CF));
phase = phase*2*pi/180; % Convert to Radians
phase = mod(phase,2*pi);
phase = unwrap(phase);
phase = phase - phase(pos);
phase = phase*180/2/pi; % Convert to Degrees
%=========================================================================%

% / Plot PR and Phase Response \
%=========================================================================%
if (nargout == 0) | strncmpi(varargin,'y',1)
    
%     figure;
%     plot(gca,psth_struct.psth_freq,pr);
%     drawnow;
%     xlabel('Centre Frequencies (Hz)')
%     title(['Power Ratio at ' num2str(CF) ' Hz'])
%     ylim([0 1]);
%     
%     figure;
%     plot(gca,freq,phase);
%     drawnow;
%     ylabel('Relative Phase Shift (deg)')
%     xlabel('Centre Frequencies (Hz)')
%     title(['Phase Response at ' num2str(CF) ' Hz'])
%     
end
%=========================================================================%