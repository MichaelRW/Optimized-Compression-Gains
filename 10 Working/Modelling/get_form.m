function [formants] = get_form(data,FS)

% get Linear prediction filter
ncoeff=2+FS/1000;           % rule of thumb for formant estimation
#a=lpc(data,ncoeff);


a = lpc_octave(data,ncoeff);
% find frequencies by root-solving
r= roots (a);                  % find roots of polynomial a
r=r(imag(r)>0.01);           % only look for roots >0Hz up to fs/2
ffreq=sort(atan2(imag(r),real(r))*FS/(2*pi));

% % plot frequency response
[h,f]=freqz(1,a,2^16,FS);
h = 20*log10(abs(h));
% plot(f,h+eps);
% legend('LP Filter');
% xlabel('Frequency (Hz)');
% ylabel('Gain (dB)');


[XMAX,IMAX,XMIN,IMIN] = extrema(h);
formants = sortrows([f(IMAX) XMAX]);

% Remove formants less than 100 Hz
if formants(1,1) < 100
    formants(1,:) = [];
end

formants = formants(1:3,1)';