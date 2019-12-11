function [f, MX]=fourier(x,Fs)

% [MX, f]=fourier(x,Fs)
%
% MX is the scaled Fourier transform, f is the array of scaled frequencies,
% x is the input waveform, and Fs is the sampling frequency.
%
% From the MathWorks Web Site, 1998
%

Fn=Fs/2;						% Nyquist frequency
% Next highest power of 2 greater than or equal to length(x):
%NFFT = 2^(nextpow2(length(x)));
NFFT = length(x); % Modified by Ian Bruce, 1998 Aug 18
% Take fft, padding with zeros, length(FFTX)==NFFT
FFTX=fft(x,NFFT);
NumUniquePts = ceil((NFFT+1)/2);
% fft is symmetric, throw away second half
FFTX=FFTX(1:NumUniquePts);
MX=abs(FFTX);            % Take magnitude of X
% Multiply by 2 to take into account the fact that we
% threw out second half of FFTX above
MX=MX*2;
MX(1)=MX(1)/2;   % Account for endpoint uniqueness
%MX(length(MX))=MX(length(MX))/2;  	% We know NFFT is even
if ~rem(NFFT,2)							% We DON'T know NFFT is even
   MX(length(MX))=MX(length(MX))/2;
end
% Scale the FFT so that it is not a function of the
% length of x.
MX=MX/length(x);
f=(0:NumUniquePts-1)*2*Fn/NFFT;
