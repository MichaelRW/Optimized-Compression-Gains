function [f, X, P] = quickfft( x, Fs )
% [ f, X, P ] = quickfft( Data, Fs );
% Calculates horizontal fft
% 
% Input:
% Data
% Fs = Sampling Frequency
% 
% Output:
% f = Frequency Axis
% X = Magnitude Spectrum
% P = Phase Spectrum in degrees
% 
% Faheem Dinath. November 28th 2007
len = size(x,2);
Pts = ceil((len+1)/2);

X = fft(x,[],2)/len;

P = angle(X);                   % Phase Response
P = unwrap(P,[],2)*180/pi; 
P=P(:,1:Pts);

X = abs(X);                     % Magnitude Response.
X=X(:,1:Pts);

f=(0:Pts-1)*Fs/len;             % Frequecy Axis
f = round(f*1000)/1000;

% if nargin == 2
%     plt = 'no';
% end
% 
% if strncmpi(plt, 'yes',1)
%     subplot(2,1,1), plot(f,X), 
%     ylabel('Abs. Magnitude'), grid on
%     subplot(2,1,2), plot(f,P)
%     ylabel('Phase [Degrees]'), grid on
%     xlabel('Frequency [Hertz]')
% end