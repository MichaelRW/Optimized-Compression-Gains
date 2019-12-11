function y = hrtffilt(x,Fs);
% Human HRTF filter based on Wiener and Ross, 1946.
%
% Note: The sampling frequency should be around 20-40 kHz
% for the bilinear transform to be accurate.  Resample the
% input to be within this range if necessary.

FH1 = 2.5e3;
FH2 = 4.2e3;
BH1 = 0.9e3;
BH2 = 2.1e3;
FHLP = 8.0e3;
Hmult = 3.0;

PH1 = 1/(2*pi*FH1)^2;
QH1 = 2*pi*BH1*PH1;

PH2 = 1/(2*pi*FH2)^2;
QH2 = 2*pi*BH2*PH2;

PHLP = 1/(2*pi*FHLP);

[b1,a1] = bilinear([1],[PH1 QH1 1],Fs);
[b2,a2] = bilinear([1],[PH2 QH2 1],Fs);
[blp, alp] = bilinear([1],[PHLP 1],Fs);

x1 = filter(b1,a1,x);
x2 = filter(b2,a2,x);
x3 = filter(blp,alp,x);
y = Hmult*(x1-x2)+x3;