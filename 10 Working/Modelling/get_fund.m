function [ Fx ] = get_fund( data, FS );

[f, X] = newfft(data, FS);
delta_f = f(end) - f(end-1);
f_low = round(50/delta_f);   % Fundamental above 50 Hz
f_high = round(300/delta_f); % Fundamental below 300 Hz

r=xcorr(X,f_high,'coeff');
r=r(f_high+1:2*f_high+1);
d = 0:(f_high);
[fmax,fx]=max(r(f_low:f_high));

if (fmax < 0.5)
    Fx = NaN;
else
    Fx = d(f_low+fx-1)*delta_f;
end