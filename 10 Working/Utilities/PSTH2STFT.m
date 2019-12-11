function [stf, f, t, clim] = PSTH2STFT(input,FS,winlen);
% [stf, f, t] = PSTH2STFT(input,16000,8192);

stf = [];
for i = 1:size(input,1)
    [A,f,t] = specgram(input(i,:),[0:(FS/2)],500000,winlen,winlen/2);
    stf(:,:,i) = abs(A);
    stf = single(stf);
end
clim = maxmin(stf);