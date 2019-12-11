function [ phoneme,phn,type,start,stop ] = import_phonemes( phonemefile, wavefile )
% [ phoneme,phn,type,start,stop ] = import_phonemes( phonemefile, wavefile )
[data,FS] = readsph(wavefile);
[start,stop,phn] = textread(phonemefile,'%f%f%s');
start = start + 1;

for i = 1:length(start)
    type{i} = phoneme_lut(phn(i));
    phoneme{i} = data(start(i):stop(i));
end