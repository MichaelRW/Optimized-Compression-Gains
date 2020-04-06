
function [ data_struct_new ] = ampl_pres( data_struct, audi_struct, pres_type, calc_details )


% AMPL_PRES applies either NAL-R or DSL amplification prescriptions
% 
% [ data_struct_out ] = ampl_pres( data_struct_in, audi_struct, pres_type )
% uses the 'audiograms' data structure to apply appropriate NAL-R or DSL
% hearing aid prescription gains to the data in data_struct_in.
% Amplificiation prescription is specified by setting 'pres_type' as either
% 'NAL' or 'DSL'. The resulting data is contained in the structure
% data_struct_out along with the original data.
%
% See Also make_data_struct


if strcmp(audi_struct.type, 'Normal')&& ~ strcmp(pres_type,'none')
    disp('Warning: You are applying amplification gains to a normal AP.');
end

%%
%=========================================================================%
%                               Variables                                 %
%=========================================================================%
data_length = length(data_struct.data);
window = 256;
shift = 128;
fft_length = 2^(ceil(log2(window))+1);
zero_padding = fft_length - window;
f = -0.5:1/fft_length:0.5; f = f*data_struct.FS; f(end) = [];
%=========================================================================%

if strncmpi(pres_type, 'DSL',3)
    %%
    %=====================================================================%
    %                       DSL Amplification Scheme                      %
    %=====================================================================%
F = [250 500 750 1000 1500 2000 3000 4000 6000 8000];
H = [0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110];
Z = [[0  3  5  7  9  12 14 17 20 22 25 29 32 36 39 43 47 51 55 59 62 66 68]'...
     [2  4  6  8  11 13 15 18 20 23 26 29 32 35 38 42 45 48 52 55 59 62 66]'...
     [3  5  7  10 12 14 17 19 22 25 28 31 34 37 40 43 47 50 54 57 61 64 68]'...
     [3  5  8  10 13 15 18 21 24 27 30 33 36 40 43 46 50 53 57 60 64 68 71]'...
     [5  8  10 13 15 18 20 23 26 29 32 35 38 42 45 48 52 55 59 62 66 70 73]'...
     [12 15 17 19 22 24 27 30 33 36 39 42 46 49 52 56 59 63 66 70 73 77 80]'...
     [16 18 20 23 25 28 30 33 36 39 42 45 48 52 55 59 62 66 69 73 76 80 83]'...
     [14 17 19 21 24 27 29 32 35 38 41 45 48 51 55 58 62 65 69 73 76 80 84]'...
     [8  11 14 17 20 23 26 29 32 36 39 43 46 50 54 58 61 65 69 73 76 81 85]'...
     [6  10 13 16 19 22 25 27 30 35 38 42 45 49 53 58 60 64 68 72 75 80 85]'];

    REAG_dB = interp2(F,H,Z,audi_struct.F,audi_struct.H, 'cubic');
    REAG_SC = 10.^(REAG_dB/20);

    % PLOTTING %
%     plot(audi_struct.F, REAG_dB)
%     title('DSL Real Ear Unadided Amplification Gains')
%     xlabel('Center Frequency (CF)')
%     ylabel('Amplification Gain (dB)')
%     xlim([0 8000]);

    % NO NEED TO APPLY THE HEAD RELATED TRANSFER FUNCTION %
    data = data_struct.data;
    gains = fftshift(interp1([0 audi_struct.F],[1 REAG_SC],abs(f),'cubic'));
    %=====================================================================%
    
elseif strncmpi(pres_type, 'NAL-R', 3)
    %%
    %=====================================================================%
    %                     NAL-R Amplification Scheme                      %
    %=====================================================================%
    
    % The NAL-R formular gives a constant for each F:
    F = [250 500 750 1000 1500 2000 3000 4000 6000 8000];
    k = [-17  -8  -3    1    0   -1   -2   -2   -2   -2];

    % The NAL-R formular averages loss at F = 500, 1000, and 2000 Hz
    H = audi_struct.orig_H;     % This should corrrespond to the F vector.
    H3FA = (H(2) + H(4) + H(6))/3;

    X = 0.15*H3FA;              % Constant
    IG_dB = X + 0.31*H + k;     % Gains
    IG_dB(IG_dB<0) = 0;         % My addition: Remove losses.    
    IG_dB = interp1(F, IG_dB, audi_struct.F, 'pchip');
    IG_SC = 10.^(IG_dB/20);
    
    % PLOTTING %
%     plot(audi_struct.F, IG_dB)
%     title('NAL-R Insertion Amplification Gains')
%     xlabel('Center Frequency (CF)')
%     ylabel('Amplification Gain (dB)')
%     xlim([0 8000]);
    
    % MUST APPLY THE HEAD RELATED TRANSFER FUNCTION %
    data = hrtffilt(data_struct.data,data_struct.FS);
    gains = fftshift(interp1([0 audi_struct.F],[1 IG_SC],abs(f),'pchip'));
    %=====================================================================%
    
else
    disp('No Amplification Prescription Applied.')
    data_struct_new = data_struct;
    return
end

%%
%=========================================================================%
%                      The Zero-Phase Filtering Loop                      %
%=========================================================================%
data = stratify(data, window, 'shift', shift);
for i = 1:size(data,1)
    if i == 1
        w = hanning(window); w(1:ceil(window/2)) = 1;
        x = data(i,:) .* w';
    elseif i == size(data,1)
        w = hanning(window); w((ceil(window/2)+1):end) = 1;
        x = data(i,:) .* w';
    else
        x = data(i,:).*hanning(window)';
    end
    x = [x zeros(1,zero_padding)];
    X = fft(x).*gains;
    x = ifft(X);
    data(i,:) = x(1:window);
end
data = overlapandadd(data,window,shift);
data = data(1:data_length);
%=========================================================================%

data_struct_new = make_data_struct( data, data_struct.FS, get_spl(data), calc_details, data_struct.SYNAPTOPATHY );

if strcmp(calc_details, 'detailed')
    %data_struct_new.ants = data_struct.formants;                % Keep the same
    %data_struct_new.approx_formants = data_struct.approx_formants;  % Keep the same
    data_struct_new.calc_details = 'detailed';
else
    data_struct_new.calc_details = 'simple';
end

% Make fundamental freq the same too?? Seems to stay the same..
data_struct_new.pres_type = pres_type;
data_struct_new.pres_gain = fftshift(gains);
data_struct_new.pres_freq = f;
% Copy over the old FFT's too??
% data_struct_new.data_orig = data_struct.data;
% data_struct_new.SPL_orig = data_struct.SPL;
data_struct_new.data_orig = data_struct;