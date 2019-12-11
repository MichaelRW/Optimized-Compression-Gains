function psth = make_psth_struct(data, FS, spl, loss, pres, CFcount, IOHC_loss, binwidth, synaptopathy)
audi = audiograms(loss, CFcount, IOHC_loss);
% added synaptopathy (050718)
data = make_data_struct( data, FS, spl, audi, pres, 'detailed',synaptopathy); % 'simple' Or 'detailed'

% PSTH may constructs the PSTH of the model fibers at given data, audiogram
% and binwidth
psth = PSTHmay(data, audi, binwidth, 'n');
end