function psth = make_psth_struct_par(data, FS, spl, loss, pres, CFcount, IOHC_loss, binwidth, synaptopathy)
audi = audiograms(loss, CFcount, IOHC_loss);
data = make_data_struct( data, FS, spl, audi, pres, 'detailed', synaptopathy); % 'simple' Or 'detailed'

psth = PSTHmay_par(data, audi, binwidth, 'n');
end