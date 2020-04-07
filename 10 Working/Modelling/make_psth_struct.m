
function [ psth ] = make_psth_struct( data, FS, spl, loss, pres, CFcount, IOHC_loss, binwidth, synaptopathy )
% function to create and fill:
%     audiogram struct audi 
%     data struct data
%     psth struct psth
%
% audiograms() creates struct audi and computes the cihc and cohc
% corresponding to the audiogram features specified by the input type of hearing
% loss (1-15), and the input IOHC_loss ('IHC only, OHC only, or Mixed')% 
audi = audiograms( loss, CFcount, IOHC_loss );

% Added synaptopathy (050718).
data = make_data_struct( data, FS, spl, audi, pres, 'detailed', synaptopathy); % 'simple' Or 'detailed'

% PSTH may constructs the PSTH of the model fibers at given data, audiogram and binwidth.
psth = PSTHmay( data, audi, binwidth, 'n' );


end


