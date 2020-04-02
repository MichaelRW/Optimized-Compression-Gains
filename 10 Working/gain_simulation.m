
function [ collector ] = gain_simulation( data_file, spl, adj, loss, pres, CFcount, IOHC_loss, binwidth, synaptopathy )

%GAIN_SIMULATION
% function to prepare input signals (phoneme chains) and receive optimal
% gains using
    % getWindow()
    % make_psth_struct()
    % psth_err_mean() 
%
% Ex:
% data_file        -> eg. 'SA1' or any other .WAV file from the TIMIT DB.
% spl = 60         -> The program will normalize the test sentence to that
%                     particular SPL.
% adj = [-40 40]   -> The program will vary gain adjustments from -40 to
%                     +40 db in steps of 5 db.
% loss             -> Type of hearing loss number. See Also audiograms
% pres             -> Either 'DSL', 'NAL-R', or 'None'
% CFcount          -> The number of CFs to use. 40 is usually good.
% IOHC_loss        -> purely outer hair cell loss ('OHCL'), inner hair cell
%                     loss ('IHCL'), or 'Mixed' hair cell loss.
% binwidth         -> Neurogram Binwidth. Fine: 10e-6, Mean: 80e-6;
%
% Example:
% gain_simulation( 'SA1', 60, [-40 40], 5, 'NAL-R', 40, 'Mixed', 10e-6 );

%% Generate MAT Filename

% save_name = [data_file '_spl_' num2str(spl) '_adj_' num2str(adj(1)) '_' ...
%     num2str(adj(end)) '_loss_' num2str(loss) '_pres_' pres '_CFcount_' ...
%     num2str(CFcount) '_IOHCimp_' IOHC_loss '_binwidth_' ...
%     num2str(binwidth*10e6)];

% data with size (length of whole (TIMIT)-wav-file x 1 (mono)), fs is the one from wav-file
% (16 kHz)
[ data, FS ]= readsph( [data_file '.WAV'] );
    
    % set spl (in dB SPL) to input parameter spl (ususally 45,65,85 dB SPL (see F. Dinath thesis))
    data = set_spl( data(:)', spl );
    
    % plot data for visual sanity check
    %figure; plot(data); grid on; xlim([1 length(data)]); title( sprintf('TIMIT Sentence %s', data_file) );

% read out start and end positions in the wav file of the phonemes from the phn-file.
% still using FS = 16kHz, positions from the textfile also provided at this
% FS
    % keep textread
[ start, stop, ~ ] = textread( [data_file '.PHN'],'%f%f%s' );


%%%% ??????
    start = start + 1;

% compute the step sizes of the whole adjustment vector, 5 dB in F.Dinath
% thesis
%adjustmentStepSizes = diff(adj);
    
%%%%%%%%% please look at
    % this operation does weird things, do not use!
    % adj = round( adj/adjustmentStepSizes(1) ) * adjustmentStepSizes(1);  % FIXME

    % this in most cases does nothing
    % adj = round(adj/5) * 5;


% force step size to 5 dB
%     adj = adj(1):5:adj(end);

sentence = [];

% make struct to later save parameters and results
collector = struct('data_file',data_file,'spl',spl,'adj',adj,'loss',loss,...
    'pres',pres,'CFcount',CFcount,'IOHC_loss',IOHC_loss,'binwidth',binwidth,...
    'spl_col',[],'adj_col',[], 'psth_err', []);

%--------------------------------------------------------------------------
% this loop contains all steps for computing the optimal gain for the 
% sentence part constisting of all phonemes up to the newly added phoneme 
% of this loop iteration

phonemeReferenceIndex = 2;

% phoneme-loop
% for l = 2:1:( length(start) - 1 )
for control_phoneme = 2: 1 : 5   %(length(start)-1)
    
     disp( ['Phoneme Number: ' num2str(control_phoneme) - 1] );
    
    % Read data vector from start of this phoneme to stop of this phoneme.
    % --> data vector of the phoneme that is going to be added in this 
    % "phoneme-loop" run
    data_next_phone = data( start(control_phoneme):1:stop(control_phoneme) );  
    
    % Append phonemes with every loop interation with no SPL adjustment, as
    % input to make_psth_struct()
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % question sentence would be changed in its SPL, but shouldn't 
    % normal-hearing and hearing impaired compared like this: 
    % NH: no amplification | HI: pres + adj ?
    pre_sentence = [ sentence adj_spl(data_next_phone, 0) ];
        %
        figure; plot(pre_sentence); grid on; axis( [control_phoneme length(data) -0.03 0.03] ); title(sprintf('Phoneme %d', control_phoneme - 1)); shg;
    
    % compute psth (normal hearing) for all CFs with data  
    psth_normal = make_psth_struct( pre_sentence, FS,  get_spl( pre_sentence ), ...
        1, 'none', CFcount, IOHC_loss, binwidth, 'healthy' );
    
    window = get_window( start(control_phoneme), stop(control_phoneme), FS, psth_normal.psth_time, psth_normal.psth_freq, start, stop, phonemeReferenceIndex );
    
    error_m = psth_err_mean( psth_normal, psth_normal, window, 'make' );  % Structure "error_m" that later contains the optimal gain, etc.
    
    
    % This loop finds optimal gain for one phoneme.
    %
    % When the gain for the previous phonemes are optimized, what is then
    % the optimal gain for the next phoneme.
    %
    % Total result is then the adjustment, where the difference of both
    % PSTHs is a minimum.
    for j = 1: length(adj)
        
        fprintf( 1, '\n\tGain adjustment: %d of %d -> %d dB\n\n', j, numel(adj), adj(j) );
        
        pre_sentence = [ sentence adj_spl( data_next_phone, adj(j)) ];        
            %
%             figure; plot(pre_sentence); grid on; title(sprintf('Phoneme %d', l - 1)); shg; keyboard;
        
        % Calculates the PSTH of the impaired system with newly adapted SPL stimulus.
        psth_impair = make_psth_struct( pre_sentence, FS,  get_spl(pre_sentence), ...
            loss, pres, CFcount, IOHC_loss, binwidth, synaptopathy );
        
        % Applies window to the working phoneme.
        window = getWindow( start(control_phoneme), stop(control_phoneme), FS, psth_impair.psth_time, psth_impair.psth_freq, start, stop, phonemeReferenceIndex );
        
        % Appends the currently tested gain adjustment step.
        error_m.ADJ = [ error_m.ADJ adj(j) ];
        
        % calculates error for the newly adjusted stim and gives back
        % difference of the normal and impaired psth, which makes smallest
        % psth error.
        error_m = psth_err_mean( psth_normal, psth_impair, window, error_m, adj(j) );
        
    end
    
    % Appends all phonemes with now optimal gain after one another.
    sentence = [  sentence  adj_spl( data_next_phone, error_m.psth_opti )  ];
    
    % Appends the optimal gain adjustment of this phoneme.
    collector.adj_col = [  collector.adj_col  error_m.psth_opti  ];
    collector.spl_col = [  collector.spl_col  get_spl(data_next_phone)  ];
    collector.psth_err = [  collector.psth_err  error_m.psth  ];
    
    
%     directoryname = 'Sentence_Simulations';    
%         mkdir(directoryname)
%     
%     save('-mat7-binary', [ '.' filesep directoryname filesep save_name '.mat'], 'collector');
%         pause(10);  % 10 sec, pause is needed, because saving process takes some time
        
end


