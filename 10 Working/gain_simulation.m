function [ collector ] = gain_simulation( data_file, spl, adj, loss, pres, CFcount, IOHC_loss, binwidth, synaptopathy )
%GAIN_SIMULATION
% function to prepare input signals (phoneme chains) and compute optimal
% gains using the functions:
%
%     make_psth_struct()
%     get_window()
%     psth_err_mean()
%
%
% Output:
% collector: struct with specification fields data_file, spl, adj, loss,
% pres, CFcount,IOHC_loss, binwidth and the result fields spl_col, adj_col
% psth_err
%
% Input:
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
%
% Example:
% collector = gain_simulation( data_file, spl, adj,    loss, pres,   CFcount,...
% IOHC_loss, binwidth, synaptopathy )

% collector = gain_simulation( 'SA1',     60, [-40 40], 5,  'NAL-R', 40,...
% 'Mixed',   10e-6,    'healthy');

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
% FS, also: keep textread against Matlab recommendation
[ start, stop, ~ ] = textread( [data_file '.PHN'],'%f%f%s' );


%%%%% question to Michael: why would one do this?
start = start + 1;

% this is the variable, that will contain the phonemes with their optimal
% gain (see line xxxx for comparison)
sentence = [];

% make struct to later save parameters and results
collector = struct('data_file',data_file,'spl',spl,'adj',adj,'loss',loss,...
    'pres',pres,'CFcount',CFcount,'IOHC_loss',IOHC_loss,'binwidth',binwidth,...
    'spl_col',[],'adj_col',[], 'psth_err', []);

%% phoneme-loop -----------------------------------------------------------
% this loop contains all steps for computing the optimal gain for the current 
% sentence part consisting of all prior phonemes (attenuated with their optimal 
% gain respectively) and the newly added phoneme of this loop iteration 
% (controlled by control_phoneme). The newly added phoneme is tested with
% all adjustment steps from vetor adj (see loop in line xxx)

% for control_phoneme = 2:1:( length(start) - 1 )
for control_phoneme = 4%2: 1 : 5   %(length(start)-1)
    
    disp( ['Phoneme Number: ' num2str(control_phoneme) - 1] );
    
    % Read data vector from start of this phoneme to stop of this phoneme.
    % --> data vector of the phoneme that is going to be added in this
    % "phoneme-loop" run
    data_next_phone = data( start(control_phoneme):1:stop(control_phoneme) );
    
    % Here the new phone is added to the analyzed data without a change of
    % SPL
    % second input argument to adj-spl() is a gain in SPL (added to current SPL)
    pre_sentence = [ sentence adj_spl(data_next_phone, 0) ];
    
    figure; plot(pre_sentence); grid on; 
    axis( [control_phoneme length(data) -0.03 0.03] );
    title(sprintf('Phoneme %d', control_phoneme - 1)); shg;
    
    % make_psth_struct( data, FS, spl, loss, pres, CFcount, IOHC_loss, binwidth, synaptopathy )
    % make psth struct and compute psth (normal hearing) for all CFs with so-far data
    psth_normal = make_psth_struct( pre_sentence, FS,  get_spl( pre_sentence ), ...
        1, 'none', CFcount, IOHC_loss, binwidth, 'healthy' );
    
    % Applies window (ramps) to the working phoneme.why?
    window = get_window( start(control_phoneme), stop(control_phoneme), FS,...
             psth_normal.psth_time, psth_normal.psth_freq);
    
    % creates the struct that will later contain the error metrics between the
    % different neurograms, 
    % sanity check note: this run of psth_error_mean just creates the
    % struct and should contain empty fields except for fields 'type' and
    % 'meth'
    error_m = psth_err_mean( psth_normal, psth_normal, window, 'make' );  % Structure "error_m" that later contains the optimal gain, etc.
    
    % sanity check for the normal psth is comparing the psth with itself,
    % error_m_sanity.psth should be zero
    % error_m_sanity = psth_err_mean( psth_normal, psth_normal, window, error_m);  % Structure "error_m" that later contains the optimal gain, etc.
    
    %% adj-loop -----------------------------------------------------------
    % This loop finds optimal gain for the newly added phoneme by looping
    % through the adjustment steps (vector adj) and computing an error
    % between the normal psth and the impaired psth of the adjusted phoneme
    %
    % Final result (collector.adj_col) is then the adjustment step for this
    % phoneme, at which the difference of both PSTHs (normal & no gain, 
    % impaired & prescription & adjustment) is at its minimum.
    for control_adj = 1: length(adj)
        
        fprintf( 1, '\n\tGain adjustment: %d of %d -> %d dB\n\n', control_adj, numel(adj), adj(control_adj) );
        
        % Here, after the old phones( with respective optimized gain), 
        % the new phone is added using the adjustment step of this loop 
        % iteration (control_adj) 
        % second input argument to adj-spl() is a gain in SPL (added to current SPL) 
        pre_sentence = [ sentence adj_spl( data_next_phone, adj(control_adj)) ];
       
        % figure; plot(pre_sentence); grid on; 
        % title(sprintf('Phoneme %d', l - 1)); shg; keyboard;
        
        % Calculates the PSTH of the impaired system with the newly
        % adjusted phoneme
        psth_impair = make_psth_struct( pre_sentence, FS,  get_spl(pre_sentence), ...
            loss, pres, CFcount, IOHC_loss, binwidth, synaptopathy );
        
        % Applies window (ramps) to the working phoneme.why?
        window = get_window( start(control_phoneme), stop(control_phoneme), FS, ...
                 psth_impair.psth_time, psth_impair.psth_freq);
        
        % Appends the currently tested gain adjustment step.
        error_m.ADJ = [ error_m.ADJ adj(control_adj) ];
        
        % calculates error for the newly adjusted pre-sentence and returns
        % the difference of the normal and impaired psth, which makes smallest
        % psth error in field psth_opti.
        error_m = psth_err_mean( psth_normal, psth_impair, window, error_m, adj(control_adj) );
        
        tttt=1;
    end
    
    % Appends the currently analyzed phoneme with the computed optimal gain
    % to the sequence of phonemes.
    sentence = [  sentence  adj_spl( data_next_phone, error_m.psth_opti )  ];
    
    % Appends the optimal gain adjustment of this phoneme.
    % collector.adj_col contains the adjustment steps for the different
    % phonemes, at which the error between normal and impaired psth were
    % minimal, final size [1 x length(start)]
    collector.adj_col = [  collector.adj_col  error_m.psth_opti  ];
    
    % collector.spl_col contains the sound pressure level of the different
    % phonemes, which can be different from the input SPL, since this is an
    % overall SPL, to which the whole sentence is normalized
    % final size [1 x length(start)]
    collector.spl_col = [  collector.spl_col  get_spl(data_next_phone)  ];
    
    % appends the psth errors (difference between normal and impair for the
    % currently analyzed phoneme for all adjustment steps
    % final size [1 x length(start) * length(adj)]
    collector.psth_err = [  collector.psth_err  error_m.psth  ];
    
    
    %     directoryname = 'Sentence_Simulations';
    %         mkdir(directoryname)
    %
    %     save('-mat7-binary', [ '.' filesep directoryname filesep save_name '.mat'], 'collector');
    %         pause(10);  % 10 sec, pause is needed, because saving process takes some time
    
end
tttt=nan;

