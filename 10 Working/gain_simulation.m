
function [ collector ] = gain_simulation( data_file, spl, adj, loss, pres, CFcount, IOHC_loss, binwidth, synaptopathy )

%GAIN_SIMULATION
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

save_name = [data_file '_spl_' num2str(spl) '_adj_' num2str(adj(1)) '_' ...
    num2str(adj(end)) '_loss_' num2str(loss) '_pres_' pres '_CFcount_' ...
    num2str(CFcount) '_IOHCimp_' IOHC_loss '_binwidth_' ...
    num2str(binwidth*10e6)];

[ data, FS ]= readsph( [data_file '.WAV'] );
    data = set_spl( data(:)', spl );
    %
    figure; plot(data); grid on; xlim([1 length(data)]); title( sprintf('TIMIT Sentence %s', data_file) );

[ start, stop, ~ ] = textread( [data_file '.PHN'],'%f%f%s' );
    start = start + 1;
    
adjustmentStepSizes = diff(adj);
    adj = round( adj/adjustmentStepSizes(1) ) * adjustmentStepSizes(1);  % FIXME
%
% adj = round(adj/5) * 5;
%     adj = adj(1):5:adj(end);

sentence = [];

collector = struct('data_file',data_file,'spl',spl,'adj',adj,'loss',loss,'pres',pres,...
    'CFcount',CFcount,'IOHC_loss',IOHC_loss,'binwidth',binwidth,'spl_col',[],'adj_col',[], 'psth_err', []);

%--------------------------------------------------------------------------
% this loop contains all steps for computing the optimal gain for the sentence
% part constisting of all phonemes up to the newly added phoneme of this
% loop iteration

phonemeReferenceIndex = 2;

% for l = 2:1:( length(start) - 1 )
for l = 2:4 %: length(start)
    
    disp( ['Phoneme Number: ' num2str(l) - 1] );
    
    next_phone = data( start(l):1:stop(l) );  % Read data vector from start of this phoneme to stop of this phoneme.
    
    pre_sentence = [ sentence adj_spl(next_phone, 0) ];  % Append phonemes with every loop interation with no SPL adjustment.
        %
        figure; plot(pre_sentence); grid on; axis( [1 length(data) -0.03 0.03] ); title(sprintf('Phoneme %d', l - 1)); shg;
    
    psth_normal = make_psth_struct( pre_sentence, FS,  get_spl( pre_sentence ), ...
        1, 'none', CFcount, IOHC_loss, binwidth, 'healthy' );
    
    window = getWindow( start(l), stop(l), FS, psth_normal.psth_time, psth_normal.psth_freq, start, stop, phonemeReferenceIndex );
    
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
        
        pre_sentence = [ sentence adj_spl( next_phone, adj(j)) ];        
            %
%             figure; plot(pre_sentence); grid on; title(sprintf('Phoneme %d', l - 1)); shg; keyboard;
        
        % Calculates the PSTH of the impaired system with newly adapted SPL stimulus.
        psth_impair = make_psth_struct( pre_sentence, FS,  get_spl(pre_sentence), ...
            loss, pres, CFcount, IOHC_loss, binwidth, synaptopathy );
        
        % Applies window to the working phoneme.
        window = getWindow( start(l), stop(l), FS, psth_impair.psth_time, psth_impair.psth_freq, start, stop, phonemeReferenceIndex );
        
        % Appends the currently tested gain adjustment step.
        error_m.ADJ = [ error_m.ADJ adj(j) ];
        
        % calculates error for the newly adjusted stim and gives back
        % difference of the normal and impaired psth, which makes smallest
        % psth error.
        error_m = psth_err_mean( psth_normal, psth_impair, window, error_m, adj(j) );
        
    end
    
    % Appends all phonemes with now optimal gain after one another.
    sentence = [  sentence  adj_spl( next_phone, error_m.psth_opti )  ];
    
    % Appends the optimal gain adjustment of this phoneme.
    collector.adj_col = [  collector.adj_col  error_m.psth_opti  ];
    collector.spl_col = [  collector.spl_col  get_spl(next_phone)  ];
    collector.psth_err = [  collector.psth_err  error_m.psth  ];
    
    
%     directoryname = 'Sentence_Simulations';    
%         mkdir(directoryname)
%     
%     save('-mat7-binary', [ '.' filesep directoryname filesep save_name '.mat'], 'collector');
%         pause(10);  % 10 sec, pause is needed, because saving process takes some time
        
end


