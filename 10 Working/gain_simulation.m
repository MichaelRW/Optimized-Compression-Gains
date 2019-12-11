function gain_simulation(data_file,spl,adj,loss,pres,CFcount,IOHC_loss,binwidth, synaptopathy)
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

% Creates a save directory and switches into it.                          %
%=========================================================================%
save_name = [data_file '_spl_' num2str(spl) '_adj_' num2str(adj(1)) '_' ...
    num2str(adj(end)) '_loss_' num2str(loss) '_pres_' pres '_CFcount_' ...
    num2str(CFcount) '_IOHCimp_' IOHC_loss '_binwidth_' ...
    num2str(binwidth*10e6)];
%=========================================================================%
[data,FS]= readsph([data_file '.WAV']);
data = set_spl(data(:)',spl);

[start,stop,~] = textread([data_file '.PHN'],'%f%f%s');
start = start + 1;

adj = round(adj/5)*5; adj = adj(1):5:adj(end);

sentence = [];

collector = struct('data_file',data_file,'spl',spl,'adj',adj,'loss',loss,'pres',pres,...
    'CFcount',CFcount,'IOHC_loss',IOHC_loss,'binwidth',binwidth,'spl_col',[],'adj_col',[], 'psth_err', []);

%--------------------------------------------------------------------------
% this loop contains all steps for computing the optimal gain for the sentence
% part constisting of all phonemes up to the newly added phoneme of this
% loop iteration

for l = 1 %: length(start)
    
    disp(['phone no ' num2str(l)])
    
    % read data vector from start of this phoneme to stop of this phoneme
    next_phone = data(start(l):stop(l));
    
    % appends phonemes with every loop interation
    %                       adj_spl gives adjusted data vector, 0 means no
    %                       adjustment
    pre_sentence = [sentence adj_spl(next_phone,0)];
    
    
    % psth =      make_psth_struct(data,         FS,  spl,                  ...
    % loss, pres, CFcount, IOHC_loss, binwidth, synaptopathy)
    
    % where get_spl simply calculates the SPL from the current stimulus
    psth_normal = make_psth_struct(pre_sentence, FS,  get_spl(pre_sentence), ...
        1, 'none', CFcount, IOHC_loss, binwidth, 'healthy');
    
    %
    window = getWindow(start(l), stop(l), FS, psth_normal.psth_time, psth_normal.psth_freq);
    
    
    % set up struct error_m, which later contains the optimal gain etc.
    error_m      = psth_err_mean( psth_normal, psth_normal, window, 'make' );
    
    %----------------------------------------------------------------------
    
    % this loop finds optimal gain for one phoneme, when the ones before
    % are already optimized; so if the one phoneme is already optimized,
    % what is then the optimal gain for the next one
    % total result is then the adjustment, where difference of both psths is minimum
    for j = 1: length(adj)
        
        
        disp(['adj no ' num2str(j)])
        
        %
        pre_sentence = [sentence adj_spl(next_phone,adj(j))];
        
        
        % calculates psth of the impaired system with newly adapted SPL stimulus
        psth_impair = make_psth_struct(pre_sentence, FS,  get_spl(pre_sentence), ...
            loss, pres, CFcount, IOHC_loss, binwidth, synaptopathy);
        
        % applies window to the phoneme 
        window = getWindow(start(l), stop(l), FS, psth_impair.psth_time, psth_impair.psth_freq);
        
        % appends the currently tested gain adjustment step
        error_m.ADJ = [error_m.ADJ adj(j)];
        
        % calculates error for the newly adjusted stim and gives back
        % difference of the normal and impaired psth, which makes smallest
        % psth error.
        error_m     = psth_err_mean( psth_normal, psth_impair, window, error_m );
        
    end
    
    % appends all phonemes with now optimal gain after one another
    sentence = [sentence adj_spl(next_phone,error_m.psth_opti)];
    
    % appends the optimal gain adjustment of this phoneme
    collector.adj_col = [collector.adj_col error_m.psth_opti];
    collector.spl_col = [collector.spl_col get_spl(next_phone)];
    collector.psth_err = [collector.psth_err error_m.psth ];
    
    directoryname = 'Sentence_Simulations';
    
    mkdir(directoryname)
    
    save('-mat7-binary', [ '.' filesep directoryname filesep save_name '.mat'], 'collector')
    
    % 10 sec, pause is needed, because saving process takes some time
    pause(10)
end


