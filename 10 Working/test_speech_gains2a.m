%% Environment
close all;
%clear all; clc;
set(0, 'DefaultFigureWindowStyle', 'docked');
addpath( genpath(pwd), '-begin' );

if(exist('OCTAVE_VERSION','builtin')~=0)
    pkg load signal  % For Octave
end
%mexANmodel


%% Gain Simulation
% gain_simulation() computes optimal gains for all phonemes of a sentence
% and saves the results in the struct collector in a .mat-file named similar to
% the following example.
%
% Example MAT filename:  "SA1_spl_45_adj_-40_40_loss_2_pres_NAL-R_CFcount_30_IOHCimp_Mixed_binwidth_100.mat"

% gain_simulation( 'SA1', 45, [-40:5:40], 2, 'NAL-R', 30, 'Mixed', 10e-6, 'healthy' );
%% fdinath_thesis
%  45, 65, or 85 dB SPL
% [-40 to +40] dB in 5 dB increments
% CF count: 30
%% 
sentence = 'SA1';% other TIMIT sentences available
spl_basic = 85;%[45, 65, 85]
adj_steps = [-40:5:40];
hl = 2;% 1 to 15 with 1 being "normal"
pres = 'DSL';%'DSL';% or 'NAL-R';
num_CF = 30;%30;
ohc_ihc = 'Mixed';% {'IHCL', 'OHCL', 'Mixed'} type of hearing loss
binwidth = 100e-6;% or for 'FINE': 10e-6 
synaptopathy = 'healthy';%{'healthy','50all','lossall_low','IHCproportional'}

% using gain_simulation
collector =    gain_simulation( sentence, spl_basic, adj_steps, hl, pres, ...
                                num_CF, ohc_ihc, binwidth, synaptopathy);

directoryname = 'Sentence_Simulations';

save_name = ['Sentence_Simulations' ...
             '_spl_' num2str(spl_basic),...
             '_adj_' num2str(adj_steps(1)) '_' num2str(adj_steps(end))...
             '_loss_' num2str(hl)...
             '_pres_' pres...
             '_CFcount_' num2str(num_CF)...
             '_IOHCimp_' ohc_ihc...
             '_binwidth_' num2str(binwidth*1e6)...
             '_synaptopathy_' synaptopathy];

save([ '.' filesep directoryname filesep save_name '.mat'], 'collector');
pause(10);  % 10 sec, pause is needed, because saving process takes some time
%%
% PLOTTING AS IN DINATH THESIS (on p.57)
% figure('Name', 'AVG'); plot(collector.spl_col,collector.adj_col, 'x', 'Linewidth', 4)
% ylim([-40,40])
% xlim([0 100])

% PLOTTING ALL ERRORS 
% reshape psth_error matrix to
% be: [number of input adjustment steps x number of phonemes]
% psth_errors_matrix = reshape(collector_ft_dsl.psth_err,length(collector.adj),length(collector.adj_col))
% phoneme_tobeplotted = 7;
% figure; plot(collector_ft_dsl.adj,psth_erros_matrix(:,phoneme_tobeplotted))

