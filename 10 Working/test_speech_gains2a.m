

%% Environment

close all;
clear all; clc;
set(0, 'DefaultFigureWindowStyle', 'docked');

addpath( genpath(pwd), '-begin' );

%pkg load signal  % For Octave
%mexANmodel


%% Gain Simulation

% "gain_simulation.m" computes optimal gains for all phonemes of a sentence
% and saves the results in the struct collector in a .mat-file named similar to
% the following example.
%
% Example MAT filename:  "SA1_spl_45_adj_-40_40_loss_2_pres_NAL-R_CFcount_30_IOHCimp_Mixed_binwidth_100.mat"

% gain_simulation( 'SA1', 45, [-40:5:40], 2, 'NAL-R', 30, 'Mixed', 10e-6, 'healthy' );


%% fdinath_thesis

%  45, 65, or 85 dB SPL

% ?40 to +40 dB in ?5 dB increments

% CF count: 30

% [ collector ] = gain_simulation( data_file, spl, adj,      loss, pres,    CFcount, IOHC_loss, binwidth, synaptopathy )
%collector_mr =    gain_simulation( 'SA1',     45, [-40:20:40], 2, 'NAL-R', 5,       'Mixed',    100e-6,   'healthy' )



% sentence = 'SA1';
% spl_basic = 45;
% adj_steps = [-40:5:40];
% hl = 2;
% pres = 'NAL-R';
% num_CF = 30;
% ohc_ihc = 'Mixed';
% binwidth = 10e-6;
% synaptopathy = 'healthy';

% collector_ft =    gain_simulation( 'SA1', 45, [-40:5:40], 2, 'NAL-R', 30,       'Mixed',    binwidth,  healthy  )

sentence = 'SA1';
spl_basic = 85;
%adj_steps = [-40:5:40];
adj_steps = [-40:5:40];
hl = 2;
pres = 'DSL';
%num_CF = 30;
num_CF = 30;
ohc_ihc = 'Mixed';
binwidth = 100e-6;
synaptopathy = 'healthy';


collector =    gain_simulation( sentence, spl_basic, adj_steps, hl, pres, num_CF,       ohc_ihc,    binwidth, synaptopathy    );


% Usage:q

directoryname = 'Sentence_Simulations';

%mkdir(directoryname)
% save_name = ['Sentence_Simulations' '_spl_' num2str(45) '_adj_' num2str(-40) '_' ...
% num2str(40) '_loss_' num2str(2) '_pres_' 'NAL-R' '_CFcount_' ...
% num2str(5) '_IOHCimp_' 'Mixed' '_binwidth_' ...
% num2str(binwidth*1e6)];

save_name = ['Sentence_Simulations' '_spl_' num2str(spl_basic) '_adj_' num2str(adj_steps(1)) '_' ...
    num2str(adj_steps(end)) '_loss_' num2str(hl) '_pres_' pres '_CFcount_' ...
    num2str(num_CF) '_IOHCimp_' ohc_ihc '_binwidth_' ...
    num2str(binwidth*1e6) '_synaptopathy_' synaptopathy];


%load(save_name)

%
save([ '.' filesep directoryname filesep save_name '.mat'], 'collector');
pause(10);  % 10 sec, pause is needed, because saving process takes some time

figure('Name', 'AVG'); plot(collector.spl_col,collector.adj_col, 'x', 'Linewidth', 4)
ylim([-40,40])
xlim([0 100])



% reshape psth_error matrix to
% be: [number of input adjustment steps x number of phonemes]
psth_errors_matrix = reshape(collector_ft_dsl.psth_err,length(collector.adj),length(collector.adj_col))


phoneme_tobeplotted = 7;
figure; plot(collector_ft_dsl.adj,psth_erros_matrix(:,phoneme_tobeplotted))




%     save('-mat7-binary', [ '.' filesep directoryname filesep save_name '.mat'], 'collector_mr');
%         pause(10);  % 10 sec, pause is needed, because saving process takes some time
%






%        collector = gain_simulation(    'SA1', 45, [-40:20:40], 2, 'NAL-R', 5, 'Mixed', 10e-6,  'healthy' )
%
%                 mkdir(directoryname)
%save_name = ['Sentence_Simulations' '_spl_' num2str(45) '_adj_' num2str(-40) '_' ...
%num2str(40) '_loss_' num2str(2) '_pres_' 'NAL-R' '_CFcount_' ...
%num2str(5) '_IOHCimp_' 'Mixed' '_binwidth_' ...
%num2str(10e-6*10e6)];
%
%
%
%
%
%
%
%figure('Name', 'FT'); plot(collector.spl_col,collector.adj_col)
%
%








% gain_simulation( data_file, spl, adj, loss, pres, CFcount, IOHC_loss, binwidth, synaptopathy )
%
% where the inputs respectively are:
%
% data_file        -> eg. 'SA1' or any other .WAV file from the TIMIT DB.
% spl = 60         -> The program will normalize the test sentence to that
%                     particular SPL.
% adj = [-40 40]   -> The program will vary gain adjustments from -40 to
%                     +40 db in steps of 5 dB.
% loss             -> Type of hearing loss number. See Also audiograms
% pres             -> Hearing aid prescription algorithm:
%                     Either 'DSL', 'NAL-R', or 'None'
% CFcount          -> The number of CFs to use. 40 is usually good.
% IOHC_loss        -> purely outer hair cell loss ('OHCL'), inner hair cell
%                     loss ('IHCL'), or 'Mixed' hair cell loss.
% binwidth         -> Neurogram Binwidth. Fine: 10e-6, Mean: 80e-6


