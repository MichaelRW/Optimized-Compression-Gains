

%% Environment

close all; clear all; clc;
set(0, 'DefaultFigureWindowStyle', 'docked');

addpath( genpath(pwd), '-begin' );

% pkg load signal  % For Octave
% mexANmodel


%% Gain Simulation

% "gain_simulation.m" computes optimal gains for all phonemes of a sentence 
% and saves the results in the struct collector in a .mat-file named similar to
% the following example.
%
% Example MAT filename:  "SA1_spl_45_adj_-40_40_loss_2_pres_NAL-R_CFcount_30_IOHCimp_Mixed_binwidth_100.mat"

% gain_simulation( 'SA1', 45, [-40:5:40], 2, 'NAL-R', 30, 'Mixed', 10e-6, 'healthy' );


collector = gain_simulation( 'SA1', 45, [-40:40:0], 2, 'NAL-R', 5, 'Mixed', 10e-6, 'healthy' )


% Usage:
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


