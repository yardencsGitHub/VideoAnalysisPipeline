function [ax,r,p] = LongRangeLockedSingleDayManualROIs_function(ax,Day,ignore_entries,join_entries,sylidx,syllabels_sequence,ROIs,locktoonset,spikes,order_flag)
%% 
% This script creates single ROI single day alignments to complex sequences
% with only one variable (duration or syllable type)
% Inputs:
%   ax - axes where to plot
%   Day - text rep. of date
%   ignore_entries - A vector of label numbers to ignore completely. 
%   join_entries - A cell of vectors, each containing a >1 number of labels
%   sylidx - index of syllable in the sequence to align to
%   syllabels_sequence - the sequence of phrase identities to lock. Insert
%       nan to keep it free
%   ROIs - the roi # 
%   locktoonset - 0/1 to lock to phrase onset or offset
%   spikes - type of variable (0-3):
%       0 - denoised Ca
%       1 - deconvolved Spikes
%       2 - Fluorescence
%       3 - HMM state (signal \ noise)
%   order_flag - the correlate. Use positive integers to indicate durations of
%   phrases. A vector of positive integers will result in summation of
%   durations. Use negative integers for type (no vectors)

%% change according to workstation
BaseDir = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/';
GithubDir = '/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/'
%%
bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation4'};
bird_params = bird1_params;
delete_frames = 1;
n_del_frames = 6;
hvc_offset = 0.04;
mulcnt = 2;
edges = [0.1 0.1];
opacity_factor = 0.5;
%%
addpath(genpath([GithubDir 'small-utils']));
addpath(genpath([GithubDir 'VideoAnalysisPipeline']));
bird_name = bird_params{1}; 
bird_folder_name = bird_params{2}; 
template_file = bird_params{3}; 
annotation_file = bird_params{4}; 
CNMFEfolder = [GithubDir 'CNMF_E'];

%% Folders that contain data
% Folders on laptop:
laptop_mov_folder = [BaseDir bird_folder_name '/movs'];
laptop_wav_folder = [BaseDir bird_folder_name '/movs/wav'];
laptop_gif_folder = [BaseDir bird_folder_name '/movs/wav/gif'];
laptop_annotated_dir = [BaseDir bird_folder_name '/movs/wav/annotated'];
laptop_annotated_images_dir = [BaseDir bird_folder_name '/movs/wav/annotated/images'];
DamagedFolder = [BaseDir bird_folder_name '/too_large_or_damaged/'];
laptop_manualROI_folder = [BaseDir bird_folder_name '/ManualROIs'];

