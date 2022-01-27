function neuron = helper_function_single_day_CNMFE_pipeline(bird_id,Day,path_to_parameters_file,varargin)
% This script runs the pipeline that: 
% >> 1. takes movies from a single day and prepares a single file to use in CNMFE, 
% >> 2. runs the CNMFE scripts, 
% >> 3. runs the manual postprocessing
% >> 4. save results

% Note: Another script will be responsible for breaking the results into the
% original files, and to matching with the numbererd, song containing,
% files for each bird.

% All parameters, shared by all days must be saved in one file
params = load(path_to_parameters_file);

GitHubFolder = '/Users/yardenc/Documents/GitHub/';
DataFolder = params.DataFolder;
output_folder = params.output_folder;
output_filename = [bird_id '_' Day '_4cnmfe.mat'];
max_frames = params.max_frames;
n_del_frames = params.n_del_frames;
resize_factor = params.resize_factor;
% these varargins allow skipping some steps 
skip_step_1 = false; 
skip_step_2 = false;
skip_step_3 = false;
skip_step_4 = false;
vid_onsets = [];
neuron = [];
nparams = length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
		case 'skip_step_1'
			skip_step_1=varargin{i+1};
        case 'skip_step_2'
			skip_step_2=varargin{i+1};
        case 'skip_step_3'
			skip_step_3=varargin{i+1};
        case 'skip_step_4'
			skip_step_4=varargin{i+1};
        case 'vid_onsets'
			vid_onsets=varargin{i+1};
        case 'githubfolder'
			GitHubFolder=varargin{i+1};
        case 'neuron'
			neuron=varargin{i+1};
         
    end
end

% 1. Prepare the concatenated file.
% This requires deleting frames that are later re-inserted.
if ~skip_step_1
    vid_onsets = helper_function_create_concatenated_vidMat(Day,'githubfolder',GitHubFolder,'datafolder',DataFolder,...
                                                            'output_folder',output_folder,'output_filename',output_filename,...
                                                            'max_frames',max_frames,'n_del_frames',n_del_frames,'resize_factor',resize_factor);
end

% 2. Run the CNMFE pipeline
if ~skip_step_2
    addpath(genpath([GitHubFolder 'CNMF_E']));
    neuron = large_data_1p_pipeline_with_parameters(fullfile(output_folder,output_filename),params);
end

% 3. manual post processing
if ~skip_step_3
    neuron.manual_prune_and_merge();
end

% 4. save the workspace for future analysis
if ~skip_step_4
    neuron.orderROIs('snr');
    cnmfe_path = neuron.save_workspace();
    % save neurons shapes
    neuron.save_neurons();
end


pars_envs = struct('memory_size_to_use', 40, ...   % GB, memory space you allow to use in MATLAB %8
    'memory_size_per_patch', 3.9, ...   % GB, space for loading data within one patch %0.6
    'patch_dims', [32, 32]);  %GB, patch size [64,64]

% -------------------------      SPATIAL      -------------------------  %
gSig = 5.5;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering %3
gSiz = 12;          % pixel, neuron diameter %13
ssub = 2;           % spatial downsampling factor
with_dendrites = false;   % with dendrites or not %false
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    updateA_search_method = 'dilate';  %#ok<UNRCH>
    updateA_bSiz = 5;
    updateA_dist = neuron.options.dist;
else
    % determine the search locations by selecting a round area
    updateA_search_method = 'ellipse'; %#ok<UNRCH>
    updateA_dist = 5;
    updateA_bSiz = neuron.options.dist;
end
spatial_constraints = struct('connected', true, 'circular', false);  % you can include following constraints: 'circular'
spatial_algorithm = 'hals_thresh';

% -------------------------      TEMPORAL     -------------------------  %
Fs = 30;             % frame rate %10
tsub = 3;           % temporal downsampling factor %1
deconv_flag = true;     % run deconvolution or not 
deconv_options = struct('type', 'ar2', ... % model of the calcium traces. {'ar1', 'ar2'} %ar1
    'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'smin', -5, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true, ...% optimize the baseline);
    'max_tau', 100);    % maximum decay time (unit: frame);

nk = 3;             % detrending the slow fluctuation. usually 1 is fine (no detrending)
% when changed, try some integers smaller than total_frame/(Fs*30)
detrend_method = 'spline';  % compute the local minimum as an estimation of trend.

% -------------------------     BACKGROUND    -------------------------  %
bg_model = 'ring';  % model of the background {'ring', 'svd'(default), 'nmf'}
nb = 1;             % number of background sources for each patch (only be used in SVD and NMF model)
ring_radius = 14;  % when the ring model used, it is the radius of the ring used in the background model. %18
%otherwise, it's just the width of the overlapping area
num_neighbors = []; % number of neighbors for each neuron
bg_ssub = 2;        % downsample background for a faster speed 

% -------------------------      MERGING      -------------------------  %
show_merge = false;  % if true, manually verify the merging step
merge_thr = 0.3;     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation] %0.65
method_dist = 'mean';   % method for computing neuron distances {'mean', 'max'}
dmin = 10;       % minimum distances between two neurons. it is used together with merge_thr %5
dmin_only = 2;  % merge neurons if their distances are smaller than dmin_only.
merge_thr_spatial = [0.8, 0.1, -inf];  % merge components with highly correlated spatial shapes (corr=0.8) and small temporal correlations (corr=0.1) %[0.8, 0.4, -inf];

% -------------------------  INITIALIZATION   -------------------------  %
K = [];             % maximum number of neurons per patch. when K=[], take as many as possible.
min_corr = 0.4;     % minimum local correlation for a seeding pixel %0.8
min_pnr = 7;       % minimum peak-to-noise ratio for a seeding pixel %8
min_pixel = gSig^2;      % minimum number of nonzero pixels for each neuron
bd = 0;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
frame_range = [];   % when [], uses all frames
save_initialization = false;    % save the initialization procedure as a video.
use_parallel = true;    % use parallel computation for parallel computing
show_init = true;   % show initialization results
choose_params = true; % manually choose parameters
center_psf = true;  % set the value as true when the background fluctuation is large (usually 1p data)
% set the value as false when the background fluctuation is small (2p)

% -------------------------  Residual   -------------------------  %
min_corr_res = 0.4;
min_pnr_res = 6;
seed_method_res = 'auto';  % method for initializing neurons from the residual
update_sn = true;

% ----------------------  WITH MANUAL INTERVENTION  --------------------  %
with_manual_intervention = true;

% -------------------------  FINAL RESULTS   -------------------------  %
save_demixed = true;    % save the demixed file or not
kt = 3;                 % frame intervals

