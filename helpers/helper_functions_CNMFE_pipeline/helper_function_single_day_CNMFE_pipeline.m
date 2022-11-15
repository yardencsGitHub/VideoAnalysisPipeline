function [neuron,vid_onsets,file_names] = helper_function_single_day_CNMFE_pipeline(bird_id,Day,path_to_parameters_file,varargin)
% This script runs the pipeline that: 
% >> 1. takes movies from a single day and prepares a single file to use in CNMFE, 
% >> 2. runs the CNMFE scripts, 
% >> 3. runs the manual postprocessing (carried in nested function)
% >> 4. save results (carried in nested function)

% Note: Another script will be responsible for breaking the results into the
% original files, and to matching with the numbererd, song containing,
% files for each bird.

% All parameters, shared by all days must be saved in one file
load(path_to_parameters_file);

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
adjust_pnr_cn_thresholds = false;
vid_onsets = [];
neuron = [];
nparams = length(varargin);
files_to_skip = {};
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
        case 'file_names'
			file_names=varargin{i+1};
        case 'githubfolder'
			GitHubFolder=varargin{i+1};
        case 'neuron'
			neuron=varargin{i+1};
        case 'adjust_pnr_cn_thresholds'
			adjust_pnr_cn_thresholds=varargin{i+1};
        case 'files_to_skip'
			files_to_skip=varargin{i+1};
         
    end
end

% 1. Prepare the concatenated file.
% This requires deleting frames that are later re-inserted.
if ~skip_step_1
    [vid_onsets,file_names] = helper_function_create_concatenated_vidMat(Day,'githubfolder',GitHubFolder,'datafolder',DataFolder,...
                                                            'output_folder',output_folder,'output_filename',output_filename,...
                                                            'max_frames',max_frames,'n_del_frames',n_del_frames,'resize_factor',resize_factor, ...
                                                            'bird_id',bird_id,'files_to_skip',files_to_skip);
    save(fullfile(output_folder,'temp_vid_onset_and_file_names.mat'),'vid_onsets','file_names');
end

% 2. Run the CNMFE pipeline
if ~skip_step_2
    addpath(genpath([GitHubFolder 'CNMF_E_CohenLab']));
    neuron = large_data_1p_pipeline_with_parameters(fullfile(output_folder,output_filename),params,...
                                                    'skip_step_3',skip_step_3,'skip_step_4',skip_step_4,'adjust_pnr_cn_thresholds',adjust_pnr_cn_thresholds);
end

source_extraction_target_folder = fullfile(output_folder,...
                               [output_filename(1:end-4) '_source_extraction']);
if ~exist(source_extraction_target_folder)
    mkdir(source_extraction_target_folder);
end                              
parsing_output_file = fullfile(source_extraction_target_folder,...
                               'parsing_record.mat');
save(parsing_output_file,'vid_onsets','file_names');




