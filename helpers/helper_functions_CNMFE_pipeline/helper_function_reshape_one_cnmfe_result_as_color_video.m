function helper_function_reshape_one_cnmfe_result_as_color_video(path_day_cnmfe_results,path_params_file,file_number,path_output,varargin)
% Prepare colored cnmfe videos
% 
% To annotate diff. cells in different colors from the cnmfe results I created this script for my ERC talk:
% 
% It needs:
% - the path to my folder of daily cnmfe results.
% - Path to the parameter file used to create those results
% - The file number in the day of recording (not in the annotation files) 
% - A path put the movie in
% 
% For example:
%  
% helper_function_reshape_one_cnmfe_result_as_color_video('/Users/yardenc/Documents/Projects_temp_data_processing/Cohen2022_CanaryHVCdynamics/temp_data_processing_folder/lrb853_2017-07-05_4cnmfe_source_extraction','/Users/yardenc/Documents/GitHub/Cohen2022_CanaryHVCdynamics/data/configs/config_cnmfe_lrb853_2022_01_27.mat',19,pwd)

% check that we're in a valid folder
d = dir(fullfile(path_day_cnmfe_results,'parsing_record.mat'));
if isempty(d)
    disp('Not a valid CNMFE results folder');
    return;
end
file_prefix = 'baseROIdata_';
Fs = 30;
nparams = numel(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
		case 'file_prefix'
			file_prefix=varargin{i+1};
        case 'fs'
			Fs=varargin{i+1};
    end
end
% remember the structure of filenames bird_id_counter_nframes_yyyy_mm_dd_hh_MM_ss
load(fullfile(path_day_cnmfe_results,'parsing_record.mat'),'vid_onsets','file_names');
load(path_params_file,'params');
n_del_frames = params.n_del_frames;

frames_dir = dir(fullfile(path_day_cnmfe_results,'frames*')); frames_dir = frames_dir.name;
log_dir = dir(fullfile(path_day_cnmfe_results,frames_dir,'LOGS*')); log_dir = log_dir(end).name;
matfiles = dir(fullfile(path_day_cnmfe_results,frames_dir,log_dir,'*.mat'));
matfile = matfiles(cellfun(@(x)strcmp(x,'intermediate_results.mat'),{d.name}) == 0).name;

load(fullfile(path_day_cnmfe_results,frames_dir,log_dir,matfile),'neuron');
vid_onsets = [vid_onsets; size(neuron.C,2)];
colors_2_use = distinguishable_colors(size(neuron.C,1),'k');
%build 3d matrix and dff for each file

deleted_frames_filler = zeros(neuron.options.d1,neuron.options.d2,3);
figure;imshow(deleted_frames_filler);set(gca,'nextplot','replacechildren');
mtframe = getframe(gcf);

curr_idx = file_number; 
start_frame_idx = vid_onsets(curr_idx);
end_frame_idx = vid_onsets(curr_idx+1)-1;

v = VideoWriter(fullfile(path_output,'test3.avi'),'Motion JPEG AVI');
open(v);
for mtframen=1:n_del_frames
    writeVideo(v,mtframe);
end

for frame_n = start_frame_idx:end_frame_idx
    currframe = zeros(neuron.options.d1,neuron.options.d2,3);
    for celln = 1:size(neuron.C,1)
        curradd = cat(3,ones(neuron.options.d1,neuron.options.d2)*colors_2_use(celln,1),ones(neuron.options.d1,neuron.options.d2)*colors_2_use(celln,2),ones(neuron.options.d1,neuron.options.d2)*colors_2_use(celln,3));
        currframe = currframe + repmat(neuron.reshape(neuron.A(:,celln)*neuron.C(celln,frame_n),2),1,1,3).*curradd*0.05;
    end
    imshow(currframe);
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)
