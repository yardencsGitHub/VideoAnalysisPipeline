function helper_function_reshape_single_day_cnmfe_spiking_as_matrices(path_day_cnmfe_results,path_params_file,path_original_results,path_output,varargin)
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
parsed_files_datetimes = {};
for file_idx = 1:numel(file_names)
    tokens = split(file_names{file_idx},'_');
    parsed_files_datetimes(file_idx) = {char(join(tokens(4:9),'_'))};
end
load(path_params_file,'params');
n_del_frames = params.n_del_frames;
if isempty(path_original_results)
    orig_file_names = file_names;
else
    orig_file_names = dir(fullfile(path_original_results,[file_prefix '*.mat']));
end

frames_dir = dir(fullfile(path_day_cnmfe_results,'frames*')); frames_dir = frames_dir.name;
log_dir = dir(fullfile(path_day_cnmfe_results,frames_dir,'LOGS*')); log_dir = log_dir(end).name;
matfiles = dir(fullfile(path_day_cnmfe_results,frames_dir,log_dir,'*.mat'));
matfile = matfiles(cellfun(@(x)strcmp(x,'intermediate_results.mat'),{d.name}) == 0).name;

load(fullfile(path_day_cnmfe_results,frames_dir,log_dir,matfile),'neuron');
vid_onsets = [vid_onsets; size(neuron.S,2)];

%build 3d matrix and dff for each file

deleted_frames_filler = zeros(neuron.options.d1,neuron.options.d2,n_del_frames);
disp('-------- Creating data files --------')
for file_num = 1:numel(orig_file_names)
    if isempty(path_original_results)
        tokens = split(orig_file_names{file_num},'_');
    else
        tokens = split(orig_file_names(file_num).name,'_');
    end
    curr_datetime = char(join(tokens(4:9),'_'));
    if strcmp(curr_datetime(end-3:end),'.mat')
        curr_datetime(end-3:end) = [];
    end
    curr_idx = find(cellfun(@(x)strcmp(x,curr_datetime),parsed_files_datetimes));
    start_frame_idx = vid_onsets(curr_idx);
    end_frame_idx = vid_onsets(curr_idx+1)-1;
    cnmfe_spikes = [zeros(size(neuron.S,1),n_del_frames) neuron.S(:,start_frame_idx:end_frame_idx)];
    cnmfe_vidTimes = [0:(size(cnmfe_spikes,2)-1)]/Fs;
    
    %save spike estimates etc
    S_file_name = fullfile(path_output,['cnmfe_spikes_' char(join(tokens(2:end),'_'))]);
    save(S_file_name,'cnmfe_vidTimes','cnmfe_spikes');
    disp(['Processed file ' num2str(file_num) ' of ' num2str(numel(orig_file_names))])
end


