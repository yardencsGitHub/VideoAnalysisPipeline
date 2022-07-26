function helper_function_reshape_single_day_cnmfe_result_as_matrices(path_day_cnmfe_results,path_params_file,path_original_results,path_output,varargin)
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

orig_file_names = dir(fullfile(path_original_results,[file_prefix '*.mat']));

frames_dir = dir(fullfile(path_day_cnmfe_results,'frames*')); frames_dir = frames_dir.name;
log_dir = dir(fullfile(path_day_cnmfe_results,frames_dir,'LOGS*')); log_dir = log_dir(end).name;
matfiles = dir(fullfile(path_day_cnmfe_results,frames_dir,log_dir,'*.mat'));
matfile = matfiles(cellfun(@(x)strcmp(x,'intermediate_results.mat'),{d.name}) == 0).name;

load(fullfile(path_day_cnmfe_results,frames_dir,log_dir,matfile),'neuron');
vid_onsets = [vid_onsets; size(neuron.C,2)];

%build 3d matrix and dff for each file

deleted_frames_filler = zeros(neuron.options.d1,neuron.options.d2,n_del_frames);
disp('-------- Creating data files --------')
for file_num = 1:numel(orig_file_names)
    tokens = split(orig_file_names(file_num).name,'_');
    curr_datetime = char(join(tokens(4:9),'_'));
    if strcmp(curr_datetime(end-3:end),'.mat')
        curr_datetime(end-3:end) = [];
    end
    curr_idx = find(cellfun(@(x)strcmp(x,curr_datetime),parsed_files_datetimes));
    start_frame_idx = vid_onsets(curr_idx);
    end_frame_idx = vid_onsets(curr_idx+1)-1;
    cnmfe_vidMat = cat(3,deleted_frames_filler,neuron.reshape(neuron.A*neuron.C(:,start_frame_idx:end_frame_idx),2));
    cnmfe_dff = [zeros(size(neuron.C,1),n_del_frames) neuron.C(:,start_frame_idx:end_frame_idx)];
    cnmfe_vidTimes = [0:(size(cnmfe_dff,2)-1)]/Fs;
    
    %save dff etc
    dff_file_name = fullfile(path_output,['cnmfe_dff_' char(join(tokens(2:end),'_'))]);
    save(dff_file_name,'cnmfe_vidTimes','cnmfe_dff');

    %save vidMat
    dff_file_name = fullfile(path_output,['cnmfe_vidMat_' char(join(tokens(2:end),'_'))]);
    save(dff_file_name,'cnmfe_vidTimes','cnmfe_vidMat');
    disp(['Processed file ' num2str(file_num) ' of ' num2str(numel(orig_file_names))])
end
disp('-------- Saving ROI files --------')
%save contour image and contours
clear ROI;
ROI.coordinates = neuron.show_contours(0.6,[],[],1);
ROI.type = 'image';
haxes = get(gcf,'Children');
ROI.reference_image = uint8(256*haxes.Children(end).CData);
ROI_file_name = fullfile(path_output,'cnmfe_ROI_struct.mat');
save(ROI_file_name,'ROI','-v7.3');
contour_fig_file_name = fullfile(path_output,'cnmfe_ROI_image.fig');
hgsave(gcf,contour_fig_file_name);
contour_image_file_name = fullfile(path_output,'cnmfe_ROI_image.png');
saveas(gcf,contour_image_file_name);


