function X = helper_function_reshape_single_day_cnmfe_result_as_footprints(path_day_cnmfe_results,varargin)
% check that we're in a valid folder
d = dir(fullfile(path_day_cnmfe_results,'parsing_record.mat'));
if isempty(d)
    disp('Not a valid CNMFE results folder');
    return;
end
xpixels = 240; ypixels = 320;
GitHubFolder = '/Users/yardenc/Documents/GitHub/';
nparams = numel(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
		case 'xpixels'
			xpixels=varargin{i+1};
        case 'ypixels'
			ypixels=varargin{i+1};
        case 'githubfolder'
			GitHubFolder=varargin{i+1};
    end
end

addpath(genpath(fullfile(GitHubFolder,'CNMF_E')))
frames_dir = dir(fullfile(path_day_cnmfe_results,'frames*')); frames_dir = frames_dir.name;
log_dir = dir(fullfile(path_day_cnmfe_results,frames_dir,'LOGS*')); log_dir = log_dir(end).name;
matfiles = dir(fullfile(path_day_cnmfe_results,frames_dir,log_dir,'*.mat'));
matfile = matfiles(cellfun(@(x)strcmp(x,'intermediate_results.mat'),{d.name}) == 0).name;
load(fullfile(path_day_cnmfe_results,frames_dir,log_dir,matfile),'neuron');

num_neurons = size(neuron.A,2);
X = zeros(num_neurons,xpixels,ypixels);

for neuron_n = 1:num_neurons
    x = reshape(neuron.A(:,neuron_n),xpixels,ypixels);
    X(neuron_n,:,:) = x/max(x(:));
end
end