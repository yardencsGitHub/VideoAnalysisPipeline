function vid_onsets = helper_function_create_concatenated_vidMat(Day,varargin)
%%
% Concatenate videos from a full day
% Make a movie that cnmfe can load.
%Day = '2017-05-16';
GitHubFolder = '/Users/yardenc/Documents/GitHub/';
DataFolder = '/Volumes/Labs/cohen/yardenc/DATA/lrb853_15/movs/';
output_filename = 'Day_full.mat';
output_folder = '/Users/yardenc/Desktop/temp';
max_frames = 1900000;
n_del_frames = 0;
resize_factor = 0.5;

nparams = length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
		case 'githubfolder'
			GitHubFolder=varargin{i+1};
        case 'datafolder'
			DataFolder=varargin{i+1};
        case 'output_filename'
			output_filename=varargin{i+1};
        case 'output_folder'
			output_folder=varargin{i+1};
        case 'max_frames'
			max_frames=varargin{i+1};
        case 'n_del_frames'
			n_del_frames=varargin{i+1};
        case 'resize_factor'
			resize_factor=varargin{i+1};
    end
end

%%
addpath(fullfile(GitHubFolder,'VideoAnalysisPipeline/helpers'));
addpath(genpath(fullfile(GitHubFolder,'/FinchScope/Analysis Pipeline')));
files = dir(fullfile(DataFolder,[Day '*.mov']));
%%
vid_onsets = [1];
cd(output_folder);
% first file
filename = fullfile(files(1).folder,files(1).name);
[vidMat, ~, ~] = Prepare_SingleFile_Raw_Video_Audio(filename); vidMat = imresize(vidMat,resize_factor);
Y = uint16(squeeze(vidMat(:,:,n_del_frames+1:end)));
%imwrite(squeeze(uint16(squeeze(vidMat(:,:,n_del_frames+1)))),output_filename);
%for i=(n_del_frames+2):size(vidMat,3)
%    imwrite(squeeze(uint16(vidMat(:,:,i))),output_filename,'WriteMode','append');
%end
%disp('done writing first file to tiff');
%num_frames = size(vidMat,3) - n_del_frames;

for filenum = 2:numel(files)
    if size(Y,3) > max_frames
        break;
    end
    vid_onsets = [vid_onsets; size(Y,3)+1];
    disp(['writing file ' num2str(filenum) ' / ' num2str(numel(files))])
    
    filename = fullfile(files(filenum).folder,files(filenum).name);
    [vidMat, ~, ~] = Prepare_SingleFile_Raw_Video_Audio(filename); vidMat = imresize(vidMat,resize_factor);
    Y = cat(3,Y,uint16(squeeze(vidMat(:,:,n_del_frames+1:end))));
%     for i=(n_del_frames+1):size(vidMat,3)
%         imwrite(squeeze(uint16(vidMat(:,:,i))),output_filename,'WriteMode','append');
%     end
%     num_frames = num_frames + size(vidMat,3) - n_del_frames;
end
disp('done making movie')
Ysiz = size(Y);

save(fullfile(output_folder,output_filename),'Y','Ysiz','-v7.3');


%%