%%
% This script allows to choose a single max max proj. images. and pick ROIs
% from it. Then, it creates ROI annotation image, and calculates all ROI 
% extracted dff data from that day to be saved in per-day folders.
% it also creates folders for the images

last_idx = 0;
init_idx = 0;
last_date = '2017_02_29';
last_time = '00_00_00';
bird_name = 'lbr3022';
bird_folder_name = 'lbr3022';

% Folders on Data desktop:
desktop_mov_folder = ['/Users/yardenc/Documents/Experiments/Imaging/CanaryData/' bird_folder_name '/movs'];
desktop_storage_folder = ['/Volumes/home/Data/Imaging/' bird_folder_name '/RawData'];
desktop_max_projections_dir = ['/Users/yardenc/Documents/Experiments/Imaging/CanaryData/' bird_folder_name '/movs/MaxProj'];



addpath(genpath('/Users/yardenc/Documents/GitHub/VideoAnalysisPipeline'));
addpath(genpath('/Users/yardenc/Documents/GitHub/FreedomScope/Analysis Pipeline'));
clear dff;
n_del_frames = 5; % # of frames to ignore in the beginning of wach file in calculations
filt_rad = 50; filt_sigma = 45; % highpass filter 
h = fspecial('gaussian',filt_rad,filt_sigma);
date_loc = [13 4]; % location of date string in file name w.r.t. end
MaxMaxProjDir = ['/Users/yardenc/Documents/Experiments/Imaging/CanaryData/' bird_name '/movs/MaxProj/max_images'];
RawDataDir = ['/Volumes/home/Data/Imaging/' bird_name '/RawData'];
ManualMaxProjROIsDir = ['/Users/yardenc/Documents/Experiments/Imaging/CanaryData/' bird_name '/ManualROIs'];
% create folders
if ~exist(ManualMaxProjROIsDir,'dir')
    mkdir(ManualMaxProjROIsDir);
end
cd(ManualMaxProjROIsDir);
if ~exist('ROIimages','dir')
    mkdir('ROIimages');
end
if ~exist('ROIdata','dir')
    mkdir('ROIdata');
end

% pick file and extract ROIs manually
[FILENAME, PATHNAME, FILTERINDEX] = uigetfile(fullfile(MaxMaxProjDir,'*.png'), 'choose day');
daynum = datenum(FILENAME(end-date_loc(1):end-date_loc(2)));
daystr = FILENAME(end-date_loc(1):end-date_loc(2));
a = 0;
if exist(fullfile(pwd,'ROIdata',daystr,['ROI_' daystr '.mat']))
    a = input('use old ROI? ');
end
cd('ROIdata');
if ~exist(daystr,'dir');
    mkdir(daystr);
end
cd(daystr);
roi_filename = ['ROI_' daystr];
if (a == 1)
    load(fullfile(pwd,['ROI_' daystr '.mat']));
else
    I = imread(fullfile(PATHNAME,FILENAME));
    I = imresize(I,2);
    [ROI fig_handle]=FS_Image_ROI(I,'save_dir','');
    hgsave(fig_handle,fullfile(ManualMaxProjROIsDir,'ROIimages',[FILENAME(1:end-3) 'mat']));
    saveas(fig_handle,fullfile(ManualMaxProjROIsDir,'ROIimages',[FILENAME(1:end-3) 'png']));
% extract dff traces from all files from that day
    
    save(roi_filename,'ROI');
end
FILES = dir(fullfile(RawDataDir,'*.mat'));
days = [];
for fnum = 1:numel(FILES)
    tokens = regexp(FILES(fnum).name,'_','split');
    days = [days; datenum([tokens{4} '-' tokens{5} '-' tokens{6}])];
end
days_to_process = find(days == daynum);

for day_cnt = 1:numel(days_to_process)
    curr_day = days_to_process(day_cnt);
    filename = FILES(curr_day).name;
    load(fullfile(RawDataDir,filename),'vidMat','vidTimes');
    [rows,columns,frames] = size(vidMat);
    
    Y = vidMat;
%     base = imfilter(Y,h,'circular','replicate');
%     Y = (Y - base);
%     Y = Y - min(min(min((Y(:,:,n_del_frames+1:end)))));
%     Y = Y / max(max(max((Y(:,:,n_del_frames+1:end)))));
    b1 = quantile(Y(:,:,n_del_frames+1:end),0.05,3);
    b1 = (b1-min(b1(:)))/max(b1(:));
    b1 = b1/median(b1(:));
    
    
    bt = median(reshape(Y,rows*columns,frames));
    b2 = bsxfun(@times,b1,reshape(bt,1,1,numel(bt)));
    c = Y-b2; c = c - min(min(min(c(:,:,n_del_frames+1:end))));
    %dff_temp = bsxfun(@rdivide,bsxfun(@minus,c,quantile(c,0.05,3)),quantile(c,0.05,3));
    
    Y = reshape(c,rows*columns,frames);
    for roi_num = 1:length(ROI.stats)
        cr = ROI.coordinates{roi_num}(:,2)+(ROI.coordinates{roi_num}(:,1)-1)*rows;
        roi_signal = sum(Y(cr,:));
        
        dff(roi_num,:) = (roi_signal - quantile(roi_signal(n_del_frames+1:end),0.05))/quantile(roi_signal(n_del_frames+1:end),0.05);
        
    end
    
    data_filename = ['ROIdata' filename(8:end)];
    save(data_filename,'dff','vidTimes');
    clear dff;
end

