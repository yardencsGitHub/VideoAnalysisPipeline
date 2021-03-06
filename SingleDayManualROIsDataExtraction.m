%%
% This script allows to choose a single max max proj. images. and pick ROIs
% from it. Then, it creates ROI annotation image, and calculates all ROI 
% extracted dff data from that day to be saved in per-day folders.
% it also creates folders for the images
bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation5_alexa'};
bird3_params = {'lbr3009' 'lbr3009' 'lbr3009_template_4TF' 'lbr3009_annotation_4TF'};
bird_params = bird3_params;
bird_name = bird_params{1}; 
bird_folder_name = bird_params{2}; 
template_file = bird_params{3}; 
annotation_file = bird_params{4}; 

last_idx = 0;
init_idx = 0;
last_date = '2017_02_29';
last_time = '00_00_00';
%bird_name = 'lrb85315'; %'lbr3009';
%bird_folder_name = 'lrb853_15'; % 'lbr3009';

% Folders on Data desktop:
%bird_name = 'lrb85315';
%bird_folder_name = 'lrb853_15';
%template_file = 'lrb85315template';
%annotation_file = 'lrb85315auto_annotation5_fix';



laptop_mov_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs'];
laptop_wav_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav'];
laptop_gif_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/gif'];
laptop_storage_folder = ['/Volumes/CanaryData/DATA/lrb853_15/movs/'];
laptop_annotated_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated'];
laptop_annotated_images_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated/images'];
DamagedFolder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/too_large_or_damaged/'];
laptop_manualROI_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/ManualROIs'];




desktop_mov_folder = ['/Users/yardenc/Documents/Experiments/Imaging/CanaryData/' bird_folder_name '/movs'];
desktop_storage_folder = ['/Volumes/home/Data/Imaging/' bird_folder_name '/RawData'];
desktop_max_projections_dir = ['/Users/yardenc/Documents/Experiments/Imaging/CanaryData/' bird_folder_name '/movs/MaxProj'];



addpath(genpath('/Users/yardenc/Documents/GitHub/VideoAnalysisPipeline'),'-end');
addpath(genpath('/Users/yardenc/Documents/GitHub/FreedomScope/Analysis Pipeline'),'-end');
clear dff;
n_del_frames = 5; % # of frames to ignore in the beginning of wach file in calculations
% filt_rad = 50; filt_sigma = 45; % highpass filter 
% h = fspecial('gaussian',filt_rad,filt_sigma);

filt_rad = 150; filt_sigma = 145; % highpass filter 
h = fspecial('gaussian',filt_rad,filt_sigma);

date_loc = [13 4]; % location of date string in file name w.r.t. end
MaxMaxProjDir = ['/Users/yardenc/Documents/Experiments/Imaging/CanaryData/' bird_folder_name '/movs/MaxProj/max_images'];
RawDataDir = ['/Volumes/home/Data/Imaging/' bird_folder_name '/RawData'];
ManualMaxProjROIsDir = ['/Users/yardenc/Documents/Experiments/Imaging/CanaryData/' bird_folder_name '/ManualROIs'];
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
    base = imfilter(Y,h,'circular','replicate');
%     base = imfilter(Y,h,'circular','replicate');
%     Y = (Y - base);
%     Y = Y - min(min(min((Y(:,:,n_del_frames+1:end)))));
%     Y = Y / max(max(max((Y(:,:,n_del_frames+1:end)))));
%     b1 = quantile(Y(:,:,n_del_frames+1:end),0.05,3);
%     b1 = (b1-min(b1(:)))/(max(b1(:))-min(b1(:)));
%     b1 = b1/median(b1(:));
%     
%     
%     bt = median(reshape(Y,rows*columns,frames));
%     b2 = bsxfun(@times,b1,reshape(bt,1,1,numel(bt)));
%     c = Y-b2; c = c - min(min(min(c(:,:,n_del_frames+1:end))));
%     %dff_temp = bsxfun(@rdivide,bsxfun(@minus,c,quantile(c,0.05,3)),quantile(c,0.05,3));
%     
%     Y = reshape(c,rows*columns,frames);
    
    c = Y-base;
    c = c-min(c(:));
    Y = reshape(c,rows*columns,frames);
    for roi_num = 1:length(ROI.stats)
        cr = ROI.coordinates{roi_num}(:,2)+(ROI.coordinates{roi_num}(:,1)-1)*rows;
        roi_signal = sum(Y(cr,:));
        
        dff(roi_num,:) = (roi_signal - quantile(roi_signal(n_del_frames+1:end),0.05))/quantile(roi_signal(n_del_frames+1:end),0.05);
        
    end
    
    data_filename = ['baseROIdata' filename(8:end)];
    save(data_filename,'dff','vidTimes');
    clear dff;
end

