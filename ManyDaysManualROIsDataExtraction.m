%%
% This script repeats the ROI extraction for already defined masks 
% in order to allow changing the bg removal and normalization

last_idx = 0;
init_idx = 0;
last_date = '2017_02_29';
last_time = '00_00_00';
bird_name = 'lrb85315'; %'lbr3022'; %'lbr3009';
bird_folder_name = 'lrb853_15'; %'lbr3022'; % 'lbr3009';

% Folders on Data desktop:
desktop_mov_folder = ['/Users/yardenc/Documents/Experiments/Imaging/CanaryData/' bird_folder_name '/movs'];
desktop_storage_folder = ['/Volumes/home/Data/Imaging/' bird_folder_name '/RawData'];
desktop_max_projections_dir = ['/Users/yardenc/Documents/Experiments/Imaging/CanaryData/' bird_folder_name '/movs/MaxProj'];



addpath(genpath('/Users/yardenc/Documents/GitHub/VideoAnalysisPipeline'));
addpath(genpath('/Users/yardenc/Documents/GitHub/FreedomScope/Analysis Pipeline'));


n_del_frames = 6; % # of frames to ignore in the beginning of wach file in calculations
filt_rad = 150; filt_sigma = 145; % highpass filter 
h = fspecial('gaussian',filt_rad,filt_sigma);
date_loc = [13 4]; % location of date string in file name w.r.t. end
MaxMaxProjDir = ['/Users/yardenc/Documents/Experiments/Imaging/CanaryData/' bird_folder_name '/movs/MaxProj/max_images'];
RawDataDir = ['/Volumes/home/Data/Imaging/' bird_folder_name '/RawData'];
%['/Volumes/CanaryData/DATA/' bird_folder_name '/RawData']; %temporary.. until guttata comes back
%
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
% [FILENAME, PATHNAME, FILTERINDEX] = uigetfile(fullfile(MaxMaxProjDir,'*.png'), 'choose day');
% daynum = datenum(FILENAME(end-date_loc(1):end-date_loc(2)));
% daystr = FILENAME(end-date_loc(1):end-date_loc(2));
% a = 0;
% if exist(fullfile(pwd,'ROIdata',daystr,['ROI_' daystr '.mat']))
%     a = input('use old ROI? ');
% end
dirs = dir(fullfile(MaxMaxProjDir,'*.png'));
for dirnum = 1:numel(dirs)
    FILENAME = dirs(dirnum).name;
    daynum = datenum(FILENAME(end-date_loc(1):end-date_loc(2)));
    daystr = FILENAME(end-date_loc(1):end-date_loc(2));

    clear dff;
    cd(ManualMaxProjROIsDir);
    cd('ROIdata');
    if ~exist(daystr,'dir')
        mkdir(daystr);
    end
    cd(daystr);
    roi_filename = ['nonoverlap_newROI_' daystr];
    try
        load(fullfile(pwd,['nonoverlap_newROI_' daystr '.mat']));
    
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
        %     Y = (Y - base);
        %     Y = Y - min(min(min((Y(:,:,n_del_frames+1:end)))));
        %     Y = Y / max(max(max((Y(:,:,n_del_frames+1:end)))));
        %%% The median based removal    
%         b1 = quantile(Y(:,:,n_del_frames+1:end),0.05,3);
%         b1 = (b1-min(b1(:)))/(max(b1(:))-min(b1(:)));
%         b1 = b1/median(b1(:));
%          bt = median(reshape(Y,rows*columns,frames));
%         b2 = bsxfun(@times,b1,reshape(bt,1,1,numel(bt)));
%         c = Y-b2; %c = c - min(min(min(c(:,:,n_del_frames+1:end))));
        %%% end of median based 
            %dff_temp = bsxfun(@rdivide,bsxfun(@minus,c,quantile(c,0.05,3)),quantile(c,0.05,3));
            c = Y-base;
            c = c-min(c(:));
            Y = reshape(c,rows*columns,frames);
            for roi_num = 1:length(ROI.stats)
                cr = ROI.coordinates{roi_num}(:,2)+(ROI.coordinates{roi_num}(:,1)-1)*rows;
                roi_signal = sum(Y(cr,:));

                dff(roi_num,:) = (roi_signal - quantile(roi_signal(n_del_frames+1:end),0.05))/quantile(roi_signal(n_del_frames+1:end),0.05);

            end

            data_filename = ['NonoverlapBaseROIdata' filename(8:end)];
            save(data_filename,'dff','vidTimes');
            clear dff;
        end
    catch em
        display(fullfile(pwd,['ROI_' daystr '.mat']));
    end
end
