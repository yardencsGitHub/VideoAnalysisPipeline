function [hndls,syl_types] = ColorOverlayMaxProjMaps_function(Day,sylidx,syllabels_sequence,order_flag,tags_to_color,base_colors,min_pixel_value,varargin)
%% 
% This script creates single ROI single day alignments to complex sequences
% with only one variable (duration or syllable type). 
% This script requires that max. projection images are prepared using
% 'CreatePhraseMaxProjections.m' helper function of the
% VideoAnalysisPipeLine repo
%
% Inputs:
%   Day - text rep. of date 
%   if use_cohen2020 is set to 1 (see below) than this is also use in
%   pointing to data folders. Otherwise it is used for figure titles only
%   (but it must be supplied .. for safety)
%   sylidx - index of syllable in the sequence to align to
%   syllabels_sequence - the sequence of phrase identities to lock. Insert
%       nan to keep it free
%   ROI - the roi # 
%   locktoonset - 0/1 to lock to phrase onset or offset
%   spikes - type of variable (0-3):
%       0 - denoised Ca
%       1 - deconvolved Spikes
%       2 - Fluorescence
%       3 - HMM state (signal \ noise)
%   order_flag - the correlate. Use positive integers to indicate durations of
%   phrases. A vector of positive integers will result in summation of
%   durations. Use negative integers for phrase type (no vectors)
%   tags_to_color - which tags to use in creating the overlay
%   base_colors - a cell array of rgb values corresponding to
%   'tags_to_color'. for example {[0 0 1], [1,1,0]} will be red and cyan
%   min_pixel_value - lower values are set to zero
%
% Additional optional inputs:
%   use_cohen2020 - legacy. set to 1 to supply parameters as in the code developed for "Hidden neural states underlie canary song syntax" (2020). 
%   ignore_entries - A vector of label numbers to ignore completely. 
%   join_entries - A cell array of vectors, each containing a >1 number of
%                  labels to be grouped
%   githubdir - where required repositories live
%   display_opt - set to 1 to print the list of participating files and
%                 phrases.
%   bird_number - if 'bird_params' is not set, set to 1-3 to use data prepared for Cohen+al 2020
%   bird_params - this cell array will be all the parameters: 
%               bird_name = bird_params{1}; 
%               bird_folder_name = bird_params{2}; can be different than
%               the birds name
%               path to template file = bird_params{3}; 
%               path to annotation file = bird_params{4}; 
%               path to ROIs data + file_prefix = bird_params{5};  (e.g.
%               /folder/folder/prefix_)
%               path to roi map = bird_params{6};
%               path to max. projection files + prefix = bird_params{7}
%
%           if use_cohen2022 is set to 1 these will be overriden by the
%           parameters in the code below.
%
%   rois - which ROIs to highlight or mark (default:0)
%   set_bg_to_level_zero - quantile of intensity level set to 0 (default 0)
%   intensity_gain - global multiplier of intensity (default 1.0)
%   nonlinear_stretch - set to a 2-vector [a,b] to strech color saturation
%   equally:
%            a = saturation_sigm_threshold
%            b = saturation_sigm_power
%            the saturation (s) stretch: s = s/(1+exp(-(s-a)*b))
%   patch_coor - where to place the scale bar
%   
%% Optional parameters
use_cohen2020 = 0;
BaseDir = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/';
GithubDir = '/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub';
display_opt = 0;
ignore_entries = [-1];
join_entries = {};
% The dataset from "Hidden neural states underlie canary song syntax" (Cohen+al 2020 Nature)
bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix' 'NonoverlapBaseROIdata_' 'nonoverlap_newROI_' 'none'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation5_alexa' 'baseROIdata_' 'ROI_' 'none'};
bird3_params = {'lbr3009' 'lbr3009' 'lbr3009_template_4TF' 'lbr3009auto_annotation1_fix' 'baseROIdata_' 'ROI_' 'none'};
bird_params = bird2_params;
max_phrase_gap = 0.5;
bird_number = 1;
ROIs = 0;
set_bg_to_level_zero = 0;
intensity_gain = 1.0;
saturation_sigm_threshold = 0.5;
saturation_sigm_power = 10;
nonlinear_stretch = 0;
patch_coor = [300 300];
show_color_bars = 0;
show_scale_bar = 0;
%% allow controlling parameters as function pair inputs
nparams=length(varargin);
for i=1:2:nparams
	switch lower(varargin{i})
        case 'use_cohen2020'
            use_cohen2020 = varargin{i+1};
        case 'bird_number'
             bird_number = varargin{i+1};
             if ~ismember('bird_params',varargin)
                
			         switch varargin{i+1}
                         case 1
                             bird_params = bird1_params;
                         case 2
                             bird_params = bird2_params;
                         case 3
                             bird_params = bird3_params;
                     end 
             end
        case 'bird_params'
            bird_params = varargin{i+1};
        case 'max_phrase_gap'
			max_phrase_gap=varargin{i+1}; 
        case 'display_opt'
			display_opt=varargin{i+1}; 
        case 'rois'
            ROIs = varargin{i+1};
        case 'set_bg_to_level_zero'
            set_bg_to_level_zero = varargin{i+1};
        case 'intensity_gain'
            intensity_gain = varargin{i+1};
        case 'nonlinear_stretch'
            nonlinear_stretch = 1;
            temp_stretch = varargin{i+1};
            saturation_sigm_threshold = temp_stretch(1);
            saturation_sigm_power = temp_stretch(2);
        case 'patch_coor'
            patch_coor = varargin{i+1};
        case 'githubdir'
            GithubDir = varargin{i+1};
        case 'ignore_entries'
            ignore_entries = varargin{i+1};
        case 'join_entries'
            join_entries = varargin{i+1};
        case 'ignore_dates'
            ignore_dates = varargin{i+1};   
        case 'show_color_bars'
            show_color_bars = varargin{i+1};      
        case 'show_scale_bar'
            show_scale_bar = varargin{i+1};
    end
end
            

%%
%addpath(genpath([GithubDir 'small-utils']),'-end');
addpath(genpath([GithubDir 'VideoAnalysisPipeline']),'-end');
addpath(genpath([GithubDir 'BirdSongBout']),'-end');
bird_name = bird_params{1}; 
bird_folder_name = bird_params{2}; 
template_file = bird_params{3}; 
annotation_file = bird_params{4}; 
file_prefix = bird_params{5}; 
roi_map_prefix = bird_params{6}; 
maxproj_file_prefix = bird_params{7};


%% Folders that contain data
if use_cohen2020
    manualROI_folder = [BaseDir bird_folder_name '/ManualROIs'];
    maxproj_phrases_folder = [BaseDir bird_folder_name '/movs/PhraseMaxProj/'];
end

%%
flag = 0;
join_entries = join_entries(:);
if ~isempty(join_entries)
    for i = 1:numel(join_entries)
        if ~isempty(intersect(join_entries{i},ignore_entries))
            flag = 1;
        end
        for j = i+1:numel(join_entries)
            if  ~isempty(intersect(join_entries{i},join_entries{j}))
                flag = 1;
            end
        end
    end
end
   
if flag == 1
    r = []; p = [];
    disp(['join or ignore lists overlap'])
    return;
end

if (any(order_flag == 0) || ...
   (any(order_flag < 0) & numel(order_flag) > 1) || ...
   any(order_flag > numel(syllabels_sequence)))
    r = []; p = [];
    disp(['Illegal order format'])
    return;
end
      

%%
AlphaNumeric = '{}ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';
target_sequence_str = '';
if use_cohen2020
    cd(manualROI_folder);
end
%load('syl_dur_gap_stat.mat');
load(template_file,'templates');
%syllables = [[templates.wavs.segType] -1 102 103];
syllables = [templates.wavs.segType];
syllables = [[-1000 1000] syllables setdiff([-1 102 103],syllables)];
n_syllables = numel(syllables)-2;
for cnt = 1:numel(syllabels_sequence)
    if isnan(syllabels_sequence(cnt))
        target_sequence_str = [target_sequence_str '.'];
    else
        target_sequence_str = [target_sequence_str AlphaNumeric(syllables == syllabels_sequence(cnt))];
    end
end

freq_min = 300; freq_max = 8000;
colors = distinguishable_colors(n_syllables,'w');
colors = [0.25 0.25 0.25; 0.25 0.25 0.25 ;colors];
load(annotation_file);
ord = [];
dates = [];
for i = 1:numel(keys)
    tokens = regexp(keys{i},'_','split');
    ord = [ord; str2num(tokens{2})];
    dates = [dates; char(join(tokens(3:5),'_'))];
end
[locs,indx] = sort(ord);
elements = elements(indx);
keys = keys(indx);
dates = dates(indx,:);
%%
if use_cohen2020
    cd([manualROI_folder '/ROIdata/' Day]);
end
FILES = dir([file_prefix bird_name '*.mat']);

FILES = {FILES.name};
hits = [];
max_syls = 0;
durations = []; % durations of selected elements only
sylidx_durations = []; % duration of all numbered phrases in sequence and the onset and offset of the target phrase
sequence_total_lengths = []; % the total duration, including gaps, of the entire sequence
sequence_onsets = []; % durations of all phrases in sequence
for fnum = 1:numel(FILES)
    fname = FILES{fnum};
    tokens = regexp(fname,'_','split');
    if use_cohen2020
        loc = find(locs == str2num(tokens{3}));
    else
        name_idx = find(cellfun(@(x)strcmp(x,bird_name),tokens));
        loc = find(locs == str2num(tokens{name_idx+1}));
    end
    
    try
        ignore_locs = find(ismember(elements{loc}.segType,ignore_entries));
    catch eem
        continue;
    end
    elements{loc}.segAbsStartTimes(ignore_locs) = [];
    elements{loc}.segFileStartTimes(ignore_locs) = [];
    elements{loc}.segFileEndTimes(ignore_locs) = [];
    elements{loc}.segType(ignore_locs) = [];  
    for i = 1:numel(join_entries)
        join_locs = find(ismember(elements{loc}.segType,join_entries{i}));
        elements{loc}.segType(join_locs) = join_entries{i}(1);
    end
    phrases = return_phrase_times(elements{loc});
    try   
        phrases = deal_with_time_gaps(phrases,max_phrase_gap);
    catch em
        'f';
    end
    phrases_str = cell2mat(arrayfun(@(x)AlphaNumeric(syllables == x),phrases.phraseType','UniformOutput',false));
    phrase_durations = phrases.phraseFileEndTimes - phrases.phraseFileStartTimes;
    gap_durations = phrases.phraseFileStartTimes(2:end) - phrases.phraseFileEndTimes(1:end-1);
    
    
    try
        endpoint = regexp(phrases_str,target_sequence_str,'end');
    catch em
        endpoint = [];
    end
    if ~isempty(endpoint)
        
        phrase_locs = endpoint - numel(target_sequence_str) + sylidx;
        
        for phrase_loc = 1:numel(phrase_locs)
            phrase_idx = endpoint(phrase_loc) - numel(target_sequence_str) + [1:numel(target_sequence_str)];
            try
            if ~any(gap_durations(phrase_idx(1:end-1)) > max_phrase_gap)
                
                phrasenum = phrase_locs(phrase_loc);
                tonset = phrases.phraseFileStartTimes(phrasenum);
                toffset = phrases.phraseFileEndTimes(phrasenum);
                syls_in_phrase = find(elements{loc}.segFileEndTimes >= tonset & ...
                elements{loc}.segFileStartTimes <= toffset);
                max_syls = max(max_syls,numel(syls_in_phrase));
                hits = [hits; loc phrasenum phrases.phraseFileStartTimes(phrase_idx(1)) phrases.phraseType(phrase_idx)'];
                sylidx_durations = [sylidx_durations; phrase_durations(phrase_idx(~isnan(syllabels_sequence))) tonset toffset]; %(sylidx)       
                durations = [durations; sum(phrase_durations(phrase_idx(sylidx)))];
                sequence_onsets = [sequence_onsets; (phrases.phraseFileStartTimes(phrase_idx)-tonset)*locktoonset + ...
                    (phrases.phraseFileEndTimes(phrase_idx)-toffset)*(1-locktoonset)];
                sequence_total_lengths = [sequence_total_lengths; phrases.phraseFileEndTimes(phrase_idx(end)) - phrases.phraseFileStartTimes(phrase_idx(1))];
                
            end
            catch em
                '4';
            end
            
        end
    end
end

try 
    [id_flags,dur_idx] = sort(hits(:,3 - order_flag));
    durations = durations(dur_idx);
    hits = hits(dur_idx,:);
    sylidx_durations = sylidx_durations(dur_idx,:);
catch em
    'f';
    p=1; r=0; gnames=''; return;
end


%%

clear MaxIm;
syl_types = unique(id_flags);
for syl_type_cnt = 1:numel(syl_types)
    MaxIm{syl_type_cnt} = [];
end
%cd(maxproj_phrases_folder);
for cnt = 1:size(hits,1)
    fnum = hits(cnt,1);
    phrasenum = hits(cnt,2);
    fname = keys{fnum};
    if use_cohen2020
        load(fullfile(maxproj_phrases_folder,['PhraseMaxProj_' fname(1:end-3) 'mat']));
    else
        load([maxproj_file_prefix fname(1:end-3) 'mat']);
    end
    phrases1 = return_phrase_times(elements{fnum});
    phrases1 = deal_with_time_gaps(phrases1,max_phrase_gap);
    try
    tonset = phrases.phraseFileStartTimes(phrasenum);
    toffset = phrases.phraseFileEndTimes(phrasenum);
    
    if (display_opt)
        display([num2str([cnt tonset]) ' ' fname]);
    end
    
    MaxIm{syl_types == id_flags(cnt)} = cat(3,MaxIm{syl_types == id_flags(cnt)},phrases.maxImages(:,:,phrasenum));
    catch eem
        'f';
        id_flags(cnt) = -1;
    end
         
end
id_flags(id_flags == -1) = [];
%% Now create images
hndls = [];
if use_cohen2020
    roi_map_filename = [manualROI_folder '/ROIdata/' Day '/' roi_map_prefix Day '.mat'];
    load(roi_map_filename);
else
    load(roi_map_prefix);
end
mask = zeros(size(ROI.reference_image));
roimask = mask;
[nrows,ncols] = size(ROI.reference_image);
textlocs = [];
roi_fg = figure('Position',[1308         682         770/2         551/2]); 
roi_ax = axes(roi_fg);
for roi_n = 1:numel(ROI.coordinates)
    temp_coor = ROI.coordinates{roi_n};
    temp_s_coor = size(temp_coor);
    if temp_s_coor(1) < temp_s_coor(2)
        temp_coor = transpose(temp_coor);
    end
    roimask(round(temp_coor(:,1)*nrows+temp_coor(:,2))) = 1*ismember(roi_n,ROIs);
    mask(round(temp_coor(:,1)*nrows+temp_coor(:,2)))=1;
    textlocs = [textlocs; mean(temp_coor(:,1)) mean(temp_coor(:,2))];
    if ~use_cohen2020
        plot(roi_ax,temp_coor(:,1),temp_coor(:,2),'Color',[0.5 0.5 0.5],'LineWidth',1);
        hold on;
    end
end
if use_cohen2020
    figure('Position',[1308         682         770         551]); 
    imagesc(mask); hold on; colormap(1-gray)
end
for roi_n = 1:numel(ROI.coordinates)
    text(roi_ax,textlocs(roi_n,1),textlocs(roi_n,2),num2str(roi_n),'FontSize',16,'HorizontalAlignment','center','Color','k');
end
if use_cohen2020
    contour(mask,[0.5 0.5],'g');
end
xticks([]); yticks([]); set(gca,'Position',[0.0842    0.0871    0.8312    0.8711]); %0.8312    0.8711
if use_cohen2020
    hndls = [hndls gca];
else
    hndls = [hndls roi_ax];
    roi_ax.XLim = [0 ncols];
    roi_ax.YLim = [0 nrows];
    roi_ax.YDir = 'reverse';
end
% create and plot the overlay image
filt_rad = 10; filt_sigma = 3; % lowpass filter (TODO: this should not be hard coded)
filt_rad = 5; filt_sigma = 1.5;
h = fspecial('gaussian',filt_rad,filt_sigma);
mn = cat(3,MaxIm{:});
mn = imfilter(mn,h,'circular','replicate');
%mn = zscore(mn,0,3); zscoring is a bad idea, IT CREATES BIAS
I = zeros(nrows,ncols,3); 
maxnorm = 0;
for layernum = 1:numel(tags_to_color)
    Iaddition = nanmean(mn(:,:,ismember(id_flags,tags_to_color{layernum})),3); % mean or median?
    if set_bg_to_level_zero ~= 0
        Iaddition = Iaddition - quantile(Iaddition(:),set_bg_to_level_zero); %Iaddition(Iaddition(:)<0)=0;
    end
    maxnorm = max(maxnorm,max(Iaddition(:)));
    I = I + repmat(Iaddition,1,1,3).*repmat(reshape(base_colors{layernum},1,1,3),nrows,ncols);    
end
%I = I-min(I(:));
%I = I - min(min(sqrt(sum(I.*I,3))));
I = I/maxnorm; I(I<0)=0;
%min_pixel_value = quantile(I(:),min_pixel_value);
min_pixel_value = quantile(reshape(sqrt(sum(I.*I,3)),1,nrows*ncols),min_pixel_value);

%I(sqrt(sum(I.*I,3)) < min_pixel_value,:) = 0;
I = I.*repmat(sqrt(sum(I.*I,3)) > min_pixel_value,1,1,3);
I = I - min_pixel_value/sqrt(3); I(I<0)=0;
% stretch colors uniformly

I = I/max(max(sqrt(sum(I.*I,3))))*sqrt(3)*intensity_gain; %
if (nonlinear_stretch == 1)
    I1 = rgb2hsv(I); I1(:,:,2) = 1./(1+exp(-(I1(:,:,2)-saturation_sigm_threshold)*saturation_sigm_power));
    I = hsv2rgb(I1);
end
im_fh = figure('Position',[1308         682         770         551]); im_ax = axes(im_fh); imshow(I,'Parent',im_ax); %/max(I(:))*1.0);
hold on;
if ~ismember(0,ROIs)
    if use_cohen2020
        contour(roimask,[0.5 0.5],'Color',[0.5 0.5 0.5]);
    end
    for roi_n = 1:numel(ROI.coordinates)
        if ismember(roi_n,ROIs)
            text(im_ax,textlocs(roi_n,1),textlocs(roi_n,2),num2str(roi_n),'FontSize',48,'HorizontalAlignment','center','Color',[0.5 0.5 0.5]);
        end
    end
end
% for roi_n = 1:numel(ROI.coordinates)
%     temp_coor = ROI.coordinates{roi_n};
%     temp_s_coor = size(temp_coor);
%     if temp_s_coor(1) < temp_s_coor(2)
%         temp_coor = transpose(temp_coor);
%     end
%     plot(temp_coor(:,1),temp_coor(:,2),'Color','w');
% end
if show_scale_bar
    patch(im_ax,patch_coor(1)+[0 50*5/6 50*5/6 0],patch_coor(2)+ [0 0 5 5],'w');
end
xticks([]); yticks([]);
if use_cohen2020
    hndls = [hndls gca];
else
    hndls = [hndls im_ax];
end
if show_color_bars
    bar_fh = figure; bar_ax = axes(bar_fh); %color bars
    I1 = I/sqrt(3);%/max(I(:));
    maxcol = max(max(sqrt(sum(I1.^2,3))));
    for layernum = 1:numel(tags_to_color)
        subplot(1,numel(tags_to_color),layernum);
        %maxcol = max(max(sum(I1.*repmat(reshape(base_colors{layernum},1,1,3),480,640),3)))/sqrt(sum(base_colors{layernum}.^2));
        %newI = base_colors{layernum}'*[0:maxcol/24:maxcol];
        newI = base_colors{layernum}'*1./(1+exp(-([0:maxcol/24:maxcol]-saturation_sigm_threshold)*saturation_sigm_power));
        imshow(reshape(newI',25,1,3));%,'Parent',bar_ax);
    end
end


