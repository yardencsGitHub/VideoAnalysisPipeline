function [results,syllables] = per_phrase_activity_mass(birdnum,ignore_dates,ignore_entries,join_entries,include_zero,varargin)
% Prerequisit: Run HMM_matlab to extract the signal vs noise segments
bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix' 'NonoverlapBaseROIdata_'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation5_alexa' 'baseROIdata_'};
bird3_params = {'lbr3009' 'lbr3009' 'lbr3009_template_4TF' 'lbr3009auto_annotation1_fix' 'baseROIdata_'};
switch birdnum
    case 1
        bird_params = bird1_params;
    case 2
        bird_params = bird2_params;
    case 3
        bird_params = bird3_params;
end

if isempty(ignore_dates)
    ignore_dates = {'2017_04_19'};
end
if isempty(ignore_entries)
   ignore_entries = [-1 100 102 101 103 202 406 408 409 402 403];
end
if isempty(join_entries)
    join_entries = {[207 307 407] [404 405] [208 209] [200 309]};
end

clear results;
zscoring_type = 0;
delete_frames = 1;
n_del_frames = 6;
hvc_offset = 0.04;
min_num_active_bins = 1;
nstates = 2;  

min_hits = 10;
g = 0.85;
warp = 0;
locktoonset = 1;
mulcnt = 0.1;
spikes = 2;
edges = [0.1 0.1]; %[0.5 0.5];
if isempty(include_zero)
    include_zero = 1;
end
opacity_factor = 0.4;
max_phrase_gap = 0.5;

%% allow controlling parameters as function pair inputs
nparams=length(varargin);
for i=1:2:nparams
	switch lower(varargin{i})
		case 'delete_frames'
			delete_frames=varargin{i+1};
        case 'bird_number'
			 switch varargin{i+1}
                 case 1
                     bird_params = bird1_params;
                 case 2
                     bird_params = bird2_params;
                 case 3
                     bird_params = bird3_params;
             end
        case 'n_del_frames'
			n_del_frames=varargin{i+1};
        case 'hvc_offset'
			hvc_offset=varargin{i+1}; 
        case 'edges'
			edges=varargin{i+1}; 
        case 'zscoring_type'
			zscoring_type=varargin{i+1}; 
        case 'max_phrase_gap'
			max_phrase_gap=varargin{i+1}; 
        case 'display_opt'
			display_opt=varargin{i+1}; 
        case 'use_residuals'
			use_residuals=varargin{i+1};
        case 'min_num_active_bins'
			min_num_active_bins=varargin{i+1};
        case 'spikes'
			spikes=varargin{i+1};        
    end
end
bird_name = bird_params{1}; % 'lrb85315'; %'lbr3022'; %
bird_folder_name = bird_params{2}; %'lrb853_15'; %'lbr3022'; %
template_file = bird_params{3}; %'lrb85315template'; %'lbr3022_template';%
annotation_file = bird_params{4};
file_prefix = bird_params{5}; 
%%
%bird_name = 'lrb85315';
%bird_folder_name = 'lrb853_15';
%template_file = 'lrb85315template';
%annotation_file = 'lrb85315auto_annotation5_fix';
CNMFEfolder = '/Users/yardenc/Documents/Experiments/Code and Hardware Dev/CNMF_E';
laptop_mov_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs'];
laptop_wav_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav'];
laptop_gif_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/gif'];
laptop_storage_folder = ['/Volumes/CanaryData/DATA/lrb853_15/movs/'];
laptop_annotated_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated'];
laptop_annotated_images_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated/images'];
DamagedFolder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/too_large_or_damaged/'];
laptop_manualROI_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/ManualROIs'];
%%
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/BirdSongBout/BirdSongBout'),'-end');

%%

%syllables = [0:9 200:209 300:309 400:409 500];

cd (laptop_manualROI_folder);
syllables = []; %[0:9 200 201 203:208 300:306 308 400 401 404 500];
load(annotation_file);
for fnum = 1:numel(keys)  
    syllables = unique([syllables unique(elements{fnum}.segType)']);
end
syllables = setdiff(syllables,ignore_entries);
if (include_zero == 0)
    syllables = setdiff(syllables,0);
end
for i = 1:numel(join_entries)
    syllables = setdiff(syllables,join_entries{i}(2:end));
end
n_syllables = numel(syllables);
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
unique_dates = datestr(setdiff(unique(datenum(dates)),datenum(ignore_dates)),'yyyy_mm_dd'); %does not include 04/19th (remove for other birds)

%%
map_states = [];
dff = [];
results = [];
for Day_num = 1: size(unique_dates,1)
    
    Day = unique_dates(Day_num,:);
    cd([laptop_manualROI_folder '/ROIdata/' Day]);
    %FILES = dir('NonoverlapBaseROIdata*.mat');
    FILES = dir([file_prefix bird_name '*.mat']);
    FILES = {FILES.name};
    load(FILES{1});
    daily_res50 = cell(size(dff,1),1);
    daily_res75 = cell(size(dff,1),1);
    daily_res90 = cell(size(dff,1),1);
    daily_res95 = cell(size(dff,1),1);
    for fnum = 1:numel(FILES)
        fname = FILES{fnum};
        tokens = regexp(fname,'_','split');
        loc = find(locs == str2num(tokens{3}));
        phrases = return_phrase_times(trim_element(elements{loc})); %elements{loc}
        phrases = deal_with_time_gaps(phrases,max_phrase_gap);
        load(fname);
        load(['GaussHmmMAP_' fname]);
        dff_tmp = dff(:,n_del_frames+1:end);
        y  = dff_tmp;% - ones(size(dff_tmp,1),1)*smooth(mean(dff_tmp),100)';
        if delete_frames == 1
            if zscoring_type == 1
                y = reshape(zscore(y(:)),size(y));                
            end
            t = vidTimes(n_del_frames+1:end)+hvc_offset;
        else
            if zscoring_type == 1
                %
                %y = [(dff(:,1:n_del_frames)-mean(y(:)))/std(y(:)) y];
                y = [(dff(:,1:n_del_frames)-mean(y(:)))/std(y(:)) reshape(zscore(y(:)),size(y))];          
                %s = reshape(zscore(y(:)),size(y));
            else
                y = dff;
            end
            t = vidTimes+hvc_offset;
        end

        if spikes ==2
            s = y; 
        end  
        
        for roi_n = 1:size(dff,1)

            if spikes < 2
                [c, s, options] = deconvolveCa(detrend(y(roi_n,:)),'ar1',g,'method','constrained-foopsi');
            end
            switch spikes
                case 1
                    signal = s;
                case 0
                  signal = c; %smooth(s(ROIs(roi_n),:),3);  
                otherwise
                    signal = smooth(s(roi_n,:),3);
            end %(n_del_frames+1:end)
            map_path = map_states(roi_n,1:numel(signal));
            map_diff = diff([0 map_path 0]);
            map_on = find(map_diff == 1);
            map_off = find(map_diff == -1)-1;
            mn = mean(signal);
            for mapcnt = 1:numel(map_on)
                if mean(signal(map_on(mapcnt):map_off(mapcnt))) < mn
                    map_path(map_on(mapcnt):map_off(mapcnt)) = 0;
                end
            end
            signal = signal.*map_path';       
            phrase_count = sum(ismember(phrases.phraseType,syllables)); 
            sig_integrals = zeros(1,phrase_count);
            phrase_type_seq = sig_integrals;
            cnt = 1;
            for phrasenum = 1:numel(phrases.phraseType)
                if ~ismember(phrases.phraseType(phrasenum),syllables)
                    continue;
                end
                tonset = phrases.phraseFileStartTimes(phrasenum);
                toffset = phrases.phraseFileEndTimes(phrasenum);
                sig_integrals(cnt) = sum(signal((t >= tonset) & ...
                        (t <= toffset)));
                phrase_type_seq(cnt) = phrases.phraseType(phrasenum);
                cnt = cnt + 1; 
            end
            [I,idx] = sort(sig_integrals,'descend');
            locs50 = min(find(cumsum(I)/sum(I) >= 0.5)); sylseq50 = unique(phrase_type_seq(idx(1:locs50)));
            locs75 = min(find(cumsum(I)/sum(I) >= 0.75)); sylseq75 = unique(phrase_type_seq(idx(1:locs75)));
            locs90 = min(find(cumsum(I)/sum(I) >= 0.9)); sylseq90 = unique(phrase_type_seq(idx(1:locs90)));
            locs95 = min(find(cumsum(I)/sum(I) >= 0.95)); sylseq95 = unique(phrase_type_seq(idx(1:locs95)));
            tmp = nan*zeros(1,numel(syllables)); 
            for tmpcnt = 1:numel(phrase_type_seq) 
                tmp(syllables == phrase_type_seq(tmpcnt)) = 1*ismember(phrase_type_seq(tmpcnt),sylseq50);
            end
            daily_res50{roi_n} = [daily_res50{roi_n}; tmp]; 
            
            tmp = nan*zeros(1,numel(syllables));
            for tmpcnt = 1:numel(phrase_type_seq) 
                tmp(syllables == phrase_type_seq(tmpcnt)) = 1*ismember(phrase_type_seq(tmpcnt),sylseq75);
            end
            daily_res75{roi_n} = [daily_res75{roi_n}; tmp]; 
            
            tmp = nan*zeros(1,numel(syllables));
            for tmpcnt = 1:numel(phrase_type_seq) 
                tmp(syllables == phrase_type_seq(tmpcnt)) = 1*ismember(phrase_type_seq(tmpcnt),sylseq90);
            end
            daily_res90{roi_n} = [daily_res90{roi_n}; tmp]; 
            
            tmp = nan*zeros(1,numel(syllables));
            for tmpcnt = 1:numel(phrase_type_seq) 
                tmp(syllables == phrase_type_seq(tmpcnt)) = 1*ismember(phrase_type_seq(tmpcnt),sylseq95);
            end
            daily_res95{roi_n} = [daily_res95{roi_n}; tmp]; 
            
        end
        
    end
    
    results(Day_num).res50 = zeros(size(dff,1),numel(syllables));
    results(Day_num).res75 = zeros(size(dff,1),numel(syllables));
    results(Day_num).res90 = zeros(size(dff,1),numel(syllables));
    results(Day_num).res95 = zeros(size(dff,1),numel(syllables));
    
    for roi_n = 1:size(dff,1)
        if size(daily_res50{roi_n},1) > 1
            results(Day_num).res50(roi_n,:) = nanmean(daily_res50{roi_n});
            results(Day_num).res75(roi_n,:) = nanmean(daily_res75{roi_n});
            results(Day_num).res90(roi_n,:) = nanmean(daily_res90{roi_n});
            results(Day_num).res95(roi_n,:) = nanmean(daily_res95{roi_n});
        else
            results(Day_num).res50(roi_n,:) = daily_res50{roi_n};
            results(Day_num).res75(roi_n,:) = daily_res75{roi_n};
            results(Day_num).res90(roi_n,:) = daily_res90{roi_n};
            results(Day_num).res95(roi_n,:) = daily_res95{roi_n};
            
        end
    end
end

function element = trim_element(old_element)
    try
        element = old_element;
        trim_locs = find(ismember(element.segType,ignore_entries));
        element.segAbsStartTimes(trim_locs) = [];
        element.segFileStartTimes(trim_locs) = [];
        element.segFileEndTimes(trim_locs) = [];
        element.segType(trim_locs) = [];  
        for i = 1:numel(join_entries)
            join_locs = find(ismember(element.segType,join_entries{i}));
            element.segType(join_locs) = join_entries{i}(1);
        end
    catch em
       'd'; 
       display(fname);
    end
end
end
    
    
    
    
    
    
    