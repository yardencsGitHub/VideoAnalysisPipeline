function [sustained,locs] = locate_sustained_activity(birdnum,ignore_dates,ignore_entries,join_entries,include_zero,edges)
on_time_thr = 2;
signal_thr = 0.05;

bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix' 'NonoverlapBaseROIdata_' 'nonoverlap_newROI_'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation5_alexa' 'baseROIdata_' 'ROI_'};
bird3_params = {'lbr3009' 'lbr3009' 'lbr3009_template_4TF' 'lbr3009auto_annotation1_fix' 'baseROIdata_' 'ROI_'};

switch birdnum
        case 1
            %load('/Users/yardenc/Documents/Projects/CohenGardner2017_CanaryHVCImaging/Figures/lrb853_colors.mat');
            bird_params = bird1_params;
            ignore_dates = {'2017_04_19'};
            ignore_entries = [-1 100 102 101 103 202 406 402 403 408 409];
            join_entries = {[207 307 407] [404 405] [208 209] [200 309]};
            if isempty(include_zero)
                include_zero = 1;
            end
        case 2
            %load('/Users/yardenc/Documents/Projects/CohenGardner2017_CanaryHVCImaging/Figures/lbr3022_colors.mat');
            bird_params = bird2_params;
            %'2017_04_05' '2017_04_06'  '2017_04_11' '2017_04_12' '2017_04_16'
            %'2017_04_20' '2017_04_21' '2017_04_23' '2017_04_25' '2017_04_26' '2017_04_30' '2017_05_03'
            ignore_dates = {'2017_04_14' '2017_04_27' };
            ignore_entries = [-1 100 102 101 103];
            join_entries = {[2 206]}; %{[207 307 407] [404 405] [208 209] [200 309]};
            if isempty(include_zero)
                include_zero = 1;
            end
        case 3
            bird_params = bird3_params;
            ignore_dates = {'2017_04_27' '2017_04_28' '2017_06_05' '2017_06_06' '2017_06_07' '2017_06_08' '2017_06_09' '2017_06_12' '2017_06_13' ... 
                '2017_06_14' '2017_06_15' '2017_06_16' '2017_06_19' '2017_06_20' '2017_06_21' '2017_06_22' '2017_06_27' ...
                '2017_06_28' '2017_06_29' '2017_06_30' '2017_07_03' '2017_07_04' '2017_07_06' '2017_07_07' '2017_07_10' ...
                '2017_07_11' '2017_07_12' '2017_07_13' '2017_07_14' '2017_07_18' '2017_07_19' '2017_07_20' '2017_07_21'};
            ignore_entries = [-1 100 102 101 103];
            join_entries = {};
            if isempty(include_zero)
                include_zero = 0;
            end
end
    
bird_name = bird_params{1}; % 'lrb85315'; %'lbr3022'; %
bird_folder_name = bird_params{2}; %'lrb853_15'; %'lbr3022'; %
template_file = bird_params{3}; %'lrb85315template'; %'lbr3022_template';%
annotation_file = bird_params{4};
file_prefix = bird_params{5};

if isempty(ignore_dates)
    ignore_dates = {'2017_04_19'};
end
% if isempty(ignore_entries)
%    ignore_entries = [-1 100 102 101 103 202 406 408 409 402 403];
% end
% if isempty(join_entries)
%     join_entries = {[207 307 407] [404 405] [208 209] [200 309]};
% end

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
if isempty(edges)
    edges = [0.1 0.1]; %[0.5 0.5];
end

opacity_factor = 0.4;
max_phrase_gap = 0.5;


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
for Day_num = 1: size(unique_dates,1)
    Day = unique_dates(Day_num,:);
    cd([laptop_manualROI_folder '/ROIdata/' Day]);
    FILES = dir([file_prefix '*.mat']);
    FILES = {FILES.name};
    load(FILES{1});
    sus = repmat({{}},size(dff,1),1);
    hits = cell(size(dff,1),1);
    %clear sus hits;
    for fnum = 1:numel(FILES)
        fname = FILES{fnum};
        tokens = regexp(fname,'_','split');
        loc = find(locs == str2num(tokens{3}));
        phrases = return_phrase_times(trim_element(elements{loc})); %elements{loc}
        phrases = deal_with_time_gaps(phrases,max_phrase_gap);
        load(fname);
        load(['GaussHmmMAP_' fname]);
        display(fname);
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
                % applay threshold on signal
                map_path = map_path.*(signal' >= signal_thr);
                map_diff = diff([0 map_path 0]);
                map_on = find(map_diff == 1);
                map_off = find(map_diff == -1)-1;
                map_durations = 1/30*[map_off-map_on];
                long_signals = find(map_durations > on_time_thr);
                for sus_cnt = 1:numel(long_signals)
                    phr = find(phrases.phraseFileStartTimes < 1/30*map_off(long_signals(sus_cnt)) & ...
                        phrases.phraseFileEndTimes > 1/30*map_on(long_signals(sus_cnt)));
                    sus{roi_n} = {sus{roi_n}{:} phr};
                    hits{roi_n} = [hits{roi_n}; [str2num(tokens{3}) 1/30*map_on(long_signals(sus_cnt)) 1/30*map_off(long_signals(sus_cnt))]];  
                    
                end
        end
    end
    sustained(Day_num).sus = sus;
    sustained(Day_num).hits = hits;
    sustained(Day_num).Day = Day;
end
function element = trim_element(old_element)

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
end
end