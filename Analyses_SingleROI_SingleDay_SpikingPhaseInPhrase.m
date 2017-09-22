%% go over all ROIs, deconvolve and relate spiking events to phrase edges
% this needs to be done per neuron.. so, 
% Define a cyclic c.o.m within phrases
% This needs to exclude long syllables
clear;
% Also build a syllable specific histogram
    zscoring_type = 1;
    delete_frames = 1;
    n_del_frames = 6;
    hvc_offset = 0.035;
    clear results;
    min_hits = 3;
    g = 0.9;
    warp = 0;
    locktoonset = 1;
    mulcnt = 0.1;
    spikes = 2;
    edges = [0.1 0.1]; %[0.5 0.5];
    opacity_factor = 0.4;
    
    thr = 0.75;
    foopsi_cnt = 0;
%%
bird_name = 'lrb85315';
bird_folder_name = 'lrb853_15';
template_file = 'lrb85315template';
annotation_file = 'lrb85315auto_annotation5';
CNMFEfolder = '/Users/yardenc/Documents/Experiments/Code and Hardware Dev/CNMF_E';
laptop_mov_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs'];
laptop_wav_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav'];
laptop_gif_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/gif'];
laptop_storage_folder = ['/Volumes/CanaryData/DATA/lrb853_15/movs/'];
laptop_annotated_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated'];
laptop_annotated_images_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated/images'];
DamagedFolder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/too_large_or_damaged/'];
laptop_manualROI_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/ManualROIs'];

laptop_manualROI_analyses_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/' bird_folder_name '/ManualROIs/PhraseLockedSpikeTiming'];
%%
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/small-utils'));
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/VideoAnalysisPipeline'));
%%

syllables = [0:9 200:209 300:309 400:409 500];
long_syllables = [3 6 7 201 205 206 207 307 401 402 403 404 405 406 407 409];
midrange_syllables = [2 200 202 204 300 301 302 303 400 500];
short_syllables = [0 1 4 5 8 9 203 208 209 304 305 306 308 309 408];

cd (laptop_manualROI_folder);

n_syllables = numel(syllables);


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
unique_dates = datestr(setdiff(unique(datenum(dates)),[736804]),'yyyy_mm_dd'); %does not include 04/19-21th (remove for other birds)

%%
clear i;
for Day_num = 1: size(unique_dates,1)
    Day = unique_dates(Day_num,:);
    display(Day);
    cd([laptop_manualROI_folder '/ROIdata/' Day]);
    FILES = dir('ROIdata*.mat');
    FILES = {FILES.name};
    load(FILES{1});
    clear per_syl_phases all_syl_phases short_syl_phases midrange_syl_phases long_syl_phases;
    clear per_syl_num all_syl_num short_syl_num midrange_syl_num long_syl_num;
    clear per_syl_time all_syl_time short_syl_time midrange_syl_time long_syl_time;
    clear num_phraseType;
    num_phraseType = zeros(n_syllables,1);
    for roi_n = 1:size(dff,1)
        all_syl_phases(roi_n).data=[];
        all_syl_num(roi_n).data=[];
        all_syl_time(roi_n).data=[];
        
        short_syl_phases(roi_n).data=[];
        short_syl_num(roi_n).data=[];
        short_syl_time(roi_n).data=[];
        
        midrange_syl_phases(roi_n).data=[];
        midrange_syl_num(roi_n).data=[];
        midrange_syl_time(roi_n).data=[];
        
        long_syl_phases(roi_n).data=[];
        long_syl_num(roi_n).data=[];
        long_syl_time(roi_n).data=[];
        
        for sylnum = 1:numel(syllables)  
            per_syl_phases(roi_n,sylnum).data=[];
            per_syl_num(roi_n,sylnum).data=[];
            per_syl_time(roi_n,sylnum).data=[];
        end
    end
    for fnum = 1:numel(FILES)
            fname = FILES{fnum};
            tokens = regexp(fname,'_','split');
            loc = find(locs == str2num(tokens{3}));
            phrases = return_phrase_times(elements{loc});
            for sylnum = 1:numel(syllables)  
                if ismember(syllables(sylnum),phrases.phraseType)
                    num_phraseType(sylnum) = num_phraseType(sylnum) + 1;
                end
            end
            if ~isempty(elements{loc}.segAbsStartTimes)
                load(fname);
                dff_tmp = dff(:,n_del_frames+1:end);
                y  = dff_tmp; % - ones(size(dff_tmp,1),1)*smooth(mean(dff_tmp),100)';
                if delete_frames == 1
                    if zscoring_type == 1
                        y = reshape(zscore(y(:)),size(y));                
                    end
                    t = vidTimes(n_del_frames+1:end)-hvc_offset;
                else
                    if zscoring_type == 1
                        %
                        %y = [(dff(:,1:n_del_frames)-mean(y(:)))/std(y(:)) y];
                        y = [(dff(:,1:n_del_frames)-mean(y(:)))/std(y(:)) reshape(zscore(y(:)),size(y))];          
                        %s = reshape(zscore(y(:)),size(y));
                    else
                        y = dff;
                    end
                    t = vidTimes-hvc_offset;
                end
                
                
                for roi_n = 1:size(dff,1)
                    try %detrend
                        %[c, s, options] = deconvolveCa(detrend(y(roi_n,:)),'ar1',g,'method','constrained-foopsi');
                        %[c, s, options] = deconvolveCa((y(roi_n,:)),'ar1','method','thresholded','optimize_b',1,'optimize_pars',1);
                        [c, s, options] = deconvolveCa(y(roi_n,:),'ar2',[1.3 -0.422],'method','thresholded','optimize_b','optimize_smin');
                    catch em
                        try
                            [c, s, options] = deconvolveCa(y(roi_n,:),'ar2',[1.5 -0.56],'method','foopsi','optimize_b','optimize_smin');
                            %[c, s, options] = deconvolveCa((y(roi_n,:)),'ar1','method','foopsi','optimize_b',1);
                        catch eem
                            [c, s, options] = deconvolveCa((y(roi_n,:)-quantile(y(roi_n,:),0.05)),'ar1',g,'method','foopsi');
                        end
                        foopsi_cnt = foopsi_cnt + 1; 
                        
                        %[c, s, options] = deconvolveCa(detrend(y(roi_n,:)),'ar1',g,'method','foopsi');
                    end
                    
                    x = 1*(s > thr);
                    spikes = find(x > 0);
                    try
                        for cnt = 1:numel(spikes)
                            if (t(spikes(cnt)) >= phrases.phraseFileStartTimes(1)) && (t(spikes(cnt)) <= phrases.phraseFileEndTimes(end))
                                phrase_loc = find(t(spikes(cnt)) >= phrases.phraseFileStartTimes & ...
                                    t(spikes(cnt)) <= phrases.phraseFileEndTimes);
                                % all syllables
                                if ismember(phrases.phraseType(phrase_loc), syllables)
                                    sylnum = find(syllables == phrases.phraseType(phrase_loc));
                                    per_syl_phases(roi_n,sylnum).data=[per_syl_phases(roi_n,sylnum).data; ...
                                        exp(2*pi*i*(t(spikes(cnt))-phrases.phraseFileStartTimes(phrase_loc))/(phrases.phraseFileEndTimes(phrase_loc)-phrases.phraseFileStartTimes(phrase_loc)))];
                                    per_syl_num(roi_n,sylnum).data=[per_syl_num(roi_n,sylnum).data; ...
                                        sum(elements{loc}.segFileStartTimes <= t(spikes(cnt)) & ...
                                        elements{loc}.segFileStartTimes >= phrases.phraseFileStartTimes(phrase_loc))];
                                    per_syl_time(roi_n,sylnum).data=[per_syl_time(roi_n,sylnum).data; ...
                                        t(spikes(cnt))-phrases.phraseFileStartTimes(phrase_loc) ...
                                        -t(spikes(cnt))+phrases.phraseFileEndTimes(phrase_loc)];

                                    all_syl_phases(roi_n).data = [all_syl_phases(roi_n).data; ...
                                        per_syl_phases(roi_n,sylnum).data(end)];
                                    all_syl_num(roi_n).data = [all_syl_num(roi_n).data; ...
                                        per_syl_num(roi_n,sylnum).data(end)];
                                    all_syl_time(roi_n).data = [all_syl_time(roi_n).data; ...
                                        per_syl_time(roi_n,sylnum).data(end,:)];

                                end

                                % long syllables
                                if ismember(phrases.phraseType(phrase_loc), long_syllables)

                                    long_syl_phases(roi_n).data=[long_syl_phases(roi_n).data; ...
                                        exp(2*pi*i*(t(spikes(cnt))-phrases.phraseFileStartTimes(phrase_loc))/(phrases.phraseFileEndTimes(phrase_loc)-phrases.phraseFileStartTimes(phrase_loc)))];
                                    long_syl_num(roi_n).data=[long_syl_num(roi_n).data; ...
                                        sum(elements{loc}.segFileStartTimes <= t(spikes(cnt)) & ...
                                        elements{loc}.segFileStartTimes >= phrases.phraseFileStartTimes(phrase_loc))];
                                    long_syl_time(roi_n).data=[long_syl_time(roi_n).data; ...
                                        t(spikes(cnt))-phrases.phraseFileStartTimes(phrase_loc) ...
                                        -t(spikes(cnt))+phrases.phraseFileEndTimes(phrase_loc)];
                                end

                                % mid range syllables
                                if ismember(phrases.phraseType(phrase_loc), midrange_syllables)

                                    midrange_syl_phases(roi_n).data=[midrange_syl_phases(roi_n).data; ...
                                        exp(2*pi*i*(t(spikes(cnt))-phrases.phraseFileStartTimes(phrase_loc))/(phrases.phraseFileEndTimes(phrase_loc)-phrases.phraseFileStartTimes(phrase_loc)))];
                                    midrange_syl_num(roi_n).data=[midrange_syl_num(roi_n).data; ...
                                        sum(elements{loc}.segFileStartTimes <= t(spikes(cnt)) & ...
                                        elements{loc}.segFileStartTimes >= phrases.phraseFileStartTimes(phrase_loc))];
                                    midrange_syl_time(roi_n).data=[midrange_syl_time(roi_n).data; ...
                                        t(spikes(cnt))-phrases.phraseFileStartTimes(phrase_loc) ...
                                        -t(spikes(cnt))+phrases.phraseFileEndTimes(phrase_loc)];
                                end

                                % short syllables
                                if ismember(phrases.phraseType(phrase_loc), short_syllables)

                                    short_syl_phases(roi_n).data=[short_syl_phases(roi_n).data; ...
                                        exp(2*pi*i*(t(spikes(cnt))-phrases.phraseFileStartTimes(phrase_loc))/(phrases.phraseFileEndTimes(phrase_loc)-phrases.phraseFileStartTimes(phrase_loc)))];
                                    short_syl_num(roi_n).data=[short_syl_num(roi_n).data; ...
                                        sum(elements{loc}.segFileStartTimes <= t(spikes(cnt)) & ...
                                        elements{loc}.segFileStartTimes >= phrases.phraseFileStartTimes(phrase_loc))];
                                    short_syl_time(roi_n).data=[short_syl_time(roi_n).data; ...
                                        t(spikes(cnt))-phrases.phraseFileStartTimes(phrase_loc) ...
                                        -t(spikes(cnt))+phrases.phraseFileEndTimes(phrase_loc)];
                                end

                            end
                        end
                    catch ee
                        '(';
                    end
                end
            end
    end
    fname = ['AlignedSpikes_' bird_name '_' Day '.mat'];
    save(fullfile(laptop_manualROI_analyses_folder,fname),'per_syl_phases','all_syl_phases','short_syl_phases', ...
        'midrange_syl_phases','long_syl_phases','per_syl_num','all_syl_num','short_syl_num','midrange_syl_num', ...
        'long_syl_num','per_syl_time','all_syl_time','short_syl_time','midrange_syl_time','long_syl_time','num_phraseType');
end
                        