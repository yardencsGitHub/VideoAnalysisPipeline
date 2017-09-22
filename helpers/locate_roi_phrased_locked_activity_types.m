%%
clear results;
min_hits = 3;
g = 0.85;
warp = 0;
locktoonset = 1;
mulcnt = 0.1;
spikes = 2;
edges = [0.1 0.1]; %[0.5 0.5];
opacity_factor = 0.4;
n_del_frames = 5;

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
%%


syllables = [0:9 200:209 300:309 400:409 500];

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
unique_dates = datestr(setdiff(unique(datenum(dates)),736804),'yyyy_mm_dd'); %does not include 04/19th (remove for other birds)

%%

for Day_num = 1: size(unique_dates,1)
    results_r = [];
    results_p = [];
    results_ratio_mean = [];
    results_ratio_max = [];
    results_ratio_std = [];
    Day = unique_dates(Day_num,:);
    cd([laptop_manualROI_folder '/ROIdata/' Day]);
    FILES = dir('ROIdata*.mat');
    FILES = {FILES.name};
    for syl_cnt = 1:numel(syllables)
        sylnum = syllables(syl_cnt);
        hits = [];
        durations = [];
        for fnum = 1:numel(FILES)
            fname = FILES{fnum};
            tokens = regexp(fname,'_','split');
            loc = find(locs == str2num(tokens{3}));
            phrases = return_phrase_times(elements{loc});
            if ismember(sylnum,phrases.phraseType)
                %load(fname);
                %s = [s; zscore(dff(ROIs,n_del_frames+1:end)')];

                phrase_locs = find(phrases.phraseType == sylnum);
                for phrase_loc = 1:numel(phrase_locs)
                    phrasenum = phrase_locs(phrase_loc);

                    tonset = phrases.phraseFileStartTimes(phrasenum);
                    toffset = phrases.phraseFileEndTimes(phrasenum);
                    hits = [hits; fnum phrasenum];
                    durations = [durations; toffset-tonset];
                end
            end
        end
        [durations,dur_idx] = sort(durations);
        hits = hits(dur_idx,:);
        
        %%%%
        if size(hits,1) >= min_hits
            load(fname);
            sig_integrals = nan*zeros(size(hits,1),size(dff,1));
            sig_std_in = sig_integrals;
            sig_std_out = sig_integrals;
            for cnt = 1:size(hits,1)
                fnum = hits(cnt,1);
                phrasenum = hits(cnt,2);
                fname = FILES{fnum};
                tokens = regexp(fname,'_','split');
                loc = find(locs == str2num(tokens{3}));
                phrases = return_phrase_times(elements{loc});
                %if ismember(sylnum,phrases.phraseType)
                load(fname);
                display(fname);
                if spikes ==2
                    s = zscore(dff(:,n_del_frames+1:end)')';
                end
                dff = dff(:,n_del_frames+1:end);
                y  = dff - ones(size(dff,1),1)*smooth(mean(dff),100)';
                t = vidTimes;   

                tonset = phrases.phraseFileStartTimes(phrasenum);
                toffset = phrases.phraseFileEndTimes(phrasenum);
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
                    end
                    sig_integrals(cnt,roi_n) = sum(signal((t(n_del_frames+1:end) >= tonset) & ...
                        (t(n_del_frames+1:end) <= toffset)));
                    sig_std_in(cnt,roi_n) = quantile(signal((t(n_del_frames+1:end) >= (tonset-edges(1))) & ...
                        (t(n_del_frames+1:end) <= (toffset+edges(2)))),0.9) - ...
                        quantile(signal((t(n_del_frames+1:end) >= (tonset-edges(1))) & ...
                        (t(n_del_frames+1:end) <= (toffset+edges(2)))),0.1);
                    sig_std_out(cnt,roi_n) = quantile(signal((t(n_del_frames+1:end) <= (tonset-edges(1))) | ...
                        (t(n_del_frames+1:end) >= (toffset+edges(2)))),0.9) - ...
                        quantile(signal((t(n_del_frames+1:end) <= (tonset-edges(1))) | ...
                        (t(n_del_frames+1:end) >= (toffset+edges(2)))),0.1);

                end
                %end

                %end
            end
            [r p] = corr(sig_integrals,durations);
            ratios = sig_std_in./(sig_std_out+1e-5);
            results_r = [results_r r];
            results_p = [results_p p];
            results_ratio_mean = [results_ratio_mean mean(ratios,1)'];
            results_ratio_std = [results_ratio_std std(ratios,[],1)']; 
            results_ratio_max = [results_ratio_max quantile(ratios,0.9,1)']; %max(ratios,[],1)'];
            if (Day_num == 65 && syl_cnt == 26)
                '9';
            end
        else
            load(fname);
            results_r = [results_r nan*ones(size(dff,1),1)];
            results_p = [results_p nan*ones(size(dff,1),1)];
            results_ratio_mean = [results_ratio_mean nan*ones(size(dff,1),1)];
            results_ratio_std = [results_ratio_std nan*ones(size(dff,1),1)]; 
            results_ratio_max = [results_ratio_max nan*ones(size(dff,1),1)];
        end
    end
    results(Day_num).results_r = results_r;
    results(Day_num).results_p = results_p;
    results(Day_num).results_ratio_mean = results_ratio_mean;
    results(Day_num).results_ratio_max = results_ratio_max;
    results(Day_num).results_ratio_std = results_ratio_std;
    results(Day_num).Date = Day;
end


