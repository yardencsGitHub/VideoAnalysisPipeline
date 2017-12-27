%%
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
edges = [0.5 0.5]; %[0.5 0.5];
opacity_factor = 0.4;


%%
bird_name = 'lrb85315';
bird_folder_name = 'lrb853_15';
template_file = 'lrb85315template';
annotation_file = 'lrb85315auto_annotation5_fix';
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
    results_max = [];
    results_min = [];
    results_std = [];
    results_std_in = [];
    results_r = [];
    results_p = [];
    results_ratio_mean = [];
    results_ratio_max = [];
    results_ratio_std = [];
    results_hmm = [];
    results_hmm_out = [];
    results_hmm_time = [];
    results_hmm_time_out = [];
    results_hist_on = {};
    results_hist_off = {};
    Day = unique_dates(Day_num,:);
    cd([laptop_manualROI_folder '/ROIdata/' Day]);
    FILES = dir('NonoverlapBaseROIdata*.mat');
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
            sig_peak = sig_integrals;
            sig_bottom = sig_integrals;
            sig_hmm = sig_integrals;
            sig_hmm_time = sig_integrals;
            sig_hmm_time_out = sig_integrals;
            sig_hmm_out = sig_integrals;
            active_times_on = cell(size(dff,1),1);
            active_times_off = cell(size(dff,1),1);
            for cnt = 1:size(hits,1)
                fnum = hits(cnt,1);
                phrasenum = hits(cnt,2);
                fname = FILES{fnum};
                tokens = regexp(fname,'_','split');
                loc = find(locs == str2num(tokens{3}));
                phrases = return_phrase_times(elements{loc});
                %if ismember(sylnum,phrases.phraseType)
                load(fname);
                load(['GaussHmmMAP_' fname]);
                display(fname);
%                 if spikes ==2
%                     s = zscore(dff(:,n_del_frames+1:end)')';
%                 end
%                 dff = dff(:,n_del_frames+1:end);
%                 y  = dff - ones(size(dff,1),1)*smooth(mean(dff),100)';
%                 t = vidTimes;  

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
%                     clear sigma;
%                     try
%                         mu = [quantile(sig,0.05) max(sig)];
%                         sigma(1,1,:) = [0.05*ones(1,nstates)];
%                         CPD = condGaussCpdCreate(mu, sigma);
%                         [model, loglikHist] = hmmFit(sig, nstates, 'gauss','pi0',[zeros(1,nstates-1) 1],'emission0',CPD);%,'transPrior',[10 1;10 3],'piPrior',[100 1]);                    
%                         map_path = hmmMap(model, sig)-1;
%                         map_path = abs(median(map_path)-map_path);
%                     catch em
%                        '-'; 
%                     end
                    
                    sig_integrals(cnt,roi_n) = sum(signal((t >= tonset) & ...
                        (t <= toffset)));
                    sig_std_in(cnt,roi_n) = quantile(signal((t >= (tonset-edges(1))) & ...
                        (t <= (toffset+edges(2)))),0.9) - ...
                        quantile(signal((t >= (tonset-edges(1))) & ...
                        (t <= (toffset+edges(2)))),0.1);
                    sig_std_out(cnt,roi_n) = quantile(signal((t <= (tonset-edges(1))) | ...
                        (t >= (toffset+edges(2)))),0.9) - ...
                        quantile(signal((t <= (tonset-edges(1))) | ...
                        (t >= (toffset+edges(2)))),0.1);
                    try
                        sig_peak(cnt,roi_n) = max(signal((t >= tonset) & (t <= toffset)));
                        sig_bottom(cnt,roi_n) = min(signal((t >= tonset) & (t <= toffset)));
                        sig_hmm(cnt,roi_n) = 1*(sum(map_path((t >= tonset) & ...
                        (t <= toffset))) >= min_num_active_bins);
                        sig_hmm_time(cnt,roi_n) = sum(map_path((t >= tonset) & ...
                        (t <= toffset)))/30; 
                        sig_hmm_time_out(cnt,roi_n) = sum(map_path((t < tonset-edges(1)) | ...
                        (t > toffset+edges(2))))/30; 
                        sig_hmm_out(cnt,roi_n) = 1*(sum(map_path((t < tonset-edges(1)) | ...
                        (t > toffset+edges(2))))/sum(((t < tonset-edges(1)) | ...
                        (t > toffset+edges(2)))));% >= min_num_active_bins/30);

                        active_times_on{roi_n} = [active_times_on{roi_n} t(map_path == 1)-tonset];
                        active_times_off{roi_n} = [active_times_off{roi_n} (t(map_path == 1)-tonset)/(toffset-tonset)];
                    catch em
                        sig_peak(cnt,roi_n) = nan;
                        sig_bottom(cnt,roi_n) = nan;
                        sig_hmm(cnt,roi_n) = nan;
                        sig_hmm_out(cnt,roi_n) = nan;
                        sig_hmm_time(cnt,roi_n) = nan;
                        sig_hmm_time_out(cnt,roi_n) = nan;
                    end

                end
                %end

                %end
            end
            [r p] = corr(sig_integrals,durations);
            ratios = sig_std_in./(sig_std_out+1e-5);
            [n,x] = cellfun(@(x)hist(x,-3:1/30:3),active_times_on,'UniformOutput',0);
            n=cellfun(@(x)x/sum(x),n,'UniformOutput',0);
            results_hist_on = {results_hist_on{:} cell2mat(n)};
            
            [n,x] = cellfun(@(x)hist(x,-3:1/30:3),active_times_off,'UniformOutput',0);
            n=cellfun(@(x)x/sum(x),n,'UniformOutput',0);
            results_hist_off = {results_hist_off{:} cell2mat(n)};
         
            results_max = [results_max quantile(sig_peak,0.9,1)'];
            results_min = [results_min quantile(sig_bottom,0.1,1)'];
            results_std = [results_std nanstd(sig_peak - sig_bottom)'];
            results_std_in = [results_std_in nanmean(sig_std_in)'];
            results_r = [results_r r];
            results_p = [results_p p];
            results_ratio_mean = [results_ratio_mean nanmean(ratios,1)'];
            results_ratio_std = [results_ratio_std std(ratios,[],1)']; 
            results_ratio_max = [results_ratio_max quantile(ratios,0.9,1)']; %max(ratios,[],1)'];
            results_hmm = [results_hmm nanmean(sig_hmm,1)'];
            results_hmm_out = [results_hmm_out nanmean(sig_hmm_out,1)'];
            results_hmm_time = [results_hmm_time nanmean(sig_hmm_time,1)'];
            results_hmm_time_out = [results_hmm_time_out nanmean(sig_hmm_time_out,1)'];
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
            results_max = [results_max nan*ones(size(dff,1),1)];
            results_min = [results_min nan*ones(size(dff,1),1)];
            results_std = [results_std nan*ones(size(dff,1),1)];
            results_std_in = [results_std_in nan*ones(size(dff,1),1)];
            results_hmm = [results_hmm nan*ones(size(dff,1),1)];
            results_hmm_out = [results_hmm_out nan*ones(size(dff,1),1)];
            results_hmm_time = [results_hmm_time nan*ones(size(dff,1),1)];
            results_hmm_time_out = [results_hmm_time_out nan*ones(size(dff,1),1)];
            
            results_hist_on = {results_hist_on{:} []};
            results_hist_off = {results_hist_off{:} []};
        end
    end
    results(Day_num).results_r = results_r;
    results(Day_num).results_p = results_p;
    results(Day_num).results_ratio_mean = results_ratio_mean;
    results(Day_num).results_ratio_max = results_ratio_max;
    results(Day_num).results_ratio_std = results_ratio_std;
    results(Day_num).Date = Day;
    results(Day_num).Max = results_max;
    results(Day_num).Min = results_min;
    results(Day_num).Std = results_std;
    results(Day_num).Std_in = results_std_in;
    results(Day_num).hmm = results_hmm;
    results(Day_num).hmm_out = results_hmm_out;
    results(Day_num).hmm_time = results_hmm_time;
    results(Day_num).hmm_time_out = results_hmm_time_out;
    results(Day_num).hist_on = results_hist_on;
    results(Day_num).hist_rel = results_hist_off;
end


