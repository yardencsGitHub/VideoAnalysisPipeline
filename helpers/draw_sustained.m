function draw_sustained(sustained,daynum)
%%
ignore_dates = [];
ignore_entries = [];
join_entries = [];
include_zero = [];
templates = [];
birdnum = 1;

bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation5'};
switch birdnum
    case 1
        bird_params = bird1_params;
    case 2
        bird_params = bird2_params;
    case 3
end
bird_name = bird_params{1}; % 'lrb85315'; %'lbr3022'; %
bird_folder_name = bird_params{2}; %'lrb853_15'; %'lbr3022'; %
template_file = bird_params{3}; %'lrb85315template'; %'lbr3022_template';%
annotation_file = bird_params{4};

if isempty(ignore_dates)
    ignore_dates = {'2017_04_19'};
end
if isempty(ignore_entries)
   ignore_entries = [-1 100 102 101 103 202 406 408 409 402 403];
end
if isempty(join_entries)
    join_entries = {[207 307 407] [404 405] [208 209] [200 309]};
end

zscoring_type = 0;
delete_frames = 1;
n_del_frames = 6;
hvc_offset = 0.04; 

locktoonset = 1;
mulcnt = 0.1;
spikes = 2;

if isempty(include_zero)
    include_zero = 1;
end
opacity_factor = 0.8;
max_phrase_gap = 0.5;

laptop_mov_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs'];
laptop_wav_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav'];
laptop_gif_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/gif'];
laptop_storage_folder = ['/Volumes/CanaryData/DATA/lrb853_15/movs/'];
laptop_annotated_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated'];
laptop_annotated_images_dir = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/movs/wav/annotated/images'];
DamagedFolder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/too_large_or_damaged/'];
laptop_manualROI_folder = ['/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/' bird_folder_name '/ManualROIs'];

cd (laptop_manualROI_folder);
syllables = []; %[0:9 200 201 203:208 300:306 308 400 401 404 500];
load(annotation_file);
load(template_file);
%syllables = [[templates.wavs.segType] -1 102 103];
syllables = [templates.wavs.segType];
% % % for fnum = 1:numel(keys)  
% % %     syllables = unique([syllables unique(elements{fnum}.segType)']);
% % % end
% % % syllables = setdiff(syllables,ignore_entries);
% % % if (include_zero == 0)
% % %     syllables = setdiff(syllables,0);
% % % end
% % % for i = 1:numel(join_entries)
% % %     syllables = setdiff(syllables,join_entries{i}(2:end));
% % % end
syllables = [[-1000 1000] syllables setdiff([-1 102 103],syllables)];
n_syllables = numel(syllables)-2;
colors = distinguishable_colors(n_syllables,'w');
colors = [0.25 0.25 0.25; 0.25 0.25 0.25; colors];
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

Day = sustained(daynum).Day;
sus_num = cellfun(@(x)size(x,1),sustained(daynum).hits);
sus_loc = find(sus_num > 0);

for roi_cnt = 1:numel(sus_loc)
    roi_n = sus_loc(roi_cnt);  
    cd([laptop_manualROI_folder '/ROIdata/' Day]);
    FILES = dir('NonoverlapBaseROIdata*.mat');
    FILES = {FILES.name};
    FILES_n = [];
    for fnum = 1:numel(FILES)
        fname = FILES{fnum};
        tokens = regexp(fname,'_','split');
        FILES_n = [FILES_n; str2num(tokens{3})];
    end
    filenums = sustained(daynum).hits{roi_n}(:,1);
    figure; ax=axes; all_phrase_types = [];
    all_betas = [];
    for fnum = 1:numel(filenums)
        fname = FILES{FILES_n == filenums(fnum)};
        tokens = regexp(fname,'_','split');
        loc = find(locs == str2num(tokens{3}));
        phrases = return_phrase_times(trim_element(elements{loc})); %elements{loc}
        phrases = deal_with_time_gaps(phrases,max_phrase_gap);
        load(fname);
        dff_tmp = dff(:,n_del_frames+1:end);
        if delete_frames == 1
            if zscoring_type == 1
                y = reshape(zscore(y(:)),size(y));                
            else
                y = dff(:,n_del_frames+1:end);
            end
            t = vidTimes(n_del_frames+1:end)+hvc_offset;
        else
            if zscoring_type == 1
                y = [(dff(:,1:n_del_frames)-mean(dff_tmp(:)))/std(dff_tmp(:)) reshape(zscore(dff_tmp(:)),size(dff_tmp))];          
            else
                y = dff;
            end
            t = vidTimes+hvc_offset;
        end
        switch spikes
            case 0
                try 
                    [c, s, options] = deconvolveCa(y(ROI,:),'ar2',[1.3 -0.422],'method','thresholded','optimize_b','optimize_smin');%,'optimize_pars');
                catch em
                    [c, s, options] = deconvolveCa((y(ROI,:)),'ar2','method','foopsi','optimize_b',1);
                end
                signal = c;
            case 1
                try 
                    [c, s, options] = deconvolveCa(y(ROI,:),'ar2',[1.3 -0.422],'method','thresholded','optimize_b','optimize_smin');
                catch em
                    [c, s, options] = deconvolveCa((y(ROI,:)),'ar2','method','foopsi','optimize_b',1);
                end
                signal = s;
            case 2
                signal = smooth(y(roi_n,:),3);
            case 3
                sig = y(roi_n,:);
                clear sigma;
                nstates = 2;
                dq = (0.5-0.05)/(nstates-2+1e-5);          
                mu = [quantile(sig,0.05+[0:nstates-2]*dq) max(sig) ];
                sigma(1,1,:) = [0.05*ones(1,nstates)];
                CPD = condGaussCpdCreate(mu, sigma);
                [model, loglikHist] = hmmFit(sig, nstates, 'gauss','pi0',[zeros(1,nstates-1) 1],'emission0',CPD);                    
                path = hmmMap(model, sig)-1;
                signal = abs(median(path)-path)';
        end
        timetag = t - sustained(daynum).hits{roi_n}(fnum,2);
        phrase_nums = sustained(daynum).sus{roi_n}{fnum};
        phrase_types = phrases.phraseType(phrase_nums);
        all_phrase_types = unique(union(all_phrase_types,phrase_types));
        
        for currphrase = phrase_nums
            try
            plot3(ax,timetag(t >= phrases.phraseFileStartTimes(currphrase) & ...
                 t <= phrases.phraseFileEndTimes(currphrase)), ...
                 fnum*mulcnt*ones(1,sum(t >= phrases.phraseFileStartTimes(currphrase) & ...
                 t <= phrases.phraseFileEndTimes(currphrase))), ...
                 signal(t >= phrases.phraseFileStartTimes(currphrase) & ...
                 t <= phrases.phraseFileEndTimes(currphrase))+0, ...
                 'LineWidth',2,'Color',[colors(find(syllables == phrases.phraseType(currphrase)),:) opacity_factor]);
            hold on;
            catch em
               'd' 
            end
        end
        subsignal = signal(t >= phrases.phraseFileStartTimes(phrase_nums(1)) & t <= phrases.phraseFileEndTimes(phrase_nums(end)));
        sub_t = timetag(t >= phrases.phraseFileStartTimes(phrase_nums(1)) & t <= phrases.phraseFileEndTimes(phrase_nums(end)));
        increases=arrayfun(@(x)max(subsignal(x:x+6)),7:numel(subsignal)-6)-arrayfun(@(x)min(subsignal(x-6:x)),7:numel(subsignal)-6);
        increases = 1*(increases > 0.1);
        sub_decreases = zeros(size(subsignal)); sub_decreases(7:numel(sub_t)-6) = 1-increases;
        map_diff = diff([0; sub_decreases; 0]);
        map_on = find(map_diff == 1);
        map_off = find(map_diff == -1)-1;
        decrease_durations = (map_off-map_on)/30;
        decrease_mags=[]; 
        for decnum = 1:numel(decrease_durations)
            decrease_mags = [decrease_mags; ...
                max(subsignal(map_on(decnum):map_off(decnum))) - min(subsignal(map_on(decnum):map_off(decnum)))];
        end
        decrease_locs = find(decrease_durations >= 0.2 & decrease_mags >= 0.1); % only work on 200mSec or larger segments
        dec_fits=zeros(size(subsignal));
        for dec_num = 1:numel(decrease_locs)
            dec_signal = subsignal(map_on(decrease_locs(dec_num)):map_off(decrease_locs(dec_num)));
            min_dec_sig = min(dec_signal);
            dec_signal = dec_signal - min_dec_sig;
            dec_t = [1:numel(dec_signal)]'/30;
            beta = fminunc(@(x)sum((dec_signal-x(1)*(dec_t<x(2))-x(1)*exp(-(dec_t-x(2))/x(3)).*(dec_t>=x(2))).^2),[0.3 0.2 1]);
            dec_fits(map_on(decrease_locs(dec_num)):map_off(decrease_locs(dec_num))) = ...
                beta(1)*(dec_t<beta(2))+beta(1)*exp(-(dec_t-beta(2))/beta(3)).*(dec_t>=beta(2))+min_dec_sig;
            all_betas = [all_betas; beta(3)];
        end
        plot3(ax,sub_t, ...
                 fnum*mulcnt*ones(1,numel(sub_t)), ...
                 dec_fits, ...
                 'LineWidth',1,'Color',[0.8 0.7 0.7]);
        
    end
    title([sustained(daynum).Day ' roi #' num2str(roi_n) ' syls: ' num2str(all_phrase_types')],'Interpreter','none');
    figure; hist(all_betas); title([sustained(daynum).Day ' roi #' num2str(roi_n) ],'Interpreter','none');
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
