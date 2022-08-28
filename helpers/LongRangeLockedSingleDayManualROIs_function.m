function [hndls,outvars,outstats] = LongRangeLockedSingleDayManualROIs_function(ax,bird_params,Day,ignore_entries,join_entries,sylidx,syllabels_sequence,ROI,locktoonset,order_flag,varargin)
%% 
% This script creates single ROI single day alignments to complex sequences
% with only one variable (duration or syllable type)

% Inputs:
%   ax - axes where to plot
%   bird_params - a struct with the following fields:
%       'bird_name' - a string with the bird's name as it appears in file names
%       'base_dir' - string, path to the folder containing bird-names subfolders, each containing date-named data folders
%       'bird_folder_name' - a string with the bird's name as it appears in the folder name
%       'template_filepath' - string, full path to the template file
%       'annotation_filepath'- string, full path to the annotation file
%       'file_prefix' - string, the initial string in all data files. This
%       must be followed, in all data files, by the string 'birdname_fileIDX_YYYY_mm_dd_hh_mm_ss.mat'
%       'brainard_features_filepath' - string, full path to the calculated Brainard acoustic features for the song files
%   Day - text rep. of date
%   ignore_entries - A vector of label numbers to ignore completely. 
%   join_entries - A cell of vectors, each containing a >1 number of labels
%   sylidx - index of syllable in the sequence to align to
%   syllabels_sequence - the sequence of phrase identities to lock. Insert
%       nan to keep it free
%   ROI - the roi # 
%   locktoonset - 0/1 to lock to phrase onset or offset
%   order_flag - the correlate. Use positive integers to indicate durations of phrases. A vector of positive integers will result in summation of
%   durations. Use negative integers for type (no vectors)

% Outputs:
%   hndls: graphics handels to the created figures.
%   outvars: struct, output variables depending on what analysis was done.
%   outstats: structs, results of statistical tests.

% Optional Field,Value Inputs:
GithubDir = '/Users/yardenc/Documents/GitHub/'; % Tell the code where the github repos are
display_opt = 0; % Set to 1 to make the code display per-file data
use_residuals = -1; % What to do with signal before carrying stat tests:
% -1: do nothing. 0: remove dep. on current phrase duration. 2: remove dep.
% on phrase surations. 1: remove dep. on phrase durations and onset in song
% of target phrase.
extra_stat = 1; % set to 0 to avoid calculating stats.
compute_raster = 0; % set to 1 to create high-res raster
use_cohen2020 = 0; % set to 1 to load using cohen2020 data structure (if available).
signal_deconv_type = 'none'; % set which algorithm to run on the signal (none is the default).
% Other options:
%       'denoisedCalcium' - run foopsi and use the Ca2+ signal proxy
%       'denoisedSpikes' - run foopsi and use the spikes signal proxy
%       'HMM' - fit a 2 state HMM (signal \ noise) and use 0,1 for the two
%       states as the signal
delete_frames = 1; % will we ignore video frames?
n_del_frames = 5; % how many frames to ignore from each video
hvc_offset = 0.04; % offset time (sec) between HVC and audio
mulcnt = 2; % spacing between trials in 3d plot.
edges = [0 0]; % how much time in seconds to consider when inegrating signal outside phrases.
% set to [before after] to integrate more signal.
opacity_factor = 0.5; %opacity of 3d line plot
zscoring_type = 0; % set to 1 to z-score the signal
max_phrase_gap = 0.5; %maximal gap (sec) between phrases to be considered one song.
raster_edges = 2; % how many seconds to add around the raster.
label_mark_loc = 0; % Where (in sec) to place the line's color bar w.r.t. 0.
flag_loc_for_durations = 1; % not to be changed .. allows changing which durations to correlate with
add_sort_crit = []; % designate which columns in the behavior data matrix 'hits' to also use when re-ordering the data
% Those will be columns 3+add_sort_crit because the first 3 are already
% used; fnum phrasenum phrases.phraseFileStartTimes(phrase_idx(1))
sort_crit_loc = 0.5; % Where (in sec) to place the line's color bar w.r.t. 0 for additional sorting criteria (only if add_sort_crit is used).
multicomp = 0; % set to 1 to perform multiple comparisons after the stats.
anova_type = 1; % 1 - ANOVA, 2 - KruskalWallis

%% allow controlling parameters as function pair inputs
nparams=length(varargin);
for i=1:2:nparams
	switch lower(varargin{i})
        case 'use_cohen2020'
			use_cohen2020=varargin{i+1};
		case 'delete_frames'
			delete_frames=varargin{i+1};   
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
        case 'raster_edges'
			raster_edges=varargin{i+1}; 
        case 'compute_raster'
			compute_raster=varargin{i+1}; 
        case 'label_mark_loc'
            label_mark_loc=varargin{i+1}; 
        case 'flag_loc_for_durations'
            flag_loc_for_durations=varargin{i+1}; 
        case 'add_sort_crit'
            add_sort_crit = varargin{i+1}; 
        case 'sort_crit_loc'
            sort_crit_loc = varargin{i+1};
        case 'multicomp'
            multicomp = varargin{i+1};
        case 'anova_type'
            anova_type = varargin{i+1};
        case 'extra_stat' % will allow working woth syllable features if =8
            extra_stat = varargin{i+1};
        case 'githubdir'
            GithubDir = varargin{i+1};  
        case 'signal_deconv_type'
            signal_deconv_type = varargin{i+1};  
    end
end

%% repos
addpath(genpath([GithubDir 'VideoAnalysisPipeline']),'-end');
CNMFEfolder = fullfile(GithubDir,'CNMF_E_CohenLab');
addpath(genpath(CNMFEfolder),'-end');

%% parameters
%BaseDir = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/';
BaseDir = bird_params.base_dir; %'/Volumes/Labs/cohen/yardenc/Laptop2016_2021_backup/Documents/Experiments/Songbirds Imaging/Data/CanaryData';
%GithubDir = '/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/';
outvars.flag = 1; hndls = []; outstats.flag = 1;

%% Bird data
bird_name = bird_params.bird_name; %{1}; 
bird_folder_name = bird_params.bird_folder_name;%{2}; 
template_file = bird_params.template_filepath;%{3}; 
annotation_file = bird_params.annotation_filepath;%{4}; 
file_prefix = bird_params.file_prefix;%{5}; 
% Folders that contain data - legacy from cohen2020
if use_cohen2020
    manualROI_folder = fullfile(BaseDir,bird_folder_name,'ManualROIs');
end

%% Check that the input sequences are as expected.
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
    disp(['join or ignore lists overlap'])
    return;
end
if (any(order_flag == 0) || ...
   (any(order_flag < 0) & numel(order_flag) > 1) || ...
   any(order_flag > numel(syllabels_sequence)))
    disp(['Illegal order format'])
    return;
end

%% prepare target sequence as alphanumeric and the list of file numbers
sort_type = 1 - any(order_flag < 0); % separate sorting by phrase type vs. duration
% create figure if needed
if isempty(ax)
    h=figure('Visible','off','Position',[77          91        640         600]);
    ax = axes;
end
hndls = ax;
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
[file_indices,indx] = sort(ord);
elements = elements(indx);
keys = keys(indx);
dates = dates(indx,:);

%% Prepare the list of all files to analyze
if use_cohen2020 % legacy - load files from the folder structure supporting cohen2020
    cd([manualROI_folder '/ROIdata/' Day]);
    FILES = dir([file_prefix bird_name '*.mat']);
else
    FILES = dir(fullfile(BaseDir,bird_folder_name,Day,[file_prefix bird_name '*.mat']));
end
FILES_folder = FILES(1).folder;
FILES = {FILES.name};
hits = []; % will collect data on all files to process
durations = []; % durations of selected elements only
sylidx_durations = []; % duration of all numbered phrases in sequence and the onset and offset of the target phrase
sequence_total_lengths = []; % the total duration, including gaps, of the entire sequence
sequence_onsets = []; % durations of all phrases in sequence
if extra_stat == 8 % if we are using Brainard acoustic features (computed separately)
    AcousticFeatures = [];
    load(bird_params.brainard_features_filepath);
end
phrase_counts_and_position_in_song = []; % will collect phrase positions: last_onset_pos first_offset_pos phrasenum
for fnum = 1:numel(FILES)
    fname = FILES{fnum};
    loc = loc_file_in_locs(fname,file_indices);    
    try % in case we wish to ignore some syllable types
        ignore_locs = find(ismember(elements{loc}.segType,ignore_entries));
    catch eem
        continue;
    end
    elements{loc}.segAbsStartTimes(ignore_locs) = [];
    elements{loc}.segFileStartTimes(ignore_locs) = [];
    elements{loc}.segFileEndTimes(ignore_locs) = [];
    elements{loc}.segType(ignore_locs) = [];  
    if extra_stat == 8 % in case we use Brainard features
        feature_elements{loc}.syllable_duration(ignore_locs) = [];
        feature_elements{loc}.FF(ignore_locs) = [];
        feature_elements{loc}.time_to_half_peak(ignore_locs) = [];
        feature_elements{loc}.FF_slope(ignore_locs) = [];
        feature_elements{loc}.Amplitude_Slope(ignore_locs) = [];
        feature_elements{loc}.Spectral_Entropy(ignore_locs) = [];
        feature_elements{loc}.Temporal_Entropy(ignore_locs) = [];
        feature_elements{loc}.SpectroTemporal_Entropy(ignore_locs) = [];
    end  
    for i = 1:numel(join_entries) % in case we join some syllable types
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
    try % try to find the target sequence in the phrase seqyence
        endpoint = regexp(phrases_str,target_sequence_str,'end');
    catch em
        endpoint = [];
    end
    if ~isempty(endpoint)
        phrase_locs = endpoint - numel(target_sequence_str) + sylidx; %indices of target phrase in phrase sequence
        for phrase_loc = 1:numel(phrase_locs)
            phrase_idx = endpoint(phrase_loc) - numel(target_sequence_str) + [1:numel(target_sequence_str)]; %indices of target sequence
            try
            if ~any(gap_durations(phrase_idx(1:end-1)) > max_phrase_gap)       
                phrasenum = phrase_locs(phrase_loc); %index of target phrase
                tonset = phrases.phraseFileStartTimes(phrasenum);
                toffset = phrases.phraseFileEndTimes(phrasenum);
                syls_in_phrase = find(elements{loc}.segFileEndTimes >= tonset & ...
                elements{loc}.segFileStartTimes <= toffset); %indices of syllables in target phrase
                %max_syls = max(max_syls,numel(syls_in_phrase));
                hits = [hits; fnum phrasenum phrases.phraseFileStartTimes(phrase_idx(1)) phrases.phraseType(phrase_idx)' sum(phrase_durations(phrase_idx(sylidx)))];
                % hits holds for each discovered phrase sequence: 1. the
                % file number (in that day) 2. The phrase number in the
                % song file. 3. The time (in sec) of the beginning of the
                % target phrase sequence. 4. The types of the phrases in
                % the target sequence (there are variables there). 5. The
                % duration of the targeted phras (or phrases).
                last_onset_pos = max(find(phrases.phraseType(1:phrasenum) == -1000));
                first_offset_pos = phrasenum - 1 + min(find(phrases.phraseType(phrasenum:end) == 1000));
                phrase_counts_and_position_in_song = [phrase_counts_and_position_in_song; last_onset_pos first_offset_pos phrasenum];
                % phrase_counts_and_position_in_song holds 1. The phrase index of the last
                % onset before the current target phrase sequence. 2. The
                % first ofset after the current sequence. 3. The index of
                % the current target phrase.
                sylidx_durations = [sylidx_durations; phrase_durations(phrase_idx(~isnan(syllabels_sequence))) phrases.phraseFileStartTimes(phrase_idx(~isnan(syllabels_sequence) & ([1:numel(syllabels_sequence)]~=sylidx)))-tonset tonset toffset];
%                     phrases.phraseFileStartTimes(phrase_idx(~isnan(syllabels_sequence) & ([1:numel(syllabels_sequence)]~=sylidx)))-tonset ...
%                     phrases.phraseFileEndTimes(phrase_idx(~isnan(syllabels_sequence)))-tonset ...
                     %(sylidx)
                if sort_type == 1
                    durations = [durations; sum(phrase_durations(phrase_idx(order_flag.*(order_flag > 0))))];
                else
                    durations = [durations; sum(phrase_durations(phrase_idx(sylidx)))];
                end
                sequence_onsets = [sequence_onsets; (phrases.phraseFileStartTimes(phrase_idx)-tonset)*locktoonset + ...
                    (phrases.phraseFileEndTimes(phrase_idx)-toffset)*(1-locktoonset)];
                sequence_total_lengths = [sequence_total_lengths; phrases.phraseFileEndTimes(phrase_idx(end)) - phrases.phraseFileStartTimes(phrase_idx(1))];
                if extra_stat == 8
                    AcousticFeatures = [AcousticFeatures; ...
                        mean(feature_elements{loc}.syllable_duration(syls_in_phrase)) ...
                        mean(feature_elements{loc}.FF(syls_in_phrase)) ...
                        mean(feature_elements{loc}.time_to_half_peak(syls_in_phrase)) ...
                        mean(feature_elements{loc}.FF_slope(syls_in_phrase)) ...
                        mean(feature_elements{loc}.Amplitude_Slope(syls_in_phrase)) ...
                        mean(feature_elements{loc}.Spectral_Entropy(syls_in_phrase)) ...
                        mean(feature_elements{loc}.Temporal_Entropy(syls_in_phrase)) ...
                        mean(feature_elements{loc}.SpectroTemporal_Entropy(syls_in_phrase)) ...
                        ];
                end
            end
            catch em
                '4'; % debug point
            end
            
        end
    end
end

try 
    if sort_type == 1
        [durations,dur_idx] = sort(durations);
        id_flags = hits(dur_idx,3 + flag_loc_for_durations);
    else
        %[id_flags,dur_idx] = sort(hits(:,3 - order_flag));
        [~,dur_idx] = sortrows(hits,[3 - order_flag 3+add_sort_crit size(hits,2)]);
        id_flags = hits(dur_idx,3 - order_flag);
        durations = durations(dur_idx);
    end
    hits = hits(dur_idx,:);
    if extra_stat == 8
        AcousticFeatures = AcousticFeatures(dur_idx,:);
    end
    sylidx_durations = sylidx_durations(dur_idx,:);
    sequence_onsets = sequence_onsets(dur_idx,:);
catch em
    'f';
    return;
end
outvars.phrase_counts_and_position_in_song = phrase_counts_and_position_in_song(dur_idx,:);

%% Collect data matrices and draw signals
sig_integrals_in = [];
sig_integrals_before = [];
%sig_com_before = [];
if (compute_raster > 0)
    max_sequence_duration = max(sequence_total_lengths);
    signal_raster = zeros(size(hits,1),ceil(max_sequence_duration*1000+200+2000*raster_edges))*nan; 
    phraseType_raster = signal_raster;
end
for cnt = 1:size(hits,1)
    fnum = hits(cnt,1);
    phrasenum = hits(cnt,2);
    fname = FILES{fnum};
    loc = loc_file_in_locs(fname,file_indices);
    phrases = return_phrase_times(elements{loc});
    phrases = deal_with_time_gaps(phrases,max_phrase_gap);
   
    if use_cohen2020
        load(fname);
    else
        [vidTimes, dff] = load_dff_file(fullfile(FILES_folder,fname));
    end    
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
    switch signal_deconv_type
        case 'denoisedCalcium'
            try 
                [c, s, options] = deconvolveCa(y(ROI,:),'ar2',[1.3 -0.422],'method','thresholded','optimize_b','optimize_smin');%,'optimize_pars');
            catch em
                try
                [c, s, options] = deconvolveCa((y(ROI,:)),'ar2','method','foopsi','optimize_b',1);
                catch emm
                    'g'
                end
            end
            signal = c + options.b*(options.b < 0);
        case 'denoisedSpikes'
            try 
                [c, s, options] = deconvolveCa(y(ROI,:),'ar2',[1.3 -0.422],'method','thresholded','optimize_b','optimize_smin');
            catch em
                [c, s, options] = deconvolveCa((y(ROI,:)),'ar2','method','foopsi','optimize_b',1);
            end
            signal = s;
        case 'none'
            signal = y(ROI,:); %smooth(y(ROI,:),3);
        case 'HMM'
            sig = y(ROI,:);
            clear sigma;
            nstates = 2;
            dq = (0.5-0.05)/(nstates-2+1e-5);          
            mu = [quantile(sig,0.05+[0:nstates-2]*dq) max(sig) ];
            sigma(1,1,:) = [0.05*ones(1,nstates)];
            CPD = condGaussCpdCreate(mu, sigma);
            [model, loglikHist] = hmmFit(sig, nstates, 'gauss','pi0',[zeros(1,nstates-1) 1],'emission0',CPD);                    
            path = hmmMap(model, sig)-1;
            signal = abs(median(path)-path)';
        %case 4
        %    signal = full(dff(ROI,n_del_frames+1:end));
    end
    tonset = phrases.phraseFileStartTimes(phrasenum);
    toffset = phrases.phraseFileEndTimes(phrasenum);
    %com_edge = (tonset - edges(1))*locktoonset + (toffset + edges(2))*(1-locktoonset);
    if (display_opt)
        display([num2str([cnt tonset]) ' ' fname]);
    end
    sig_integrals_in = [sig_integrals_in; ...
                sum(signal((t >= tonset-edges(1)) & (t <= toffset+edges(2))))];
    %sig_integrals_before = [sig_integrals_before; ...
    %    sum(signal((t <= tonset-edges(1)) & (t >= hits(cnt,3))))];
    timetag = (t-tonset*locktoonset-(1-locktoonset)*toffset);
    t_phr_on = (phrases.phraseFileStartTimes-tonset*locktoonset-(1-locktoonset)*toffset);
    t_phr_off = (phrases.phraseFileEndTimes-tonset*locktoonset-(1-locktoonset)*toffset);
    if (compute_raster >= 1)
        interpolated_signal = interp1(timetag,signal,timetag(1):0.001:timetag(end)+1/30);
        interpolated_timetag = interp1(timetag,timetag,timetag(1):0.001:timetag(end)+1/30);
        interpolated_phrase_id = zeros(size(interpolated_timetag));
        for phnum = 1:numel(phrases.phraseType)
            interpolated_phrase_id((interpolated_timetag >= t_phr_on(phnum)) & (interpolated_timetag <= t_phr_off(phnum))) = phrases.phraseType(phnum);         
        end
        idxmap = [1:numel(interpolated_timetag)] - min(find(abs(interpolated_timetag) == min(abs(interpolated_timetag))))+(1000*raster_edges+100)+round((1-locktoonset)*max_sequence_duration*1000);
        t_on = (phrases.phraseFileStartTimes(1)-tonset*locktoonset-(1-locktoonset)*toffset);
        t_off = (phrases.phraseFileEndTimes(end)-tonset*locktoonset-(1-locktoonset)*toffset);
        idxs = find((interpolated_timetag >= t_on) & (interpolated_timetag <= t_off) & ...
                    (interpolated_timetag >= -raster_edges*locktoonset  + (-max_sequence_duration-raster_edges)*(1-locktoonset)) & ...
                    (interpolated_timetag <= raster_edges*(1-locktoonset)+locktoonset*(raster_edges+max_sequence_duration)));
        signal_raster(cnt,idxmap(idxs)) = interpolated_signal(idxs);
        phraseType_raster(cnt,idxmap(idxs)) = interpolated_phrase_id(idxs);
    end
        
    for currphrase = 1:numel(phrases.phraseType)
        try
        plot3(ax,timetag(t >= phrases.phraseFileStartTimes(currphrase) & ...
             t <= phrases.phraseFileEndTimes(currphrase)), ...
             cnt*mulcnt*ones(1,sum(t >= phrases.phraseFileStartTimes(currphrase) & ...
             t <= phrases.phraseFileEndTimes(currphrase))), ...
             signal(t >= phrases.phraseFileStartTimes(currphrase) & ...
             t <= phrases.phraseFileEndTimes(currphrase))+0, ...
             'LineWidth',2,'Color',[colors(find(syllables == phrases.phraseType(currphrase)),:) opacity_factor]);
         hold on;
        catch em1
            'd'; % debug point
        end
         if (currphrase < numel(phrases.phraseType))
             startidx = max(find(t <= phrases.phraseFileEndTimes(currphrase)));
             stopidx = min(find(t >= phrases.phraseFileStartTimes(currphrase+1)));
             plot3(ax,timetag(startidx:stopidx), ...
             cnt*mulcnt*ones(1,numel(startidx:stopidx)), ...
             signal(startidx:stopidx)+0, ...
             'LineWidth',2,'Color',[0 0 0 0.4],'LineStyle','-');
         end
    end
    plot3(ax,timetag(t >= phrases.phraseFileEndTimes(currphrase) & ...
             t <= phrases.phraseFileEndTimes(currphrase)+2), ...
             cnt*mulcnt*ones(1,sum(t >= phrases.phraseFileEndTimes(currphrase) & ...
             t <= phrases.phraseFileEndTimes(currphrase)+2)), ...
             signal(t >= phrases.phraseFileEndTimes(currphrase) & ...
             t <= phrases.phraseFileEndTimes(currphrase)+2)+0, ...
             'LineWidth',2,'Color',[0 0 0 0.4],'LineStyle','--');
    
    set(gca,'color','none');
    title(['Phrase #' 'TBD' ' locked Ca signals from ' datestr(Day,'yyyy-mm-dd')])
    set(gca,'XTick',[0 1]*locktoonset+(1-locktoonset)*[-1 0]);
    set(gca,'YTick',[]);
    ylabel(['ROI# ' num2str(ROI)]);
    xlabel('Real Time');
    set(gca,'FontSize',16);
    axis tight;
    xlim([-1 3]-(1-locktoonset));         
end
set(ax,'CameraPosition', [-0.3040 -360.4786 3.7984]);
%%
%%%%%%%%%%% Perform all the statistical tests
if ~isempty(durations) & (extra_stat == 1) % run only if needed
    outstats.use_residuals = use_residuals;
    if use_residuals > 0 %remove dependencies on covariates if needed
        try
            [~,~,sig_integrals_in,~,~] = mvregress([ones(size(hits,1),1) sylidx_durations(:,1:end-use_residuals)],sig_integrals_in);
        catch em1
           sig_integrals_in = linear_res(sig_integrals_in,durations);
        end
    elseif use_residuals == 0
        [~,~,sig_integrals_in,~,~] = mvregress([ones(size(hits,1),1) durations],sig_integrals_in);
    else
        'do nothing';
    end
    outvars.sig_integrals_in = sig_integrals_in;
    if sort_type ~= 1 
        if ismember(anova_type,1:2)
            if anova_type == 1
                [p,ANOVATAB,STATS] = anova1(sig_integrals_in,id_flags); %one-way ANOVA
                outstats.type = '1-way ANOVA';
            elseif anova_type == 2
                [p,ANOVATAB,STATS] = kruskalwallis(sig_integrals_in,id_flags); %Kruskal-Wallis
                outstats.type = '1-way KruskalWallis';
            end
            boxplot(sig_integrals_in,id_flags);
            hndls = [hndls; gca];
            outstats.p = p;
            r = ANOVATAB{2,5}; outstats.r = r;
            gnames = STATS.gnames; outvars.gnames = gnames;
            outvars.id_flags = id_flags;
            if (multicomp == 1 && p < 0.05)
                figure;
                try
                    outstats.COMPARISON = multcompare(STATS);
                catch em
                    outstats.COMPARISON = nan;
                end
            else
                outstats.COMPARISON = nan;
            end
        elseif ismember(anova_type,3:4)
            if anova_type == 3 % 2-way ANOVA
                g1 = id_flags;
                g2 = hits(:,3+add_sort_crit);
                [p,ANOVATAB,STATS] = anovan(sig_integrals_in,{g2,g1},'display','off','model','linear');
                outstats.type = '2-way ANOVA';
                outstats.p = p;
                hndls = [hndls; gca];
                %figure;
                try
                    outvars.COMPARISON2 = multcompare(STATS,'display','off');
                catch em
                    outvars.COMPARISON2 = nan;
                end
                [p,ANOVATAB,STATS] = anovan(sig_integrals_in,{g1,g2},'display','off','model','linear');
                hndls = [hndls; gca];
                outstats.p = p;
                r = ANOVATAB{2,5}; outstats.r = r;
                gnames = STATS.grpnames; outvars.gnames = gnames;
                outvars.id_flags = [g1 g2];
                %figure;
                try
                    outstats.COMPARISON1 = multcompare(STATS,'display','off');
                catch em
                    outstats.COMPARISON1 = nan;
                end
                %pause;
                
                %pause;
            elseif anova_type == 4 %Kruskal Wallis
                [p,ANOVATAB,STATS] = kruskalwallis(sig_integrals_in,id_flags);
                r = ANOVATAB{2,5}; outstats.r = r;
                outstats.p = p;
                outstats.type = 'KruskalWallis';
            end
        end
        %[pp,ANOVATAB,STATS] = anova1(durations,id_flags);
        
    else    
        [r, p] = corr(durations,sig_integrals_in); % Pearson correlation
        outstats.type = 'Pearson Corr.';
        outstats.r = r; outstats.p = p;
        outvars.durations = durations;
        figure('Visible','on'); plot(durations,sig_integrals_in,'bo','MarkerSize',12,'MarkerFaceColor','none','MarkerEdgeColor','k');
        hndls = [hndls; gca];
        hold on;
        segtypes = unique(id_flags);
        per_seg_stats = [];
        for segt = 1:numel(segtypes)
            ntrials = find(id_flags == segtypes(segt));
            plot(durations(ntrials),sig_integrals_in(ntrials),'o','MarkerSize',10,'MarkerFaceColor',colors(find(syllables == segtypes(segt)),:),'MarkerEdgeColor','none');
            [r_tag, p_tag] = corr(durations(ntrials),sig_integrals_in(ntrials)); 
            disp([segtypes(segt) r_tag, p_tag]);
            per_seg_stats = [per_seg_stats; segtypes(segt) r_tag, p_tag];
        end
        set(gca,'FontSize',16); xlabel('Durations'); ylabel('Signal Integral'); title(['(r,p) = ' num2str([r,p])]);
        outstats.per_seg_stats = per_seg_stats;       
%         [r1, p1] = corr(durations,sylidx_durations(:,1)); 
%         figure('Visible','on'); plot(durations,sylidx_durations(:,1),'ro','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','none');
%         set(gca,'FontSize',16); xlabel('Durations'); ylabel('idx durations'); title(['(r,p) = ' num2str([r1,p1])]);
%         [r1, p1] = corr(sig_integrals_in,sylidx_durations(:,2)); 
%         figure('Visible','on'); plot(durations,sylidx_durations(:,2),'mo','MarkerSize',10,'MarkerFaceColor','m','MarkerEdgeColor','none');
%         set(gca,'FontSize',16); xlabel('sig.integral'); ylabel('onset'); title(['(r,p) = ' num2str([r1,p1])]);
        
    end
    
else
    p = 1; r = 0; gnames = numel(durations);
end

if extra_stat == 8
    outstats.use_residuals = 'Acoustic';
    try
        [~,~,sig_integrals_in,~,~] = mvregress([ones(size(hits,1),1) AcousticFeatures],sig_integrals_in);
        outvars.sig_integrals_in = sig_integrals_in;

        if anova_type == 1
            [p,ANOVATAB,STATS] = anova1(sig_integrals_in,id_flags);
            outstats.type = '1-way ANOVA';
        elseif anova_type == 2
            [p,ANOVATAB,STATS] = kruskalwallis(sig_integrals_in,id_flags);
            outstats.type = '1-way KruskalWallis';      
        end
        outstats.p = p;
        boxplot(sig_integrals_in,id_flags);
        hndls = [hndls; gca];
        r = ANOVATAB{2,5}; outstats.r = r;
        gnames = STATS.gnames; outvars.gnames = gnmaes;
        outvars.id_flags = id_flags;
        if (multicomp == 1 && p < 0.05)
            figure;
            try
                outvars.COMPARISON = multcompare(STATS);
            catch em
                outvars.COMPARISON = nan;
            end
        else
            outvars.COMPARISON = nan;
        end
    catch em
        p = 1; r = 0; gnames = numel(durations);
    end
end

%% high-res raster plots 
if compute_raster == 1
    %line([2100 2100],[0 size(hits,1)+0.5],'Color','w');
    zeroidx = (100+raster_edges*1000)*locktoonset + (1-locktoonset)*(size(signal_raster,2)-(100+raster_edges*1000));
    figure; imagesc(signal_raster); colormap(copper); yticks([]); hold on;
    for i = 1:size(hits,1)
        for j = 1:size(sequence_onsets,2)
            line([zeroidx+1000*sequence_onsets(i,j) zeroidx+1000*sequence_onsets(i,j)],[i-0.5 i+0.5],'Color','w'); %sylidx_
        end
    end
    if sort_type ~= -1 
        if sort_type == -1 
            segtypes = unique(id_flags);
            for segt = 1:numel(segtypes)
                ntrials = find(id_flags == segtypes(segt));% ntrials(1) = ntrials(1)-0.5; ntrials(end) = ntrials(end)+0.5;
                plot(label_mark_loc*1000*ones(2,1)+zeroidx,[ntrials(1)-0.4; ntrials(end)+0.4],'Color',... %numel(ntrials)
                    colors(find(syllables == segtypes(segt)),:),'LineWidth',5);
            end
        else
            segtypes = id_flags;
            for segt = 1:numel(segtypes)
                %ntrials = find(id_flags == segtypes(segt));% ntrials(1) = ntrials(1)-0.5; ntrials(end) = ntrials(end)+0.5;
                plot(label_mark_loc*1000*ones(2,1)+zeroidx,[segt-0.4; segt+0.4],'Color',... %numel(ntrials)
                    colors(find(syllables == segtypes(segt)),:),'LineWidth',5);
            end
        end
        if ~isempty(add_sort_crit)
           segtypes = hits(:,[3+add_sort_crit]);
           outvars.segtypes = segtypes;
            for segt = 1:size(segtypes,1)
                for itrnum = 1:size(segtypes,2)
                    %ntrials = find(id_flags == segtypes(segt));% ntrials(1) = ntrials(1)-0.5; ntrials(end) = ntrials(end)+0.5;
                    plot(sort_crit_loc*1000*ones(2,1)+zeroidx+(itrnum-1)*50,[segt-0.4; segt+0.4],'Color',... %numel(ntrials)
                        colors(find(syllables == segtypes(segt,itrnum)),:),'LineWidth',3);
                end
            end 
        end
    end
    %
    if locktoonset == 1
        xticks(raster_edges*1000+[100 1100]); xticklabels([0 1]);
    else
        xticks([zeroidx-1000 zeroidx]); xticklabels([-1 0]);
    end
    hndls = [hndls; gca];
end

if compute_raster == 2 % Overlay signals as semi-transparent fills
    figure; 
    signal_raster(isnan(signal_raster)) = 0;
    phraseType_raster(isnan(phraseType_raster)) = 0;
    for row_num = 1:size(signal_raster,1)
        line_onsets = find(diff([-500 phraseType_raster(row_num,:)]) ~= 0); 
        line_offsets = find(diff([phraseType_raster(row_num,:) -500]) ~=0);
        for line_num = 1:numel(line_onsets)
            xs = [line_onsets(line_num):line_offsets(line_num) line_offsets(line_num):-1:line_onsets(line_num)];
            ys = row_num*0.03+[0.01+signal_raster(row_num,line_onsets(line_num):line_offsets(line_num)) zeros(1,numel(line_onsets(line_num):line_offsets(line_num)))-0.01]; %-fliplr(signal_raster(row_num,line_onsets(line_num):line_offsets(line_num)))
            segt = median(phraseType_raster(row_num,line_onsets(line_num):line_offsets(line_num)));
            if (segt == 0)
                col = [0.5 0.5 0.5];
            else
                col = colors(find(syllables == segt),:);
            end
            fill(xs,ys,col,'LineStyle','none','FaceAlpha',0.5);
            hold on;
        end
    end        
    yticks([]); hold on;   
    %line([2100 2100],[0 size(hits,1)+0.5],'Color','w');
    zeroidx = (100+raster_edges*1000)*locktoonset + (1-locktoonset)*(size(signal_raster,2)-(100+raster_edges*1000));
    if locktoonset == 1
        xticks(raster_edges*1000+[100 1100]); xticklabels([0 1]);
    else
        xticks([zeroidx-1000 zeroidx]); xticklabels([-1 0]);
    end
    hndls = [hndls; gca];
end
%%
%set(ax,'CameraTarget', [1 36 0.2605]);
% helper funcitons
function loc = loc_file_in_locs(fname,locs)
    tokens = regexp(fname,'_','split');
    if use_cohen2020
        loc = find(locs == str2num(tokens{3}));
    else
        name_idx = find(cellfun(@(x)strcmp(x,bird_name),tokens));
        loc = find(locs == str2num(tokens{name_idx+1}));
    end
end
function resid = linear_res(vec_a,vec_b)
    % returns the residuals of vec_a after removing the linear component of
    % vec_b
    tmpcov = cov(vec_a,vec_b,1);
    beta = tmpcov(1,2)/var(vec_b,1);
    alpha = mean(vec_a) - beta*mean(vec_b);
    resid = vec_a - alpha - beta*vec_b;
end

function [v_times, dff_mat] = load_dff_file(file_path)
% script assumes file has only matrix variables that are 3d or 1d (with
% size [1, l])
% except the 3d video matrix the file can have 1 audio times vector and one
% video times vector and that's it. They will be identified by their length
    S = whos('-file',file_path);
    tempfile = load(file_path);
    v_times=[]; dff_mat=[];
    for s_ind = 1:numel(S)
        if S(s_ind).size(1) == 1
            eval(['v_times = tempfile.' S(s_ind).name ';']);
        else
            eval(['dff_mat = tempfile.' S(s_ind).name ';']);
        end
    end
end
end



