function [ax,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,Day,ignore_entries,join_entries,sylidx,syllabels_sequence,ROI,locktoonset,spikes,order_flag)
%% 
% This script creates single ROI single day alignments to complex sequences
% with only one variable (duration or syllable type)
% Inputs:
%   ax - axes where to plot
%   Day - text rep. of date
%   ignore_entries - A vector of label numbers to ignore completely. 
%   join_entries - A cell of vectors, each containing a >1 number of labels
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
%   durations. Use negative integers for type (no vectors)

%% change according to workstation
BaseDir = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/';
GithubDir = '/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/';
display_opt = 0;
%%
bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation4'};
bird_params = bird1_params;
delete_frames = 1;
n_del_frames = 6;
hvc_offset = 0.04;
mulcnt = 2;
edges = [0.25 1];
opacity_factor = 0.5;
zscoring_type = 0;
max_phrase_gap = 0.5;
%%
%addpath(genpath([GithubDir 'small-utils']),'-end');
%addpath(genpath([GithubDir 'VideoAnalysisPipeline']),'-end');
bird_name = bird_params{1}; 
bird_folder_name = bird_params{2}; 
template_file = bird_params{3}; 
annotation_file = bird_params{4}; 
CNMFEfolder = [GithubDir 'CNMF_E'];
%addpath(genpath(CNMFEfolder),'-end');

%% Folders that contain data
% Folders on laptop:
laptop_mov_folder = [BaseDir bird_folder_name '/movs'];
laptop_wav_folder = [BaseDir bird_folder_name '/movs/wav'];
laptop_gif_folder = [BaseDir bird_folder_name '/movs/wav/gif'];
laptop_annotated_dir = [BaseDir bird_folder_name '/movs/wav/annotated'];
laptop_annotated_images_dir = [BaseDir bird_folder_name '/movs/wav/annotated/images'];
DamagedFolder = [BaseDir bird_folder_name '/too_large_or_damaged/'];
laptop_manualROI_folder = [BaseDir bird_folder_name '/ManualROIs'];

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
sort_type = 1 - any(order_flag < 0);
      
if isempty(ax)
    h=figure('Visible','off','Position',[77          91        640         600]);
    ax = axes;
end
%%
AlphaNumeric = '{}ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';
target_sequence_str = '';
cd (laptop_manualROI_folder);
%load('syl_dur_gap_stat.mat');
load(template_file);
syllables = [[templates.wavs.segType] -1 102 103];
%syllables = [templates.wavs.segType];
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
cd([laptop_manualROI_folder '/ROIdata/' Day]);
FILES = dir('NonoverlapBaseROIdata*.mat');
FILES = {FILES.name};
hits = [];
max_syls = 0;
durations = [];
sequence_durations = [];
for fnum = 1:numel(FILES)
    fname = FILES{fnum};
    tokens = regexp(fname,'_','split');
    loc = find(locs == str2num(tokens{3}));
    ignore_locs = find(ismember(elements{loc}.segType,ignore_entries));
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
                hits = [hits; fnum phrasenum phrases.phraseFileStartTimes(phrase_idx(1)) phrases.phraseType(phrase_idx)'];
                if sort_type == 1
                    durations = [durations; sum(phrase_durations(phrase_idx(order_flag.*(order_flag > 0))))];
                end
                sequence_durations = [sequence_durations; phrase_durations(phrase_idx)];
            end
            catch em
                '4';
            end
            
        end
    end
end

try 
    if sort_type == 1
        [durations,dur_idx] = sort(durations);
    else
        [durations,dur_idx] = sort(hits(:,3 - order_flag));
    end
    hits = hits(dur_idx,:);
catch em
    'f';
    p=1; r=0; gnames=''; return;
end


%%

sig_integrals_in = [];
sig_integrals_before = [];
sig_com_before = [];
max_signal = 0;
for cnt = 1:size(hits,1)
    fnum = hits(cnt,1);
    phrasenum = hits(cnt,2);
    fname = FILES{fnum};
    tokens = regexp(fname,'_','split');
    loc = find(locs == str2num(tokens{3}));
    phrases = return_phrase_times(elements{loc});
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
            signal = smooth(y(ROI,:),3);
        case 3
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
    end
    tonset = phrases.phraseFileStartTimes(phrasenum);
    toffset = phrases.phraseFileEndTimes(phrasenum);
    com_edge = (tonset - edges(1))*locktoonset + (toffset + edges(2))*(1-locktoonset);
    if (display_opt)
        display([num2str([cnt tonset]) ' ' fname]);
    end
    sig_integrals_in = [sig_integrals_in; ...
                sum(signal((t >= tonset-edges(1)) & (t <= toffset+edges(2))))];
    sig_integrals_before = [sig_integrals_before; ...
        sum(signal((t <= tonset-edges(1)) & (t >= hits(cnt,3))))];
    timetag = (t-tonset*locktoonset-(1-locktoonset)*toffset);
    sig_com_before = [sig_com_before; ...
        sum(signal((t <= com_edge) & (t >= hits(cnt,3))).*timetag(((t <= com_edge) & (t >= hits(cnt,3))))')/ ...
        sum(signal((t <= com_edge) & (t >= hits(cnt,3))))];
    
    for currphrase = 1:numel(phrases.phraseType)
        plot3(ax,timetag(t >= phrases.phraseFileStartTimes(currphrase) & ...
             t <= phrases.phraseFileEndTimes(currphrase)), ...
             cnt*mulcnt*ones(1,sum(t >= phrases.phraseFileStartTimes(currphrase) & ...
             t <= phrases.phraseFileEndTimes(currphrase))), ...
             signal(t >= phrases.phraseFileStartTimes(currphrase) & ...
             t <= phrases.phraseFileEndTimes(currphrase))+0, ...
             'LineWidth',2,'Color',[colors(find(syllables == phrases.phraseType(currphrase)),:) opacity_factor]);
         hold on;
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
r=[]; p=[];
set(ax,'CameraPosition', [-0.3040 -360.4786 3.7984]);
if ~isempty(durations)
    if sort_type ~= 1 
        [p,ANOVATAB,STATS] = anova1(sig_integrals_in,durations);
        r = ANOVATAB{2,5};
        gnames = STATS.gnames;
    else    
        [r, p] = corr(durations,sig_integrals_in); 
        gnames = figure('Visible','off'); plot(durations,sig_integrals_in,'bo','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','none');
        set(gca,'FontSize',16); xlabel('Durations'); ylabel('Signal Integral')
    end
else
    p = 1; r = 0; gnames = 0;
end
%set(ax,'CameraTarget', [1 36 0.2605]);
end

% function new_phrases = deal_with_time_gaps(phrases,max_phrase_gap)
%     if isempty(phrases.phraseType)
%         new_phrases = phrases;
%         return;
%     end
%     new_phrases.phraseType = [-1000; phrases.phraseType(1)];
%     new_phrases.phraseFileStartTimes = [phrases.phraseFileStartTimes(1) - 1, ...
%         phrases.phraseFileStartTimes(1)];
%     new_phrases.phraseFileEndTimes = [phrases.phraseFileStartTimes(1) - 0.001, ...
%         phrases.phraseFileEndTimes(1)];
%     
%     for temp_cnt = 2:numel(phrases.phraseType)
%         curr_gap = (phrases.phraseFileStartTimes(temp_cnt) - phrases.phraseFileEndTimes(temp_cnt-1));
%         if (curr_gap > max_phrase_gap) 
%             new_phrases.phraseType = [new_phrases.phraseType; 1000; -1000];
%             new_phrases.phraseFileStartTimes = [new_phrases.phraseFileStartTimes phrases.phraseFileEndTimes(temp_cnt - 1)+0.001 ...
%                 phrases.phraseFileStartTimes(temp_cnt) - min(0.005,curr_gap/10)];
%             new_phrases.phraseFileEndTimes = [new_phrases.phraseFileEndTimes phrases.phraseFileEndTimes(temp_cnt - 1)+min(0.005,curr_gap/10) ...
%                 phrases.phraseFileStartTimes(temp_cnt) - 0.001];
%         end
%         new_phrases.phraseType = [new_phrases.phraseType; phrases.phraseType(temp_cnt)];
%         new_phrases.phraseFileStartTimes = [new_phrases.phraseFileStartTimes phrases.phraseFileStartTimes(temp_cnt)];
%         new_phrases.phraseFileEndTimes = [new_phrases.phraseFileEndTimes phrases.phraseFileEndTimes(temp_cnt)];
%                
%     end
%     new_phrases.phraseType = [new_phrases.phraseType; 1000];
%     new_phrases.phraseFileStartTimes = [new_phrases.phraseFileStartTimes phrases.phraseFileEndTimes(end) + 0.001];
%     new_phrases.phraseFileEndTimes = [new_phrases.phraseFileEndTimes phrases.phraseFileEndTimes(end) + 1];
% end 

