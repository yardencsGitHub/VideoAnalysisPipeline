function [ax,r,p,gnames] = LongRangeLockedSingleDayManualROIs_function(ax,Day,ignore_entries,join_entries,sylidx,syllabels_sequence,ROI,locktoonset,spikes,order_flag,varargin)
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
use_residuals = 0;
extra_stat = 1;
compute_raster = 0;
file_prefix = 'baseROIdata_'; %'NonoverlapBaseROIdata_'; %
%%
bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix' 'NonoverlapBaseROIdata_'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation5_alexa' 'baseROIdata_'};
bird_params = bird2_params;
delete_frames = 1;
n_del_frames = 6;
hvc_offset = 0.04;
mulcnt = 2;
edges = [0 0];
opacity_factor = 0.5;
zscoring_type = 0;
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
        case 'compute_raster'
			compute_raster=varargin{i+1}; 
    end
end
            

%%
addpath(genpath([GithubDir 'small-utils']),'-end');
%addpath(genpath([GithubDir 'VideoAnalysisPipeline']),'-end');
bird_name = bird_params{1}; 
bird_folder_name = bird_params{2}; 
template_file = bird_params{3}; 
annotation_file = bird_params{4}; 
file_prefix = bird_params{5}; 
CNMFEfolder = [GithubDir 'CNMF_E'];
addpath(genpath(CNMFEfolder),'-end');

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
cd([laptop_manualROI_folder '/ROIdata/' Day]);
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
   
    loc = find(locs == str2num(tokens{3}));
    
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
                hits = [hits; fnum phrasenum phrases.phraseFileStartTimes(phrase_idx(1)) phrases.phraseType(phrase_idx)'];
                sylidx_durations = [sylidx_durations; phrase_durations(phrase_idx(~isnan(syllabels_sequence))) tonset toffset]; %(sylidx)
                if sort_type == 1
                    durations = [durations; sum(phrase_durations(phrase_idx(order_flag.*(order_flag > 0))))];
                else
                    durations = [durations; sum(phrase_durations(phrase_idx(sylidx)))];
                end
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
    if sort_type == 1
        [durations,dur_idx] = sort(durations);
    else
        [id_flags,dur_idx] = sort(hits(:,3 - order_flag));
        durations = durations(dur_idx);
    end
    hits = hits(dur_idx,:);
    sylidx_durations = sylidx_durations(dur_idx,:);
catch em
    'f';
    p=1; r=0; gnames=''; return;
end


%%

sig_integrals_in = [];
sig_integrals_before = [];
sig_com_before = [];
max_signal = 0;
if (compute_raster == 1)
    max_sequence_duration = max(sequence_total_lengths);
    signal_raster = zeros(size(hits,1),ceil(max_sequence_duration*1000+4200))*nan; 
end
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
            signal = c + options.b*(options.b < 0);
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
    if (compute_raster == 1)
        interpolated_signal = interp1(timetag,signal,timetag(1):0.001:timetag(end)+1/30);
        interpolated_timetag = interp1(timetag,timetag,timetag(1):0.001:timetag(end)+1/30);
        idxmap = [1:numel(interpolated_timetag)] - min(find(abs(interpolated_timetag) == min(abs(interpolated_timetag))))+2100+round((1-locktoonset)*max_sequence_duration*1000);
        t_on = (phrases.phraseFileStartTimes(1)-tonset*locktoonset-(1-locktoonset)*toffset);
        t_off = (phrases.phraseFileEndTimes(end)-tonset*locktoonset-(1-locktoonset)*toffset);
        idxs = find((interpolated_timetag >= t_on) & (interpolated_timetag <= t_off) & ...
                    (interpolated_timetag >= -2*locktoonset  + (-max_sequence_duration-2)*(1-locktoonset)) & ...
                    (interpolated_timetag <= 2*(1-locktoonset)+locktoonset*(2+max_sequence_duration)));
        signal_raster(cnt,idxmap(idxs)) = interpolated_signal(idxs);
    end
    
    sig_com_before = [sig_com_before; ...
        sum(signal((t <= com_edge) & (t >= hits(cnt,3))).*timetag(((t <= com_edge) & (t >= hits(cnt,3))))')/ ...
        sum(signal((t <= com_edge) & (t >= hits(cnt,3))))];
    
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
            'd';
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
r=[]; p=[]; 
set(ax,'CameraPosition', [-0.3040 -360.4786 3.7984]);
if ~isempty(durations) & (extra_stat == 1)
    if use_residuals == 1
        try
            [~,~,sig_integrals_in,~,~] = mvregress([ones(size(hits,1),1) sylidx_durations(:,1:end-1)],sig_integrals_in);
        catch em1
           'd'; 
           sig_integrals_in = linear_res(sig_integrals_in,durations);
        end
        %sig_integrals_in =
        %linear_res(sig_integrals_in,sylidx_durations(:,1)); this is set to
        %1 because the original sylidx_durations kept only the duration of
        %the target phrase .. change it it uncommenting
    end
    if sort_type ~= 1 
        [p,ANOVATAB,STATS] = anova1(sig_integrals_in,id_flags);
        r = ANOVATAB{2,5};
        gnames = STATS.gnames;
        %[pp,ANOVATAB,STATS] = anova1(durations,id_flags);
        
    else    
        [r, p] = corr(durations,sig_integrals_in); 
        gnames = figure('Visible','off'); plot(durations,sig_integrals_in,'bo','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','none');
        set(gca,'FontSize',16); xlabel('Durations'); ylabel('Signal Integral'); title(['(r,p) = ' num2str([r,p])]);
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
if compute_raster == 1
    figure; imagesc(signal_raster); colormap(copper); yticks([]); hold on;   
    %line([2100 2100],[0 size(hits,1)+0.5],'Color','w');
    zeroidx = 2100*locktoonset + (1-locktoonset)*(size(signal_raster,2)-2100);
    for i = 1:size(hits,1)
        for j = 1:size(sequence_onsets,2)
            line([zeroidx+1000*sequence_onsets(i,j) zeroidx+1000*sequence_onsets(i,j)],[i-0.5 i+0.5],'Color','w'); %sylidx_
        end
    end
    if sort_type ~= 1 
        segtypes = unique(id_flags);
        for segt = 1:numel(segtypes)
            ntrials = find(id_flags == segtypes(segt)); ntrials(1) = ntrials(1)-0.5; ntrials(end) = ntrials(end)+0.5;
            plot(zeros(numel(ntrials),1),ntrials,'Color',...
                colors(find(syllables == segtypes(segt)),:),'LineWidth',10);
        end
    end
    %
    if locktoonset == 1
        xticks([2100 3100]); xticklabels([0 1]);
    else
        xticks([zeroidx-1000 zeroidx]); xticklabels([-1 0]);
    end
end
%set(ax,'CameraTarget', [1 36 0.2605]);

end
function resid = linear_res(vec_a,vec_b)
    % returns the residuals of vec_a after removing the linear component of
    % vec_b
    tmpcov = cov(vec_a,vec_b,1);
    beta = tmpcov(1,2)/var(vec_b,1);
    alpha = mean(vec_a) - beta*mean(vec_b);
    resid = vec_a - alpha - beta*vec_b;
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

