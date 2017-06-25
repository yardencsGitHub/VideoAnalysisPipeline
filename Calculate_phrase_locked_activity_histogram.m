%% Takes ROI signals from CNMF_E results and performs calculations
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/small-utils'));
n_del_frames = 5;
event_thr = 6;
RawDIR = '/Volumes/CanaryData/DATA/lrb853_15/RawData';
DIR = '/Users/yardenc/Documents/Experiments/Imaging/Analyses/CanaryData/lrb853_15/cnmfe';
prefix = 'pruned_';
cd ('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav');
load('lrb85315auto_annotation3.mat');
keys_index = [];
for i=1:numel(keys)
    tokens = regexp(keys{i},'_','split');
    keys_index = [keys_index; str2num(tokens{2})];
end
    
cd(DIR);
file_list = dir([prefix '*.mat']);
Syllables  = [0:9,200:209, 300:309 400:409];
for sylnum = 1:numel(Syllables)
    timing_hists(sylnum).data = [];
end
for fnum = 1:numel(file_list)
    tokens = regexp(file_list(fnum).name,'_','split');
    RawName = fullfile(RawDIR,['RawData_' file_list(fnum).name(15:end)]);
    load(RawName,'vidTimes');   
    v_dt = mean(diff(vidTimes));
    vid_shift = n_del_frames*v_dt;
    
    phrases = return_phrase_times(elements{keys_index == str2num(tokens{4})});
    n_phrases = numel(phrases.phraseType);
    phrases.phraseFileStartTimes = phrases.phraseFileStartTimes - vid_shift;
    phrases.phraseFileEndTimes = phrases.phraseFileEndTimes - vid_shift;
    
    cnmfefile = fullfile(DIR,file_list(fnum).name);
    load(cnmfefile,'S','C','C_raw');
    n_neurons = size(S,1);
    
    for ncell = 1:n_neurons
        events = find(S(ncell,:) >= event_thr);
        for ev_n = 1:numel(events)
         
            loc = max(find(phrases.phraseFileStartTimes <= events(ev_n)*v_dt));
            if ~isempty(loc) && ismember(phrases.phraseType(loc),Syllables) && (events(ev_n)*v_dt <= phrases.phraseFileEndTimes(loc))
                timing_hists(Syllables == phrases.phraseType(loc)).data = ...
                    [timing_hists(Syllables == phrases.phraseType(loc)).data; ...
                    (events(ev_n)*v_dt - phrases.phraseFileStartTimes(loc))]; %/(phrases.phraseFileEndTimes(loc)- phrases.phraseFileStartTimes(loc))
            end
        end
    end
    
    
end
            

%%        
for i=1:40
figure; hist(timing_hists(i).data,20); %0:0.05:1
title(Syllables(i));
end
    
    
    % raster?
    
    % build time histogram
    