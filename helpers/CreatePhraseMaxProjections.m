bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation5_alexa'};
bird3_params = {'lbr3009' 'lbr3009' 'lbr3009_template_4TF' 'lbr3009_annotation_4TF'};
switch birdnum
    case 1
        bird_params = bird1_params;
    case 2
        bird_params = bird2_params;
    case 3
        bird_params = bird3_params;
end
ignore_entries,join_entries

bird_name = bird_params{1}; 
bird_folder_name = bird_params{2}; 
template_file = bird_params{3}; 
annotation_file = bird_params{4}; 

SourceDir = ['/Volumes/home/Data/Imaging/' bird_folder_name '/BaselineSubtractedRawData'];
TargetDir = ['/Users/yardenc/Documents/Experiments/Imaging/CanaryData/' bird_folder_name '/movs/MaxProj/PhraseMaxProj'];

load(fullfile('/Users/yardenc/Documents/Experiments/Imaging/CanaryData',bird_folder_name, 'movs','wavs',annotation_file));

for fnum = 1:numel(keys)
    fname = keys{fnum};
    ignore_locs = find(ismember(elements{fnum}.segType,ignore_entries));
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
end