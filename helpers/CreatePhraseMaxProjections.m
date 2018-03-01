bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation5_alexa'};
bird3_params = {'lbr3009' 'lbr3009' 'lbr3009_template_4TF' 'lbr3009_annotation_4TF'};
hvc_offset = 0.04;
max_phrase_gap = 0.5;

birdnum = 1;
switch birdnum
    case 1
        bird_params = bird1_params;
        ignore_dates = {'2017_04_19'};
        ignore_entries = [-1 100 102 101 103 202 406 408 409 402 403];
        join_entries = {[207 307 407] [404 405] [208 209] [200 309]};
    case 2
        bird_params = bird2_params;
        %'2017_04_05' '2017_04_06'  '2017_04_11' '2017_04_12' '2017_04_16'
        %'2017_04_20' '2017_04_21' '2017_04_23' '2017_04_25' '2017_04_26' '2017_04_30' '2017_05_03'
        ignore_dates = {'2017_04_14' '2017_04_27' };
        ignore_entries = [-1 100 102 101 103];
        join_entries = {[2 206]}; %{[207 307 407] [404 405] [208 209] [200 309]};
    case 3
        bird_params = bird3_params;
        ignore_dates = {'2017_04_27' '2017_04_28' '2017_06_05' '2017_06_06' '2017_06_07' '2017_06_08' '2017_06_09' '2017_06_12' '2017_06_13' ... 
            '2017_06_14' '2017_06_15' '2017_06_16' '2017_06_19' '2017_06_20' '2017_06_21' '2017_06_22' '2017_06_27' ...
            '2017_06_28' '2017_06_29' '2017_06_30' '2017_07_03' '2017_07_04' '2017_07_06' '2017_07_07' '2017_07_10' ...
            '2017_07_11' '2017_07_12' '2017_07_13' '2017_07_14' '2017_07_18' '2017_07_19' '2017_07_20' '2017_07_21'};
        ignore_entries = [-1 100 102 101 103];
        join_entries = {};
end


bird_name = bird_params{1}; 
bird_folder_name = bird_params{2}; 
template_file = bird_params{3}; 
annotation_file = bird_params{4}; 

SourceDir = ['/Volumes/home/Data/Imaging/' bird_folder_name '/BaselineSubtractedRawData'];
TargetDir = ['/Users/yardenc/Documents/Experiments/Imaging/CanaryData/' bird_folder_name '/movs/MaxProj/PhraseMaxProj'];

load(fullfile('/Users/yardenc/Documents/Experiments/Imaging/CanaryData/',bird_folder_name, 'movs','wav',annotation_file));

for fnum = 1:numel(keys)
    fname = keys{fnum};
    tokens = regexp(keys{fnum},'_','split');
    Daystr =  char(join(tokens(3:5),'_'));
    if ismember(Daystr,ignore_dates)
        continue;
    end
    ignore_locs = find(ismember(elements{fnum}.segType,ignore_entries));
    elements{fnum}.segAbsStartTimes(ignore_locs) = [];
    elements{fnum}.segFileStartTimes(ignore_locs) = [];
    elements{fnum}.segFileEndTimes(ignore_locs) = [];
    elements{fnum}.segType(ignore_locs) = [];  
    for i = 1:numel(join_entries)
        join_locs = find(ismember(elements{fnum}.segType,join_entries{i}));
        elements{fnum}.segType(join_locs) = join_entries{i}(1);
    end
    phrases = return_phrase_times(elements{fnum});
    try   
        phrases = deal_with_time_gaps(phrases,max_phrase_gap);
    catch em
        'f';
    end
    source_mat = fullfile(SourceDir,['BGless_RawData_' fname(1:end-3) 'mat']);
    load(source_mat);
    t = vidTimes + hvc_offset;
    phrases.maxImages = [];
    for phrasenum = 1:numel(phrases.phraseType)
        frame_numbers = find(t >= phrases.phraseFileStartTimes(phrasenum) & t <= phrases.phraseFileEndTimes(phrasenum));
        if numel(frame_numbers) == 0
            phrases.maxImages = cat(3,phrases.maxImages, nan*zeros(480,640));
        elseif numel(frame_numbers) == 1
            phrases.maxImages = cat(3,phrases.maxImages, BGlessVidMat(:,:,frame_numbers));
        else
            phrases.maxImages = cat(3,phrases.maxImages, max(BGlessVidMat(:,:,frame_numbers),[],3));
        end
    end
    target_fname = ['PhraseMaxProj_' fname(1:end-3) 'mat'];
    save(fullfile(TargetDir,target_fname),'phrases');
    display(['done: ' fullfile(TargetDir,target_fname)]);
end