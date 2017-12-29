% set of scripts to clean an annotation file
MinPhraseDuration = 0.05;

DIR = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav';
annotation_file = 'lrb85315auto_annotation5_fix.mat';
template_file = 'lrb85315template.mat';

cd(DIR);
load(annotation_file);
load(template_file);

syllables = [templates.wavs.segType];

for fnum = 1:numel(keys)
    syllables = unique([syllables unique(elements{fnum}.segType)']);
end

addpath(genpath('/Users/yardenc/Documents/GitHub/VideoAnalysisPipeline'),'-end');
%% 1. Find unlabeled segments
clc;
display('unlabeled:');
for fnum = 1:numel(keys)
    if ismember(-1,elements{fnum}.segType)
        display(keys{fnum})
    end
end

%% 2. Find short phrases
clc;
display(['Short phrases (MinPhraseDuration =  ' num2str(MinPhraseDuration)]);
for fnum = 1:numel(keys)
    phrases = return_phrase_times(elements{fnum});
    durations = phrases.phraseFileEndTimes - phrases.phraseFileStartTimes;
    if any(durations < MinPhraseDuration)
        display([keys{fnum} ' | ' num2str(find(durations < MinPhraseDuration)') ' ||| ' ...
            num2str(phrases.phraseType(durations < MinPhraseDuration)')]);
    end
end

%% 3. Find phrase breaks by potentially mislabelled syllables
clc;
display(['Enclaved:']);
for fnum = 1:numel(keys)
    phrases = return_phrase_times(elements{fnum});
    for phrasenum = 2:numel(phrases.phraseType)-1
        if (phrases.phraseType(phrasenum - 1) == phrases.phraseType(phrasenum + 1))
            display([keys{fnum} ' | ' num2str(phrasenum) ' ||| ' ...
                num2str(phrases.phraseType(phrasenum))]);
        end
    end
end
