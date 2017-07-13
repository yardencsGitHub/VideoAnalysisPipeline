function [rasters,numreps] = SingleDayPhraseAlignment(DIR,date_string,annotation_file,template_file)
addpath(genpath('/Users/yardenc/Documents/GitHub/VideoAnalysisPipeline'));
load(fullfile(DIR,template_file));
syllables = [templates.wavs.segType];
load(annotation_file);

days = [];
for fnum = 1:numel(keys)
    tokens = regexp(keys{fnum},'_','split');
    days = [days; datenum([tokens{3} '-' tokens{4} '-' tokens{5}])];
end
locs_to_process = find(days == datenum(date_string));
cd ROIdata;
cd(date_string);
load(['ROI_' date_string '.mat']);
numROIs = numel(ROI.stats);
for i = 1:numel(syllables)   
    for j = 1:numROIs
        rasters(i,j).data = {};
    end
end
numreps = zeros(numel(syllables),numROIs);

all_phrase_types = [];

for fnum = 1:numel(locs_to_process)
    filename = ['ROIdata_' keys{locs_to_process(fnum)}(1:end-3) 'mat'];
    load(filename);
    phrases = return_phrase_times(elements{locs_to_process(fnum)});
    for phrasenum = 1:numel(phrases.phraseType)
        if ~ismember(phrases.phraseType(phrasenum),[-1 100 101 102 103])
            loc = find(syllables == phrases.phraseType(phrasenum));
            idx = find((vidTimes >= phrases.phraseFileStartTimes(phrasenum)) & ...
            (vidTimes <= phrases.phraseFileEndTimes(phrasenum)));
            for roi = 1:numROIs
                trace = smooth(zscore(dff(roi,:)),3);
                rasters(loc,roi).data = {rasters(loc,roi).data{:} trace(idx)};
                numreps(loc,roi) = numreps(loc,roi) + 1;
            end
        end
        
    end
end
    
    