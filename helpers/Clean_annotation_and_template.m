%%
% This script removes labels, joins labels and sets all tags in preparation for the ML annotation
cd('/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lbr3009/movs/wav');
load('lbr3009_template');
load('lbr3009_annotation');
tags_to_ignore = [];
tags_to_join = {[9 23] [11 17] [14 25] [8 26]};

orig_tags = sort([templates.wavs.segType]);
target_tags = orig_tags;
for tag_cnt = 1:numel(orig_tags)
    a=find(cell2mat(cellfun(@(x)ismember(target_tags(tag_cnt),x),tags_to_join,'UniformOutput',0)));
    if ~isempty(a)
         target_tags(tag_cnt) = tags_to_join{a}(1);
    end
end
unique_tags = unique(target_tags);

new_keys = {};
new_elements = {};
cnt = 1;
for fnum = 1:numel(keys)
    if ~isempty(elements{fnum}.segType)
        element = elements{fnum};
        for sylidx = 1:numel(element.segType)
            if ismember(element.segType(sylidx),tags_to_ignore)
                element.segType(sylidx) = [];
                element.segAbsStartTimes(sylidx) = [];
                element.segFileStartTimes(sylidx) = [];
                element.segFileEndTimes(sylidx) = [];
            else
                a=find(cell2mat(cellfun(@(x)ismember(element.segType(sylidx),x),tags_to_join,'UniformOutput',0))); 
                if ~isempty(a)
                    element.segType(sylidx) = tags_to_join{a}(1);
                end
            end
            element.segType(sylidx) = find(unique_tags == target_tags(orig_tags == element.segType(sylidx)));
        end
        new_keys{cnt} = keys{fnum};
        new_elements{cnt} = element;
        cnt = cnt + 1;
    end
end
%%
new_templates.wavs = [];
for idx = 1:numel(tags_to_join)
    tags_to_ignore = [tags_to_ignore tags_to_join{idx}(2:end)];
end
for idx = 1:numel(templates.wavs)
    if ~ismember(templates.wavs(idx).segType,tags_to_ignore)
        new_templates.wavs = [new_templates.wavs templates.wavs(idx)];
        new_templates.wavs(end).segType = find(unique_tags == target_tags(orig_tags == templates.wavs(idx).segType));
    end
end
%%
templates = new_templates;
elements = new_elements;
keys = new_keys;

save lbr3009_template_4TF templates;
save lbr3009_annotation_4TF keys elements;