function [elements, keys] = update_annotation_from_ML_estimates(path_annotation,path_templates,path_estimates)
MinSylDuration = 0.005; % minimal syllable duration = 8 mSec
%addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/TimOs'));
params = load(path_annotation);
expr = params.elements{1}.exper;
base_struct = struct('exper',expr, ...
                     'filenum',0, ...
                     'segAbsStartTimes',[], ...
                     'segFileStartTimes',[], ...
                     'segFileEndTimes',[], ...
                     'segType',[], ...
                     'fs',48000, ...
                     'drugstatus', 'No Drug', ...
                     'directstatus', 'Undirected');

load(path_estimates);
load(path_templates);
syllables = [templates.wavs.segType];
num_files = numel(keys);
elements = params.elements;
newelements = {};
tempkeys = {};
dt = 1/3.692307692307692e+02;
dt = 1/370.5882;

cnt = 1;
for fnum = 1:num_files  
    oldcnt = find(cellfun(@(x)strcmp(x,[keys{fnum}(1:end-3) 'wav']),params.keys));
    temp = regexp(keys{fnum},'_','split');
    x = estimates{fnum};
    x = [0 x 0];
    syl_onset = find(x(1:end-1) == 0 & x(2:end) ~=0);
    syl_offset = find(x(1:end-1) ~= 0 & x(2:end) ==0);
    if numel(syl_onset) > 0 % if we have any syllables at all

        time = getFileTime(keys{cnt});
        syl_durations = (syl_offset - syl_onset) * dt;
        % remove too short syllables <10mSec
        syl_onset(syl_durations < MinSylDuration) = [];
        syl_offset(syl_durations < MinSylDuration) = [];
        y = zeros(numel(syl_onset),1);
        for sylnum = 1:numel(y)
            y(sylnum) = mode(estimates{fnum}(syl_onset(sylnum):syl_offset(sylnum)-1));
        end
        if ~isempty(oldcnt)
            elements{oldcnt}.segFileStartTimes = (syl_onset - 1) * dt;
            elements{oldcnt}.segAbsStartTimes = time + elements{cnt}.segFileStartTimes/(24*60*60);
            elements{oldcnt}.segFileEndTimes = (syl_offset - 1) * dt;
            elements{oldcnt}.segType = syllables(y)';
        else
            newelements{cnt} = base_struct;
            newelements{cnt}.filenum = temp{2};
            newelements{cnt}.segFileStartTimes = (syl_onset - 1) * dt;
            newelements{cnt}.segAbsStartTimes = time + elements{cnt}.segFileStartTimes/(24*60*60);
            newelements{cnt}.segFileEndTimes = (syl_offset - 1) * dt;
            newelements{cnt}.segType = syllables(y)';

            tempkeys{cnt} = [keys{fnum}(1:end-3) 'wav'];
            cnt = cnt + 1;  
        end
    else
        if isempty(oldcnt)
            newelements{cnt} = base_struct;
            newelements{cnt}.filenum = temp{2};
            newelements{cnt}.segFileStartTimes = [];
            newelements{cnt}.segAbsStartTimes = [];
            newelements{cnt}.segFileEndTimes = [];
            newelements{cnt}.segType = [];

            tempkeys{cnt} = [keys{fnum}(1:end-3) 'wav'];
            cnt = cnt + 1;
        end
    end
end

keys = [params.keys; tempkeys];
elements = [elements; newelements];
end

function time = getFileTime(filename)
    if strcmp(filename(end-2:end),'mat')
        filename = filename(1:end-4);
    end
    strparts = regexp(filename,'_', 'split');

    y1 = str2double(strparts{3});
    m = str2double(strparts{4});
    d = str2double(strparts{5});
    th = str2double(strparts{6});
    tm = str2double(strparts{7});
    try
        ts = strparts{8};
        ts =  str2double(ts(1:end-4));
    catch em
        ts = 0;
    end

    time = datenum(y1,m,d,th,tm,ts);

end
    
    
    