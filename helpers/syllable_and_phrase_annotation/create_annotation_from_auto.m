
function [keys, elements] = create_annotation_from_auto(DIR,old_annotation_file,template_file,auto_file_name)
    cd(DIR);
    MinSylDuration = 0.005; % minimal syllable duration = 8 mSec
    params = load(old_annotation_file);
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
    load(fullfile(DIR,auto_file_name));
    load(template_file);
    syllables = [templates.wavs.segType];
    num_files = numel(keys);
    elements = {};
    tempkeys = {};
    dt = 1/3.692307692307692e+02;
    % specific for this bird
    % trill_syllables = [0:2 4 5 8 9 200 203 208 209 300:306 308 309];
    cnt = 1;
    for fnum = 1:num_files  
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
            temp_segType = syllables(y);  
            elements{cnt} = base_struct;
            elements{cnt}.filenum = temp{2};
            elements{cnt}.segFileStartTimes = (syl_onset - 1) * dt;
            elements{cnt}.segAbsStartTimes = time + elements{cnt}.segFileStartTimes/(24*60*60);
            elements{cnt}.segFileEndTimes = (syl_offset - 1) * dt;
            elements{cnt}.segType = syllables(y)';         
            tempkeys{cnt} = [keys{fnum}(1:end-3) 'wav'];
            cnt = cnt + 1;

        end
        

    end
    keys = tempkeys;
end
    
    
    