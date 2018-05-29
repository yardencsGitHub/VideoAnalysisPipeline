function [seg_dur, seg_phr_num] = collect_sustained_transient_lengths(sustained,daynums)
seg_dur=[];
seg_phr_num = [];
for daynum_cnt = 1:numel(daynums)
    daynum = daynums(daynum_cnt);
    disp(daynum);
    Day = sustained(daynum).Day;
    sus_num = cellfun(@(x)size(x,1),sustained(daynum).hits);
    sus_loc = find(sus_num > 0);

    for roi_cnt = 1:numel(sus_loc)
        roi_n = sus_loc(roi_cnt);
        try
            seg_phr_num =  [seg_phr_num; arrayfun(@(a)a{1}(1),sustained(daynum).sus{roi_n})'];
        catch em
            'emm';
        end
        seg_dur = [seg_dur diff(sustained(daynum).hits{roi_n}(:,2:3)')];
    end
end