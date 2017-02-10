function RemoveInactiveFrames(mov_list,prefix,thr)
% must be run at the folder with the files in mov_list
if exist('reduced') ~= 7
    mkdir('reduced')
end
if isempty(thr)
    thr = 25;
end
for fnum = 1:numel(mov_list)
    params = load(mov_list{fnum});
    idx = inactive_frames_idx(params.vidMat,thr);
    if (idx(end) >= size(params.vidMat,3))
        error('All dark frames.. change threshold or check what is wrong');
    else
        maxidx = max(idx);
        delta_T = params.vidTimes(maxidx + 1);
        aud_idx = ceil(params.Aud.rate * delta_T);
        params.Aud.nrFrames = params.Aud.nrFrames - aud_idx;
        params.Aud.data(1:aud_idx) = [];
        params.Aud.TotalDuration = params.Aud.TotalDuration - delta_T;
        params.vidMat(:,:,1:maxidx) = [];
        params.vidTimes(idx) = [];
        params.vidTimes = params.vidTimes - params.vidTimes(1);
    end
    outname = [prefix mov_list{fnum}];
    save(fullfile(pwd,'reduced',outname),'-struct','params','-v7.3');
end
    