function idx = inactive_frames_idx(Y,thr)
    if ~isfloat(Y)
        Y = double(Y);
    end
    f = squeeze(mean(min(Y(:,:,:),[],1)));
    idx = find(f < thr);
end