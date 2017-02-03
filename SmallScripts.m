%% convert a bunch of mat files to video 
d = dir('aligned_reduced*.mat');
for fnum = 2:42
    load(d(fnum).name);
    write_mat_2_moviefile(double(Y),[d(fnum).name(1:end-3) '.mp4'],'MPEG-4',0);
end
%% create a video from many max images of different trials 
input_prefix = '';
%[FileName,PathName,FilterIndex] = uigetfile(fullfile(pwd,[input_prefix '*.mat']),'MultiSelect','on');

res = zeros(240,320,numel(FileName));
for fnum = 1:numel(FileName)
    load(fullfile(PathName,FileName{fnum}));
    res(:,:,fnum) = max(double(Y),[],3);
end
write_mat_2_moviefile(res,['NoAlign.mp4'],'MPEG-4',0);
%% apply transformations from sbxalign2 to all frames
%load('SingleTrial - 2016-06-03 12 38 29.mat')
Y = double(Y);
X = Y(:,:,:);
for i = 1:size(Y,3)
    tform = affine2d(r.T{i});
    X(:,:,i) = imwarp(X(:,:,i),tform,'OutputView',imref2d(size(Y(:,:,1))));
end
%% use inactive_frames_idx to remove initial dark frames and save results
d = dir('SingleTrial*.mat');
idxs=[];
for fnum = 1:numel(d)
    load(d(fnum).name);
    idx = inactive_frames_idx(Y,25);
    Y(:,:,idx) = [];
    outname = ['reduced_' d(fnum).name];
    Ysiz = size(Y)';
    save(outname,'Y','Ysiz','-v7.3');
    idxs = [idxs; idx];
end
%% use sbxalign2 to align frames in a bunch of videos seperately
d = dir('reduced_SingleTrial*.mat');
for fnum = 1:numel(d)
    load(d(fnum).name);
    Y = double(Y);
    r = sbxalign2(Y-imgaussfilt(Y,25),1:size(Y,3));
    for i = 1:size(Y,3)
        tform = affine2d(r.T{i});
        Y(:,:,i) = imwarp(Y(:,:,i),tform,'OutputView',imref2d(size(Y(:,:,i))));
    end
    save(['aligned_' d(fnum).name], 'Y','Ysiz','-v7.3');
end
%%
d = dir('SingleTrial*.mat');
res = zeros(240,320,42);
for i = 1:numel(d)
    load(d(i).name);
    res(:,:,i) = max(double(Y),[],3);
end
r = sbxalign2(res-imgaussfilt(res,100),1:size(res,3));
res = zeros(240,320,42);
for fnum = 1:numel(d)
    load(d(fnum).name);
    tform = affine2d(r.T{fnum});
    Y = imwarp(double(Y),tform,'OutputView',imref2d(size(Y(:,:,1))));
    res(:,:,fnum) = max(Y,[],3);
    save(['aligned_' d(fnum).name], 'Y','Ysiz','-v7.3');
end
write_mat_2_moviefile(res,['reduced_aligned' '.mp4'],'MPEG-4',0);
%%
d = dir('reduced_Single*.mat');
res = zeros(240,320,42);
for i = 1:numel(d)
    load(d(i).name);
    res(:,:,i) = std(double(Y),1,3);
end
[Ytmp,d1]  = AlignSingleTrial_sequential(res);
%%
d = dir('reduced_Single*.mat');
for fnum = 1:41
    load(d(fnum).name);
    Y = uint8(imwarp(double(Y),d1{fnum},'OutputView',imref2d(size(Y(:,:,1)))));
    save(['aligned_' d(fnum).name], 'Y','Ysiz','-v7.3');
end